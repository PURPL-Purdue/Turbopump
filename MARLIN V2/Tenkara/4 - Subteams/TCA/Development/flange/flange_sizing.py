## Gasketed Flange Sizing code for Integral Flanges
# Author: Louis DeSano
# Date: 08/29/2026

# Sources: ASME BPVC.VIII-1-2021 Mandatory Appendix 2
"""
Integral Type Flanges. This type covers designs
where the flange is cast or forged integrally with the noz-
zle neck, vessel or pipe wall, butt welded thereto, or at-
tached by other forms of welding of such a nature that
the flange and nozzle neck, vessel or pipe wall is consid-
ered to be the equivalent of an integral structure. In
welded construction, the nozzle neck, vessel, or pipe wall
is considered to act as a hub.

CHANGE LOG:
  1. fixed radius/diameter mixing bug in gasket/bolt-circle geometry
  2. fixed gasket-seating design load: W = (Am + Ab)*Sa/2
  3. safetyFactor now applied to Sa/Sb/Sf before sizing, not just reported after
  4. added guard: F/V/f hardcode only valid for g1 == g0
  5. added Table 2-14 rigidity index (J) check
  6. b0 now derived from physical gasket width N per Table 2-5.2, not input directly
  7. legend moved beside plot; added G/C legend entries
"""

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.patches import Circle
from matplotlib.lines import Line2D
from pathlib import Path
import shutil

from refBolt import bolt_lookup
from margins import margins

# Unit conversions
IN2M = 0.0254
PSI2PA = 6894.76
N2LBF = 0.224809
boltNames = ["#0", "#2", "#4", "#5", "#6", "#8", "#10", "1/4", "5/16", "3/8", "7/16", "1/2",
                "9/16", "5/8", "3/4", "7/8", "1", "1-1/8", "1-1/4", "1-3/8", "1-1/2"]
STANDOFF_M = 0
DIAM_FACTOR = 1.25

def basic_seating_width(N, facing_sketch="1a", column="II", w=None, T=None):
    """Table 2-5.2: basic gasket seating width b0, derived from physical gasket width N."""
    column = column.upper()
    if column not in ("I", "II"):
        raise ValueError('column must be "I" or "II"')

    if facing_sketch in ("1a", "1b"):
        return N / 2
    elif facing_sketch in ("1c", "1d"):
        if w is None or T is None:
            raise ValueError(f"facing_sketch {facing_sketch!r} requires w and T")
        return min((w + T) / 2, (w + N) / 4)
    elif facing_sketch == "2":
        if w is None:
            raise ValueError('facing_sketch "2" requires w (nubbin land width)')
        return (w + N) / 4 if column == "I" else (w + 3 * N) / 8
    elif facing_sketch == "3":
        return N / 4 if column == "I" else 3 * N / 8
    elif facing_sketch == "4":
        return 3 * N / 8 if column == "I" else 7 * N / 16
    elif facing_sketch == "5":
        return N / 4 if column == "I" else 3 * N / 8
    elif facing_sketch == "6":
        if w is None:
            raise ValueError('facing_sketch "6" (ring joint) requires w (ring width)')
        if column != "I":
            raise ValueError("Column II is not defined for facing_sketch \"6\" (ring joint) -- use column I")
        return w / 8
    else:
        raise ValueError(f"Unknown facing_sketch {facing_sketch!r}")


class FlangeSizer:
    def __init__(self, B, P, m, y, N, Sa, Sb, t, g0, g1, Sf, E, safetyFactor,
                 facing_sketch="1a", column="II", w=None, T=None):
        path = Path("Development/flange/figures")
        if path.exists() and path.is_dir():
            shutil.rmtree(path)
        path.mkdir(parents=True, exist_ok=True)
        self.path = path

        self.safetyFactor = safetyFactor

        ### DEFINITIONS ###
        #operating conditions
        self.P = P # internal design pressure [Pa]

        # gasket related
        self.m = m # gasket factor []
        self.y = y # gasket unit seating load [Pa]

        # physical gasket contact width [m]
        self.N = N
        self.facing_sketch = facing_sketch
        self.bo = basic_seating_width(N, facing_sketch, column, w, T) # basic gasket seating width [m]
        bo_mm = self.bo * 1000 #[mm]

        # gasket ID = hub OD, gasket OD = gasket ID + 2*N
        self.gasket_ID = B + (1/8 * IN2M)
        self.gasket_OD = self.gasket_ID + 2 * N

        if bo_mm <= 6:
            b_mm = bo_mm
            self.G = self.gasket_ID + N # diameter at location of gasket load [m]
        else:
            Cb = 2.5 #SI conversion factor (for mm)
            b_mm = Cb * np.sqrt(bo_mm)
            self.G = self.gasket_OD - 2 * (b_mm / 1000) # diameter at location of gasket load [m]

        self.b = b_mm / 1000

        #bolt related
        # raw yield values, kept for margin reporting
        self.Sa_yield = Sa
        self.Sb_yield = Sb
        self.Sf_yield = Sf

        self.Sa = Sa / safetyFactor # Allowable stress for bolt at atmospheric temp [Pa]
        self.Sb = Sb / safetyFactor # Allowable stress for bolt at design temp [Pa]
        self.Sf = Sf / safetyFactor # Allowable stress for flange material [Pa]

        # flange related
        self.t = t # flange thickness [m]
        self.g0 = g0 # thickness of hub at small end [m]
        self.g1 = g1 # thickness of hub at back of flange [m]
        self.B = B # Flange ID
        self.E = E # modulus of elasticity of flange material [Pa]

        # F, V, f hardcoded to g0/g1 = 1 values (Table 2-7.1); only valid for g1 == g0
        if not np.isclose(self.g0, self.g1):
            raise NotImplementedError(
                "g1 != g0: hardcoded F=0.908920, V=0.550103, f=1 are only "
                "valid for uniform-thickness hubs (g1 == g0). Implement the "
                "full Table 2-7.1 polynomial factors to support tapered hubs."
            )

        return

    def solve(self):
        # determine gasket operating load
        self.Wm1 = self.operating_load()

        # determine gasket seating load
        self.Wm2 = self.seating_load()

        #determine bolt area needed for worst case condition
        self.Am = self.bolt_area()
        print(f"Operating Load: {self.Wm1*N2LBF:0.3f} lbf")
        print(f"Seating Load: {self.Wm2*N2LBF:0.3f} lbf")
        print(f"Total Bolt Area: {self.Am/(IN2M**2):0.3f} in^2")

        # find valid bolt configurations
        validBolts = self.bolt_configs()

        # calculate moments and stresses on flange
        anyValid = False
        for config in validBolts:
            Mo = self.flange_moments(config)
            C = config["bolt_circle_diam"]
            A = C + 2 * config["diameter"] * DIAM_FACTOR # flange OD [m]
            [SH, SR, ST], J, validFlag = self.flange_stress(Mo, A)

            # Outputs
            if validFlag:
                anyValid = True
                print(f'{config["count"]}x {config["name"]} ')
                print(f'   (spacing = {config["spacing"] / IN2M:.3f} in)')
                print(f'   Flange OD: {A / IN2M:.3f}')
                print(f'   Bolt Circle: {C / IN2M:.3f}')
                print(int(Mo), int(SH/PSI2PA), int(SR/PSI2PA), int(ST/PSI2PA))
                print(f'   Rigidity index J: {J:.3f} (must be <= 1.0)')
                maxStress = max(SH, SR, ST)
                # SF already baked into Sf; compare vs yield with factor 1
                SFmargin = margins(1, self.Sf_yield, maxStress)
                self.plot_geometry(A, self.B, C, config["count"], config["diameter"], config["name"], SFmargin, J)

        if not anyValid:
            print("No bolt/flange combination satisfies the 2-8(a) stress limits at this thickness/safetyFactor.")
        return

    def operating_load(self):
        """
        Calculates operating load of a gasket
            Arguments:
                G (float) diameter at location of gasket load [m]
                P (float) internal design pressure [Pa]
                b (float) effective gasket seating width [m]
                m (float) gasket factor []

            Return:
                Wm1 (float) Gasket operating load [N]
        """
        self.H = 0.785 * self.G**2 * self.P # Hydrostatic End Force
        self.Hp = 2 * self.b * np.pi * self.G * self.m * self.P # Joint Contact Force
        Wm1 = self.H + self.Hp
        return Wm1

    def seating_load(self):
        """
        Calculates seating load of a gasket
            Arguments:
                G (float) diameter at location of gasket load [m]
                b (float) effective gasket seating width [m]
                y (float) gasket unit seating load [Pa]

            Return:
                Wm2 (float) Gasket seating load [N]
        """
        Wm2 = np.pi * self.b * self.G * self.y
        return Wm2

    def bolt_area(self):
        """
        Calculates required bolt area for gasket operation
            Arguments:
                Wm1 (float) Gasket operating load [N]
                Sb (float) Allowable stress for bolt at design temp [Pa]
                Wm2 (float) Gasket seating load [N]
                Sa (float) Allowable stress for bolt at atmospheric temp [Pa]


            Return:
                Am (float) required bolt area [m^2]
        """
        Am1 = self.Wm1 / self.Sb # Bolt area required for operating stress
        Am2 = self.Wm2 / self.Sa # Bolt area required for seating stress
        Am = max(Am1, Am2)
        if Am1 > Am2:
            print("Bolt Operating Stress is driving.")
        else:
            print("Bolt Seating Stress is driving.")

        return Am

    def bolt_spacing(self, a):
        """
        Calculates maximum bolt spacing for "lethal service" -> Ensures gasket seats evenly
        NOTE: per 2-6, this correction is technically only mandatory for
        lethal service or when specified by the user/designated agent; it
        is applied unconditionally here as a conservative default.
            Arguments:
                a (float) nominal bolt diameter [m]
                t (float) flange thickness [m]
                m (float) gasket factor []

            Return:
                Bs (float) maximum bolt spacing [m]
        """
        Bs = 2 * a + (6 * self.t) / (self.m + 0.5)
        return Bs

    def bolt_configs(self):
        """Find valid bolt size/count configurations."""

        validBolts = []

        print("Possible Bolt Configurations:")

        for boltName in boltNames:

            # Get UNF properties
            boltArea = bolt_lookup("area", boltName, "name", "UNF") * IN2M**2

            boltDiam = bolt_lookup("diameter", boltName, "name", "UNF") * IN2M

            try:
                bolt_clearance_diam = bolt_lookup("clearance_close", boltName, "name", "UNF") * IN2M
            except Exception:
                bolt_clearance_diam = boltDiam + 0.0156 * IN2M

            # Maximum allowable spacing
            maxSpacing = self.bolt_spacing(boltDiam)

            # bolt circle: gasket OD + clearance margin
            C = self.gasket_OD + 2*bolt_clearance_diam #2x since we are adding to the diameter
            circum = np.pi * C

            # Number required by bolt area
            count_area = int(np.ceil(self.Am / boltArea))

            # Number required by maximum bolt spacing
            count_spacing = int(np.ceil(circum / maxSpacing))

            # Must satisfy both requirements
            count = max(count_area, count_spacing)

            # Actual spacing with selected number of bolts
            actualSpacing = circum / count


            if actualSpacing > maxSpacing:
                continue
            if actualSpacing < (1.5*boltDiam):
                continue
            validBolts.append({
                "name": boltName,
                "count": count,
                "area": boltArea,
                "diameter": boltDiam,
                "spacing": actualSpacing,
                "max_spacing": maxSpacing,
                "bolt_circle_diam": C
            })

        return validBolts

    def flange_moments(self, config):
        a = config["diameter"]
        C = config["bolt_circle_diam"]
        Bs = config["spacing"]
        boltArea = config["area"]
        count = config["count"]

        R = (C - self.B) / 2 - self.g1 # radial distance from bolt circle to point of inter-section of hub and back of flange

        HD = 0.785 * (self.B**2) * self.P # hydrostatic end force on area inside of flange
        hD = R + 0.5 * self.g1 # radial distance from the bolt circle, to the circle on which HD acts, as prescribed in Table 2-6
        MD = HD * hD

        HG = self.Wm1 - self.H # Operating condition gasket load
        hG = (C - self.G) / 2 #radial distance from gasket load reaction to the bolt circle
        MG = HG * hG

        HT = self.H - HD
        hT = (R + self.g1 + hG) / 2 #radial distance from the bolt circle to the circle on which HT acts as prescribed in Table 2-6
        MT = HT * hT

        # seating moment: eq. (5) -- W = (Am + Ab)*Sa / 2, NOT Ab*Sa.
        Ab = boltArea * count # total area of bolts actually being used [m^2]
        W = (self.Am + Ab) * self.Sa / 2 # flange design bolt load for gasket seating [N]
        Mseating = W * (C - self.G) / 2

        # operating condition moment
        Moperating = MD + MT + MG
        Mo = max(Mseating, Moperating)

        # correction factor when bolt spacing, Bs, is > 2a + t
        if Bs > (2 * a + self.t):
            BSC = (Bs / (2 * a + self.t))**0.5
        else:
            BSC = 1

        Mo *= BSC
        return Mo

    def flange_stress(self, Mo, A):
        """
            Calculates stresses on an integral type flange based on
            Arguments:
                B (float) flange ID [m]
                t (float) flange thickness [m]
                f (float) hub stress correction factor for integral flanges []
                Mo (float) total moment acting on flange (greater of seating and operating moments) [N*m]
                g0 (float) thickness of hub at small end [m]
                g1 (float) thickness of hub at back of flange [m]

                Y (float) factor involving K (Fig 2-7.1)
                Z (float) factor involving K (Fig 2-7.1)
                T (float) factor involving K (Fig 2-7.1)
                U (float) factor involving K (Fig 2-7.1)
                F (float) factor for integral type flanges (Fig 2-7.2)
                V (float) factor for integral type flanges (Fig 2-7.3)

            Return:
                [SH, SR, ST] (list of float) flange stresses [Pa]
                J (float) Table 2-14 rigidity index []
                validFlag (bool) whether stresses satisfy 2-8(a)
        """
        K = A / self.B
        T = (K**2 * (1 + 8.55246 * np.log10(K)) - 1) / ((1.04720 + 1.9448 * K**2) * (K - 1))
        U = (K**2 * (1 + 8.55246 * np.log10(K)) - 1) / (1.36136 * (K**2 - 1) * (K - 1))
        Y = 1 / (K - 1) * (0.66845 + 5.71690 * (K**2 * np.log10(K)) / (K**2 - 1))
        Z = (K**2 + 1) / (K**2 - 1)

        # values for g0/g1 = 1 (validity enforced in __init__)
        F = 0.908920
        V = 0.550103
        f = 1

        ho = (self.B * self.g0) ** 0.5
        e = F / ho
        d = (U / V) * ho * self.g0**2
        L = (self.t * e + 1) / T + (self.t**3) / d

        SH = (f * Mo) / (L * self.g1**2 * self.B)
        SR = (1.33 * self.t * e + 1) * Mo / (L * self.t**2 * self.B)
        ST = (Y * Mo) / (self.t**2 * self.B) - Z * SR

        # Table 2-14 rigidity index for integral-type flanges: KI = 0.3
        KI = 0.3
        J = (52.14 * V * Mo) / (L * self.E * self.g0**2 * KI * ho)

        #check if stresses allow for valid configuration
        validFlag = True
        for stress in [SR, ST]:
            if stress > self.Sf:
                validFlag = False
                break
        if SH > 1.5 * self.Sf:
            validFlag = False
        if (SH + SR)/2 > self.Sf:
            validFlag = False
        if (SH + ST)/2 > self.Sf:
            validFlag = False

        return [SH, SR, ST], J, validFlag

    def plot_geometry(self, A, B, C, count, boltDiam, boltName, SFmargin, J):
        fig, ax = plt.subplots(figsize=(8, 8))

        colorB = "black"
        colorA = "dimgray"
        colorC = "tab:purple"
        colorG = "tab:green"
        colorGasket = "tab:orange"
        colorBolt = "tab:blue"

        # Bolts
        for i in range(count):
            theta = 2 * np.pi * i / count

            x = (C / 2) * np.cos(theta)
            y = (C / 2) * np.sin(theta)
            ax.add_patch(Circle((x, y), boltDiam / 2, fill=False, color=colorBolt, linewidth=1.5))

        ax.add_patch(Circle((0, 0), B / 2, fill=False, color=colorB, linewidth=2))
        ax.add_patch(Circle((0, 0), self.G / 2, fill=False, color=colorG, linewidth=2, linestyle="--"))
        ax.add_patch(Circle((0, 0), C / 2, fill=False, color=colorC, linewidth=2, linestyle="-."))
        ax.add_patch(Circle((0, 0), A / 2, fill=False, color=colorA, linewidth=2))

        # Gasket ID/OD
        ax.add_patch(Circle((0, 0), self.gasket_ID / 2, fill=False, color=colorGasket, linewidth=2))
        ax.add_patch(Circle((0, 0), self.gasket_OD / 2, fill=False, color=colorGasket, linewidth=2))

        ax.set_aspect("equal")
        ax.set_xlabel("x [m]")
        ax.set_ylabel("y [m]")
        ax.set_title(f"{count}x Bolt Flange Geometry")

        margin = 0.1 * A
        ax.set_xlim(-(A / 2 + margin), A / 2 + margin)
        ax.set_ylim(-(A / 2 + margin), A / 2 + margin)

        ax.grid(True, alpha=0.25)

        legend_handles = [
            Line2D([0], [0], color=colorB, lw=2, label=f"B: flange ID ({B / IN2M:.3f} in)"),
            Line2D([0], [0], color=colorGasket, lw=2, label=f"Gasket ID/OD (width N = {self.N / IN2M:.3f} in)"),
            Line2D([0], [0], color=colorG, lw=2, ls="--", label=f"G: gasket reaction diameter ({self.G / IN2M:.3f} in)"),
            Line2D([0], [0], color=colorC, lw=2, ls="-.", label=f"C: bolt circle ({C / IN2M:.3f} in)"),
            Line2D([0], [0], color=colorA, lw=2, label=f"A: flange OD ({A / IN2M:.3f} in)"),
            Line2D([0], [0], color=colorBolt, lw=2, label=f"Bolt: {boltName} x {count}"),
            Line2D([0], [0], color="none", label=f"Basic seating width b0 (facing {self.facing_sketch}): {self.bo / IN2M:.3f} in"),
            Line2D([0], [0], color="none", label=f"Flange thickness: {self.t/IN2M:.2f} in"),
            Line2D([0], [0], color="none", label=f"Stress margin vs yield (SF={self.safetyFactor}): {SFmargin:.2f}%"),
            Line2D([0], [0], color="none", label=f"Rigidity index J: {J:.3f}"),
        ]

        ax.legend(handles=legend_handles, loc="upper left", bbox_to_anchor=(1.02, 1.0),
                  borderaxespad=0.0, frameon=True)
        safeBoltName = boltName.replace("/", "_")
        plt.savefig(self.path / f"{safeBoltName}x{count}.png", bbox_inches="tight")
        plt.close(fig)
        return

def main():

    # Flange Geometry
    B = 3.826 * IN2M # Flange ID (Chamber ID) [m]
    # Flange Thicknesses
    g0 = 0.337 * IN2M # [m]
    g1 = g0 #equal for straight integral flange

    # working chamber pressure
    P = 440 * PSI2PA # [Pa]

    # gasket properties: vermiculite vermiculite with SS insert
    m = 2.0 # gasket factor []
    y = 2500  * PSI2PA # design seating stress [Pa]

    # physical gasket contact width, facing sketch (1a) per Table 2-5.2
    N = 3/8 * IN2M # [m]
    facing_sketch = "1a"
    column = "II"

    # bolt yield stress
    Sa = 30000 * PSI2PA # Yield stress for bolt at atmospheric temp [Pa]
    Sb = Sa # Allowable stress for bolt at design temp [Pa]

    t = 0.5 * IN2M # flange thickness [m]
    Sf = 36000 * PSI2PA # flange yield at temp [Pa]

    # modulus of elasticity of flange material at design temperature
    E = 29e6 * PSI2PA # [Pa]

    flange = FlangeSizer(B, P, m, y, N, Sa, Sb, t, g0, g1, Sf, E, safetyFactor=1.5,
                          facing_sketch=facing_sketch, column=column)
    flange.solve()
    return

if __name__ == "__main__":
    main()
