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


class FlangeSizer:
    def __init__(self, B, P, m, y, bo, Sa, Sb, t, g0, g1, Sf, safetyFactor):
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

        # effective gasket seating width [m]
        self.bo = bo
        bo_mm = bo * 1000 #[mm]

        if bo_mm <= 6:
            b_mm = bo_mm
            #calculate gasket load location as mean diameter of gasket contact face
            self.G = B + (bo / 2) # diameter at location of gasket load [m]
        else:
            Cb = 2.5 #SI conversion factor (for mm)
            b_mm = Cb * np.sqrt(bo_mm)
            #gasket load location as gasket OD - 2b
            self.G = (B + bo) -  2 * (b_mm/1000)# diameter at location of gasket load [m]

        self.b = b_mm / 1000
        
        #bolt related
        self.Sa = Sa / safetyFactor # Allowable stress for bolt at atmospheric temp [Pa]
        self.Sb = Sb / safetyFactor# Allowable stress for bolt at design temp [Pa]

        # flange related
        self.t = t # flange thickness [m]
        self.g0 = g0 # thickness of hub at small end [m]
        self.g1 = g1 # thickness of hub at back of flange [m]
        self.B = B # Flange ID
        self.Sf = Sf
        
        return
    
    def solve(self):
        # determine gasket operating load
        self.Wm1 = self.operating_load()

        # determine gasket seating load
        self.Wm2 = self.seating_load()
        
        #determine bolt area needed for worst case condition
        self.Am = self.bolt_area()
        print(f"Operating Stress: {self.Wm1*N2LBF:0.3f} lbf")
        print(f"Seating Stress: {self.Wm2*N2LBF:0.3f} lbf")
        print(f"Total Bolt Area: {self.Am/(IN2M**2):0.3f} in^2")

        # find valid bolt configurations
        validBolts = self.bolt_configs()

        # calculate moments and stresses on flange
        for config in validBolts:
            Mo = self.flange_moments(config)
            C = config["bolt_circle_diam"]
            A = C + config["diameter"]*1.5 # flange OD [m]
            [SH, SR, ST], validFlag = self.flange_stress(Mo, A)

            # Outputs
            if validFlag:
                print(f"{config["count"]}x {config["name"]} ")
                print(f"   (spacing = {config["spacing"] / IN2M:.3f} in)")
                print(f"   Flange OD: {A / IN2M:.3f}")
                print(f"   Bolt Circle: {C / IN2M:.3f}")
                print(int(Mo), int(SH/PSI2PA), int(SR/PSI2PA), int(ST/PSI2PA))
                maxStress = max(SH, SR, ST)
                SFmargin = margins(1.5, self.Sf, maxStress)
                self.plot_geometry(A, self.B, C, config["count"], config["diameter"], config["name"], self.g0, self.bo, SFmargin)
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
            except:
                bolt_clearance_diam = boltDiam + 0.0156 * IN2M
            
            # Maximum allowable spacing
            maxSpacing = self.bolt_spacing(boltDiam)

            # Determine bolt circle circumference
            C = self.B + self.g0 + self.bo + bolt_clearance_diam*1.5
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
        
        #seating moment
        Ab = boltArea * count # total area of bolts actually being used [m^2]
        W = Ab * self.Sa   # full allowable bolt load at assembly temperature [N]
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
        """
        K = A / self.B
        T = (K**2 * (1 + 8.55246 * np.log10(K)) - 1) / ((1.04720 + 1.9448 * K**2) * (K - 1))
        U = (K**2 * (1 + 8.55246 * np.log10(K)) - 1) / (1.36136 * (K**2 - 1) * (K - 1))
        Y = 1 / (K - 1) * (0.66845 + 5.71690 * (K**2 * np.log10(K)) / (K**2 - 1))
        Z = (K**2 + 1) / (K**2 - 1)
        
        # values for g0/g1 = 1
        F = 0.908920 # g0/g1 = 1 
        V = 0.550103
        f = 1

        ho = (self.B * self.g0) ** 0.5
        e = F / ho
        d = (U / V) * ho * self.g0**2
        L = (self.t * e + 1) / T + (self.t**3) / d

        SH = (f * Mo) / (L * self.g1**2 * self.B)
        SR = (1.33 * self.t * e + 1) * Mo / (L * self.t**2 * self.B)
        ST = (Y * Mo) / (self.t**2 * self.B) - Z * SR

        #check if stresses allow for valid configuration 
        validFlag = True
        for stress in [SH, SR, ST]:
            if stress > self.Sf:
                validFlag = False
                break
        if (SH + SR)/2 > self.Sf:
            validFlag = False
        if (SH + ST)/2 > self.Sf:
            validFlag = False
        
        return [SH, SR, ST], validFlag

    def plot_geometry(self, A, B, C, count, boltDiam, boltName, g0, bo, SFmargin):
        fig, ax = plt.subplots(figsize=(8, 8))

        # Bolts
        for i in range(count):
            theta = 2 * np.pi * i / count

            x = (C / 2) * np.cos(theta)
            y = (C / 2) * np.sin(theta)
            # Bolts
            ax.add_patch(Circle((x, y), boltDiam / 2, fill=False, color="tab:blue", linewidth=1.5))

        ax.add_patch(Circle((0, 0), B / 2, fill=False, color="black", linewidth=2))
        ax.add_patch(Circle((0, 0), A / 2, fill=False, color="black", linewidth=2))
        ax.add_patch(Circle((0, 0), C / 2, fill=False, color="gray", linewidth=1.5))


        # Gasket
        ax.add_patch(Circle((0, 0), (B + g0) / 2, fill=False, color="tab:orange", linewidth=2))
        ax.add_patch(Circle((0, 0), (B + g0 + bo) / 2, fill=False, color="tab:orange", linewidth=2))

        ax.set_aspect("equal")
        ax.set_xlabel("x [m]")
        ax.set_ylabel("y [m]")
        ax.set_title(f"{count}x Bolt Flange Geometry")

        margin = 0.1 * A
        ax.set_xlim(-(A / 2 + margin), A / 2 + margin)
        ax.set_ylim(-(A / 2 + margin), A / 2 + margin)

        ax.grid(True, alpha=0.25)

        legend_handles = [
            Line2D([0], [0], color="black", lw=2, label=f"Flange ID: {B / IN2M:.3f} in"),
            Line2D([0], [0], color="black", lw=2, label=f"Flange OD: {A / IN2M:.3f} in"),
            Line2D([0], [0], color="tab:blue", lw=2, label=f"Bolt: {boltName} × {count}"),
            Line2D([0], [0], color="tab:orange", lw=2, label=f"Gasket width: {bo / IN2M:.3f} in"),
            Line2D([0], [0], color="none", label=f"Flange thickness: {self.t/IN2M:.2f} in"),
            Line2D([0], [0], color="none", label=f"Stress margin (SF={self.safetyFactor}): {SFmargin:.2f}%"),
        ]
                
        ax.legend(handles=legend_handles, loc="best")
        plt.savefig(self.path / f"{boltName.replace("/", "_")}x{count}.png")
        return

def main():
    
    #Flange Geometry
    B = 3.826 * IN2M # Flange ID (Chamber ID) [m]
    #Flange Thicknesses  
    g0 = 0.337 * IN2M # [m]
    g1 = g0 #equal for straight integral flange

    # proof pressure: 1.5x chamber pressure 
    P = 440 * PSI2PA * 1.5 # [Pa]
    
    # gasket properties: vermiculite (Thermiculute 715 - coreless)
    m = 3.0 # gasket factor []
    y = 1500  * PSI2PA # design seating stress [Pa]

    # gasket seating width
    bo = 1 * IN2M # [m]

    # bolt yield stress
    Sa = 30000 * PSI2PA # Yield stress for bolt at atmospheric temp [Pa]
    Sb = Sa * 0.65 # Allowable stress for bolt at design temp [Pa]

    t = 2 * IN2M # flange thickness [m]
    Sf = 15000 * PSI2PA # flange yield at temp [Pa]
    
    flange = FlangeSizer(B, P, m, y, bo, Sa, Sb, t, g0, g1, Sf, safetyFactor=2)
    flange.solve()
    return 

if __name__ == "__main__":
    main() 