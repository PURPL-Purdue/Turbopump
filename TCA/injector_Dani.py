"""
Injector manifold and impinging doublet sizing tool.

This script computes:
- Propellant mass-flow split (fuel / oxidizer)
- Per-hole mass flow, injection velocities, and injector orifice diameters
- Impinging-jet geometry (fuel / LOX angles, sheet angle, required spacing)
- Injector ring radii/diameters at exit plane
- Injector ring radii/diameters at a back plane one plate thickness upstream
- Simple annular manifold dimensions using a Darcy–Weisbach + K model

Units:
- Pressures: Pa internally (psi converted once).
- Lengths: m internally; prints key diameters in mm.
"""

import numpy as np
import math
import matplotlib.pyplot as plt
from scipy import optimize

# ---------------------------------------------------------------------------
# Unit conversions
# ---------------------------------------------------------------------------
PSI_TO_PA = 6894.76         # 1 psi = 6894.76 Pa
INCH_TO_M = 0.0254          # 1 inch = 0.0254 m
DEG_TO_RAD = math.pi / 180.0
RAD_TO_DEG = 180.0 / math.pi

# ---------------------------------------------------------------------------
# Design parameters (EDIT THESE AS NEEDED)
# ---------------------------------------------------------------------------

# Total mass flow and mixture ratio
mdot_total = 9.0            # total propellant mass flow [kg/s]
OF_RATIO   = 2.1            # O/F = m_dot_ox / m_dot_fuel

# Propellant densities at injector conditions
rho_rp1 = 810.0             # RP-1 density [kg/m^3]
rho_lox = 1141.0            # LOX density [kg/m^3]

# Injector hydraulic parameters
Cd = 0.7                    # discharge coefficient (assumed same for both)

Pc_psi         = 500.0      # chamber static pressure [psi]
Pin_psi        = 850.0      # injector inlet pressure [psi]
Pc             = Pc_psi * PSI_TO_PA
Pin            = Pin_psi * PSI_TO_PA
stiffness      = Pc / Pin   # informative only

# Single, consistent injector pressure drop used everywhere:
deltaP_inj_psi = 350.0                        # injector ΔP [psi]
deltaP_inj     = deltaP_inj_psi * PSI_TO_PA   # injector ΔP [Pa]

# Combustion chamber / injector geometry
Length_chamber_in = 9.928   # chamber length [inch]
Length_chamber    = Length_chamber_in * INCH_TO_M  # [m]

CombDiam_in = 6.336         # chamber inner diameter [inch]
CombDiam    = CombDiam_in * INCH_TO_M              # [m]

# Axial location of impingement plane as fraction of chamber length
impinge_fraction = 0.10     # jets impinge at 10% of chamber length (dimensionless)

# Injector layout parameters (for ring placement)
marginWall  = 0.007          # clearance from outer RP1 jet to wall [m]
pairSpacing = 0.026         # radial spacing between ring mid-radii [m]

# Manifold design parameters
h_manifold = 0.0115         # manifold height (radial depth) [m]
dpFracMani = 0.08           # allowed manifold DP fraction of injector DP

fOx  = 0.02                 # Darcy friction factor for LOX manifold
fRP1 = 0.02                 # Darcy friction factor for RP-1 manifold
KOx  = 1.0                  # lumped local loss K for LOX manifold
KRP1 = 1.0                  # lumped local loss K for RP-1 manifold

NinletsOx  = 1              # number of LOX manifold feed inlets
NinletsRP1 = 1              # number of RP-1 manifold feed inlets

# Number of rings (fuel and ox share two rings in this simple model)
Nrings = 2

# ---------------------------------------------------------------------------
# Helper functions
# ---------------------------------------------------------------------------

def mdot_prop(mdot_propellant: float,
              number_holes: int,
              total_number_of_holes: int):
    """
    Distribute total propellant mass flow among a specific subset of holes.
    """
    mdot_subset    = mdot_propellant * number_holes / total_number_of_holes
    mdot_per_hole  = mdot_subset / number_holes
    return mdot_subset, mdot_per_hole


def diameter(mdot_propellant: float,
             rho: float,
             deltaP_psi: float,
             number_holes: int):
    """
    Compute injector-hole diameter for a given propellant.

    Orifice equation:
        m_dot_total = C_d * A_total * sqrt(2 * rho * ΔP)
    """
    deltaP = deltaP_psi * PSI_TO_PA  # convert once
    A_total = mdot_propellant / (Cd * math.sqrt(2.0 * rho * deltaP))
    A_hole  = A_total / number_holes
    d_hole  = math.sqrt(4.0 * A_hole / math.pi)
    return d_hole, A_hole, A_total


def velocity(Cd_local: float, deltaP_pa: float, rho: float) -> float:
    """Injection velocity through an orifice from ΔP and density."""
    return Cd_local * math.sqrt(2.0 * deltaP_pa / rho)


def max_mixing_ratio(rho_fuel: float,
                     rho_ox: float,
                     mdot_fuel: float,
                     mdot_ox: float) -> float:
    """Simple impinging-jet mixing metric (just informative)."""
    M = 1.0  # for 1-on-1 unlike impingement
    return M * (rho_ox / rho_fuel * (mdot_ox / mdot_fuel) ** 2.0) ** 0.7


def solve_thetas(mdot_fuel: float,
                 mdot_ox: float,
                 v_fuel: float,
                 v_ox: float,
                 theta_fuel_deg: float,
                 beta_target_deg: float,
                 theta_ox_min_deg: float = 5.0,
                 theta_ox_max_deg: float = 80.0,
                 dtheta_deg: float       = 0.05):
    """
    For given fuel angle θ_f and desired spray sheet angle β_target, find the
    LOX angle θ_ox that gives a resultant sheet angle closest to β_target.
    """
    theta_f = theta_fuel_deg * DEG_TO_RAD
    beta_target = beta_target_deg * DEG_TO_RAD

    p_f  = mdot_fuel * v_fuel
    p_ox = mdot_ox * v_ox

    best_err = 1e30
    best_theta_ox_deg = None
    best_beta_deg = None

    theta_ox_deg = theta_ox_min_deg
    while theta_ox_deg <= theta_ox_max_deg + 1e-9:
        theta_ox = theta_ox_deg * DEG_TO_RAD

        # Components (chamber axis = y, radial = x)
        px = p_ox * math.sin(theta_ox) - p_f * math.sin(theta_f)
        py = p_ox * math.cos(theta_ox) + p_f * math.cos(theta_f)

        beta_res = math.atan2(px, py)
        err = abs(beta_res - beta_target)

        if err < best_err:
            best_err = err
            best_theta_ox_deg = theta_ox_deg
            best_beta_deg = beta_res * RAD_TO_DEG

        theta_ox_deg += dtheta_deg

    return best_beta_deg, best_theta_ox_deg


def compute_spacing_doublet(Lc_m: float,
                            impingement_fraction: float,
                            theta_fuel_deg: float,
                            theta_ox_deg: float):
    """
    Compute impingement location and required F–O spacing from angles.
    """
    z_imp = Lc_m * impingement_fraction
    theta_f  = theta_fuel_deg * DEG_TO_RAD
    theta_ox = theta_ox_deg   * DEG_TO_RAD

    d_imp_f  = z_imp * math.tan(theta_f)
    d_imp_ox = z_imp * math.tan(theta_ox)
    d_fo     = d_imp_f + d_imp_ox

    return z_imp, d_imp_f, d_imp_ox, d_fo


# -------- NEW: back-side ring geometry from plate thickness + angle --------

def back_radius_from_exit(R_exit: float,
                          thickness: float,
                          theta_deg: float):
    """
    Given the radius of a jet at the exit plane (chamber side) and the plate
    thickness, compute the radius of the same jet at a back plane that is
    'thickness' upstream, assuming a straight hole at angle theta.

    Geometry (in cross-section):
      Δr = thickness * tan(theta)
      R_back = R_exit - Δr   (since going upstream towards the center)

    Inputs:
      R_exit    -- radius at exit plane [m]
      thickness -- plate thickness [m]
      theta_deg -- injection angle from chamber axis [deg]

    Returns:
      (R_back [m], Δr [m])
    """
    theta = theta_deg * DEG_TO_RAD
    dr = thickness * math.tan(theta)   # radial change over thickness
    R_back = R_exit - dr               # going upstream, radius decreases
    return R_back, dr


def angle_from_rings(R_exit: float,
                     R_back: float,
                     thickness: float):
    """
    Given the radius at exit and at a back plane, plus the axial separation
    (thickness), compute the implied injection angle.

      Δr = R_exit - R_back
      theta = atan(Δr / thickness)

    All in meters, returns angle in degrees.
    """
    dr = R_exit - R_back
    theta = math.atan2(dr, thickness)
    return theta * RAD_TO_DEG


# ---------------------------------------------------------------------------

def design_manifold(mdot: float,
                    rho: float,
                    Lpath: float,
                    dp_inj: float,
                    dp_frac: float,
                    f: float,
                    Ktot: float,
                    h: float,
                    w_init: float):
    """
    Size rectangular manifold cross-section (height fixed, solve for width).

    Darcy–Weisbach:
        ΔP_f = f * (L / Dh) * (ρ v^2 / 2)

    Local losses:
        ΔP_loc = Ktot * (ρ v^2 / 2)

    Choose w such that:
        ΔP_total = ΔP_f + ΔP_loc = dp_frac * dp_inj
    """
    dp_allow = dp_frac * dp_inj

    def dp_residual(w):
        if w <= 0.0:
            return 1e9
        A = w * h
        v = mdot / (rho * A)
        Dh = 2.0 * w * h / (w + h)
        dyn = rho * v * v / 2.0
        dpF = f * (Lpath / Dh) * dyn
        dpLoc = Ktot * dyn
        dpTot = dpF + dpLoc
        return dpTot - dp_allow

    # Bracket
    w_min = w_init / 10.0
    w_max = w_init * 10.0
    f_min = dp_residual(w_min)
    f_max = dp_residual(w_max)

    if f_min * f_max > 0.0:
        for _ in range(6):
            w_min /= 2.0
            w_max *= 2.0
            f_min = dp_residual(w_min)
            f_max = dp_residual(w_max)
            if f_min * f_max <= 0.0:
                break

    if f_min * f_max > 0.0:
        raise RuntimeError("Could not bracket root for manifold width; check inputs.")

    sol = optimize.root_scalar(dp_residual, bracket=[w_min, w_max], method="brentq")
    w = sol.root

    # Final properties
    A  = w * h
    v  = mdot / (rho * A)
    Dh = 2.0 * w * h / (w + h)
    dyn = rho * v * v / 2.0
    dpF  = f * (Lpath / Dh) * dyn
    dpLoc = Ktot * dyn
    dpTot = dpF + dpLoc
    ratio = dpTot / dp_inj

    return w, Dh, A, v, dpF, dpLoc, dpTot, ratio


def print_mani(label: str,
               mani,
               L_path: float,
               dp_inj: float,
               h_man: float):
    """Pretty-print manifold sizing results."""
    w, Dh, A, v, dpF, dpLoc, dpTot, ratio = mani
    print(f"{label}:")
    print(f"  width w             = {w:.6f} m")
    print(f"  height h            = {h_man:.6f} m")
    print(f"  hydraulic diam Dh   = {Dh:.6f} m")
    print(f"  L_path              = {L_path * 1e3:.2f} mm")
    print(f"  flow area A         = {A * 1e6:.2f} mm^2")
    print(f"  mean velocity v     = {v:.3f} m/s")
    print(f"  Δp_friction         = {dpF:.1f} Pa")
    print(f"  Δp_local            = {dpLoc:.1f} Pa")
    print(f"  Δp_total            = {dpTot:.1f} Pa")
    print(f"  (Δp_total / Δp_inj) = {ratio:.3f}")
    print()


# --- Plot helpers ----------------------------------------------------------

def ring_points(radius_mm: float, N: int):
    ang = np.linspace(0.0, 2.0 * math.pi, N, endpoint=False)
    x = radius_mm * np.cos(ang)
    y = radius_mm * np.sin(ang)
    return x, y


def plot_injector_layout(D_c_mm: float,
                         R_ring_f_mm,
                         R_ring_ox_mm,
                         N_hole_ring: int,
                         margin_extra_mm: float = 5.0):
    """Visualize injector layout: chamber wall + rings + hole locations."""
    combRad = D_c_mm / 2.0

    x_f_in,  y_f_in  = ring_points(R_ring_f_mm[0],  N_hole_ring)
    x_f_out, y_f_out = ring_points(R_ring_f_mm[1],  N_hole_ring)
    x_ox_in,  y_ox_in  = ring_points(R_ring_ox_mm[0], N_hole_ring)
    x_ox_out, y_ox_out = ring_points(R_ring_ox_mm[1], N_hole_ring)

    fig, ax = plt.subplots(figsize=(8, 8))

    chamber = plt.Circle((0.0, 0.0), combRad, fill=False, linewidth=2)
    ax.add_artist(chamber)

    for R in R_ring_f_mm:
        ax.add_artist(plt.Circle((0.0, 0.0), R, color="red", fill=False,
                                 linewidth=0.7, alpha=0.4))
    for R in R_ring_ox_mm:
        ax.add_artist(plt.Circle((0.0, 0.0), R, color="blue", fill=False,
                                 linewidth=0.7, alpha=0.4))

    ax.scatter(x_f_in,  y_f_in,  s=20, color="red",     label="RP-1 inner")
    ax.scatter(x_f_out, y_f_out, s=20, color="darkred", label="RP-1 outer")
    ax.scatter(x_ox_in,  y_ox_in,  s=20, color="blue",  label="LOX inner")
    ax.scatter(x_ox_out, y_ox_out, s=20, color="navy",  label="LOX outer")

    ax.set_aspect("equal", "box")
    lim = combRad + margin_extra_mm
    ax.set_xlim(-lim, lim)
    ax.set_ylim(-lim, lim)
    ax.set_xlabel("x [mm]")
    ax.set_ylabel("y [mm]")
    ax.set_title("Injector layout: rings and hole locations")
    ax.legend()
    ax.grid(True)
    plt.tight_layout()
    plt.show()


# ---------------------------------------------------------------------------
# Main script
# ---------------------------------------------------------------------------
if __name__ == "__main__":

    print("=== Injector / Manifold Design Tool ===")
    print(f"Total mass flow mdot_total  = {mdot_total:.3f} kg/s")
    print(f"O/F ratio                   = {OF_RATIO:.3f}")
    print(f"Injector ΔP (all uses)      = {deltaP_inj_psi:.1f} psi ({deltaP_inj:.0f} Pa)")
    print()

    # Split total mass flow
    mdot_kero = mdot_total / (1.0 + OF_RATIO)
    mdot_lox  = mdot_total * OF_RATIO / (1.0 + OF_RATIO)

    print(f"Fuel (RP-1) mass flow       = {mdot_kero:.3f} kg/s")
    print(f"Oxidizer (LOX) mass flow    = {mdot_lox:.3f} kg/s")
    print()

    # Number of holes
    num_holes_rp1_inj = int(input("Number of kerosene (RP-1) holes: "))
    num_holes_lox_inj = int(input("Number of LOX holes: "))

    # Plate thickness for CAD geometry
    plate_thickness_mm = float(input("Injector plate thickness [mm]: "))
    plate_thickness = plate_thickness_mm * 1e-3  # [m]

    # Mass flow per hole
    mdot_rp1_inj, mdot_rp1_per_hole = mdot_prop(mdot_kero,
                                                num_holes_rp1_inj,
                                                num_holes_rp1_inj)
    mdot_lox_inj, mdot_lox_per_hole = mdot_prop(mdot_lox,
                                                num_holes_lox_inj,
                                                num_holes_lox_inj)

    print("=== Mass flow per hole ===")
    print(f"RP-1 per hole ({num_holes_rp1_inj} holes): {mdot_rp1_per_hole:.5f} kg/s")
    print(f"LOX  per hole ({num_holes_lox_inj} holes): {mdot_lox_per_hole:.5f} kg/s")
    print()

    # Injection velocities
    v_rp1 = velocity(Cd, deltaP_inj, rho_rp1)
    v_lox = velocity(Cd, deltaP_inj, rho_lox)

    print("=== Injection velocity ===")
    print(f"RP-1 velocity: {v_rp1:.2f} m/s")
    print(f"LOX  velocity: {v_lox:.2f} m/s")
    print()

    # Injector-hole diameters
    d_rp1, A_rp1_hole, A_rp1_total = diameter(mdot_rp1_inj, rho_rp1, deltaP_inj_psi, num_holes_rp1_inj)
    d_lox,  A_lox_hole,  A_lox_total  = diameter(mdot_lox_inj,  rho_lox,  deltaP_inj_psi, num_holes_lox_inj)

    print("=== Injector-hole diameters ===")
    print(f"RP-1 hole diameter: {d_rp1*1e3:.3f} mm")
    print(f"LOX  hole diameter: {d_lox*1e3:.3f} mm")
    print()

    # Mixing metric
    MME = max_mixing_ratio(rho_rp1, rho_lox, mdot_rp1_per_hole, mdot_lox_per_hole)
    print("=== Max 'mixing metric' (MME) ===")
    print(f"MME ≈ {MME:.3f}")
    print()

    # Angles
    theta_rp1_deg   = float(input("Fuel angle θ_f (deg from chamber axis): "))
    beta_target_deg = float(input("Desired spray sheet angle β (deg): "))

    beta_result_deg, theta_lox_deg = solve_thetas(
        mdot_rp1_per_hole,
        mdot_lox_per_hole,
        v_rp1,
        v_lox,
        theta_rp1_deg,
        beta_target_deg
    )

    print("=== Injection angles (impinging doublet) ===")
    print(f"Fuel angle   θ_f   = {theta_rp1_deg:.3f} deg")
    print(f"LOX angle    θ_ox  = {theta_lox_deg:.3f} deg")
    print(f"Result sheet β_res = {beta_result_deg:.3f} deg")
    print()

    # Impingement geometry
    z_imp, d_imp_f, d_imp_ox, d_fo_req = compute_spacing_doublet(
        Length_chamber,
        impinge_fraction,
        theta_rp1_deg,
        theta_lox_deg
    )

    print("=== Impingement geometry ===")
    print(f"Impingement axial location z_imp = {z_imp*1e3:.2f} mm")
    print(f"Fuel lateral offset at impingement  = {d_imp_f*1e3:.2f} mm")
    print(f"LOX  lateral offset at impingement  = {d_imp_ox*1e3:.2f} mm")
    print(f"Required hole spacing d_FO          = {d_fo_req*1e3:.2f} mm")
    print()

    # -----------------------------------------------------------------------
    # Injector ring radii at exit plane
    # -----------------------------------------------------------------------
    combRad = CombDiam / 2.0
    radius_fuel = d_rp1 / 2.0
    radius_lox  = d_lox  / 2.0

    # FO pair mid-radii (fuel inside, LOX outside)
    Rmid_outer = combRad - marginWall - radius_fuel - d_fo_req / 2.0
    Rmid_inner = Rmid_outer - pairSpacing

    # Actual ring radii at exit plane
    Rf_inner   = Rmid_inner - d_fo_req / 2.0
    Rf_outer   = Rmid_outer - d_fo_req / 2.0
    Rlox_inner = Rmid_inner + d_fo_req / 2.0
    Rlox_outer = Rmid_outer + d_fo_req / 2.0

    Rring_rp1_exit = np.array([Rf_inner, Rf_outer])
    Rring_lox_exit = np.array([Rlox_inner, Rlox_outer])

    print("=== Injector ring diameters at EXIT plane (chamber side) ===")
    print(f"Chamber inner diameter D_c         = {CombDiam*1e3:.2f} mm")
    print(f"RP-1 inner ring D_f_in_exit        = {2*Rf_inner*1e3:.2f} mm")
    print(f"RP-1 outer ring D_f_out_exit       = {2*Rf_outer*1e3:.2f} mm")
    print(f"LOX inner ring D_ox_in_exit        = {2*Rlox_inner*1e3:.2f} mm")
    print(f"LOX outer ring D_ox_out_exit       = {2*Rlox_outer*1e3:.2f} mm")
    print()

    # -----------------------------------------------------------------------
    # NEW: ring radii at BACK plane (one plate thickness upstream)
    # -----------------------------------------------------------------------
    Rf_inner_back,  dr_f_in  = back_radius_from_exit(Rf_inner,  plate_thickness, theta_rp1_deg)
    Rf_outer_back,  dr_f_out = back_radius_from_exit(Rf_outer,  plate_thickness, theta_rp1_deg)
    Rlox_inner_back, dr_ox_in  = back_radius_from_exit(Rlox_inner, plate_thickness, theta_lox_deg)
    Rlox_outer_back, dr_ox_out = back_radius_from_exit(Rlox_outer, plate_thickness, theta_lox_deg)

    Rring_rp1_back = np.array([Rf_inner_back, Rf_outer_back])
    Rring_lox_back = np.array([Rlox_inner_back, Rlox_outer_back])

    print("=== Injector ring diameters at BACK plane (upstream by thickness) ===")
    print(f"Plate thickness t                  = {plate_thickness_mm:.2f} mm")
    print(f"RP-1 inner ring D_f_in_back        = {2*Rf_inner_back*1e3:.2f} mm")
    print(f"RP-1 outer ring D_f_out_back       = {2*Rf_outer_back*1e3:.2f} mm")
    print(f"LOX inner ring D_ox_in_back        = {2*Rlox_inner_back*1e3:.2f} mm")
    print(f"LOX outer ring D_ox_out_back       = {2*Rlox_outer_back*1e3:.2f} mm")
    print()

    print("=== Injector plate distances ===")
    print(f"RP-1  distance       = {2*(Rf_inner-Rf_inner_back)*1e3:.2f} mm")
    print(f"Ox distance       = {2*(Rlox_inner-Rlox_inner_back)*1e3:.2f} mm")

    # Quick consistency check: recover angles from rings + thickness
    theta_f_from_geom_in  = angle_from_rings(Rf_inner,  Rf_inner_back,  plate_thickness)
    theta_f_from_geom_out = angle_from_rings(Rf_outer,  Rf_outer_back,  plate_thickness)
    theta_ox_from_geom_in = angle_from_rings(Rlox_inner, Rlox_inner_back, plate_thickness)
    theta_ox_from_geom_out= angle_from_rings(Rlox_outer, Rlox_outer_back, plate_thickness)

    print("=== Angle check from CAD-style geometry ===")
    print(f"Fuel angle from inner ring pair    ≈ {theta_f_from_geom_in:.3f} deg")
    print(f"Fuel angle from outer ring pair    ≈ {theta_f_from_geom_out:.3f} deg")
    print(f"LOX angle from inner ring pair     ≈ {theta_ox_from_geom_in:.3f} deg")
    print(f"LOX angle from outer ring pair     ≈ {theta_ox_from_geom_out:.3f} deg")
    print("(All should match the design angles within small numerical error.)")
    print()

    # -----------------------------------------------------------------------
    # Manifold path lengths and sizing (using EXIT-plane ring radii)
    # -----------------------------------------------------------------------
    circ_ox_exit  = 2.0 * math.pi * Rring_lox_exit
    circ_rp1_exit = 2.0 * math.pi * Rring_rp1_exit

    Lpath_ox   = circ_ox_exit  / (2.0 * NinletsOx)
    Lpath_rp1  = circ_rp1_exit / (2.0 * NinletsRP1)
    Lpath_ox_avg = 0.5 * (Lpath_ox[0] + Lpath_ox[1])

    w_guess = 0.02   # 2 cm initial width guess

    mdot_lox_per_ring  = mdot_lox_inj  / Nrings
    mdot_rp1_per_ring  = mdot_rp1_inj  / Nrings

    mani_lox_inner = design_manifold(mdot_lox_per_ring, rho_lox,
                                     Lpath_ox[0], deltaP_inj, dpFracMani,
                                     fOx, KOx, h_manifold, w_guess)
    mani_lox_outer = design_manifold(mdot_lox_per_ring, rho_lox,
                                     Lpath_ox[1], deltaP_inj, dpFracMani,
                                     fOx, KOx, h_manifold, w_guess)
    mani_lox_main = design_manifold(mdot_lox, rho_lox,
                                    Lpath_ox_avg, deltaP_inj, dpFracMani,
                                    fOx, KOx, h_manifold, w_guess)

    mani_rp1_inner = design_manifold(mdot_rp1_per_ring, rho_rp1,
                                     Lpath_rp1[0], deltaP_inj, dpFracMani,
                                     fRP1, KRP1, h_manifold, w_guess)
    mani_rp1_outer = design_manifold(mdot_rp1_per_ring, rho_rp1,
                                     Lpath_rp1[1], deltaP_inj, dpFracMani,
                                     fRP1, KRP1, h_manifold, w_guess)

    print("=== Manifold design results ===")
    print_mani("RP-1 inner ring manifold", mani_rp1_inner, Lpath_rp1[0], deltaP_inj, h_manifold)
    print_mani("RP-1 outer ring manifold", mani_rp1_outer, Lpath_rp1[1], deltaP_inj, h_manifold)
    print_mani("LOX inner ring manifold",  mani_lox_inner, Lpath_ox[0],  deltaP_inj, h_manifold)
    print_mani("LOX outer ring manifold",  mani_lox_outer, Lpath_ox[1],  deltaP_inj, h_manifold)
    print_mani("LOX main manifold (avg)",  mani_lox_main,  Lpath_ox_avg, deltaP_inj, h_manifold)

    # -----------------------------------------------------------------------
    # Injector layout plot (EXIT-plane rings)
    # -----------------------------------------------------------------------
    D_c_mm = CombDiam * 1e3
    R_ring_f_mm   = Rring_rp1_exit * 1e3
    R_ring_lox_mm = Rring_lox_exit * 1e3
    N_hole_ring   = int(num_holes_lox_inj / Nrings)

    plot_injector_layout(D_c_mm, R_ring_f_mm, R_ring_lox_mm, N_hole_ring)
