import sys
from pathlib import Path
import numpy as np
from pint import Quantity as Q_

# ---------------- Import Path Setup ----------------

_here = Path(__file__).resolve().parent

sys.path[:] = [
    p for p in sys.path
    if Path(p or ".").resolve() != _here
]

sys.path.insert(0, str(_here.parent))

# ---------------- Project Imports ----------------

import Main
from Inducer.RadialEquilibrium import InducerRadialEquilibrium


def main():

    # ---------------- Pump and Inducer ----------------

    pump = Main.CreatePump(Main.Propellant.LOX)
    pump.AddInducer(Q_(65, "m"))

    inducer = pump.Inducer

    # ---------------- Radial Equilibrium ----------------

    equilibrium = InducerRadialEquilibrium(
        inducer,
        incidence_ratio=0,
        blade_thickness=Q_(0.5, "mm")
    )

    # ---------------- Necessary Parameters ----------------

    R_tip = inducer.R_tip
    z_inlet = Q_(0, "m")
    z_outlet = inducer.L_ax

    omega = inducer.DP.N_shaft.to("rad/s").magnitude
    U_tip = omega * R_tip.to("m").magnitude
   # ---------------- Blade Geometry ----------------

    T_Diameter = (2 * inducer.R_tip).to("m").magnitude
    Blade_Num = inducer.Z_blade

    Blade_Spacing = (np.pi * T_Diameter) / Blade_Num

    # ---------------- Relative Tip Velocities ----------------

    def relative_tip_velocity(z):
        C_m = equilibrium.c_axial(
            R_tip, z
        ).to("m/s").magnitude

        C_theta = equilibrium.c_theta(
            R_tip, z
        ).to("m/s").magnitude

        return np.hypot(C_m, C_theta - U_tip)

    W1_tip = relative_tip_velocity(z_inlet)
    W2_tip = relative_tip_velocity(z_outlet)

    # ---------------- Inlet Angles ----------------

    Blade_Angle_rad = float(equilibrium.beta1)
    Incidence_Angle_rad = float(equilibrium.alpha1)

    # ---------------- Maximum Blade Thickness ----------------

    # Original thickness equation, retained for reference.
    # Requires validated local angles and a consistent flow solution.
    #
    Blade_Max_Thickness = Blade_Spacing * (np.sin(Blade_Angle_rad) - (W1_tip / W2_tip) * np.sin(Blade_Angle_rad - Incidence_Angle_rad))

    

    # ---------------- Display ----------------

    print("\n========== BLADE PARAMETERS ==========")
    print(f"Tip Diameter       = {T_Diameter * 1000:.3f} mm")
    print(f"Number of Blades   = {Blade_Num}")
    print(f"Blade Spacing      = {Blade_Spacing * 1000:.3f} mm")
    print(f"Inlet Blade Angle  = {np.degrees(Blade_Angle_rad):.3f} deg")
    print(f"Inlet Incidence    = {np.degrees(Incidence_Angle_rad):.3f} deg")
    print(f"Inlet Relative W   = {W1_tip:.3f} m/s")
    print(f"Outlet Relative W  = {W2_tip:.3f} m/s")
    print(f"Maximum Thickness  = {Blade_Max_Thickness * 1000:.3f} mm")

    print("\n========== INLET ANGLES ==========")
    print(
        f"Blade Angle = "
        f"{np.degrees(Blade_Angle_rad):.3f} deg"
    )
    print(
        f"Incidence Angle = "
        f"{np.degrees(Incidence_Angle_rad):.3f} deg"
    )

    print("\n========== RELATIVE TIP VELOCITIES ==========")
    print(f"W1_tip = {W1_tip:.3f} m/s")
    print(f"W2_tip = {W2_tip:.3f} m/s")

    return {
        "W1_tip": W1_tip,
        "W2_tip": W2_tip,
        "Blade_Angle_rad": Blade_Angle_rad,
        "Incidence_Angle_rad": Incidence_Angle_rad,
        "Blade_Max_Thickness": Blade_Max_Thickness,
    }


if __name__ == "__main__":
    main()