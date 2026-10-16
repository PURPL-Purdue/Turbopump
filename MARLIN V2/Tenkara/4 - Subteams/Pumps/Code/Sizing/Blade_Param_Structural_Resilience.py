import numpy as np
import matplotlib.pyplot as plt
 
 
def get_float(prompt: str) -> float:
    """Keep asking until the user enters a valid number."""
    while True:
        try:
            return float(input(prompt))
        except ValueError:
            print("Please enter a valid number.")
 
 
def main():
    # ---------------- Inputs ----------------
    H_Diameter = get_float("Enter hub diameter (H_Diameter): ")
    T_Diameter = get_float("Enter tip diameter (T_Diameter): ")
    Blade_Angle = get_float("Enter blade angle in degrees (Blade_Angle): ")
 
    Blade_Angle_rad = np.deg2rad(Blade_Angle)
 
    # ---------------- Parameters ----------------
    H_Radius = H_Diameter / 2
    T_Radius = T_Diameter / 2
 
    r = np.linspace(H_Radius, T_Radius, 500)  
 
    Blade_Lead = 2 * np.pi * r * np.tan(Blade_Angle_rad) #SP 8052
    Cascade_Solidity = 2.5 #SP 8052
    Blade_Num = 3 #Investigation of Inducer Geometry Variations on Centrifugal Pump Hydraulic Performance: Numerical Approach
    Chord_Length = Cascade_Solidity * Blade_Spacing #SP 8052
    Blade_Spacing = Chord_Length / Cascade_Solidity #SP 8052
    Blade_Axial_Length = Cascade_Solidity * Blade_Lead #SP 8052
    Axial_Clearance = Blade_Axial_Length / Cascade_Solidity #SP 8052
 
    # ---------------- Display ----------------
    print("\n=========== Parameters ===========")
    print(f"H_Diameter         = {H_Diameter}")
    print(f"T_Diameter         = {T_Diameter}")
    print(f"Blade_Angle        = {Blade_Angle} deg ({Blade_Angle_rad:.4f} rad)")
    print(f"H_Radius           = {H_Radius}")
    print(f"T_Radius           = {T_Radius}")
    print(f"Cascade_Solidity   = {Cascade_Solidity}")
    print(f"Blade_Num          = {Blade_Num}")
    print(f"Blade_Spacing      = {Blade_Spacing:.4f}")
    print(f"Chord_Length       = {Chord_Length:.4f}")
    print(f"r (radius)         : {r.size} points, {r[0]:.4f} to {r[-1]:.4f}")
    print(f"Blade_Lead         : min {Blade_Lead.min():.4f}, max {Blade_Lead.max():.4f}")
    print(f"Blade_Axial_Length : min {Blade_Axial_Length.min():.4f}, max {Blade_Axial_Length.max():.4f}")
    print(f"Axial_Clearance    : min {Axial_Clearance.min():.4f}, max {Axial_Clearance.max():.4f}")
    print("==================================\n")
 
    # ---------------- Plots ----------------
    fig, axes = plt.subplots(3, 1, figsize=(8, 10), sharex=True)
 
    # Chord length is a single value, so draw it as a constant line across the radius
    axes[0].plot(r, np.full_like(r, Chord_Length), color="tab:blue")
    axes[0].set_ylabel("Chord Length")
    axes[0].set_title("Chord Length vs Radius")
 
    axes[1].plot(r, Blade_Lead, color="tab:orange")
    axes[1].set_ylabel("Blade Lead")
    axes[1].set_title("Blade Lead vs Radius")
 
    axes[2].plot(r, Blade_Axial_Length, color="tab:green")
    axes[2].set_ylabel("Blade Axial Length")
    axes[2].set_xlabel("Radius, r")
    axes[2].set_title("Blade Axial Length vs Radius")
 
    for ax in axes:
        ax.grid(True)
 
    fig.tight_layout()
    plt.show()
 
 
if __name__ == "__main__":
    main()
 
