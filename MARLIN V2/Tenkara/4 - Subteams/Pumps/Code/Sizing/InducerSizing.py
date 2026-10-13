import numpy as np
import matplotlib.pyplot as plt
from pint import Quantity as Q_
from Inducer import Inducer, InducerDesignPoint, InducerInletDesignConstraints

g = Q_(9.80, 'm/s^2')

# 1. Nominal Design Point
design_point = InducerDesignPoint(
    Q=Q_(1.69, 'L/s'),
    N=Q_(40000, 'rpm'),
    NPSH_a=Q_(105.46, 'm'),
    NPSH_r=Q_(150.0, 'm')
)

# 2. Sweep Ranges (Defined directly as variables)
N_range = Q_(np.linspace(35000, 40000, 6), 'rpm')
NPSH_a_sweep = Q_(np.linspace(5, 150, 100), 'm')

# 3. Geometry Constraints
constraints = InducerInletDesignConstraints(
    hub_tip_ratio=0.3
)

# 4. Instantiate and Calculate Nominal Point
inducer = Inducer(design_point, constraints, g=g)

# 5. Print Sizing Results
print("\n--- Inputs ---")
print(f"Volumetric Flow (Q)      = {inducer.DP.Q.to('L/s'):.2f} ({inducer.DP.Q.to('gal/min'):.2f})")
print(f"Suction Head Available   = {inducer.DP.NPSH_a.to('m'):.1f} ({inducer.DP.NPSH_a.to('ft'):.2f})")
print(f"Shaft Speed (n)          = {inducer.DP.N.to('rpm'):.0f}")
print(f"Suction Head Required    = {inducer.DP.NPSH_r.to('m'):.1f}")
print(f"Head Rise (H)            = {inducer.Delta_H.to('m'):.1f}")
print(f"Hub-to-Tip Ratio (rho)   = {inducer.C.hub_tip_ratio:.2f}")

print("\n--- Inlet Sizing ---")
print(f"Suction Specific Speed   = {inducer.N_ss:.2f}")
print(f"Corrected Nss (Nss')     = {inducer.N_ss_corrected:.2f}")
print(f"Flow Coefficient (phi)   = {inducer.phi:.6f}")
print(f"Tip Diameter (D_t)       = {inducer.D_tip.to('in'):.4f} ({inducer.D_tip.to('ft'):.6f})")
print(f"Hub Diameter (D_h)       = {inducer.D_hub.to('in'):.4f} ({inducer.D_hub.to('ft'):.6f})")

# 6. Parametric Plotting
plt.figure(figsize=(9, 5.5))

for N_val in N_range:
    d_tip_in = []
    Nss_sweep = []
    for npsh_val in NPSH_a_sweep:
        dp = InducerDesignPoint(
            Q=design_point.Q,
            N=N_val,
            NPSH_a=npsh_val,
            NPSH_r=design_point.NPSH_r,
        )
        inducer_stage = Inducer(dp, constraints, g=g)
        Nss_sweep.append(inducer_stage.N_ss_corrected)
        d_tip_in.append(inducer_stage.D_tip.to('in').magnitude)

    plt.plot(Nss_sweep, d_tip_in, label=f"{N_val.magnitude:.0f} RPM", lw=1.8)

# Mark the nominal design point on the plot
plt.plot(
    inducer.N_ss_corrected,
    inducer.D_tip.to('in').magnitude,
    'ro',
    markersize=6,
    label=f"Nominal ({inducer.DP.N.magnitude:.0f} RPM, {inducer.D_tip.to('in'):.3f})"
)

plt.xlabel('Corrected Suction Specific Speed, $Nss$', fontsize=11)
plt.ylabel('Inducer Tip Diameter, $D_t$ [in]', fontsize=11)
plt.title('Inducer Tip Diameter vs. $Nss$ (35,000 – 40,000 RPM)', fontsize=12, pad=10)
plt.grid(True, linestyle='--', alpha=0.5)
plt.legend(frameon=True, loc='upper right')
plt.tight_layout()
plt.show()