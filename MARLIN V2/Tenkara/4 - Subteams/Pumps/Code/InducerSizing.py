import numpy as np
import matplotlib.pyplot as plt
from pint import Quantity as Q_
from dataclasses import replace
from Inducer import Inducer, InducerDesignPoint, InducerInletDesignConstraints

g = Q_(9.81, 'm/s^2')

# 1. Nominal Design Points
design_points = {
    "LOx": InducerDesignPoint(
        Q=Q_(1.69, 'L/s'),
        N=Q_(35000, 'rpm'),
        NPSH_a=Q_(105.46, 'm'),
        NPSH_r=Q_(150.0, 'm')
    ),
    "IPA": InducerDesignPoint(
        Q=Q_(2.739, 'L/s'),
        N=Q_(35000, 'rpm'),
        NPSH_a=Q_(105.46, 'm'),
        NPSH_r=Q_(150.0, 'm')
    ),
}

# 2. Sweep Ranges (Defined directly as variables)
hub_tip_sweep = np.linspace(0.2, 0.6, 100)

# 3. Geometry Constraints
constraints = InducerInletDesignConstraints(
    hub_tip_ratio=0.5,
    sigma=0.03,
    K=1.1,
    psi_1=0.08,
    axial_length_ratio=1.0
)

# 4. Instantiate and Calculate Nominal Points
inducers = {name: Inducer(dp, constraints, g=g) for name, dp in design_points.items()}

# 5. Print Sizing Results
for name, inducer in inducers.items():
    print(f"\n===== {name} Inducer =====")
    print("\n--- Inputs ---")
    print(f"Volumetric Flow (Q)      = {inducer.DP.Q.to('L/s'):.2f} ({inducer.DP.Q.to('gal/min'):.2f})")
    print(f"Suction Head Available   = {inducer.DP.NPSH_a.to('m'):.1f} ({inducer.DP.NPSH_a.to('ft'):.2f})")
    print(f"Shaft Speed (n)          = {inducer.DP.N.to('rpm'):.0f}")
    print(f"Suction Head Required    = {inducer.DP.NPSH_r.to('m'):.1f}")
    print(f"Head Rise (H)            = {inducer.Delta_H.to('m'):.1f}")
    print(f"Hub-to-Tip Ratio (rho)   = {inducer.C.hub_tip_ratio:.2f}")
    print(f"Cavitation Param (sigma) = {inducer.C.sigma:.3f}")
    print(f"Constant (K)             = {inducer.C.K:.2f}")
    print(f"Baseline Head Coeff psi1 = {inducer.C.psi_1:.4f}")
    print(f"Axial Length Ratio L/D_t = {inducer.C.axial_length_ratio:.2f}")

    print("\n--- Inlet Sizing ---")
    print(f"Suction Specific Speed   = {inducer.N_ss:.2f}")
    print(f"Flow Coefficient (phi)   = {inducer.phi:.6f}")
    print(f"Estimated NPSH3          = {inducer.NPSH_3.to('m'):.2f} ({inducer.NPSH_3.to('ft'):.2f})")
    print(f"Tip Diameter (D_t)       = {inducer.D_tip.to('in'):.4f} ({inducer.D_tip.to('ft'):.6f})")
    print(f"Hub Diameter (D_h)       = {inducer.D_hub.to('in'):.4f} ({inducer.D_hub.to('ft'):.6f})")

    print("\n--- Hub Taper ---")
    print(f"Tip Speed (U_t)          = {inducer.U_tip.to('m/s'):.2f} ({inducer.U_tip.to('ft/s'):.2f})")
    print(f"Target Head Coeff (psi2) = {inducer.psi_2:.4f}")
    print(f"Inlet Hub Radius (R_h1)  = {inducer.R_hub_1.to('mm'):.3f}")
    print(f"Outlet Hub Radius (R_h2) = {inducer.R_hub_2.to('mm'):.3f}")
    print(f"Axial Length (L_a)       = {inducer.L_ax.to('mm'):.3f}")
    print(f"Hub Taper Angle (alpha)  = {inducer.alpha_taper.to('deg'):.2f}")

# 6. Parametric Plotting
fig, axes = plt.subplots(1, len(inducers), figsize=(13, 5.5), sharey=True)

for ax, (name, inducer) in zip(axes, inducers.items()):
    d_tip_in = []
    d_hub_in = []
    for ht_val in hub_tip_sweep:
        c = replace(constraints, hub_tip_ratio=ht_val)
        inducer_stage = Inducer(inducer.DP, c, g=g)
        d_tip_in.append(inducer_stage.D_tip.to('in').magnitude)
        d_hub_in.append(inducer_stage.D_hub.to('in').magnitude)

    ax.plot(hub_tip_sweep, d_tip_in, label="Tip Diameter $D_t$", lw=1.8)
    ax.plot(hub_tip_sweep, d_hub_in, label="Hub Diameter $D_h$", lw=1.8)

    # Mark the nominal design point on the plot
    ax.plot(
        inducer.C.hub_tip_ratio,
        inducer.D_tip.to('in').magnitude,
        'ro',
        markersize=6,
        label=f"Nominal $D_t$ ({inducer.D_tip.to('in'):.3f})"
    )
    ax.plot(
        inducer.C.hub_tip_ratio,
        inducer.D_hub.to('in').magnitude,
        'ks',
        markersize=6,
        label=f"Nominal $D_h$ ({inducer.D_hub.to('in'):.3f})"
    )

    ax.set_xlabel('Hub-to-Tip Ratio, $\\rho$', fontsize=11)
    ax.set_title(f'{name} ({inducer.DP.Q.to("L/s").magnitude:.2f} L/s, {inducer.DP.N.magnitude:.0f} RPM, $\\sigma$ = {inducer.C.sigma:.2f})', fontsize=12, pad=10)
    ax.grid(True, linestyle='--', alpha=0.5)
    ax.legend(frameon=True, loc='center left')

axes[0].set_ylabel('Inducer Diameter [in]', fontsize=11)
fig.suptitle('Inducer Tip and Hub Diameter vs. Hub-to-Tip Ratio', fontsize=13)
plt.tight_layout()
plt.show()