import numpy as np
from Impeller import DesignPoint, InputGeometry, Impeller
from pint import Quantity as Q_
from matplotlib import pyplot as plt

g = Q_(9.81, 'm/s^2')           # gravitational acceleration

opt_m_dot = Q_(1.94, 'kg/s')    # mass flow rate at BEP
rho = Q_(1141, 'kg/m^3')        # fluid density
opt_dP = Q_(33, 'bar')          # total pressure rise at BEP

design_point: DesignPoint = DesignPoint(
    Q = opt_m_dot/rho,
    H = (opt_dP / (g * rho)),    # developed head at BEP
    N_shaft = Q_(40000, 'rpm'),      # shaft speed
    n_hyd_BEP = 0.5
)
geometry: InputGeometry = InputGeometry(
    Z_blade = 6,
    Beta2B = Q_(10,'deg').to('rad'),
    D2 = Q_(2, 'in'),               # impeller outlet diameter
    b2 = Q_(0.15, 'in'),             # impeller outlet height
    thk2 = Q_(0.04, 'in'),           # blade thickness at exit
)

imp = Impeller(geometry, design_point)

print("\n--- Design point parameters ---")
print(f"Flowrate (Q)                = {imp.DP.Q.to("L/s"):.2f}")
print(f"Headrise (ΔH)               = {imp.DP.H.to('m'):.1f}")
print(f"Shaft Speed (N)             = {imp.DP.N_shaft:.0f}")
print(f"Specific speed (imperial)   = {imp.SpecificSpeed:.0f}")

# ------------------------------------------------------------------
# Derived impeller characteristics

print("\n--- Derived impeller characteristics ---")
print(f"Slip factor (σ)             = {imp.GEOM.WiesnerSlip:.4f}")
print(f"Outlet tip speed (U₂)       = {imp.U_2_design.to('m/s'):.3f}")
print(f"Outlet area (A₂)            = {imp.GEOM.Area2.to('in^2'):.3f}")
print(f"Meridional velocity (Cm₂)   = {imp.C_m2_design.to('m/s'):.3f}")
print(f"Head coefficient (ψ)        = {imp.HeadCoeff:.4f}")
print(f"Flow coefficient (ϕ)        = {imp.FlowCoeff:.4f}")

imp.PlotPerformanceHQ(Q_([20000, 25000, 30000, 35000, 40000], 'rpm'))
vel, _ = imp.GetOutletVelocities()

vel.plot(unit='m/s', station=2)

plt.show()