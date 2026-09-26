import numpy as np
from Impeller import DesignPoint, InputGeometry, Impeller
from pint import Quantity as Q_
from matplotlib import pyplot as plt
import CoolProp.CoolProp as CP

g = Q_(9.81, 'm/s^2')				# gravitational acceleration

opt_m_dot = Q_(2.15, 'kg/s')		# mass flow rate at BEP
rho = Q_(786, 'kg/m^3')				# fluid density
opt_dP = Q_(44.63, 'bar')			# total pressure rise at BEP

design_point: DesignPoint = DesignPoint(
	Q = opt_m_dot/rho,				# volumetric flow rate at BEP
	H = (opt_dP / (g * rho)),    	# developed head at BEP
	N_shaft = Q_(35000, 'rpm'),   	# shaft speed
	n_hyd_BEP = 0.5
)
ipa_geometry: InputGeometry = InputGeometry(
	Z_blade = 6,
	Beta2B = Q_(20,'deg').to('rad'),# blade angle at exit, relative to tangent
    d_1 = Q_(1.5, 'in'),			# impeller inlet diameter
	d_2 = Q_(2.5, 'in'),			# impeller outlet diameter
	b_2 = Q_(0.15, 'in'),			# impeller outlet height
	thk2 = Q_(0.04, 'in'),			# blade thickness at exit
    d_hub=Q_(1.2, 'in')				# hub diameter
)

ipa_imp = Impeller(ipa_geometry, design_point)

print("\n--- Design point parameters ---")
print(f"Flowrate (Q)                = {ipa_imp.DP.Q.to("L/s"):.2f}")
print(f"Headrise (ΔH)               = {ipa_imp.DP.H.to('m'):.1f}")
print(f"Shaft Speed (N)             = {ipa_imp.DP.N_shaft:.0f}")
print(f"Specific speed (imperial)   = {ipa_imp.SpecificSpeed:.0f}")
print(f"Specific speed (metric)     = {ipa_imp.n_q:.1f}")

# ------------------------------------------------------------------
# Derived impeller characteristics

print("\n--- Derived impeller characteristics ---")
print(f"Slip factor (σ)             = {ipa_imp.WiesnerSlip:.4f}")
print(f"Outlet tip speed (U₂)       = {ipa_imp.U_2_design.to('m/s'):.3f}")
print(f"Outlet area (A₂)            = {ipa_imp.Area2.to('in^2'):.3f}")
print(f"Meridional velocity (Cm₂)   = {ipa_imp.C_m2_design.to('m/s'):.3f}")
print(f"Head coefficient (ψ)        = {ipa_imp.HeadCoeff:.4f}")
print(f"Flow coefficient (ϕ)        = {ipa_imp.FlowCoeff:.4f}")

# ------------------------------------------------------------------
# Performance
specific_work, _, _ = ipa_imp.GetSpecificWork()
power_consumption = (opt_m_dot * specific_work)


print("\n--- Performance characteristics ---")
print(f"Power consumption, nominal  = {power_consumption.to('kW'):.2f}")

feed_pressure = Q_(150, 'psi')		# inlet feed pressure, from tank pressure
atm_press = Q_(1, 'atm')			# atmospheric pressure
IPA_temp = Q_(25, 'C')				# assume IPA temperature is 25 degrees C

p_vapor_IPA = Q_(6060, 'Pa')		# vapor pressure at inlet (pressurized) - inputting manually bc coolprop has no ipa

NPSH_a = (feed_pressure - p_vapor_IPA) / rho / g

print("\n--- Inlet conditions ---")
print(f"NPSH inception              = {ipa_imp.NPSH_i.to('m'):.0f}")
print(f"NPSH available              = {NPSH_a.to('m'):.0f}")
print(f"Inlet flow velocity         = {(ipa_imp.DP.Q / ipa_imp.Area1).to('m/s'):.1f}")

ipa_imp.PlotPerformanceHQ(Q_([20000, 25000, 30000, 35000], 'rpm'))
vel, _ = ipa_imp.GetOutletVelocities()

vel.Plot(unit='m/s', station=2)

ipa_imp.SweepInletDiam()

plt.show()