import numpy as np
from pint import Quantity as Q_
from _IPA_impeller import ipa_imp
from Inducer import Inducer, design_points, constraints
from Pump import Pump

ipa_pump = Pump(Inducer(design_points['LOx'], constraints), ipa_imp)	# inducer design point and constraints from Inducer.py

g = Q_(9.81, 'm/s^2')				# gravitational acceleration

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
print(f"Inlet area (A₁)             = {ipa_imp.Area1.to('in^2'):.3f}")
print(f"Outlet area (A₂)            = {ipa_imp.Area2.to('in^2'):.3f}")
print(f"Meridional velocity (Cm₂)   = {ipa_imp.C_m2_design.to('m/s'):.3f}")
print(f"Head coefficient (ψ)        = {ipa_imp.HeadCoeff:.4f}")
print(f"Flow coefficient (ϕ)        = {ipa_imp.FlowCoeff:.4f}")

# ------------------------------------------------------------------
# Performance
specific_work, _, _ = ipa_imp.GetSpecificWork()
power_consumption = (ipa_imp.DP.Q * ipa_imp.DP.rho * specific_work)


print("\n--- Performance characteristics ---")
print(f"Power consumption, nominal  = {power_consumption.to('kW'):.2f}")

feed_pressure = Q_(150, 'psi')		# inlet feed pressure, from tank pressure
atm_press = Q_(1, 'atm')			# atmospheric pressure
IPA_temp = Q_(25, 'C')				# assume IPA temperature is 25 degrees C

p_vapor_IPA = Q_(6060, 'Pa')		# vapor pressure at inlet (pressurized) - inputting manually bc coolprop has no ipa

NPSH_a = (feed_pressure - p_vapor_IPA) / ipa_imp.DP.rho / g

print("\n--- Inlet conditions ---")
print(f"NPSH inception              = {ipa_imp.NPSH_i.to('m'):.0f}")
print(f"NPSH available              = {NPSH_a.to('m'):.0f}")
print(f"Inlet flow velocity         = {(ipa_imp.DP.Q / ipa_imp.Area1).to('m/s'):.1f}")

ipa_imp.PlotPerformanceHQ(Q_([20000, 25000, 30000, 35000], 'rpm'))

ipa_imp.SweepInletDiam()

ipa_imp.PlotMeridional(color='green')
ipa_pump.PlotMeridional(color='green')
ipa_imp.PotentialFlow.Plot()
ipa_imp.ShroudHubRevolve()

ipa_imp.SweepVelocities().show()

# torque = power_consumption / ipa_imp.DP.N_shaft
# mu = 0.2 # static coefficient of friction, aluminum on hardened steel
# shaft_diameter = Q_(0.4, 'in')
# mean_radius: Q_[float] = (ipa_imp.d_hub + shaft_diameter)/4 #np.sqrt((ipa_imp.d_hub/2)**2 + (shaft_diameter/2)**2)
# axial_force = (torque / mean_radius / mu).to('kN')
# area = 1/4 * (ipa_imp.d_hub**2 - shaft_diameter**2)
# stress = (axial_force / area).to('MPa')
# print(mean_radius)
# print(torque.to('kN*m'))
# print(axial_force)
# print(stress)