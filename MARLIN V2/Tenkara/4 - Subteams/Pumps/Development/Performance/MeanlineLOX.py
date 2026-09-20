import numpy as np
from Impeller import DesignPoint, InputGeometry, Impeller
from pint import Quantity as Q_
from matplotlib import pyplot as plt
import CoolProp.CoolProp as CP

g = Q_(9.81, 'm/s^2')				# gravitational acceleration

opt_m_dot = Q_(1.94, 'kg/s')		# mass flow rate at BEP
rho = Q_(1141, 'kg/m^3')			# fluid density
opt_dP = Q_(33, 'bar')				# total pressure rise at BEP

design_point: DesignPoint = DesignPoint(
	Q = opt_m_dot/rho,
	H = (opt_dP / (g * rho)),    	# developed head at BEP
	N_shaft = Q_(35000, 'rpm'),   	# shaft speed
	n_hyd_BEP = 0.5
)
lox_geometry: InputGeometry = InputGeometry(
	Z_blade = 6,
	Beta2B = Q_(20,'deg').to('rad'),
	D2 = Q_(2, 'in'),				# impeller outlet diameter
	b2 = Q_(0.15, 'in'),			# impeller outlet height
	thk2 = Q_(0.04, 'in'),			# blade thickness at exit
)

lox_imp = Impeller(lox_geometry, design_point)

print("\n--- Design point parameters ---")
print(f"Flowrate (Q)                = {lox_imp.DP.Q.to("L/s"):.2f}")
print(f"Headrise (ΔH)               = {lox_imp.DP.H.to('m'):.1f}")
print(f"Shaft Speed (N)             = {lox_imp.DP.N_shaft:.0f}")
print(f"Specific speed (imperial)   = {lox_imp.SpecificSpeed:.0f}")

# ------------------------------------------------------------------
# Derived impeller characteristics

print("\n--- Derived impeller characteristics ---")
print(f"Slip factor (σ)             = {lox_imp.GEOM.WiesnerSlip:.4f}")
print(f"Outlet tip speed (U₂)       = {lox_imp.U_2_design.to('m/s'):.3f}")
print(f"Outlet area (A₂)            = {lox_imp.GEOM.Area2.to('in^2'):.3f}")
print(f"Meridional velocity (Cm₂)   = {lox_imp.C_m2_design.to('m/s'):.3f}")
print(f"Head coefficient (ψ)        = {lox_imp.HeadCoeff:.4f}")
print(f"Flow coefficient (ϕ)        = {lox_imp.FlowCoeff:.4f}")

# ------------------------------------------------------------------
# Performance
specific_work, _, _ = lox_imp.GetSpecificWork()
power_consumption = (opt_m_dot * specific_work)

# cavitation criteria: p_inlet <= p_vapor
# although in reality, vaporization can be delayed 
SAFE_CAVITATION_NUMBER = 1
feed_pressure = Q_(150, 'psi')		# inlet feed pressure, from tank pressure
atm_press = Q_(1, 'atm')			# atmospheric pressure
LO2_temp = Q_(						# assume LO2 temperature is saturation temperature at atmospheric pressure
    CP.PropsSI('T', 'P', atm_press.to('Pa').magnitude, 'Q', 0, 'Oxygen'), 'K')

p_vapor_LO2 = Q_(					# vapor pressure at inlet (pressurized)
    CP.PropsSI('P', 'T', LO2_temp.magnitude, 'Q', 0, 'Oxygen'), 'Pa')

inlet_diameter = Q_(1, 'in')		# impeller inlet diam TODO: make it an attriubte of impeller class

# cavitation_num = (p_i - p_vapor) / [1/2 * rho * U^2]
u_i = inlet_diameter/2 * lox_imp.DP.N_shaft # TODO: It is more accurate to use the relative velocity w_i rather than tip velocity if there is significant preswhirl (there is with an inducer)
p_inlet = Q_(
    SAFE_CAVITATION_NUMBER * 1/2 * rho * u_i**2 + p_vapor_LO2,
	'Pa')

NPSH_i = (p_inlet - p_vapor_LO2) / rho / g
NPSH_a = (feed_pressure - p_vapor_LO2) / rho / g # TODO: add penalty due to dynamic pressure using inlet velocity (continuity)

print("\n--- Performance characteristics ---")
print(f"Power consumption, nominal  = {power_consumption.to('kW'):.2f}")
print(f"NPSH @ inception           = {NPSH_i.to('m'):.0f}")
print(f"NPSH available              = {NPSH_a.to('m'):.0f}")

lox_imp.PlotPerformanceHQ(Q_([20000, 25000, 30000, 35000], 'rpm'))
vel, _ = lox_imp.GetOutletVelocities()

vel.Plot(unit='m/s', station=2)

plt.show()