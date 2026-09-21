import numpy as np
from Impeller import DesignPoint, InputGeometry, Impeller
from pint import Quantity as Q_
from matplotlib import pyplot as plt
import CoolProp.CoolProp as CP

g = Q_(9.81, 'm/s^2')				# gravitational acceleration

opt_m_dot = Q_(2.03, 'kg/s')		# mass flow rate at BEP
rho = Q_(1141, 'kg/m^3')			# fluid density
opt_dP = Q_(33, 'bar')				# total pressure rise at BEP

design_point: DesignPoint = DesignPoint(
	Q = opt_m_dot/rho,				# volumetric flow rate at BEP
	H = (opt_dP / (g * rho)),    	# developed head at BEP
	N_shaft = Q_(35000, 'rpm'),   	# shaft speed
	n_hyd_BEP = 0.5
)
lox_geometry: InputGeometry = InputGeometry(
	Z_blade = 6,
	Beta2B = Q_(20,'deg').to('rad'),# blade angle at exit, relative to tangent
	d_2 = Q_(2.2, 'in'),			# impeller outlet diameter
	b_2 = Q_(0.07, 'in'),			# impeller outlet height
	thk2 = Q_(0.04, 'in'),			# blade thickness at exit
    d_hub=Q_(0.5, 'in')
)

lox_imp = Impeller(lox_geometry, design_point)

print("\n--- Design point parameters ---")
print(f"Flowrate (Q)                = {lox_imp.DP.Q.to("L/s"):.2f}")
print(f"Headrise (ΔH)               = {lox_imp.DP.H.to('m'):.1f}")
print(f"Shaft Speed (N)             = {lox_imp.DP.N_shaft:.0f}")
print(f"Specific speed (imperial)   = {lox_imp.SpecificSpeed:.0f}")
print(f"Specific speed (metric)     = {lox_imp.n_q:.1f}")

# ------------------------------------------------------------------
# Derived impeller characteristics

print("\n--- Derived impeller characteristics ---")
print(f"Slip factor (σ)             = {lox_imp.WiesnerSlip:.4f}")
print(f"Outlet tip speed (U₂)       = {lox_imp.U_2_design.to('m/s'):.3f}")
print(f"Outlet area (A₂)            = {lox_imp.Area2.to('in^2'):.3f}")
print(f"Meridional velocity (Cm₂)   = {lox_imp.C_m2_design.to('m/s'):.3f}")
print(f"Head coefficient (ψ)        = {lox_imp.HeadCoeff:.4f}")
print(f"Flow coefficient (ϕ)        = {lox_imp.FlowCoeff:.4f}")

# ------------------------------------------------------------------
# Performance
specific_work, _, _ = lox_imp.GetSpecificWork()
power_consumption = (opt_m_dot * specific_work)


print("\n--- Performance characteristics ---")
print(f"Power consumption, nominal  = {power_consumption.to('kW'):.2f}")

feed_pressure = Q_(150, 'psi')		# inlet feed pressure, from tank pressure
atm_press = Q_(1, 'atm')			# atmospheric pressure
LO2_temp = Q_(						# assume LO2 temperature is saturation temperature at atmospheric pressure
    CP.PropsSI('T', 'P', atm_press.to('Pa').magnitude, 'Q', 0, 'Oxygen'), 'K')

p_vapor_LO2 = Q_(					# vapor pressure at inlet (pressurized)
    CP.PropsSI('P', 'T', LO2_temp.magnitude, 'Q', 0, 'Oxygen'), 'Pa')

NPSH_a = (feed_pressure - p_vapor_LO2) / rho / g

print("\n--- Inlet conditions ---")
print(f"NPSH inception              = {lox_imp.NPSH_i.to('m'):.0f}")
print(f"NPSH available              = {NPSH_a.to('m'):.0f}")
print(f"Inlet flow velocity         = {(lox_imp.DP.Q / lox_imp.Area1).to('m/s'):.1f}")

lox_imp.PlotPerformanceHQ(Q_([20000, 25000, 30000, 35000, 40000], 'rpm'))
vel, _ = lox_imp.GetOutletVelocities()

vel.Plot(unit='m/s', station=2)

vel1 = lox_imp.GetInletVelocities()
vel1.Plot(unit='m/s', station=1)

# TODO: could probably move into the impeller class and clean it up
npsh = []
diam = []
min_NPSH = Q_(np.inf, 'm')
opt_diam = Q_(-1, 'in')

for d_1 in np.linspace(lox_imp.d_hub, lox_imp.d_2, 100):
    lox_imp.d_1 = d_1
    NPSH_i = lox_imp.NPSH_i
    
    if NPSH_i < min_NPSH:
        min_NPSH = NPSH_i.to('m')
        opt_diam = d_1.to('in')
    
    npsh.append(NPSH_i.magnitude)
    diam.append(d_1.to('in').magnitude)

plt.figure()
plt.plot(diam, npsh)
plt.plot([opt_diam.magnitude], [min_NPSH.magnitude], 'o', label=f'optimal {opt_diam:.2f}, {min_NPSH:.0f}')
plt.xlabel('Inlet diameter (in.)')
plt.ylabel('NPSH_i (m)')
plt.title("NPSH inception vs inlet diameter")
plt.grid()
plt.legend()

plt.show()