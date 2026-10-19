import numpy as np
from pint import Quantity as Q_
import CoolProp.CoolProp as CP
import Pump
from Inducer import print_parameters

lox_pump = Pump.LOX()
lox_imp = lox_pump.Impeller

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
print(f"Inlet area (A₁)             = {lox_imp.Area1.to('in^2'):.3f}")
print(f"Outlet area (A₂)            = {lox_imp.Area2.to('in^2'):.3f}")
print(f"Axial extension (Z_e)       = {lox_imp.Z_e:.4f}")
print(f"Meridional velocity (Cm₂)   = {lox_imp.C_m2_design.to('m/s'):.3f}")
print(f"Head coefficient (ψ)        = {lox_imp.HeadCoeff:.4f}")
print(f"Flow coefficient (ϕ)        = {lox_imp.FlowCoeff:.4f}")

# ------------------------------------------------------------------
# Performance
specific_work, _, _ = lox_imp.GetSpecificWork()
power_consumption = (lox_imp.DP.Q * lox_imp.DP.rho * specific_work)
shaft_torque = power_consumption / lox_imp.DP.N_shaft
shaft_diam = Q_(0.5, 'in')
shear_stress = shaft_torque * shaft_diam / 2 / (np.pi/32 * shaft_diam**4)

print("\n--- Performance characteristics ---")
print(f"Power consumption, nominal  = {power_consumption.to('kW'):.2f}")
print(f"Shaft torque, nominal       = {shaft_torque.to('N*m'):.1f}")
print(f"Shear stress                = {shear_stress.to('MPa'):.1f}")

atm_press = Q_(1, 'atm')			# atmospheric pressure
LO2_temp = Q_(						# assume LO2 temperature is saturation temperature at atmospheric pressure
    CP.PropsSI('T', 'P', atm_press.to('Pa').magnitude, 'Q', 0, 'Oxygen'), 'K')

p_vapor_LO2 = Q_(					# vapor pressure at inlet (pressurized)
    CP.PropsSI('P', 'T', LO2_temp.magnitude, 'Q', 0, 'Oxygen'), 'Pa')

feed_pressure = Q_(150, 'psi')
rho = lox_imp.DP.rho
g = Q_(9.81, 'm/s^2')
NPSH_a = (feed_pressure - p_vapor_LO2) / rho / g

print("\n--- Inlet conditions ---")
print(f"NPSH inception              = {lox_imp.NPSH_i.to('m'):.0f}")
print(f"NPSH available              = {NPSH_a.to('m'):.0f}")
print(f"Inlet flow velocity         = {(lox_imp.DP.Q / lox_imp.Area1).to('m/s'):.1f}")

# ------------------------------------------------------------------
# Inducer

lox_inducer = lox_pump.Inducer
print_parameters(lox_inducer, 'LOX')
lox_inducer.PlotMeridional(name='LOX Inducer')
lox_inducer.PlotSwirl(name='LOX Inducer')

def dothething():
    lox_imp.PlotPerformanceHQ(Q_([20000, 25000, 30000, 35000], 'rpm'))

    #lox_imp.SweepInletDiam()

    lox_imp.PlotMeridional(color='blue')
    lox_pump.PlotMeridional(color='blue')
    lox_imp.PotentialFlow.Plot()
    lox_imp.ShroudHubRevolve()
    #lox_imp.GetVelocities(station=1).Plot().show()
    lox_imp.SweepVelocities().show()