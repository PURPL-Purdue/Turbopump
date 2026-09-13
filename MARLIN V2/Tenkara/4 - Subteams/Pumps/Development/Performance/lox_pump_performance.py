from pint import Quantity as Q_
import numpy as np
import matplotlib.pyplot as plt
import EmpiricalRelations as emp

# ------------------------------------------------------------------
# Design point parameters
g = Q_(9.81, 'm/s^2')                   # gravitational acceleration
opt_m_dot = Q_(1.94, 'kg/s')            # mass flow rate at BEP
rho = Q_(1141, 'kg/m^3')                # fluid density
opt_flow_rate = opt_m_dot / rho         # volumetric flow rate at BEP
opt_pressure_rise = Q_(488, 'psi')      # total pressure rise at BEP
opt_H = (opt_pressure_rise / (g * rho)) # developed head at BEP
N_shaft = Q_(40000, 'rpm')              # shaft speed
omega = N_shaft.to('rad/s')

# ------------------------------------------------------------------
# Basic impeller parameters
Z = 6                           # number of impeller blades
beta2b = Q_(25,'deg').to('rad') # blade backsweep angle at outlet (from tangent)
head_coeff = 0.45               # head coefficient  gH/U2^2  at BEP
flow_coeff = 0.1                # flow coefficient  Cm2/U2   at BEP (~0.08-0.13 typ.)
n_hyd_BEP = 0.5                 # Hydraulic efficiency at BEP

# Wiesner slip factor
# https://manual.cfturbo.com/en/bl_te_wiesner.html
sigma = 1 - np.sqrt(np.sin(beta2b)) / Z**0.7
#print(sigma)
# Solve for outlet tip speed U2 and diameter D2 from the design point
U2 = np.sqrt(g * opt_H / head_coeff).to('m/s')
D2 = (2 * U2 / omega).to('in')

# Outlet meridional velocity & blade width from continuity at BEP
Cm2_design = Q_(flow_coeff * U2, 'm/s')           
Area2 = opt_flow_rate / Cm2_design      # Q / C_m2
b2 = (Area2 / (np.pi * D2)).to('in')    # blade height at exit

print("=== Derived impeller geometry / design parameters ===")
print(f"Slip factor sigma   = {sigma:.3f}")
print(f"Outlet tip speed U2 = {U2:.2f}")
print(f"Impeller OD  D2     = {D2:.2f}")
print(f"Outlet width b2     = {b2:.3f}")


Q_sweep = Q_(np.linspace(0, 2.0 * opt_flow_rate.to('L/s').magnitude, 200), 'L/s')

# https://ntrs.nasa.gov/api/citations/19950013379/downloads/19950013379.pdf
# page 4 and 5

# Theoretical Euler head as a function of Q
#   H_theoretical(Q) = U2^2/g - [U2 / (g*Area2*tan(beta2b))] * Q
#       Or
#   H_theoretical(Q) = U2^2/g - [omega*cot(beta2b)/(2*pi*b2*g)] * Q
#   https://youtu.be/H7XfYO_-cEg?t=2705

def H_euler(Q):
    return sigma*U2**2/g - omega/np.tan(beta2b)/(2*np.pi*b2*g) * Q

print(f"Theoretical shutoff head H0 = {H_euler(Q_(0, 'L/s')):.2f}")

H_theoretical = H_euler(Q_sweep)

H_predict = []

for q in Q_sweep:

    flowrate = Q_(q, 'm^3/s')

    f = emp.FlowSpeedRatio(
        flowrate, N_shaft, opt_flow_rate.to('m^3/s'), N_shaft
        )

    C_m2 = flowrate / Area2
    slip = U2 * (1 - sigma)
    W_u2 = C_m2 / np.tan(beta2b) + slip
    C_u2 = (U2 - W_u2)

    h_euler = (C_u2 * U2 / g).to('m')
    n_hydraulic = emp.HydraulicEfficiency(f) * n_hyd_BEP

    H_predict.append((h_euler * n_hydraulic).to('m').magnitude)

    

fig, ax1 = plt.subplots(figsize=(6, 4.5))
 
ax1.plot(Q_sweep.to('L/s').magnitude, H_theoretical.to('m').magnitude,
         '--', color='gray', label='Theoretical Euler')

ax1.plot(Q_sweep.to('L/s').magnitude, H_predict,
         '-', color='red', label=f"Empirical prediction ({100*n_hyd_BEP:.0f}% @BEP)")

ax1.set_xlabel('Volumetric flow rate [L/s]')
ax1.set_ylabel('Head [m]')
ax1.set_title('H-Q Curve')
ax1.set_ylim(ymin=0)
ax1.legend()
ax1.grid(True, alpha=0.3)
 
plt.tight_layout()
plt.show()