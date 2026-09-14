from pint import Quantity as Q_
import numpy as np
import matplotlib.pyplot as plt
import EmpiricalRelations as emp

# ------------------------------------------------------------------
# Physical constants
g = Q_(9.81, 'm/s^2')           # gravitational acceleration

# ------------------------------------------------------------------
# Design point parameters
opt_m_dot = Q_(1.94, 'kg/s')    # mass flow rate at BEP
rho = Q_(1141, 'kg/m^3')        # fluid density
opt_Q = opt_m_dot / rho         # volumetric flow rate at BEP
opt_dP = Q_(33, 'bar')          # total pressure rise at BEP
opt_H = (opt_dP / (g * rho))    # developed head at BEP
N_shaft = Q_(40000, 'rpm')      # shaft speed
Ns = (                          # imperial specific speed       
        N_shaft * np.sqrt(opt_Q.to('gal/min')) / opt_H.to('ft')**(3/4)
    ).magnitude

print("\n--- Design point parameters ---")
print(f"Flowrate (Q)             = {opt_Q.to("L/s"):.2f}")
print(f"Headrise (ΔH)            = {opt_H.to('m'):.1f}")
print(f"Shaft Speed (N)         =  {N_shaft:.0f}")
print(f"Specific speed (imperial)   = {Ns:.0f}")

# ------------------------------------------------------------------
# Basic impeller parameters
Z = 6                           # number of impeller blades
beta2b = Q_(10,'deg').to('rad') # blade backsweep angle at outlet (from tangent)
D2 = Q_(2.2, 'in')              # impeller outlet diameter
b2 = Q_(0.1, 'in')              # impeller outlet height
n_hyd_BEP = 0.5                 # Hydraulic efficiency at BEP, empirically chosen prediction

# ------------------------------------------------------------------
# Derived impeller characteristics

sigma = 1 - np.sqrt(np.sin(beta2b)) / Z**0.7 # Wiesner slip factor: https://manual.cfturbo.com/en/bl_te_wiesner.html
U2_opt = (N_shaft * D2/2).to('m/s')   # outlet tip speed
Area2 = (np.pi * D2 * b2).to('in^2')    # Outlet area, continuity TODO: account for metal blockage
C_m2_design = (opt_Q / Area2)   # Meridional flow velocity at design point

head_coeff = (g*opt_H/n_hyd_BEP/U2_opt**2).to('dimensionless')  # head coefficient / stage loading
flow_coeff = (C_m2_design / U2_opt).to('dimensionless')         # flow coefficient / flow factor

print("\n--- Derived impeller characteristics ---")
print(f"Slip factor (σ)             = {sigma:.4f}")
print(f"Outlet tip speed (U₂)       = {U2_opt.to('m/s'):.3f}")
print(f"Outlet area (A₂)            = {Area2.to('in^2'):.3f}")
print(f"Meridional velocity (Cm₂)   = {C_m2_design.to('m/s'):.3f}")
print(f"Head coefficient (ψ)        = {head_coeff:.4f}")
print(f"Flow coefficient (ϕ)        = {flow_coeff:.4f}")

# ------------------------------------------------------------------
# Ideal geometry parameters

print("\n--- Ideal geometry parameters ---")
# Compute impeller diameter at derived head coefficient
D2_opt = 2 * (np.sqrt(g*opt_H / head_coeff) / N_shaft).to('in')
print(f"Ideal outlet diameter:      = {D2_opt:.3f}")

Q_sweep = Q_(np.linspace(0, 2.0 * opt_Q.to('L/s').magnitude, 50), 'L/s')

# https://ntrs.nasa.gov/api/citations/19950013379/downloads/19950013379.pdf
# page 4 and 5

# Theoretical Euler head as a function of Q
#   H_theoretical(Q) = U2^2/g - [U2 / (g*Area2*tan(beta2b))] * Q
#       Or
#   H_theoretical(Q) = U2^2/g - [N_shaft*cot(beta2b)/(2*pi*b2*g)] * Q
#   https://youtu.be/H7XfYO_-cEg?t=2705

def H_euler(Q):
    return sigma*U2_opt**2/g - N_shaft/np.tan(beta2b)/(2*np.pi*b2*g) * Q

H_theoretical = H_euler(Q_sweep)

N_sweep = Q_(np.array([
    20000, 25000, 30000, 35000, 40000
]), 'rpm')

H_predict = []
Q_predict = []

for shaft_speed in N_sweep:

    head_trace = []
    flow_trace = []
    
    for flowrate in Q_sweep:

        f = emp.FlowSpeedRatio(
            flowrate.to('m^3/s'), shaft_speed,
            opt_Q.to('m^3/s'), N_shaft
            )
        if f > 2: continue

        U2 = shaft_speed * D2/2
        C_m2 = flowrate / Area2
        slip = U2 * (1 - sigma)
        W_u2 = C_m2 / np.tan(beta2b) + slip
        C_u2 = (U2 - W_u2)

        h_euler = (C_u2 * U2 / g).to('m')

        # empirical prediction, valid for 0 < f < 2
        n_hydraulic = emp.HydraulicEfficiency(f) * n_hyd_BEP

        head_trace.append((h_euler * n_hydraulic).to('m').magnitude)
        flow_trace.append(flowrate.to('L/s').magnitude)

    H_predict.append(head_trace)
    Q_predict.append(flow_trace)
    

fig, ax1 = plt.subplots(figsize=(6, 4.5))

# Plot theoretical Euler (linear)
ax1.plot(Q_sweep.to('L/s').magnitude, H_theoretical.to('m').magnitude,
         '--', color='gray', label='Theoretical Euler')

# Plot RPM swept traces
for i in range(len(H_predict)):

    ax1.plot(Q_predict[i], H_predict[i],
         '-', color='red', label = 'Empircal prediction' if i == 0 else None)
    ax1.annotate(f"{N_sweep[i].magnitude:.0f} RPM",
                     xy=(Q_predict[i][-1], H_predict[i][-1]),
                     xytext=(4, 0), textcoords='offset points',
                     fontsize=8, va='center')

# Plot design point head and flowrate
ax1.plot(opt_Q.to('L/s').magnitude, opt_H.to('m').magnitude,
         'o', color='black', markersize=4, label=f'Design point',
         zorder=5)

ax1.set_xlabel('Volumetric flow rate [L/s]')
ax1.set_ylabel('Head [m]')
ax1.set_title('H-Q Curve')
ax1.set_ylim(ymin=0)
ax1.legend()
ax1.grid(True, alpha=0.3)
 
plt.tight_layout()
plt.show()