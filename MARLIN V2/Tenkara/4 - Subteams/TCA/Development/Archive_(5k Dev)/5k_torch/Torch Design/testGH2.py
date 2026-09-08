import numpy as np
from pyfluids import Fluid, FluidsList, Input

# Given
p1_psi = 464.844
p2_psi = 150
T0_C = 25
mdot_lbm_s = 0.004
d_in = 0.04264192

# Constants
psi2Pa = 6894.76
lbm2kg = 0.453592
in2m = 0.0254

# Convert
p1 = p1_psi * psi2Pa         # [Pa]
p2 = p2_psi * psi2Pa         # [Pa]
T0 = T0_C + 273.15           # [K]
mdot = mdot_lbm_s * lbm2kg   # [kg/s]
d = d_in * in2m              # [m]
A = np.pi * (d / 2)**2          # [m^2]

# GH2 properties
gamma = 1.41                 # for GH2
R = 4124                     # J/kg·K
Tt = 2 * T0 / (gamma + 1)

pt = p1 * (2 / (gamma + 1)) ** (gamma / (gamma - 1))
v = np.sqrt((2 * gamma / (gamma + 1)) * R * T0)

# Step 3: Report
print(f"GH2 velocity: {v:.2f} m/s")


fuel = Fluid(FluidsList.Hydrogen).with_state(
    Input.pressure(p1),
    Input.temperature(T0_C)
)

rho0 = fuel.density  # [kg/m^3]
gamma = 1.41
rho_star = rho0 * (2 / (gamma + 1)) ** (1 / (gamma - 1))

v_cont = mdot / (rho_star * A)

print(v_cont)
print(A)
print(fuel.sound_speed)