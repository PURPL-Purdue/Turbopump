import pint 
import cea 
import numpy as np

from cea_setup import cea_setup

ureg = pint.UnitRegistry()

Dt = 12 * ureg.cm
Pc = 65 * ureg.bar

reac_names = ["H2", "O2(L)"]
reac_temps = np.array([298, 90.17]) * ureg.K
fuel_weights = np.array([1.0, 0.0])
oxidant_weights = np.array([0.0, 1.0])
cont_ratio = 1.8
of_ratio = 5.5
Prat = 1.7384

solution = cea_setup(reac_names, oxidant_weights, fuel_weights, of_ratio, reac_temps.magnitude, Pc.to(ureg.bar).magnitude, Prat, cont_ratio)
num_pts = solution.num_pts
T = solution.T
P = solution.P
rho = solution.density
enthalpy = solution.enthalpy
energy = solution.energy
gibbs = solution.gibbs_energy
entropy = solution.entropy
M_1n = solution.M
MW = solution.MW
cp_eq = solution.cp_eq
cp_fr = solution.cp_fr
cv_eq = solution.cv_eq
cv_fr = solution.cv_fr
Mach = solution.Mach
gamma_s = solution.gamma_s
v_sonic = solution.sonic_velocity
ae_at = solution.ae_at
c_star = solution.c_star
Cf = solution.coefficient_of_thrust
Isp = solution.Isp
Isp_vac = solution.Isp_vacuum
visc = solution.viscosity

# bartz equation for hg
Mu0 = visc[0]

#hg = (0.026 / Dt**0.2) * ((Mu0**0.2 * Cp0) / Pr0**0.6) * (Rhoinf * Vinf) ** 0.8 *(Rhoam / Rhoinf)**0.8 * (Muam/Mu0)**0.2