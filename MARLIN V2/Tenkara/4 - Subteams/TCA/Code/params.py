### Similar script to setpoint.py but determines the engine propreties for a given enine ###

from matplotlib import contour
import numpy as np
import cea
from pint import UnitRegistry
import yaml
import pandas as pd

from contour_script import contour_length_mm
from contour_script import bell_nozzle

ureg = UnitRegistry()

# --------------------------------------------------------------------------- #
# 1. Load engine setpoints
# --------------------------------------------------------------------------- #

with open('Inputs/TCA_params.yaml') as f:
    p = yaml.safe_load(f)

F = p['Thrust_target'] * ureg.lbf ## Target thrust
pc = p['chamber_pressure'] * ureg.psi ## Target chamber pressure
Dc = p['chamber_diameter'] * ureg.inch ## Chamber diameter
tc = p['chamber_thickness'] * ureg.inch ## Chamber thickness

# Fixed design parameters
pe = p['exit_pressure'] * ureg.psi

eta_cstar = p['cstar_efficiency']      # C* efficiency factor
eta_cf    = p['cf_efficiency']         # Thrust coefficient efficiency factor
L_star    = p['L_star'] * ureg.m       # Characteristic chamber length
alpha     = p['alpha_divergence'] * ureg.deg  # Conical chamber half-angle
of_ratio  = p['of_ratio']              # Oxidizer/fuel mixture ratio

# --------------------------------------------------------------------------- #
# 2. CEA (combustion equilibrium) setup
# --------------------------------------------------------------------------- #

### CEA Setup ###

reac_names      = [p['propellants']['fuel'], p['propellants']['oxidizer']]
T_reactant      = np.array([p['reactant_T']['fuel'], p['reactant_T']['oxidizer']])
fuel_weights    = np.array([1.0, 0.0])
oxidant_weights = np.array([0.0, 1.0])

reac = cea.Mixture(reac_names)
prod = cea.Mixture(reac_names, products_from_reactants=True)
solver   = cea.RocketSolver(prod, reactants=reac, transport=True)
solution = cea.RocketSolution(solver)

# Mass fractions for the target O/F ratio, and reactant enthalpy (mixed,
# non-dimensionalized by the gas constant) needed for the CEA solve calls.
weights = reac.of_ratio_to_weights(oxidant_weights, fuel_weights, of_ratio)
hc      = reac.calc_property(cea.ENTHALPY, weights, T_reactant) / cea.R  # plain float

# --------------------------------------------------------------------------- #
# 3. Run CEA 
# --------------------------------------------------------------------------- #

# Pressure ratio (chamber/exit) needed for the CEA IAC (infinite-area
# combustor) run that gives us C* and Cf at this chamber pressure.
pi_p = [pc.to(ureg.psi).magnitude / pe.to(ureg.psi).magnitude]
pc = pc.to(ureg.bar)

# --- IAC run: get C* and Cf for this chamber pressure ---
solver.solve(solution, weights, pc.magnitude, pi_p, iac=True, hc=hc)

Cf    = solution.coefficient_of_thrust[-1] * eta_cf
cstar = solution.c_star[-1] * eta_cstar * ureg.m / ureg.s  # m/s
Isp   = cstar * Cf / ureg.g0

# --- Throat sizing from thrust, C*, and Cf ---
mdot = F.to(ureg.N) / (Cf * cstar)                    # kg/s
At   = (mdot * cstar) / pc.to(ureg.Pa)        # m^2 (throat area)

A_c   = np.pi / 4 * Dc.to(ureg.m) ** 2   # chamber cross-sectional area, m^2
ac_at = A_c / At                          # contraction ratio

# --- Finite-area-combustor run with this contraction ratio ---
# (needed to get the correct expansion ratio ae_at / T for this ac_at)
solver.solve(solution, weights, pc.magnitude, pi_p, ac_at=ac_at, iac=False, hc=hc)

# --- Nozzle geometry ---
Dt = 2 * np.sqrt(At / np.pi)                # throat diameter, m
Ae = solution.ae_at[-1] * At.to(ureg.m ** 2)  # exit area, m^2
De = 2 * np.sqrt(Ae / np.pi)                # exit diameter, m

# --- Chamber geometry: conical converging section + cylindrical
#     section sized to hit the target L* (characteristic length) ---
L_cone = (Dc.to(ureg.m) / 2 - Dt.to(ureg.m) / 2) / np.tan(alpha.to(ureg.rad))
V_cone = (np.pi / 3) * ((Dc.to(ureg.m) / 2) ** 2+ Dc.to(ureg.m) / 2 * Dt.to(ureg.m) / 2 + (Dt.to(ureg.m) / 2) ** 2) * L_cone

Vc_total = L_star.to(ureg.m) * At.to(ureg.m ** 2)  # total required chamber volume
V_cyl    = Vc_total - V_cone                       # remaining volume -> cylindrical section

Lc     = V_cyl / A_c.to(ureg.m ** 2)   # cylindrical section length, m
Ltotal = Lc + L_cone                   # total chamber length, m

values = {
    'F_lbf': F.to(ureg.lbf).magnitude,
    'pc_psi': pc.to(ureg.psi).magnitude,
    'D_chamber_in': Dc.to(ureg.inch).magnitude,
    't_chamber_in': tc.to(ureg.inch).magnitude,
    'Lc_m': Lc.to(ureg.m).magnitude,
    'Ltotal_m': Ltotal.to(ureg.m).magnitude,
    'ac_at': ac_at.magnitude,
    'ae_at': solution.ae_at[-1],
    'mdot_kg/s': mdot.to(ureg.kg / ureg.s).magnitude,
    'Dt_in': Dt.to(ureg.inch).magnitude,
    'De_in': De.to(ureg.inch).magnitude,
    'C_star_ms': cstar.to(ureg.m / ureg.s).magnitude,
    'Cf': Cf,
    'Isp_s': Isp.to(ureg.s).magnitude,
    'T_chamber_K': solution.T[0],
}

for name, value in values.items():
    print(f"{name:20s}: {value}")