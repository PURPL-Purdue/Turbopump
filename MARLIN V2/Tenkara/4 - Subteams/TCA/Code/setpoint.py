"""
For a grid of target thrust levels, chamber pressures, and chamber
diameters/wall thicknesses, this script:
  1. Runs NASA CEA
  2. Sizes the throat/exit/chamber geometry from those results.
  3. Estimates chamber wall stress (thin-wall pressure vessel, von Mises).
  4. Generates a bell nozzle contour and its arc length.
  5. Filters the resulting design matrix against length/diameter limits.
  6. Writes the surviving designs to Outputs/designs.csv.
"""

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
# 1. Load engine inputs
# --------------------------------------------------------------------------- #

with open('Inputs/setpoint_params.yaml') as f:
    p = yaml.safe_load(f)

### Engine Parameters/Inputs ###

# Sweep ranges
F = np.arange(
    p['Thrust_array']['Min'],
    p['Thrust_array']['Max'],
    p['Thrust_array']['Step'],
) * ureg.lbf

pc_array = np.arange(
    p['Chamber_pressure_array']['Min'],
    p['Chamber_pressure_array']['Max'],
    p['Chamber_pressure_array']['Step'],
) * ureg.psi

# Candidate chamber geometries (paired diameter/thickness)
chamber_diameters = np.array(p['Chamber_diameters_in']) * ureg.inch
chamber_thickness = np.array(p['Chamber_thicknesses_in']) * ureg.inch

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
# 3. Build the design matrix: sweep pc -> thrust -> chamber geometry
# --------------------------------------------------------------------------- #

### Build design matrix ###
designs = []

for pc in pc_array:
    # Pressure ratio (chamber/exit) needed for the CEA IAC (infinite-area
    # combustor) run that gives us C* and Cf at this chamber pressure.
    pi_p = [pc.to(ureg.psi).magnitude / pe.to(ureg.psi).magnitude]
    pc = pc.to(ureg.bar)

    # --- IAC run: get C* and Cf for this chamber pressure ---
    solver.solve(solution, weights, pc.magnitude, pi_p, iac=True, hc=hc)

    Cf    = solution.coefficient_of_thrust[-1] * eta_cf
    cstar = solution.c_star[-1] * eta_cstar * ureg.m / ureg.s  # m/s
    Isp   = cstar * Cf / ureg.g0

    for F_target in F:
        # --- Throat sizing from thrust, C*, and Cf ---
        F_N  = F_target.to(ureg.N)
        mdot = F_N / (Cf * cstar)                    # kg/s
        At   = (mdot * cstar) / pc.to(ureg.Pa)        # m^2 (throat area)

        for D_c, t_c in zip(chamber_diameters, chamber_thickness):
            A_c   = np.pi / 4 * D_c.to(ureg.m) ** 2   # chamber cross-sectional area, m^2
            ac_at = A_c / At                          # contraction ratio

            # --- Finite-area-combustor run with this contraction ratio ---
            # (needed to get the correct expansion ratio ae_at / T for this ac_at)
            solver.solve(solution, weights, pc.magnitude, pi_p,
                         ac_at=ac_at, iac=False, hc=hc)

            # --- Nozzle geometry ---
            Dt = 2 * np.sqrt(At / np.pi)                # throat diameter, m
            Ae = solution.ae_at[-1] * At.to(ureg.m ** 2)  # exit area, m^2
            De = 2 * np.sqrt(Ae / np.pi)                # exit diameter, m

            # --- Chamber geometry: conical converging section + cylindrical
            #     section sized to hit the target L* (characteristic length) ---
            L_cone = (D_c.to(ureg.m) / 2 - Dt.to(ureg.m) / 2) / np.tan(alpha.to(ureg.rad))
            V_cone = (np.pi / 3) * (
                (D_c.to(ureg.m) / 2) ** 2
                + D_c.to(ureg.m) / 2 * Dt.to(ureg.m) / 2
                + (Dt.to(ureg.m) / 2) ** 2
            ) * L_cone

            Vc_total = L_star.to(ureg.m) * At.to(ureg.m ** 2)  # total required chamber volume
            V_cyl    = Vc_total - V_cone                       # remaining volume -> cylindrical section

            Lc     = V_cyl / A_c.to(ureg.m ** 2)   # cylindrical section length, m
            Ltotal = Lc + L_cone                   # total chamber length, m

            # --- Thin-wall pressure vessel stress check (von Mises) ---
            sigma_th = (pc.to(ureg.Pa) * D_c.to(ureg.m)) / (2 * t_c.to(ureg.m))          # hoop stress (seamless pipe)
            sigma_ax = (pc.to(ureg.Pa) * D_c.to(ureg.m)) / (4 * t_c.to(ureg.m) * 0.6)    # axial stress (weld coeff. 0.6)
            sigma_vM = np.sqrt(sigma_th ** 2 + sigma_ax ** 2 - sigma_th * sigma_ax)      # von Mises stress

            # --- Bell nozzle contour + total engine length ---
            angles, contour, R2 = bell_nozzle(
                solution.ae_at[-1],
                Dt.to(ureg.mm).magnitude / 2,
                80,
                ac_at,
                alpha.to(ureg.deg).magnitude,
                Lc.to(ureg.mm).magnitude,
            )
            total_length_mm = contour_length_mm(contour)

            designs.append({
                'F_lbf'             : F_target.to(ureg.lbf).magnitude,
                'pc_psi'           : pc.to(ureg.psi).magnitude,
                'D_chamber_in'     : D_c.to(ureg.inch).magnitude,
                't_chamber_in'     : t_c.to(ureg.inch).magnitude,
                'Lc_m'             : Lc.to(ureg.m).magnitude,
                'ac_at'            : ac_at.magnitude,
                'ae_at'            : solution.ae_at[-1],
                'mdot_kg/s'        : mdot.to(ureg.kg / ureg.s).magnitude,  # kg/s
                'Dt_in'            : Dt.to(ureg.inch).magnitude,          # in
                'De_in'            : De.to(ureg.inch).magnitude,          # in
                'C_star_ms'        : cstar.to(ureg.m / ureg.s).magnitude, # m/s
                'Cf'               : Cf,
                'Isp_s'            : Isp.to(ureg.s).magnitude,            # s
                'T_chamber_K'      : solution.T[0],
                'sigma_vM_MPa'     : sigma_vM.to(ureg.MPa).magnitude,     # MPa
                'total_length_mm'  : total_length_mm,
            })

df = pd.DataFrame(designs)

# --------------------------------------------------------------------------- #
# 4. Filter the design matrix against size constraints and save
# --------------------------------------------------------------------------- #

# Size thresholds
max_L_m = p['max_L'] * ureg.m   # meters
min_L_m = p['min_L'] * ureg.m   # meters
max_De  = p['max_D'] * ureg.m   # meters (compared against De in inches below)

# Filter out designs exceeding thresholds
df = df[
    (df['total_length_mm'] <= max_L_m.to(ureg.mm).magnitude)
    & (df['total_length_mm'] >= min_L_m.to(ureg.mm).magnitude)
    & (df['De_in'] <= max_De.to(ureg.inch).magnitude)
]

print(f"Designs after filtering: {len(df)}")
print(df.to_string())
df.to_csv(p['output_directory'], index=False)