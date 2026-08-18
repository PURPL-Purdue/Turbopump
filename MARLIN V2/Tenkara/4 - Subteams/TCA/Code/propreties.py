import pandas as pd
import numpy as np
from pint import UnitRegistry
import cea
import yaml


#######################################################################################
# FUNCTIONS
# =====================================================================================

# =====================================================================================
# Fit a degree-N polynomial (property vs. temperature)
# =====================================================================================

def generate_poly_coeffs_csv(properties_df, of_ratio, pc_psia,
                              output_path, degree):
    property_map = {
        'Gamma':                            ('gamma', 1.0),
        'Cp_KJ_kgK':                        ('Cp [(KJ/kg-K)]', 1.0),   
        'Thermal_Conductivity_W_mK':        ('k [W/m-K]', 1.0),
        'Viscosity_Pa_s':                   ('Viscosity [Pa*s]', 1.0),     
        'Molecular_Weight_kg_kmol':         ('MW [kg/kmol]', 1.0),
    }

    T_all = properties_df['T_chamber [K]'].astype(float).values
    T_unique, unique_idx = np.unique(T_all, return_index=True)

    n_pts = len(T_unique)
    if n_pts < degree + 1:
        raise ValueError(
            f"Only {n_pts} unique temperature points available but degree={degree} "
            f"needs at least {degree + 1}. Lower the degree or add more CEA solve points."
        )

    T_points = [float(t) for t in T_unique]

    csv_rows = []
    for prop_name, (col, factor) in property_map.items():
        if col not in properties_df.columns:
            print(f"Warning: column '{col}' not found in properties_df, skipping '{prop_name}'")
            continue

        Y_all = properties_df[col].astype(float).values * factor
        Y_unique = Y_all[unique_idx]

        coeffs = np.polyfit(T_unique, Y_unique, degree)  # high -> low
        coeffs = [float(c) for c in coeffs]              # strip np.float64

        Y_points = [float(y) for y in Y_unique]

        csv_rows.append({
            'O/F': of_ratio,
            'Pc_psia': pc_psia,
            'Property': prop_name,
            'Degree': degree,
            'Coefficients_high->low': coeffs,
            'T_points_K': T_points,
            'Y_points': Y_points,
        })

    out_df = pd.DataFrame(csv_rows)
    out_df.to_csv(output_path, index=False)
    return out_df


# =====================================================================================
# Function to get positional areas from CSV files
# =====================================================================================

def get_positional_areas(csv_files):
    if not csv_files:
        print("No CSV files provided.")
        return

    for path in csv_files:
        df = pd.read_csv(path)
        df = df.sort_values("x_mm")
        x = df["x_mm"].values * ureg.mm
        y = df["y_mm"].values * ureg.mm  
    return x.to(ureg.m).magnitude, y.to(ureg.m).magnitude

with open('Inputs/TCA_params.yaml') as f:
    p = yaml.safe_load(f)
ureg = UnitRegistry()


#######################################################################################
# ENGINE PAREMETERS/INPUTS
# =====================================================================================

F     = p['Thrust_target']  * ureg.lbf      # Target thrust             [lbf]
pc    = p['chamber_pressure'] * ureg.psi    # Chamber pressure          [psia]
pe    = p['exit_pressure'] * ureg.psi       # Exit pressure             [psia]
ac_at = p['contraction_ratio']             # Contraction ratio
of_ratio = p['of_ratio']

# =====================================================================================
# Get expansion ratio from contour
# =====================================================================================

x, y = get_positional_areas(["Outputs/contour.csv"])
rt = np.interp(0.0, x, y) * ureg.m
ae_at = np.array((y / rt.magnitude)**2)

# =====================================================================================
# CEA Setup
# =====================================================================================

reac_names = [p['propellants']['fuel'], p['propellants']['oxidizer']]
T_reactant = np.array([298.15, 90.17]) * ureg.K
fuel_weights = np.array([1.0, 0.0])
oxidant_weights = np.array([0.0, 1.0])

reac = cea.Mixture(reac_names)
prod = cea.Mixture(reac_names, products_from_reactants=True)

solver = cea.RocketSolver(prod, reactants=reac, transport=True)
solution = cea.RocketSolution(solver)

weights = reac.of_ratio_to_weights(oxidant_weights, fuel_weights, of_ratio) #Convert OF to weights.
hc = reac.calc_property(cea.ENTHALPY, weights, T_reactant.to(ureg.kelvin).magnitude)/cea.R

# =====================================================================================
# Set propreties array, and sort area ratios and positions into subsonic and supersonic
# =====================================================================================

properties = []

throat_i = np.where(ae_at == ae_at.min())[0]

x_before = x[:throat_i[0]]
areas_before_throat = ae_at[:throat_i[0]]

x_after = x[throat_i[-1] + 1:]
areas_after_throat = ae_at[throat_i[-1] + 1:]

# =====================================================================================
# Run CEA solver for each area ratio, subsonic and supersonic solutions
# =====================================================================================

for xpos, subar in zip(x_before, areas_before_throat): 
    solver.solve(solution, weights, pc.to(ureg.bar).magnitude, subar=[subar], ac_at=ac_at, iac=False, hc=hc)
    print(solution.P)
    properties.append({
        'x_m': xpos,
        'gamma' : solution.gamma_s[-1],
        'Cp [(KJ/kg-K)]' : solution.cp[-1],
        'k [W/m-K]' : (solution.conductivity_eq[-1] * ureg.watt / (ureg.centimeter * ureg.kelvin)).to(ureg.watt / (ureg.meter * ureg.kelvin)).magnitude,
        'MW [kg/kmol]' : solution.MW[-1],
        'R [kJ/kg-K]' : cea.R / (solution.MW[-1]),
        'T_chamber [K]'  : solution.T[-1],
        'P_chamber [bar]' : solution.P[-1],
        'ae_at' : solution.ae_at[-1],
        'Viscosity [Pa*s]': (solution.viscosity[-1] * ureg.millipoise).to(ureg.pascal * ureg.second).magnitude
        })

for xpos, supar in zip(x_after, areas_after_throat):
    solver.solve(solution, weights, pc.to(ureg.bar).magnitude, supar=[supar], ac_at=ac_at, iac=False, hc=hc)

    properties.append({
        'x_m': xpos,
        'gamma' : solution.gamma_s[-1],
        'Cp [(KJ/kg-K)]' : solution.cp[-1],
        'k [W/m-K]' : (solution.conductivity_eq[-1] * ureg.watt / (ureg.centimeter * ureg.kelvin)).to(ureg.watt / (ureg.meter * ureg.kelvin)).magnitude,
        'MW [kg/kmol]' : solution.MW[-1],
        'R [kJ/kg-K]' : cea.R / (solution.MW[-1]),
        'T_chamber [K]'  : solution.T[-1],
        'P_chamber [bar]' : solution.P[-1],
        'ae_at' : solution.ae_at[-1],
        'Viscosity [Pa*s]': (solution.viscosity[-1] * ureg.millipoise).to(ureg.pascal * ureg.second).magnitude
        })

# =====================================================================================
# Save all propreties into a csv file and generate polynomial coefficients for each property
# =====================================================================================

df = pd.DataFrame(properties)
df.to_csv('Outputs/properties.csv', index=False)

generate_poly_coeffs_csv(
    properties_df=df,
    of_ratio=of_ratio,
    pc_psia=pc.magnitude,
    output_path='Outputs/turbopump_poly_coeffs.csv',
    degree=2
)