import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import os
from pint import UnitRegistry
import cea
import yaml

#### Function to get positional areas from CSV files ####
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

with open('TCA_params_test.yaml') as f:
    p = yaml.safe_load(f)
ureg = UnitRegistry()

### Engine Parameters/Inputs

F     = p['Thrust_target']  * ureg.lbf      # Target thrust             [lbf]
pc    = p['Chamber_pressure'] * ureg.psi    # Chamber pressure          [psia]
pe    = p['Exit_pressure'] * ureg.psi       # Exit pressure             [psia]
ac_at = p['contraction_ratio']             # Contraction ratio
of_ratio = p['of_ratio']

### Get expansion ratio from contour ###
x, y = get_positional_areas(["CSV_DXF_OUTPUTS/nozzle_contour2600-600psi.csv"])
rt = np.interp(0.0, x, y) * ureg.m
ae_at = np.array((y / rt.magnitude)**2)

### CEA Setup ###
reac_names = ["C2H5OH(L)", "O2(L)"]
T_reactant = np.array([298.15, 90.17]) * ureg.K
fuel_weights = np.array([1.0, 0.0])
oxidant_weights = np.array([0.0, 1.0])

reac = cea.Mixture(reac_names)
prod = cea.Mixture(reac_names, products_from_reactants=True)

solver = cea.RocketSolver(prod, reactants=reac, transport=True)
solution = cea.RocketSolution(solver)

weights = reac.of_ratio_to_weights(oxidant_weights, fuel_weights, of_ratio) #Convert OF to weights.
hc = reac.calc_property(cea.ENTHALPY, weights, T_reactant.to(ureg.kelvin).magnitude)/cea.R

### Set propreties array, and sort area ratios and positions into subsonic and supersonic ###
properties = []

throat_i = np.where(ae_at == ae_at.min())[0]

x_before = x[:throat_i[0]]
areas_before_throat = ae_at[:throat_i[0]]

x_after = x[throat_i[-1] + 1:]
areas_after_throat = ae_at[throat_i[-1] + 1:]

### Run CEA solver for each area ratio, subsonic and supersonic solutions ###
for xpos, subar in zip(x_before, areas_before_throat): 
    solver.solve(solution, weights, pc.to(ureg.bar).magnitude, subar=[subar], ac_at=ac_at, iac=False, hc=hc)
    print(solution.P)
    properties.append({
        'x_m': xpos,
        'gamma' : solution.gamma_s[-1],
        'Cp [(KJ/kg-K)]' : solution.cp[-1],
        'k [uW/cm-K]' : solution.conductivity_eq[-1],
        'MW [kg/kmol]' : solution.MW[-1],
        'R [kJ/kg-K]' : cea.R / (solution.MW[-1]),
        'T_chamber [K]'  : solution.T[-1],
        'P_chamber [bar]' : solution.P[-1],
        'ae_at' : solution.ae_at[-1]})

for xpos, supar in zip(x_after, areas_after_throat):
    solver.solve(solution, weights, pc.to(ureg.bar).magnitude, supar=[supar], ac_at=ac_at, iac=False, hc=hc)
    #print(solution.P)

    properties.append({
        'x_m': xpos,
        'gamma' : solution.gamma_s[-1],
        'Cp [(KJ/kg-K)]' : solution.cp[-1],
        'k [uW/cm-K]' : solution.conductivity_eq[-1],
        'MW [kg/kmol]' : solution.MW[-1],
        'R [kJ/kg-K]' : cea.R / (solution.MW[-1]),
        'T_chamber [K]'  : solution.T[-1],
        'P_chamber [bar]' : solution.P[-1],
        'ae_at' : solution.ae_at[-1]})
    

### Save all propreties into an excel ###
df = pd.DataFrame(properties)
df.to_excel('properties.xlsx', index=False)
