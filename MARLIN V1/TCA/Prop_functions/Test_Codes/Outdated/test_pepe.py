import numpy as np
import cea
from pint import UnitRegistry

ureg = UnitRegistry()

reac_names = ["C3H8O,2propanol", "O2(L)"]
T_reactant = np.array([298.15, 90.17]) * ureg.K
fuel_weights = np.array([1.0, 0.0])
oxidant_weights = np.array([0.0, 1.0])
of_ratio = 1.6

pc = np.array([15, 60]) * ureg.bar  # Chamber pressure (psi)
patm = 1.01325 * ureg.bar  # Atmospheric pressure (bar)
pi_p = [pc / patm]   # Pressure ratio chamber to exit - Perfectly expanded engine
T = 5000 * ureg.lbf  # Target thrust (lbf)

reac = cea.Mixture(reac_names)
prod = cea.Mixture(reac_names, products_from_reactants=True)

solver = cea.RocketSolver(prod, reactants=reac)
solution = cea.RocketSolution(solver)

weights = reac.of_ratio_to_weights(oxidant_weights, fuel_weights, of_ratio) #Convert OF to weights.

#Compute chamber enthalpy. Normalized.
hc = reac.calc_property(cea.ENTHALPY, weights, T_reactant.to(ureg.kelvin).magnitude)/cea.R

solver.solve(solution, weights, pc.to(ureg.bar).magnitude, pi_p, iac=False, hc=hc)

