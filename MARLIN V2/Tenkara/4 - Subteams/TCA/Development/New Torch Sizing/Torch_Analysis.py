# Torch_Analysis
import numpy as np
from pyfluids import Fluid, FluidsList, Input
from New_TCA_Igniter_Sizing import diameter_to_area
from New_TCA_Igniter_Sizing import area_to_diameter
import cea as cea


R = 8.314 # Gas constant
CM_TO_IN = 1/2.54 # Centimeters to Inches
IN_TO_CM = 2.54 # Inches to Centimeters 
T_AMB_CELSIUS = 20 # [deg C]
T_AMB_KELVIN = 293.15 # [deg K]

## Fixed torch parameters
A_t = 1    # Throat area
Df = 1 * IN_TO_CM     # Fuel injector diameter (input in inches)
Dox = 1 * IN_TO_CM    # Oxidizer injector diameter
Tc = 300     # Chamber temperature
Tox = 22    # Oxidizer line temperature
Tf = 22    # Fuel line temperature
k = 1.2     # gamma
k_ox = 1      # Oxidizer gamma
k_fuel = 1     # Fuel gamma
mdot = 1 # Mass flow rate
of_ratio = 1 # O/F Ratio


reac_names = ["CH4", "O2"]
T_reactant = np.array([T_AMB_KELVIN, T_AMB_KELVIN])  # Reactant temperatures (K)
fuel_weights = np.array([1.0, 0.0])
ox_weights = np.array([0.0, 1.0])

reac_names = ["CH4", "O2"]
reac = cea.Mixture(reac_names)
prod = cea.Mixture(reac_names, products_from_reactants=True)
solver = cea.RocketSolver(prod, reactants=reac)
solution = cea.RocketSolution(solver)
weights = reac.of_ratio_to_weights(ox_weights, fuel_weights, of_ratio=2.5)
hc = reac.calc_property(cea.ENTHALPY, weights, T_reactant)/cea.R

solver.solve(solution, weights, p_c, hc=hc, iac=True)

isp = solution.Isp[1] / G0
rho_c = solution.density[0]      #unknown units
k_c = solution.gamma_s[0]

def chamber_pressure(throatArea, k, massflow, chamberTemp):
#   Input:
    # Throat area
    # Heat capacity ratio
    # Mass flow rate
    # Chamber temperature
#   Output:
    # Chamber Pressure
    return (massflow/(throatArea*k))*((np.sqrt(k*R*chamberTemp))/(np.sqrt((2/(k+1))**((k+1)/(k-1)))))

def critical_pressure(k):
#   Input:
    # gamma
#   Output:
    # Critical Pressure
    return (((k+1)/(2))**(k/(k-1)))

mdot_ox = mdot * (of_ratio)/(of_ratio + 1)
mdot_f = mdot - mdot_ox

Pc = chamber_pressure(A_t, k, mdot, Tc)
Pox = chamber_pressure(diameter_to_area(Dox), k, mdot_ox, Tox)
Pf = chamber_pressure(area_to_diameter(Df), k, mdot_f, Tf)

Pcrit_ox = critical_pressure(k_ox)
Pcrit_f = critical_pressure(k_fuel)
Pcrit_chamber = critical_pressure(k)

