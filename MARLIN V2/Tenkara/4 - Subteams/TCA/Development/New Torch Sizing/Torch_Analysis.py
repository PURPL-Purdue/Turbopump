# Torch_Analysis
import numpy as np
from pyfluids import Fluid, FluidsList, Input
from New_TCA_Igniter_Sizing import diameter_to_area
from New_TCA_Igniter_Sizing import area_to_diameter
from New_TCA_Igniter_Sizing import choked_backpressure
import cea as cea


R = 8.314 # Gas constant
CM_TO_IN = 1/2.54 # Centimeters to Inches
IN_TO_CM = 2.54 # Inches to Centimeters
BAR_TO_PA = 1e5 # Bar to Pascals 
T_AMB_CELSIUS = 20 # [deg C]
T_AMB_KELVIN = 293.15 # [deg K]



def chamber_pressure(throatArea, k, massflow, chamberTemp):
#   Input:
    # Throat area
    # Heat capacity ratio
    # Mass flow rate
    # Chamber temperature
#   Output:
    # Chamber Pressure
    return (massflow/(throatArea*k))*((np.sqrt(k*R*chamberTemp))/(np.sqrt((2/(k+1))**((k+1)/(k-1)))))

def get_k(p_c, OF):

    reac_names = ["CH4", "O2"]
    T_reactant = np.array([T_AMB_KELVIN, T_AMB_KELVIN])
    fuel_weights = np.array([1.0, 0.0])
    ox_weights = np.array([0.0, 1.0])
    p_c = p_c / BAR_TO_PA
    reac = cea.Mixture(reac_names)
    prod = cea.Mixture(reac_names, products_from_reactants=True)
    solver = cea.RocketSolver(prod, reactants=reac)
    solution = cea.RocketSolution(solver)
    weights = reac.of_ratio_to_weights(ox_weights, fuel_weights, of_ratio=OF)
    hc = reac.calc_property(cea.ENTHALPY, weights, T_reactant) / cea.R
    solver.solve( solution, weights, p_c, hc=hc, iac=True)
    k_c = solution.gamma_s[0]

    return k_c

def critical_pressure(k):
#   Input:
    # gamma
#   Output:
    # Critical Pressure
    return (((k+1)/(2))**(k/(k-1)))

def solve_chamber_pressure(throatArea, massflow, chamberTemp, OF, k):

    for i in range(100):
        p_c = chamber_pressure(throatArea, k, massflow, chamberTemp)
        k_new = get_k(p_c, OF)

        if abs(k_new - k) < 1e-5:
            k = k_new
            break

        k = k_new
        
    return p_c, k

def main():
    ## Fixed torch parameters
    A_t = 1    # Throat area
    Df = 1 * IN_TO_CM     # Fuel injector diameter (input in inches)
    Dox = 1 * IN_TO_CM    # Oxidizer injector diameter
    Tc = 300     # Chamber temperature
    Tox = 22    # Oxidizer line temperature
    Tf = 22    # Fuel line temperature
    k_c = 1.2     # gamma
    k_ox = 1.40      # Oxidizer gamma
    k_fuel = 1.31     # Fuel gamma
    mdot = 1 # Mass flow rate
    of_ratio = 1 # O/F Ratio
    stiffness = 0.3 # Stiffness of the system

    mdot_ox = mdot * (of_ratio)/(of_ratio + 1)
    mdot_f = mdot - mdot_ox

    k_c, Pc = solve_chamber_pressure(A_t, mdot, Tc, of_ratio, k_c)
    Pox = chamber_pressure(diameter_to_area(Dox), k_ox, mdot_ox, Tox)
    Pf = chamber_pressure(diameter_to_area(Df), k_fuel, mdot_f, Tf)

    Pcrit_ox = critical_pressure(k_ox)
    Pcrit_f = critical_pressure(k_fuel)
            
    if Pox < Pcrit_ox:
        print("Oxidizer is choked")
        pox_line = choked_backpressure(k_ox, stiffness, Pc)
    else:
        print("Oxidizer is not choked")
        pox_ratio = 1 / (1 + stiffness)   # pc/pline [--]
        pox_line = Pc / pox_ratio          # line pressure [Pa]

    if Pf < Pcrit_f:
        print("Fuel is choked")
        pf_line = choked_backpressure(k_fuel, stiffness, Pc)
    else:
        print("Fuel is not choked")
        pf_ratio = 1 / (1 + stiffness)   # pc/pline [--]
        pf_line = Pc / pf_ratio          # line pressure [Pa]

    print("Chamber Pressure: ", Pc / BAR_TO_PA, " bar")
    print("Oxidizer Line Pressure: ", pox_line / BAR_TO_PA, " bar")
    print("Fuel Line Pressure: ", pf_line / BAR_TO_PA, " bar")



