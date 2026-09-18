# Torch_Analysis
import numpy as np
from pyfluids import Fluid, FluidsList, Input
from New_TCA_Igniter_Sizing import diameter_to_area
from New_TCA_Igniter_Sizing import area_to_diameter

R = 8.314 # Gas constant
CM_TO_IN = 1/2.54 # Centimeters to Inches
IN_TO_CM = 2.54 # Inches to Centimeters 

## Fixed torch parameters
A_t = 1    # Throat area
Df = 1 * IN_TO_CM     # Fuel injector diameter (input in inches)
Dox = 1 * IN_TO_CM    # Oxidizer injector diameter
Tc = 1     # Chamber temperature
Tox = 22    # Oxidizer line temperature
Tf = 22    # Fuel line temperature
k = 1     # gamma
k_ox = 1      # Oxidizer gamma
k_fuel = 1     # Fuel gamma
mdot = 1 # Mass flow rate
of_ratio = 1 # O/F Ratio

def chamber_pressure(throatArea, k, massflow, chamberTemp):
#   Input:
    # Throat area
    # Heat capacity ratio
    # Mass flow rate
    # Chamber temperature
#   Output:
    # Chamber Pressure
    return (massflow/(throatArea*k))*((np.sqrt(k*R*chamberTemp))/(np.sqrt((2/(k+1))**((k+1)/(k-1)))))



mdot_ox = mdot * (of_ratio)/(of_ratio + 1)
mdot_f = mdot - mdot_ox

Pc = chamber_pressure(A_t, k, mdot, Tc)
Pox = chamber_pressure(diameter_to_area(Dox), k, mdot_ox, Tox)
Pf = chamber_pressure(area_to_diameter(Df), k, mdot_f, Tf)