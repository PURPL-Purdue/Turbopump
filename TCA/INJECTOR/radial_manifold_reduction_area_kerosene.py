import numpy as np
import CoolProp.CoolProp as CP
from matplotlib import pyplot as plt

kg_into_lbm = 2.20462 # Conversion factor from kg to lbm
psi_into_pa = 6894.76 # Conversion factor from psi to Pa
meters_into_inches = 39.37 # Conversion factor from meters to inches
meters_into_foot = 3.28084 # Conversion factor from meters to feet
metersSquare_into_inchesSquare = 1550 # Conversion factor from m^2 to in^2

mdot = 9 #[kg/s]
OF_Ratio = 2.1 # O/F ratio 
mdot_kero = mdot/(1+OF_Ratio) #Fuel mass flow (RP-1), from global O/F
P_in = 710.12 # [psi]
rho_rp1 = CP.PropsSI('D','T',298,'P',P_in*psi_into_pa,'Dodecane') # Density of RP-1 at inlet conditions
number_of_elements_type1 = 8 # Number of two mini-manifolds (type1: 6 holes)
number_of_elements_type2 = 8 # Number of two mini-manifolds (type2: 4 holes)
number_of_elements = number_of_elements_type1 + number_of_elements_type2 # Total number of elements in the manifold
number_of_holes = number_of_elements_type1*6 + number_of_elements_type2*4
mass_flow_rate_per_element = mdot_kero/number_of_holes # Mass flow rate through

# Dimensions of the manifold
length_manifold = 16.70e-3 # Length of the manifold [m]
length_inch = length_manifold * meters_into_inches # Length of the manifold [inches]    

def calculate_area_step(initial_area, total_mass_flow_rate, mass_flow_rate_per_element, density, number_steps):
    area_steps = [initial_area] # creating an array to hold area values
    mass_flow_rate_steps = [total_mass_flow_rate] # creating an array to hold mass flow rate values
    print(f"initial_mass_flow_rate {total_mass_flow_rate*kg_into_lbm:.4f} lbm/s\n")
    mass_flow_rate_steps.append(total_mass_flow_rate/2) # first step mass flow rate. The flow divides into two branches
    print(f"initial_area {initial_area*metersSquare_into_inchesSquare:.4f} in^2\n")
    height_manifold = []
    for _ in range(number_steps):
        current_mass_flow_rate = mass_flow_rate_steps[-1]
        new_area = current_mass_flow_rate / (density*velocity) # new area calculation
        area_steps.append(new_area)
        mass_flow_rate_steps.append(current_mass_flow_rate - 4*mass_flow_rate_per_element)
        print(f"new_area: {new_area*metersSquare_into_inchesSquare:.4f} in^2\n")
        print(f"current_mass_flow_rate: {current_mass_flow_rate*kg_into_lbm:.4f} lbm/s\n")
        height = new_area/length_manifold # height of the manifold at this step
        height_manifold.append(height)
        print(f"height of the manifold at this step: {height*meters_into_inches:.4f} inches\n")
        current_mass_flow_rate = mass_flow_rate_steps[-1]
        new_area = current_mass_flow_rate / (density*velocity) # new area calculation
        area_steps.append(new_area)
        mass_flow_rate_steps.append(current_mass_flow_rate - 6*mass_flow_rate_per_element)
        print(f"new_area: {new_area*metersSquare_into_inchesSquare:.4f} in^2\n")
        print(f"current_mass_flow_rate: {current_mass_flow_rate*kg_into_lbm:.4f} lbm/s\n")
        print("----\n")
        height = new_area/length_manifold # height of the manifold at this step
        height_manifold.append(height)
        print(f"height of the manifold at this step: {height*meters_into_inches:.4f} inches\n")
    return area_steps, mass_flow_rate_steps

dynamics_pressure = 0.01 * P_in*psi_into_pa # assuming 1% pressure as dynamic pressure
velocity = np.sqrt(2*dynamics_pressure/rho_rp1) # velocity from dynamic pressure [m/s] - must be constant throughout the manifold
initial_area = mdot_kero/(rho_rp1*velocity) # initial area assuming all flow at the start [m^2]
print(f"The velocity would be {velocity*meters_into_foot:.4f} ft/s thoughout all the manifold.\n")

number_steps = number_of_elements/4 # number of steps in the manifold
area_list, mass_flow_rate_list = calculate_area_step(initial_area, mdot_kero, mass_flow_rate_per_element, rho_rp1, int(number_steps))

# PRINT RESUTLS
print(f"Initial area: {initial_area*metersSquare_into_inchesSquare:.4f} in^2\n")
print(f"Final area: {area_list[-1]*metersSquare_into_inchesSquare:.4f} in^2\n")
print(f"Number of steps: {number_steps}\n")
print(f"Mass flow rate per element: {mass_flow_rate_per_element*kg_into_lbm:.4f} lbm/s\n")