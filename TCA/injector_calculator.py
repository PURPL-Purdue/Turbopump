"""
This program computes the geometry and flow characteristics required to design a liquid
rocket engine injector element consisting of a fuel jet and an oxidizer jet that collide
at a prescribed inpingement angle and location.

Author: Joaquin Alarcon
"""

import numpy as np
import math
import matplotlib.pyplot as plt

# == Conversion between units ==
psi_into_pa = 6894.76         # Convert psi -> Pascal
meters_into_inches = 39.37    # Convert meters -> inches
degrees_into_rad = np.pi/180  # Convert degrees -> radians

# == Desing Parameters (Change as needed) == 
mdot = 10 #Total propellant mass flow rate [kg/s] (CHANGE AS NEEDED)
OF_Ratio = 2.1 #Mixture ratio O/F 
rho_rp1 = 810 #RP1 Density [kg/m^3] at injector conditions
rho_lox = 1141 #LOX Density [kg/m^3] at injector conditions
inj_press_drop = 0.6 #Fraction of chamber pressure allocated to injector
P0 = 500*psi_into_pa #Chamber stagnation pressure [Pa]
Cd = 0.7 #Discharge coefficient
delta_P = 300 #Injector pressure drop [psi] (converted inside functions)
theta_1 = 40 #Initial guess / design choice for injector angle [deg]
Longitude_chamber = 319.05 #Combustion chamber length [mm]
impinge_fraction = 0.10 #streams will impinge at 10% of chamber length
#ratio_inj_cooling = 1 #Fraction of fuel mass flow sent to cooling vs. injection 

mdot_kero = mdot/(1+OF_Ratio) #Fuel mass flow (RP-1), from global O/F
mdot_lox = mdot*OF_Ratio/(1+OF_Ratio) #Oxidizer mass flow (LOX), from global O/F

delta_pressure_injector = P0 * inj_press_drop #Injector pressure drop [Pa]

# === Injector element layout (design description) ===  
#The injector uses a doublet impinging pair at the center of the faceplate.
#Around the center, three concentric rings of holes distribute additional fuel/oxidizer 
#The central two holes form the primary impinging element (RP1 + LOX)
#The outer ring of fuel holes is directed toward the chamber wall to provide film cooling.

#defining how many holes of each type are included in the design.
num_holes_rp1_cooling = int(input("Number of holes designed for cooling: "))
num_holes_rp1_inj = int(input("Number of kerosene holes for the injector: "))
num_holes_lox_inj = int(input("Number of Liquid oxygen holes for the injector: "))

#This function distribute the total propellant mass flow among a specific set of injector holes
    #mdot_propellant --- the total mass flow rate of that propelant [kg/s]
    #number_holes --- number of holes of this specific type (e.g., cooling, main injection)
    #total_number_of_holes --- total number of holes that share the propellant flow
def mdot_prop(mdot_propellant, number_holes, total_number_of_holes):
    mdot_propel = mdot_propellant * number_holes/(total_number_of_holes)
    mdot_per_hole = mdot_propel/number_holes
    return mdot_propel, mdot_per_hole

#This function computes the required injector orifice diameter for a given mass flow 
    #mdot_per_prop --- mass flow rate through this specific parth [kg/s]
    #rho --- propellant density at injector conditions
    #delta_pressure --- pressure drop across the orifice [psi]
    #number_holes --- number of identical holes sharing this mass flow
def diameter(mdot_per_prop, rho, delta_pressure, number_holes): 
    #following the orifice flow model, and solving for the Area
    Total_area = mdot_per_prop/(Cd*np.sqrt(2*rho*delta_pressure*psi_into_pa))
    Area_per_hole = Total_area/number_holes
    d_hole = 2*np.sqrt(Area_per_hole/np.pi) 
    return d_hole

#This function computes the ideal jet exit velocity through an injector orifice using Bernoulli's 
#incompressible flow relation with a pressure drop
    #delta_pressure_injector --- pressure drop across the injector 
    #rho --- propellant density at injector conditions [kg/m^3]
def velocity(delta_pressure_injector, rho):
    velocity = np.sqrt((2*delta_pressure_injector)/rho)
    return velocity

#This function solves for the oxidizer jet angle (theta_ox) required so that the resulting 
#momentum vector of the fuel and oxidizer streams 
#forms a desired impingement angle beta relative to the vertical
    """
    Parameters:
    theta_fuel_deg: injection angle of the fuel jet (degrees from vertical)
    beta_deg: Desired impingement angle of the resultant momentum vector
    mdot_fuel: Mass flow rate of the fuel jet [kg/s]
    v_fuel: Exit velocity of the fuel jet [m/s]
    mdot_ox: mass flow rate of the oxidizer jet [m/s]
    v_ox: exit velocity of the oxidizer jet [m/s]
    theta_ox_min_deg, theta_ox_max_deg: search range fro possible oxidizer angles (degrees)
    dtheta_deg: step size for brute-force scan (degrees)
    
    Method:
    - Momentum of each jet is decomposed into X and Y components. Fuel momentum is taken
      with negative X component (pointing inward). Oxidizer momentum uses positive X
      direction.
    - The resultant vector angle is computed via:
            beta_res = atan2(Mx, My)
    - The angle that minimize |beta_res - beta_target| is selected.

    Notes:
    This method assumes 2-D planar impingement.
    Relies only in conservation of momentum, no geometry included. 
    """
def solve_theta_ox_from_beta(theta_fuel_deg, beta_deg, mdot_fuel, v_fuel, mdot_ox, v_ox, 
                             theta_ox_min_deg=5.0, theta_ox_max_deg=80.0, dtheta_deg=0.05):
    theta_f = theta_fuel_deg * degrees_into_rad #converting input into radians
    beta = beta_deg * degrees_into_rad #converting input into radians
    #Momenutm magnittudes for each jet
    M_fuel = mdot_fuel * v_fuel
    M_ox = mdot_ox * v_ox 
    best_theta_ox = None #best angle found so far
    best_err = 1e9 # initialize with an extremely large error
    theta_ox_deg = theta_ox_min_deg 
    #sweep oxidizer angle from min to max
    while theta_ox_deg <= theta_ox_max_deg:
        theta_ox = theta_ox_deg * degrees_into_rad
        #momentum decomposition (fuel will be consider negative direction)
        Mx = -M_fuel * np.sin(theta_f) + M_ox * np.sin(theta_ox)
        My =  M_fuel * np.cos(theta_f) + M_ox * np.cos(theta_ox)
        beta_res = math.atan2(Mx, My) #resultant impingement angle
        err = abs(beta_res - beta) #calculating error relative with the desired data
        #keep the best match
        if err < best_err:
            best_err = err
            best_theta_ox = theta_ox
        theta_ox_deg += dtheta_deg #increment sweep angle
    return best_theta_ox/degrees_into_rad #returning the data in rads

#This function computes the impingement distance (along the chamber axis) and the required spacing
#between the fuel and the oxidizer orifices for a double injector.
"""
Parameters:
Lc: combustion chamber length [mm]
impingement_fraction: Fraction of the Lc where the jets are intended to collide 
theta_fuel_deg: Fuel jet angle relative to the vertical (degrees)
theta_ox_deg: oxidizer jet angle relative to the vertical (degrees)
"""
def compute_spacing_doublet(Lc, impingement_fraction, theta_fuel_deg, theta_ox_deg):
    #Convert the chamber length from mm into meters and multiply by the chosen fraction
    #This sets where, along the chamber axis, the two jets are expected to collide
    impinge_distance = Lc*1e-3*impingement_fraction
    #converting jet angles into radians
    theta_f = theta_fuel_deg*degrees_into_rad
    theta_ox = theta_ox_deg*degrees_into_rad
    #Compute the required spacing between fuel + oxidizer holes
    #for small-angle geometry, lateral spacing = z * (tan(theta_fuel) + tan(theta_ox))
    #where z is the distance at which the jets meet
    d_fo = impinge_distance * (np.tan(theta_f) + np.tan(theta_ox))
    return impinge_distance, d_fo

mdot_rp1_cooling, mdot_rp1_cooling_per_hole = mdot_prop(mdot_kero,num_holes_rp1_cooling,num_holes_rp1_cooling + num_holes_lox_inj)
mdot_rp1_inj, mdot_rp1_inj_per_hole = mdot_prop(mdot_kero, num_holes_rp1_inj, num_holes_rp1_inj + num_holes_rp1_cooling)
mdot_lox_inj, mdot_lox_inj_per_hole = mdot_prop(mdot_lox, num_holes_lox_inj, num_holes_lox_inj)

#Printing the values
print("=== Mass flow rate per hole ===")
print(f"RP1 colling mass flow rate per hole ({num_holes_rp1_cooling}): {mdot_rp1_cooling_per_hole:.5f} kg/s")
print(f"RP1 injection mass flow rate per hole ({num_holes_rp1_inj}): {mdot_rp1_inj_per_hole:.5f} kg/s")
print(f"LOX injection mass flow rate per hole ({num_holes_lox_inj}): {mdot_lox_inj_per_hole:.5f} kg/s")

# Velocity through the holes
v_lox = velocity(delta_pressure_injector, rho_lox)
v_rp1 = velocity(delta_pressure_injector, rho_rp1)
print("=== Velocity through injector holes ===")
print(f"Velocity through LOX injector holes: {v_lox:.5f} m/s")
print(f"Velocity through RP1 injector holes: {v_rp1:.5f} m/s")

#Diameter of the holes
diameter_cooling_rp1 = diameter(mdot_rp1_cooling, rho_rp1, delta_P, num_holes_rp1_cooling)
diameter_inj_rp1 = diameter(mdot_rp1_inj, rho_rp1, delta_P, num_holes_rp1_inj)
diameter_inj_lox = diameter(mdot_lox_inj, rho_lox, delta_P, num_holes_lox_inj)
print("=== Diameter of each hole ===")
print(f"The diameter of each RP1 cooling hole is {diameter_cooling_rp1:.5f} m")
print(f"The diameter of each RP1 hole into the injector is {diameter_inj_rp1:.5f} m")
print(f"The diameter of each liquid oxygen hole into the injector is {diameter_inj_lox:.5f} m")

#Angle of the oxidizer given a angle of the fuel and a angle of resultant spray after impingement (beta)
print("=== Angle of lox ===")
beta = float(input("Enter an angle for the impingement (relative to vertical): "))
theta_rp1 = float(input("Enter an angle for the rp1 fuel orifice (relative to vertical): "))
theta_lox = solve_theta_ox_from_beta(theta_rp1, beta, mdot_rp1_inj_per_hole, v_rp1, mdot_lox_inj_per_hole, v_lox)
print(f"The angle of the oxidizer, liquid oxygen, respective to the vertical is: {theta_lox:.2f}°")

#Distance of the two orifices that will allow both streams to meet at a single point for a given
#impingement distance (in this case 10% of the chamber length). 
print("=== Distances between orifices and the impingement distance ===")
impinge_d, d_oo = compute_spacing_doublet(Longitude_chamber, impinge_fraction, theta_rp1, theta_lox)
print(f"Distance between RP1 - lox orifices: {d_oo} m")

print(f"Impingement distance: {impinge_d} m")
