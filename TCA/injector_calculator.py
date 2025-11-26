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
psi_into_pa = 6894.76 # Convert psi -> Pascal
meters_into_inches = 39.37 # Convert meters -> inches
degrees_into_rad = np.pi/180 # Convert degrees -> radians

# == Desing Parameters (Change as needed) == 
mdot = 9 #Total propellant mass flow rate [kg/s] (CHANGE AS NEEDED)
OF_Ratio = 2.1 #Mixture ratio O/F 
rho_rp1 = 810 #RP1 Density [kg/m^3] at injector conditions
rho_lox = 1141 #LOX Density [kg/m^3] at injector conditions
inj_press_drop = 5/7 #Fraction of chamber pressure allocated to injector
P0 = 700*psi_into_pa #Chamber stagnation pressure [Pa]
Cd = 0.7 #Discharge coefficient
delta_P = 200 #Injector pressure drop [psi] (converted inside functions)
theta_1 = 40 #Initial guess / design choice for injector angle [deg]
Length_chamber = 400 #Combustion chamber length [mm]
impinge_fraction = 0.10 #streams will impinge at 10% of chamber length
distance_between_holes = 0.05 #[m]
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
def velocity(Cd, delta_pressure_injector, rho):
    velocity = Cd*np.sqrt((2*delta_pressure_injector)/rho)
    return velocity

#This function determines the oxidizer injection angle (theta_ox) for a doublet impinging injector
#such that the post-impingement spray sheet forms a desired angle beta with respect to the chamber 
#axis, for a given fuel injection angle theta_fuel. 

    """
    Parameters:
    mdot_fuel: Mass flow rate of the fuel jet [kg/s]
    mdot_ox: Mass flow rate of the oxidizer jet [kg/s]
    v_fuel: Exit velocity of the fuel jet [m/s]
    v_ox: Exit velocity of the oxidizer jet [m/s]
    theta_fuel: Fuel injection angle measured from the chamber axis [deg]
    beta_target_deg: Desired resultant spray angle after the impingement [deg]
    theta_ox_min_deg, theta_ox_max_deg: Search range for oxidizer injection angle [deg]
    dtheta_deg: angular resolution of the brute-force sweep [deg]
    
    Method:
    - Momentum of each jet is decomposed into X and Y components. Fuel momentum is taken
      with negative X component (pointing inward). Oxidizer momentum uses positive X
      direction.
    - For each candidate oxidizer angle theta_ox in the search range:
        * Convert theta_fuel and theta_ox to radians.
        * Resolve each jet momentum into X/Y components (X lateral, Y axial):
              p_f = ( -P_f * sin(theta_f),  P_f * cos(theta_f) )
              p_ox = ( P_ox * sin(theta_ox), P_ox * cos(theta_ox) )
          (fuel is taken pointing inward).
        * Sum both vectors to obtain the resultant momentum:
              M = p_f + p_ox = (Mx, My)
        * Compute the sheet angle:
              beta_res = atan2(Mx, My)
        * Evaluate the error |beta_res − beta_target| and keep the angle
          that minimizes this error.
    - The function returns the best beta_res (in degrees) and the
      corresponding theta_ox (in degrees).

    Notes:
    This model assumes planar (2-D) impingement and neglects droplet breakup physics    
    No viscous losses or jet spreading are included: pure momentum geometry only
    The penalty factor is unrealistic injector designs with extreme asymmetry. 
    """
def solve_thetas(mdot_fuel, mdot_ox, v_fuel, v_ox, theta_fuel, beta_target_deg, 
                 theta_ox_min_deg=5.0, theta_ox_max_deg=80.0, dtheta_deg=0.05):
    theta_f_rad = theta_fuel*degrees_into_rad
    beta_target = beta_target_deg*degrees_into_rad
    p_f = mdot_fuel * v_fuel
    p_ox = mdot_ox * v_ox
    best_err = 1e12
    best_theta_ox_deg = None
    best_beta_deg = None 
    theta_ox_deg = theta_ox_min_deg
    while theta_ox_deg <= theta_ox_max_deg:
        theta_ox_rad = theta_ox_deg*degrees_into_rad
        px = p_ox*np.sin(theta_ox_rad) - p_f*np.sin(theta_f_rad)
        py = p_ox*np.cos(theta_ox_rad) + p_f*np.cos(theta_f_rad)
        R = px/py
        beta_res = np.arctan(R)
        err = abs(beta_res - beta_target)
        if err < best_err:
            best_err = err
            best_theta_ox_deg = theta_ox_deg
            best_beta_deg = beta_res/degrees_into_rad
        theta_ox_deg += dtheta_deg
    return best_beta_deg, best_theta_ox_deg

#This function computes the geometric relationships of a doublet impinging injection:
"""
Parameters:
Lc: combustion chamber length [mm]
impingement_fraction: Fraction of the Lc where the jets are intended to collide 
theta_fuel_deg: Fuel jet angle relative to the vertical [deg]
theta_ox_deg: oxidizer jet angle relative to the vertical [deg]
Notes:
- Assumes straight, ballistic jets with small lateral deflection.
- Suitable for preliminary injector geometry sizing.
"""
def compute_spacing_doublet(Lc, impingement_fraction, theta_fuel_deg, theta_ox_deg):
    #convert chamber length from mm to m and compute axial impingement location
    z_imp = Lc * 1e-3 *impingement_fraction
    #converting jet angles into radians
    theta_f = theta_fuel_deg*degrees_into_rad
    theta_ox = theta_ox_deg*degrees_into_rad
    #Lateral displacement of each jet at the impingement plane:
    # x = z * tan(theta)
    d_imp_f = z_imp * (np.tan(theta_f))
    d_imp_ox = z_imp * (np.tan(theta_ox))
    #required injector-to-injector spacing for the jets to meet at z_imp
    d_fo = z_imp * (np.tan(theta_f) + np.tan(theta_ox))
    return z_imp, d_imp_f, d_imp_ox, d_fo

#compute total and per-hole mass flow for RP1 and LOX injection holes
#mdot_prop(...) distributes the propellant mass flow across the specified number of holes. 
mdot_rp1_inj, mdot_rp1_inj_per_hole = mdot_prop(mdot_kero, num_holes_rp1_inj, num_holes_rp1_inj)
mdot_lox_inj, mdot_lox_inj_per_hole = mdot_prop(mdot_lox, num_holes_lox_inj, num_holes_lox_inj) 
print("=== Mass flow rate per hole ===")
print(f"RP1 injection mass flow rate per hole ({num_holes_rp1_inj}): {mdot_rp1_inj_per_hole:.5f} kg/s")
print(f"LOX injection mass flow rate per hole ({num_holes_lox_inj}): {mdot_lox_inj_per_hole:.5f} kg/s")

#computing injection velocity for each propellant
v_lox = velocity(Cd, delta_pressure_injector, rho_lox)
v_rp1 = velocity(Cd, delta_pressure_injector, rho_rp1)
print("=== Velocity through injector holes ===")
print(f"Velocity through LOX injector holes: {v_lox:.5f} m/s")
print(f"Velocity through RP1 injector holes: {v_rp1:.5f} m/s")

#compute the required orifice diameter for RP-1 and LOX jets using the orifice flow model.
#each diameter is sized to pass the per-hole mass flow at the specified ΔP and fluid density
diameter_inj_rp1 = diameter(mdot_rp1_inj, rho_rp1, delta_P, num_holes_rp1_inj)
diameter_inj_lox = diameter(mdot_lox_inj, rho_lox, delta_P, num_holes_lox_inj)
print("=== Diameter of each hole ===")
print(f"The diameter of each RP1 hole into the injector is {diameter_inj_rp1:.5f} m")
print(f"The diameter of each liquid oxygen hole into the injector is {diameter_inj_lox:.5f} m")

#computes the LOX injection angle rquired so that the resultant 
#post-impingement sheet angle matches the desired beta.
print("=== Injection Angle Results ===")
#user provides fuel angles and target sheet angle beta.
theta_rp1_deg = float(input("Fuel angle (deg): "))
beta_des = float(input("Enter an angle for the impingement (relative to vertical): "))
beta_res_deg, theta_lox_deg = solve_thetas(mdot_rp1_inj_per_hole, mdot_lox_inj_per_hole, v_rp1, v_lox, theta_rp1_deg, beta_des)
print(f"RP-1 angle θ_f: {theta_rp1_deg:.3f}°")
print(f"LOX angle θ_f: {theta_lox_deg:.3f}°")
print(f"Resultant sheet angle β: {beta_res_deg:.3f}°")

# Impinging geometry check
#calculate impingement geometry using the final RP-1 and LOX angles:
z_imp, d_rp1, d_lox, d_fo_req = compute_spacing_doublet(Length_chamber, impinge_fraction, theta_rp1_deg, theta_lox_deg)
print("=== Impingement geometry ===")
#Axial distance from the faceplate where both jets collide
print(f"Axial impingement distance from the faceplate: {z_imp:.4f} m")
#Required spacing between the RP-1 and LOX holes so the jets meet at z_imp
print(f"Required fuel-LOX spacing d_fo: {d_fo_req*1e3:.4f} mm")
#how far the RP-1 jet travels sideways before hitting the LOX jets.
print(f"RP1 lateral displacement at z_imp: {d_rp1*1e3:.2f} mm")
#how far the LOX jet travels sideways before hitting the RP-1 jets.
print(f"LOX lateral displacement at z_imp: {d_lox*1e3:.2f} mm")

# Momentum Ratio
J = (rho_lox*(v_lox)**2)/(rho_rp1*(v_rp1)**2)
print("=== Momentum Ratio ===")
print(f"The momentum-flux ratio is: {J:.2f}")

# == PLOTS == 
#Plotting theta_fuel vs theta_oxidizer
theta_f_range = np.linspace(0, 80, 500)
beta_list = [0, 5, 10, 15]
plt.figure(figsize = (9,6))
for beta in beta_list:
    theta_ox_curve = []
    for theta_f in theta_f_range:   
        beta_ins, theta_ox = solve_thetas(mdot_rp1_inj_per_hole, mdot_lox_inj_per_hole, v_rp1, v_lox, theta_f, beta)
        theta_ox_curve.append(theta_ox)
    plt.plot(theta_f_range, theta_ox_curve, label=f"β = {beta}°")
plt.xlabel("Fuel angle θ_f [deg]")
plt.ylabel("Required LOX angle θ_ox [deg]")
plt.title("Relationship between θ_f and θ_ox for different desired β")
plt.legend()
plt.grid(True)
plt.tight_layout()
plt.show()

#Plotting theta_fuel vs distance between orifices
theta_f_range = np.linspace(0, 80, 500)
beta_list = [0, 5, 10, 15]
z_ipm = Length_chamber * 1e-3 * impinge_fraction
plt.figure(figsize=(9,6))
for beta in beta_list:
    d_between_orifices_curve = []
    for theta_f in theta_f_range:
        beta_ins, theta_ox_deg = solve_thetas(mdot_rp1_inj_per_hole, mdot_lox_inj_per_hole, v_rp1, v_lox, theta_f, beta)
        theta_f_rad = theta_f*degrees_into_rad
        theta_ox_rad = theta_ox_deg*degrees_into_rad
        d_between_orifices = z_imp*(np.tan(theta_f_rad)+np.tan(theta_ox_rad))
        d_between_orifices_curve.append(d_between_orifices)
    plt.plot(theta_f_range, d_between_orifices_curve, label=f"β = {beta}°")
plt.xlabel("Fuel angle θ_f [deg]")
plt.ylabel("Required injector spacing d_between [m]")
plt.title("Fuel angle vs required spacing between fuel and LOX orifices")
plt.grid(True)
plt.legend()
plt.tight_layout()
plt.show()

