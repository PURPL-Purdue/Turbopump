import numpy as np
import matplotlib.pyplot as plt
from spicy import optimize
import csv
from CoolProp.CoolProp import PropsSI

# == Conversion between units ==
psi_into_pa = 6894.76 # Convert psi -> Pascal
meters_into_inches = 39.37 # Convert meters -> inches
degrees_into_rad = np.pi/180 # Convert degrees -> radians

# == Desing Parameters (Change as needed) == 
mdot = 9 #Total propellant mass flow rate [kg/s] (CHANGE AS NEEDED)
OF_Ratio = 2.1 #Mixture ratio O/F Core
rho_rp1 = 810 #RP1 Density [kg/m^3] at injector conditions
rho_lox = 1141 #LOX Density [kg/m^3] at injector conditions
ratio_film_cooling = 0.05 #Fraction of fuel mass flow sent to film cooling vs. injection
Cd = 0.7 #Discharge coefficient

Pc = 500*psi_into_pa #Chamber stagnation pressure [Pa]
Pin = 850*psi_into_pa #injector inlet pressure [Pa]
delta_P = 350 #Injector pressure drop [psi] (converted inside functions)
inj_press_drop = Pc/Pin #Fraction of chamber pressure allocated to injector
stifness = delta_P*psi_into_pa/Pc

sigma_rp1 = PropsSI('SURFACE_TENSION','T',298,'Q',0,'Dodecane') #Surface tension of RP-1 at injector conditions [N/m]
sigma_lox = PropsSI('SURFACE_TENSION','T',90,'Q',0,'Oxygen') #Surface tension of LOX at injector conditions [N/m]

Length_chamber = 9.928/meters_into_inches * 1e3 #Combustion chamber length [mm]
impinge_fraction = 0.022 #streams will impinge at 2.2% of chamber length
distance_between_holes = 0.05 #[m]
CombDiam = 6.336/meters_into_inches #chamber inner diameter [m]
marginWall = 0.0125 #clearance from outer RP1 jet to wall [m]
pairSpacing = 0.0287 #spacing between mid-radii of FO pairs [m]
thickness = 0.014 #[m]

#MANIFOLD DESING
h_manifold = 0.0155 #m
dpFracMani = 0.08 #dp manifold ~9% of injector dp
fOx = 0.01 #darcy friction factor LOX
fRP1 = 0.01 #darcy friction factor RP1
KOx = 1 #lumped local loss K LOX
KRP1 = 1 #lumped local loss K RP1
NinletsOx = 1 #number of LOX manifold inlets
NinletsRP1 = 1 #number of RP1 manifold inlets

#ratio_inj_cooling = 1 #Fraction of fuel mass flow sent to cooling vs. injection 

mdot_kero = mdot/(1+OF_Ratio) #Fuel mass flow (RP-1), from global O/F
mdot_lox = mdot*OF_Ratio/(1+OF_Ratio) #Oxidizer mass flow (LOX), from global O/F

delta_pressure_injector = Pin * inj_press_drop #Injector pressure drop [Pa]

# === Injector element layout (design description) === 

#defining how many holes of each type are included in the design
num_holes_rp1_inj = int(input("Number of kerosene holes for the injector: "))
num_holes_lox_inj = int(input("Number of Liquid oxygen holes for the injector: "))
N_film_Holes = 80 #Number of film cooling holes

#Number of quadlets
num_quadlets = num_holes_lox_inj // 2
print(f"Number of quadlets: {num_quadlets}")

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
    velocity = Cd*np.sqrt((2*delta_pressure_injector)/rho)
    return velocity

#This functions determines the oxidizer (theta_ox) injection angle for a quadlet impinging injector
#such that the post-impingement spray sheet forms a desired angle beta with respect to the chamber 
#axis (in this case beta will be 0 degrees), for a given fuel injection angle theta_fuel

def find_optimal_theta_ox(theta_fuel_deg, fuel_mass_flow, ox_mass_flow, fuel_velocity, ox_velocity):
    theta_f_rad = theta_fuel_deg*degrees_into_rad
    #transverse component factor for fuel (toward center)
    sin_theta_f = np.sin(theta_f_rad)
    required_sin_theta_ox = (fuel_mass_flow * fuel_velocity * sin_theta_f) / (ox_mass_flow * ox_velocity)
    theta_ox_rad = np.arcsin(required_sin_theta_ox)
    theta_ox_deg = theta_ox_rad/degrees_into_rad
    return theta_ox_deg

def spacing_from_center(Length_chamber, impinge_fraction, theta):
    dist_impingement = Length_chamber*impinge_fraction
    r_from_center = dist_impingement*np.tan(theta*degrees_into_rad)
    return r_from_center

def new_global_OF_ratio(RP1_mass_flow_core):
    mdot_fuel_total = RP1_mass_flow_core/(1 - ratio_film_cooling)
    mdot_fuel_film = mdot_fuel_total - RP1_mass_flow_core
    OF_global_new = mdot_lox/(mdot_fuel_total)
    mdot_total_new = mdot_lox + mdot_fuel_total
    return mdot_fuel_film, mdot_fuel_total, OF_global_new, mdot_total_new


#compute total and per-hole mass flow for RP1 and LOX injection holes
#mdot_prop(...) distributes the propellant mass flow across the specified number of holes. 
mdot_rp1_inj, mdot_rp1_inj_per_hole = mdot_prop(mdot_kero, num_holes_rp1_inj, num_holes_rp1_inj)
mdot_lox_inj, mdot_lox_inj_per_hole = mdot_prop(mdot_lox, num_holes_lox_inj, num_holes_lox_inj) 
mdot_fuel_film, mdot_fuel_total, OF_global_new, mdot_total_new = new_global_OF_ratio(mdot_rp1_inj_per_hole*num_holes_rp1_inj)
print("\n=== Mass flow rate per hole ===")
print(f"RP1 injection mass flow rate per hole ({num_holes_rp1_inj}): {mdot_rp1_inj_per_hole:.5f} kg/s")
print(f"LOX injection mass flow rate per hole ({num_holes_lox_inj}): {mdot_lox_inj_per_hole:.5f} kg/s")
print(f"RP1 film cooling mass flow rate ({N_film_Holes} holes): {mdot_fuel_film:.5f} kg/s")
print(f"New global O/F ratio considering film cooling: {OF_global_new:.3f}")
print(f"New total mass flow considering film cooling: {mdot_total_new:.3f} kg/s")

#computing injection velocity for each propellant
v_lox = velocity(delta_P*psi_into_pa, rho_lox)
v_rp1 = velocity(delta_P*psi_into_pa, rho_rp1)
print("\n=== Velocity through injector holes ===")
print(f"Velocity through LOX injector holes: {v_lox:.5f} m/s")
print(f"Velocity through RP1 injector holes: {v_rp1:.5f} m/s")

#compute the required orifice diameter for RP-1 and LOX jets using the orifice flow model.
#each diameter is sized to pass the per-hole mass flow at the specified ΔP and fluid density
diameter_inj_rp1 = diameter(mdot_rp1_inj, rho_rp1, delta_P, num_holes_rp1_inj)
diameter_inj_lox = diameter(mdot_lox_inj, rho_lox, delta_P, num_holes_lox_inj)
diameter_film_cooling = diameter(mdot_fuel_film, rho_rp1, delta_P, N_film_Holes)
print("\n=== Diameter of each hole ===")
print(f"The diameter of each RP1 hole into the injector is {diameter_inj_rp1*1e3:.8f} mm")
print(f"The diameter of each liquid oxygen hole into the injector is {diameter_inj_lox*1e3:.8f} mm")
print(f"The diameter of each film cooling hole into the injector is {diameter_film_cooling*1e3:.8f} mm")

#computes the LOX injection angle rquired so that the resultant 
#post-impingement sheet angle matches the desired beta.
print("\n=== Injection Angle Results ===")
#user provides fuel angles and target sheet angle beta.
theta_rp1_deg = float(input("Fuel angle (deg): "))
theta_lox_deg = find_optimal_theta_ox(theta_rp1_deg, mdot_rp1_inj, mdot_lox_inj, v_rp1, v_lox)
print(f"Optimal LOX injection angle (deg): {theta_lox_deg:.3f}")

# Impinging geometry check
#calculate impingement geometry using the final RP-1 and LOX angles:
r_fuel = spacing_from_center(Length_chamber, impinge_fraction, theta_rp1_deg)
r_ox = spacing_from_center(Length_chamber, impinge_fraction, theta_lox_deg)
print("\n=== Impingement Geometry ===")
print(f"Radius from centerline for RP1 jets at impingement: {r_fuel:.5f} mm")
print(f"Radius from centerline for LOX jets at impingement: {r_ox:.5f} mm")

#Weber Number check
We_lox = (rho_lox * (v_lox*np.sin(theta_lox_deg*degrees_into_rad))**2 * diameter_inj_lox) / sigma_lox
We_rp1 = (rho_rp1 * (v_rp1*np.sin(theta_rp1_deg*degrees_into_rad))**2 * diameter_inj_rp1) / sigma_rp1
print("\n=== Weber Number Check ===")
print(f"Weber Number for LOX jets: {We_lox:.2f}")
print(f"Weber Number for RP1 jets: {We_rp1:.2f}")


