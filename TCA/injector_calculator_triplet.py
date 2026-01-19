import numpy as np
import math
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
OF_Ratio = 2.1 #Mixture ratio O/F 
rho_rp1 = 810 #RP1 Density [kg/m^3] at injector conditions
rho_lox = 1141 #LOX Density [kg/m^3] at injector conditions
Cd = 0.8 #Discharge coefficient

Pc = 500*psi_into_pa #Chamber stagnation pressure [Pa]
Pin = 850*psi_into_pa #injector inlet pressure [Pa]
delta_P = 350 #Injector pressure drop [psi] (converted inside functions)
inj_press_drop = Pc/Pin #Fraction of chamber pressure allocated to injector
stifness = delta_P*psi_into_pa/Pc

Length_chamber = 9.928/meters_into_inches * 1e3 #Combustion chamber length [mm]
impinge_fraction = 0.032 #streams will impinge at 3.2% of chamber length
distance_between_holes = 0.05 #[m]
CombDiam = 6.336/meters_into_inches #chamber inner diameter [m]
marginWall = 0.0125 #clearance from outer RP1 jet to wall [m]
pairSpacing = 0.0287 #spacing between mid-radii of FO pairs [m]
thickness_lox = 0.014 #[m]
thickness_rp1 = 0.014 #[m]

sigma_rp1 = PropsSI('SURFACE_TENSION','T',298,'Q',0,'Dodecane') #Surface tension of RP-1 at injector conditions [N/m]
sigma_lox = PropsSI('SURFACE_TENSION','T',90,'Q',0,'Oxygen') #Surface tension of LOX at injector conditions [N/m]
mu_lox = PropsSI('VISCOSITY','T',90,'Q',0,'Oxygen')  # LOX viscosity [Pa.s] at injector conditions
mu_rp1 = PropsSI('VISCOSITY','T',298,'Q',0,'Dodecane')  # RP1 viscosity [Pa.s] at injector conditions

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

delta_pressure_injector = delta_P*psi_into_pa #Injector pressure drop [Pa]

#defining how many holes of each type are included in the design.
num_elements = int(input("Number of injector elements (triplets): "))
num_holes_rp1_inj = num_elements
num_holes_lox_inj = num_elements*2

#defining the angle of the rp1 injector holes
#the angle is with respect to the vertical line of the impingement jet.
angle_lox_inj = float(input("Angle of kerosene holes for the injector (degrees): ")) 

#This function computes the required injector orifice diameter for a given mass flow 
    #mdot_per_prop --- mass flow rate through this specific parth [kg/s]
    #rho --- propellant density at injector conditions
    #delta_pressure --- pressure drop across the orifice [Pa]
    #number_holes --- number of identical holes sharing this mass flow
def diameter(mdot_per_prop, rho, delta_pressure, number_holes): 
    #following the orifice flow model, and solving for the Area
    Total_area = mdot_per_prop/(Cd*np.sqrt(2*rho*delta_pressure))
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

def distance_from_centerline(angle_deg):
    angle_rad = angle_deg * degrees_into_rad
    z_impinge = Length_chamber * impinge_fraction
    o_f_distance = z_impinge * math.tan(angle_rad)
    return o_f_distance, z_impinge

def counterbore(diameter, theta):
    L_d_ratio = 3.5
    thickness = L_d_ratio * (diameter) * math.tan(theta * degrees_into_rad)
    return thickness

print(f"Delta Pressure across injector: {delta_pressure_injector/psi_into_pa:.2f} psi")

diameter_rp1 = diameter(mdot_kero, rho_rp1, delta_pressure_injector, num_holes_rp1_inj)
diameter_lox = diameter(mdot_lox, rho_lox, delta_pressure_injector, num_holes_lox_inj)
print(f"RP1 injector orifice diameter: {diameter_rp1:.7f} m")
print(f"LOX injector orifice diameter: {diameter_lox:.7f} m")

velocity_rp1 = velocity(delta_pressure_injector, rho_rp1)
velocity_lox = velocity(delta_pressure_injector, rho_lox)
print(f"RP1 injector exit velocity: {velocity_rp1:.2f} m/s")
print(f"LOX injector exit velocity: {velocity_lox:.2f} m/s")

o_f_distance, z_impinge = distance_from_centerline(angle_lox_inj)
print(f"The impingement point is located at: {z_impinge:.4f} mm from the injector face")
print(f"The distance from centerline for RP1 jets is: {o_f_distance:.4f} mm")

thickness_lox = counterbore(diameter_lox, angle_lox_inj)
print(f"Recommended thickness for a stable L/D ratio: {thickness_lox*1e3:.4f} mm")

rho_mix = (OF_Ratio+1)/((OF_Ratio/rho_lox)+(1/rho_rp1))
We_c = rho_mix*(velocity_lox*np.sin(angle_lox_inj*degrees_into_rad))**2*diameter_rp1/sigma_rp1
print(f"Calculated Weber number at chamber centerline: {We_c:.2f}")

def darcy_friction_factor_colebrook(density, velocity, diameter, viscosity, ks):
    Re = (density * velocity * diameter) / viscosity
    if Re < 2300:
        return 64/Re #This is since it is laminar
    f = 0.25/(np.log10(ks/(3.7*diameter) + 5.74/(Re**0.9)))**2 # A good guess would be the Swamee-Jain equation
    for _ in range(50):
        f_old = f 
        f_inverse_root = -2.0*np.log10(ks/(3.7*diameter) + (2.51/(Re*np.sqrt(f)))) # Colebrook-White equation
        f = 1.0/(f_inverse_root**2) # Update friction factor
        if abs(f - f_old) < 1e-15: # Convergence check
            break
    return f

def pressure_drop(density, velocity, diameter, length, viscosity, ks):
    f = darcy_friction_factor_colebrook(density, velocity, diameter, viscosity, ks) # Darcy friction factor
    delta_p = f * (length/diameter) * (0.5 * density * velocity**2) # Pressure drop calculation
    return f, delta_p

# == LOX Pressure Drop Calculation ==
length_pipe_lox = thickness_lox/np.cos(angle_lox_inj*degrees_into_rad)  # LOX pipe length [m] (CHANGE AS NEEDED)
print(f"LOX pipe length for pressure drop calc: {length_pipe_lox:.4f} m")
f_lox, Delta_P_lox = pressure_drop(rho_lox, velocity_lox, diameter_lox, length_pipe_lox, mu_lox, ks=0.045e-3)  # Pressure drop in Pascals
print(f"LOX Pressure Drop: {Delta_P_lox/psi_into_pa:.2f} psi")
print(f"LOX Darcy Friction Factor: {f_lox:.5f}")