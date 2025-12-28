"""
This program computes the main geometric and flow parameters needed to design
a liquid rocket injector plate. It calculates the injector element layout,
orifice locations and angles, and the fuel and oxidizer manifold dimensions
so that the required flow rates and pressure drops are satisfied.

Author: Joaquin Alarcon
"""

import numpy as np
import math
import matplotlib.pyplot as plt
from spicy import optimize
import csv

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
impinge_fraction = 0.022 #streams will impinge at 2.2% of chamber length
distance_between_holes = 0.05 #[m]
CombDiam = 6.336/meters_into_inches #chamber inner diameter [m]
marginWall = 0.0168 #clearance from outer RP1 jet to wall [m]
pairSpacing = 0.0212 #spacing between mid-radii of FO pairs [m]
thickness = 0.013 #[m]

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
#The injector uses a doublet impinging pair at the center of the faceplate.
#Around the center, three concentric rings of holes distribute additional fuel/oxidizer 
#The central two holes form the primary impinging element (RP1 + LOX)
#The outer ring of fuel holes is directed toward the chamber wall to provide film cooling.

#defining how many holes of each type are included in the design.
num_holes_rp1_inj = int(input("Number of kerosene holes for the injector: "))
num_holes_lox_inj = int(input("Number of Liquid oxygen holes for the injector: "))

#Number of rings per propellant
Nrings = 2 #RP1 - LOX = LOX = RP1   

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

def dist(thickness, theta):
    theta_rad = theta*degrees_into_rad
    dist = thickness*(np.tan(theta_rad))
    return dist

def max_mixing_ratio(rho_fuel, rho_ox, mdot_fuel, mdot_ox):
    M = 1 #M=1 for 1-on-1 unlike impingement
    MME = M*(rho_ox/rho_fuel * (mdot_ox/mdot_fuel)**2)**0.7
    return MME

def design_manifold(d_orifice, n_orifices, h, rho, m_dot):
    #Total exit area (orifices):
    A_exit = n_orifices * (np.pi*(d_orifice**2)/4)
    #required manifold cross-sectional area:
    #based on the 4:1 rule for pressure uniformity
    req_manifold_area = 3 * A_exit
    #geometry with rounded corners (r = h/4)
    r = h/6
    #area loss at the 4 corners compared to a perfect rectangle
    corner_area_loss = (r**2)*(4 - np.pi)
    w = (req_manifold_area + corner_area_loss)/h
    #manifold internal velocity
    v = m_dot/(rho*req_manifold_area)
    return A_exit, req_manifold_area, h, w, r, v

def calculate_manifold_pressure_drop(w, h, r, length, rho, mu, mdot, epsilon = 0.015):
    area = w*h - (4 - np.pi)*(r**2)
    perimeter = 2*(w + h - 4*r) + 2*np.pi*r
    #hydraulic diameter
    Dh = (4*area/perimeter)
    velocity = mdot/(rho*area)
    #reynolds number
    Re = (rho*velocity*Dh)/mu
    #friction factor (using the Haaland equation)
    if Re > 2300:
        rel_roughness = (epsilon/1000)/Dh
        f_inv_sqrt = -1.8 * np.log10(((rel_roughness/3.7)**1.11) + (6.9/Re))
        f = (1/f_inv_sqrt)**2
    else:
        f = 64/Re
    delta_P = f * (length/Dh) * (rho*velocity**2)/2
    delta_P_psi = delta_P / psi_into_pa
    return delta_P_psi

#compute total and per-hole mass flow for RP1 and LOX injection holes
#mdot_prop(...) distributes the propellant mass flow across the specified number of holes. 
mdot_rp1_inj, mdot_rp1_inj_per_hole = mdot_prop(mdot_kero, num_holes_rp1_inj, num_holes_rp1_inj)
mdot_lox_inj, mdot_lox_inj_per_hole = mdot_prop(mdot_lox, num_holes_lox_inj, num_holes_lox_inj) 
print("\n=== Mass flow rate per hole ===")
print(f"RP1 injection mass flow rate per hole ({num_holes_rp1_inj}): {mdot_rp1_inj_per_hole:.5f} kg/s")
print(f"LOX injection mass flow rate per hole ({num_holes_lox_inj}): {mdot_lox_inj_per_hole:.5f} kg/s")

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
print("\n=== Diameter of each hole ===")
print(f"The diameter of each RP1 hole into the injector is {diameter_inj_rp1:.8f} m")
print(f"The diameter of each liquid oxygen hole into the injector is {diameter_inj_lox:.8f} m")

MME = max_mixing_ratio(rho_rp1, rho_lox, mdot_rp1_inj_per_hole, mdot_lox_inj_per_hole)
print("\n=== Max Mixing Ratio ===")
print(f"The maximum mixing ratio is {MME}")

#computes the LOX injection angle rquired so that the resultant 
#post-impingement sheet angle matches the desired beta.
print("\n=== Injection Angle Results ===")
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
print("\n=== Impingement geometry ===")
#Axial distance from the faceplate where both jets collide
print(f"Axial impingement distance from the faceplate: {z_imp:.4f} m")
#Required spacing between the RP-1 and LOX holes so the jets meet at z_imp
print(f"Required fuel-LOX spacing d_fo: {d_fo_req*1e3:.4f} mm")
#how far the RP-1 jet travels sideways before hitting the LOX jets.
print(f"RP1 lateral displacement at z_imp: {d_rp1*1e3:.2f} mm")
#how far the LOX jet travels sideways before hitting the RP-1 jets.
print(f"LOX lateral displacement at z_imp: {d_lox*1e3:.2f} mm")

#distance between inner and outer orifices
dist_fuel = dist(thickness, theta_rp1_deg)
dist_ox = dist(thickness, theta_lox_deg)
print(f"\nThe distance between the inner and the outer orifice for fuel is {dist_fuel:.8f}")
print(f"The distance between the inner and the outer orifice for oxidizer is {dist_ox:.8f}")

# Momentum Ratio
J = (rho_lox*(v_lox)**2)/(rho_rp1*(v_rp1)**2)
print("\n=== Momentum Ratio ===")
print(f"The momentum-flux ratio is: {J:.2f}")

#MANIFOLD OUTPUTS
combRad = CombDiam / 2.0
radius_Ox = diameter_inj_lox / 2.0
radius_fuel = diameter_inj_rp1 / 2.0
#mid-radii of FO pairs (inner, outer):
RmidPairOuter = combRad - marginWall - radius_fuel - d_fo_req/2
RmidPairInner = RmidPairOuter - pairSpacing
Rmid_list = np.array([RmidPairInner, RmidPairOuter])

p_list = 2*np.pi*Rmid_list / (num_holes_lox_inj/Nrings)
circ_elements_list = 2*np.pi*Rmid_list

Rf_inner = RmidPairInner - d_fo_req/2
ROx_inner = RmidPairInner + d_fo_req/2
Rf_outer = RmidPairOuter + d_fo_req/2
ROx_outer = RmidPairOuter - d_fo_req/2

Rring_rp1 = np.array([Rf_inner, Rf_outer])
Rring_lox = np.array([ROx_inner, ROx_outer])

circ_ox = 2*np.pi*Rring_lox
circ_rp1 = 2*np.pi*Rring_rp1

Lpath_ox = circ_ox / (2*NinletsOx)
Lpath_rp1 = circ_rp1 / (2*NinletsRP1)
Lpath_ox_avg = (Lpath_ox[0] + Lpath_ox[1])/2

dpInjOx = delta_P*psi_into_pa
dpInjRP1 = delta_P*psi_into_pa
DguessOx = 0.3 #m
DguessRP1 = 0.3 #m
w_guess = 0.20
mdot_lox_per_ring = mdot_lox_inj / Nrings
mdot_rp1_per_ring = mdot_rp1_inj / Nrings
#CAD Design Parameters
print("\n=== CAD Design Parameters ===")
print("=== Injector ring diameters at EXIT plane (chamber side) ===")
print("Inner Ring Diameters:")
print(f"Inner RP1 Ring (Diameter): {Rf_inner*1e3*2} mm")
print(f"Inner LOX Ring (Diameter): {ROx_inner*1e3*2} mm")
print("Outer Ring Diameters:")
print(f"Outer LOX Ring (Diameter): {ROx_outer*1e3*2} mm")
print(f"Outer RP1 Ring (Diameter): {Rf_outer*1e3*2} mm\n")
print("=== Injector ring diameters at BACK plane (upstream by thickness) ===")
print(f"RP-1  inner ring D_f_in_back = {(-1)*dist_fuel*1e3:.2f} mm")
print(f"RP-1 outer ring D_f_out_back = {dist_fuel*1e3:.2f} mm")
print(f"LOX  inner ring D_ox_in_back = {(-1)*dist_ox*1e3:.2f} mm")
print(f"LOX outer ring D_ox_out_back = {dist_ox*1e3:.2f} mm")
print("\n=== Manifold Design Results ===")
A_exit_ox, A_mani_ox, h_ox, w_ox, r_ox, v_mani_ox = design_manifold(diameter_inj_lox, num_holes_lox_inj, h_manifold, rho_lox, mdot_lox_per_ring*2)
A_exit_f_in, A_mani_f_in, h_f_in, w_f_in, r_f_in, v_mani_f_in = design_manifold(diameter_inj_rp1, num_holes_rp1_inj/2, h_manifold, rho_rp1, mdot_rp1_per_ring)
A_exit_f_out, A_mani_f_out, h_f_out, w_f_out, r_f_out, v_mani_f_out = design_manifold(diameter_inj_rp1, num_holes_rp1_inj/2, h_manifold, rho_rp1, mdot_rp1_per_ring)
print("=== LOX Manifold ===")
print(f"LOX manifold required cross-sectional area: {A_mani_ox*1e6:.2f} mm^2")
print(f"LOX manifold width w: {w_ox*1e3:.2f} mm")
print(f"LOX manifold height h: {h_ox*1e3:.2f} mm")
print(f"LOX manifold corner radius r: {r_ox*1e3:.2f} mm")
print(f"LOX manifold internal velocity: {v_mani_ox:.2f} m/s\n")
print("=== RP-1 Inner Manifold ===")
print(f"RP-1 inner manifold required cross-sectional area: {A_mani_f_in*1e6:.2f} mm^2")
print(f"RP-1 inner manifold width w: {w_f_in*1e3:.2f} mm")
print(f"RP-1 inner manifold height h: {h_f_in*1e3:.2f} mm")
print(f"RP-1 inner manifold corner radius r: {r_f_in*1e3:.2f} mm")
print(f"RP-1 inner manifold internal velocity: {v_mani_f_in:.2f} m/s\n")
print("=== RP-1 Outer Manifold ===")
print(f"RP-1 outer manifold required cross-sectional area: {A_mani_f_out*1e6:.2f} mm^2")
print(f"RP-1 outer manifold width w: {w_f_out*1e3:.2f} mm")
print(f"RP-1 outer manifold height h: {h_f_out*1e3:.2f} mm")
print(f"RP-1 outer manifold corner radius r: {r_f_out*1e3:.2f} mm")
print(f"RP-1 outer manifold internal velocity: {v_mani_f_out:.2f} m/s\n")

#compute manifold pressure drops
dp_mani_ox = calculate_manifold_pressure_drop(w_ox, h_ox, r_ox, Lpath_ox_avg, rho_lox, 0.00019, mdot_lox_per_ring*2)
dp_mani_rp1_in = calculate_manifold_pressure_drop(w_f_in, h_f_in, r_f_in, Lpath_rp1[0], rho_rp1, 0.0002, mdot_rp1_per_ring)
dp_mani_rp1_out = calculate_manifold_pressure_drop(w_f_out, h_f_out, r_f_out, Lpath_rp1[1], rho_rp1, 0.0002, mdot_rp1_per_ring)
print("=== Manifold Pressure Drops ===")
print(f"LOX manifold pressure drop: {dp_mani_ox:.2f} psi")
print(f"RP-1 inner manifold pressure drop: {dp_mani_rp1_in:.2f} psi")
print(f"RP-1 outer manifold pressure drop: {dp_mani_rp1_out:.2f} psi")

# == EXCEL ==

def export_fusion_parameters(filename, parameters):
    with open(filename, "w", newline = "", encoding = "utf-8") as fo:
        writer = csv.writer(fo)
        writer.writerow(["Name", "Unit", "Expression", "Value", "Comments", "Favorite"])
        for name, unit, value in parameters:
            if unit == "ul":
                expression = f"{value}"
            else:
                expression = f"{value} {unit}"
            writer.writerow([name, unit, expression, value, "", "False"])
            
w_avg_lox = (ROx_outer*1e3*2 - dist_ox*1e3 + ROx_inner*1e3*2 + dist_ox*1e3)/2

params = [
    ("rp1_in", "mm", Rf_inner*2*1e3), ("rp1_out", "mm", Rf_outer*2*1e3), ("ox_in", "mm", ROx_inner*2*1e3), ("ox_out", "mm", ROx_outer*2*1e3),
    ("rp1_hole_diameter", "mm", diameter_inj_rp1*1e3), ("rp1_hole_diameter2", "mm", diameter_inj_rp1*1e3), ("ox_hole_diameter", "mm", diameter_inj_lox*1e3),
    ("ox_hole_diameter2", "mm", diameter_inj_lox*1e3), ("negative_rp1_d_angle", "mm", -dist_fuel*1e3), ("positive_rp1_d_angle", "mm", dist_fuel*1e3),
    ("positive_ox_d_angle", "mm", dist_ox*1e3), ("negative_ox_d_angle", "mm", -dist_ox*1e3), ("rp1_hole_diameter3", "mm", diameter_inj_rp1*1e3), ("rp1_hole_diameter4", "mm", diameter_inj_rp1*1e3),
    ("ox_hole_diameter3", "mm", diameter_inj_lox*1e3), ("ox_hole_diameter4", "mm", diameter_inj_lox*1e3), ("holes_per_ring", "", int(num_holes_lox_inj/Nrings)),
    ("positive_rp1_in_w_half", "mm", w_f_in*1e3/2), ("negative_rp1_in_w_half", "mm", -w_f_in*1e3/2), ("positive_rp1_out_w_half", "mm", w_f_out*1e3/2),
    ("negative_rp1_out_w_half", "mm", -w_f_out*1e3/2), ("manifold_height", "mm", h_manifold*1e3), ("inj_plate_thickness", "mm", thickness*1e3),
    ("w_avg_lox", "mm", w_avg_lox), ("positive_ox_w_half", "mm", w_ox*1e3/2)]

export_fusion_parameters("injector_parameters.csv", params)

# == PLOTS == 

def ring_points(radius_mm, N):
    ang = np.linspace(0, 2*np.pi, N, endpoint=False)
    x = radius_mm * np.cos(ang)
    y = radius_mm * np.sin(ang)
    return x, y

def plot_injector_layout(D_c_mm, R_ring_f_mm, R_ring_ox_mm, N_hole_ring, margin_extra_mm=5):
    combRad = D_c_mm / 2.0 #chamber radius
    x_f_in,  y_f_in  = ring_points(R_ring_f_mm[0], N_hole_ring)
    x_f_out, y_f_out = ring_points(R_ring_f_mm[1], N_hole_ring)
    x_ox_in,  y_ox_in  = ring_points(R_ring_ox_mm[0], N_hole_ring)
    x_ox_out, y_ox_out = ring_points(R_ring_ox_mm[1], N_hole_ring)
    fig, ax = plt.subplots(figsize=(10, 10))
    chamber = plt.Circle((0, 0), combRad, color='gray',
                         fill=False, linewidth=2)
    ax.add_artist(chamber)
    for R in R_ring_f_mm:
        ax.add_artist(plt.Circle((0, 0), R, color='red',
                                 fill=False, linewidth=0.7, alpha=0.4))
    for R in R_ring_ox_mm:
        ax.add_artist(plt.Circle((0, 0), R, color='blue',
                                 fill=False, linewidth=0.7, alpha=0.4))
    ax.scatter(x_f_in,  y_f_in, s=18, color='red', label='RP-1 inner')
    ax.scatter(x_f_out, y_f_out, s=18, color='darkred', label='RP-1 outer')
    ax.scatter(x_ox_in,  y_ox_in, s=18, color='royalblue', label='LOX inner')
    ax.scatter(x_ox_out, y_ox_out, s=18, color='navy', label='LOX outer')
    ax.set_aspect('equal', 'box')
    ax.set_xlabel('x [mm]')
    ax.set_ylabel('y [mm]')
    ax.set_title('4-ring impinging injector layout')
    ax.legend(loc='upper right')
    ax.grid(True, linestyle=':', linewidth=0.5)
    margin = combRad + margin_extra_mm
    ax.set_xlim(-margin, margin)
    ax.set_ylim(-margin, margin)
    plt.tight_layout()
    plt.show()

holes_per_ring = int(num_holes_lox_inj/Nrings)
plot_injector_layout(2*combRad * 1e3, Rring_rp1 * 1e3, Rring_lox * 1e3, holes_per_ring)

#Plotting manifold geometry. h vs w
w_init_ref = 0.2
h_vals = np.linspace(0.005, 0.025, 50)
def manifold_plot(diameter, num_holes, rho, mdot):
    w_list = []
    for h in h_vals:
        Area_exit, A_mani, h_mani, w_mani, r_mani, v_mani = design_manifold(diameter, num_holes, h, rho, mdot)
        w_list.append(w_mani)
    return np.array(w_list)
w_RP1_inner = manifold_plot(diameter_inj_rp1, num_holes_rp1_inj/2, rho_rp1, mdot_rp1_per_ring)
w_RP1_outer = manifold_plot(diameter_inj_rp1, num_holes_rp1_inj/2, rho_rp1, mdot_rp1_per_ring)
w_LOX_main = manifold_plot(diameter_inj_lox, num_holes_lox_inj, rho_lox, mdot_lox_inj)
fig, axs = plt.subplots(2,2,figsize=(10,8))
axs[0,0].plot(h_vals*1e3, w_RP1_inner*1e3, marker='o')
axs[0,0].set_title("RP-1 Inner Manifold")
axs[0,0].set_xlabel("h [mm]")
axs[0,0].set_ylabel("w [mm]")
axs[0,0].grid(True)
axs[0,1].plot(h_vals*1e3, w_RP1_outer*1e3, marker='o', color='coral')
axs[0,1].set_title("RP-1 Outer Manifold")
axs[0,1].set_xlabel("h [mm]")
axs[0,1].set_ylabel("w [mm]")
axs[0,1].grid(True)
axs[1,0].plot(h_vals*1e3, w_LOX_main*1e3, marker='o', color='navy')
axs[1,0].set_title("LOX Main Manifold")
axs[1,0].set_xlabel("h [mm]")
axs[1,0].set_ylabel("w [mm]")
axs[1,0].grid(True)
plt.tight_layout()
plt.show()

#Plotting theta_fuel vs theta_oxidizer
theta_f_range = np.linspace(0, 80, 500)
beta_list = [beta_res_deg - 5, beta_res_deg, beta_res_deg + 5]
plt.figure(figsize = (9,6))
for beta in beta_list:
    theta_ox_curve = []
    for theta_f in theta_f_range:   
        beta_ins, theta_ox = solve_thetas(mdot_rp1_inj_per_hole, mdot_lox_inj_per_hole, v_rp1, v_lox, theta_f, beta)
        theta_ox_curve.append(theta_ox)
    plt.plot(theta_f_range, theta_ox_curve, label=f"β = {beta:.3f}°")
plt.xlabel("Fuel angle θ_f [deg]")
plt.ylabel("Required LOX angle θ_ox [deg]")
plt.title("Relationship between θ_f and θ_ox for different desired β")
plt.legend()
plt.grid(True)
plt.tight_layout()
plt.show()

#Plotting theta_fuel vs distance between orifices
theta_f_range = np.linspace(0, 80, 500)
beta_list = [beta_res_deg - 5, beta_res_deg, beta_res_deg + 5]
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
    plt.plot(theta_f_range, d_between_orifices_curve, label=f"β = {beta:.3f}°")
plt.xlabel("Fuel angle θ_f [deg]")
plt.ylabel("Required injector spacing d_between [m]")
plt.title("Fuel angle vs required spacing between fuel and LOX orifices")
plt.grid(True)
plt.legend()
plt.tight_layout()
plt.show()
