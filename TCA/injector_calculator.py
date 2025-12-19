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

# == Conversion between units ==
psi_into_pa = 6894.76 # Convert psi -> Pascal
meters_into_inches = 39.37 # Convert meters -> inches
degrees_into_rad = np.pi/180 # Convert degrees -> radians

# == Desing Parameters (Change as needed) == 
mdot = 9 #Total propellant mass flow rate [kg/s] (CHANGE AS NEEDED)
OF_Ratio = 2.1 #Mixture ratio O/F 
rho_rp1 = 810 #RP1 Density [kg/m^3] at injector conditions
rho_lox = 1141 #LOX Density [kg/m^3] at injector conditions
Cd = 0.7 #Discharge coefficient

inj_press_drop = 2/7 #Fraction of chamber pressure allocated to injector
Pc = 500*psi_into_pa #Chamber stagnation pressure [Pa]
Pin = 850*psi_into_pa #injector inlet pressure [Pa]
stitfness = Pc/Pin
delta_P = 350 #Injector pressure drop [psi] (converted inside functions)
delta_pressure_injector = delta_P*psi_into_pa #Injector pressure drop [Pa]

Length_chamber = 9.928/meters_into_inches #Combustion chamber length [mm]
impinge_fraction = Length_chamber*0.1 #10% of chamber length
distance_between_holes = 0.05 #[m]
CombDiam = 6.336/meters_into_inches #chamber inner diameter [m]
marginWall = 0.01 #clearance from outer RP1 jet to wall [m]
pairSpacing = 0.019 #spacing between mid-radii of FO pairs [m]

#MANIFOLD DESING
h_manifold = 0.0115 #m
dpFracMani = 0.08 #dp manifold ~10% of injector dp
fOx = 0.02 #darcy friction factor LOX
fRP1 = 0.02 #darcy friction factor RP1
KOx = 1 #lumped local loss K LOX
KRP1 = 1 #lumped local loss K RP1
NinletsOx = 1 #number of LOX manifold inlets
NinletsRP1 = 1 #number of RP1 manifold inlets

#ratio_inj_cooling = 1 #Fraction of fuel mass flow sent to cooling vs. injection 

mdot_kero = mdot/(1+OF_Ratio) #Fuel mass flow (RP-1), from global O/F
mdot_lox = mdot*OF_Ratio/(1+OF_Ratio) #Oxidizer mass flow (LOX), from global O/F



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

def max_mixing_ratio(rho_fuel, rho_ox, mdot_fuel, mdot_ox):
    M = 1 #M=1 for 1-on-1 unlike impingement
    MME = M*(rho_ox/rho_fuel * (mdot_ox/mdot_fuel)**2)**0.7
    return MME

def design_manifold(mdot, rho, Lpath, dp_inj, dp_frac, f, Ktot, h, w_init):
    dp_allow = dp_frac*dp_inj
    def func(w):
        A = w*h
        v = mdot / (rho * A)
        Dh = 2*w*h/(w+h)
        dpF = f * (Lpath / Dh) * (rho * v**2 / 2.0)
        dpLoc = Ktot * (rho * v**2 / 2.0)
        dpTot = dpF + dpLoc
        return dpTot - dp_allow
    w_min = w_init / 10.0
    w_max = 10.0 * w_init
    f_min = func(w_min)
    f_max = func(w_max)
    if f_min * f_max > 0:
        factor = 10.0
        for _ in range(5):
            w_min /= factor
            w_max *= factor
            f_min = func(w_min)
            f_max = func(w_max)
            if f_min * f_max <= 0:
                break
    sol = optimize.root_scalar(func, bracket=[w_min, w_max], method="brentq")
    w = sol.root
    A = w*h
    P = 2*(w+h)
    Dh = 4*A/P
    v = mdot/(rho*A)
    dpF = f*(Lpath / Dh)*(rho*v**2/2.0)
    dpLoc = Ktot*(rho*v**2/2.0)
    dpTot = dpF+dpLoc
    ratio = dpTot/dp_inj
    return w, Dh, A, v, dpF, dpLoc, dpTot, ratio

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

MME = max_mixing_ratio(rho_rp1, rho_lox, mdot_rp1_inj_per_hole, mdot_lox_inj_per_hole)
print("=== Max Mixing Ratio ===")
print(f"The maximum mixing ratio is {MME}")

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
mani_lox_inner = design_manifold(mdot_lox_inj, rho_lox, Lpath_ox[0], dpInjOx, dpFracMani, fOx, KOx, h_manifold, w_guess)
mani_lox_outer = design_manifold(mdot_lox_inj, rho_lox, Lpath_ox[1], dpInjOx, dpFracMani, fOx, KOx, h_manifold, w_guess)
mani_lox = design_manifold(mdot_lox, rho_lox, Lpath_ox_avg, dpInjOx, dpFracMani, fOx, KOx, h_manifold, w_guess)
mani_rp1_inner = design_manifold(mdot_rp1_per_ring, rho_rp1, Lpath_rp1[0], dpInjRP1, dpFracMani, fRP1, KRP1, h_manifold, w_guess)
mani_rp1_outer = design_manifold(mdot_rp1_per_ring, rho_rp1, Lpath_rp1[1], dpInjRP1, dpFracMani, fRP1, KRP1, h_manifold, w_guess)
def print_mani(label, mani, L_path, dp_inj, h_man):
    w, Dh, A, v, dpF, dpLoc, dpTot, ratio = mani
    print(f"{label}:")
    print(f"The width is: {w:.7f} m")
    print(f"The height of the manifold is {h_man:.5f} m")
    print(f"The hydraulic diameter is {Dh:.2f} m")
    print(f"L_path (manifold length): {L_path*1e3:.2f} mm")
    print(f"Flow area A: {A*1e6:.2f} mm^2")
    print(f"Mean velocity v: {v:.3f} m/s")
    print(f"Δp_friction (dpF): {dpF:.1f} Pa")
    print(f"Δp_local (dpLoc): {dpLoc:.1f} Pa")
    print(f"Δp_total manifold (dpTot): {dpTot:.1f} Pa")
    print(f"dpTot / dp_inj (target ratio): {dpTot/dp_inj:.3f}")
    print()
print("=== Manifold desing results ===")
print_mani("RP-1 inner ring", mani_rp1_inner, Lpath_rp1[0], dpInjRP1, h_manifold)
print_mani("RP-1 outer ring", mani_rp1_outer, Lpath_rp1[1], dpInjRP1, h_manifold)
print_mani("LOX Main Manifold", mani_lox, Lpath_ox_avg, dpInjOx, h_manifold)

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
    fig, ax = plt.subplots(figsize=(20, 20))
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
def manifold_plot(mdot, rho, Lpath, dp_inj, dp_frac, f, Ktot, w_init_ref):
    w_list = []
    for h in h_vals:
        w, Dh, A, v, dpF, dpLoc, dpTot, ratio = design_manifold(mdot, rho, Lpath, dp_inj, dp_frac, f, Ktot, h, w_init_ref)
        w_list.append(w)
    return np.array(w_list)
w_RP1_inner = manifold_plot(mdot_rp1_per_ring, rho_rp1, Lpath_rp1[0], dpInjRP1, dpFracMani, fRP1, KRP1, w_init_ref)
w_RP1_outer = manifold_plot(mdot_rp1_per_ring, rho_rp1, Lpath_rp1[1], dpInjRP1, dpFracMani, fRP1, KRP1, w_init_ref)
w_LOX_inner = manifold_plot(mdot_lox_per_ring, rho_lox, Lpath_ox[0], dpInjOx, dpFracMani, fOx, KOx, w_init_ref)
w_LOX_outer = manifold_plot(mdot_lox_per_ring, rho_lox, Lpath_ox[1], dpInjOx, dpFracMani, fOx, KOx, w_init_ref)
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
axs[1,0].plot(h_vals*1e3, w_LOX_inner*1e3, marker='o', color='navy')
axs[1,0].set_title("LOX Inner Manifold")
axs[1,0].set_xlabel("h [mm]")
axs[1,0].set_ylabel("w [mm]")
axs[1,0].grid(True)
axs[1,1].plot(h_vals*1e3, w_LOX_outer*1e3, marker='o', color='red')
axs[1,1].set_title("LOX Outer Manifold")
axs[1,1].set_xlabel("h [mm]")
axs[1,1].set_ylabel("w [mm]")
axs[1,1].grid(True)
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
