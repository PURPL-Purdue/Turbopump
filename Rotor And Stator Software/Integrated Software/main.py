"""
Core Turbine Generation Script
Author: Amanjyoti Mridha
"""

import numpy as np
import yaml
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import pandas as pd
import os
import pickle
import json

import turbine_design.structure_functions as sf
import turbine_design.stator_functions as stf
import turbine_design.analysis_functions as af
import turbine_design.turbine_blade_gen as tbg
import turbine_design.turbine_functions as tf
import turbine_design.isentropic_relations as isentropic

import turbine_design.tabbed_figures as tabbed_figures
import turbine_design.helpers as helpers

## OVVERALL CONFIGURATION

yaml_path = "../../params.yaml"
output_table_path = "output_tables/"

with open(yaml_path, 'r') as file:
    params = yaml.safe_load(file)

if os.path.exists(output_table_path) == False:
    os.makedirs(output_table_path)

PLOT_AT_ONCE = True # plot all together (True) or sequentially as individual windows
if PLOT_AT_ONCE:
    pw = tabbed_figures.plotWindow("Turbine Design Plots")
    window_parent = pw.MainWindow
else:
    pw = None

USE_INPUT_GUI = True # use terminal or qt for input questions

## HOT GAS CONSTANTS

P_0 = params['S_inlet_pressure'] * params['psi_to_pa'] # [Pa]
P_e = params['S_design_exit_pressure'] * params['psi_to_pa'] # [Pa]
P_A = P_e # [Pa]
m_dot = params['S_mass_flow_rate'] * params['lbm_to_kg'] # [kg/s]
T_0 = params['S_inlet_temperature'] # [K]
gamma = params['S_inlet_specific_heat'] # [unitless]
R = params['S_inlet_gas_constant'] # [J/(kg*K)]
m_m = params['S_inlet_molar_mass'] # [kg/mol]
num_nozzle = params['S_number_nozzles'] # [unitless]

## GAS CONSTANTS TABLE

# TODO: SAVE GAS CONSTANTS TO A TABLE

## NOZZLE CALCULATIONS

R_S = stf.calc_R_S(R, m_m) # [J/(kg*K)]
rho_0 = stf.calc_rho_0(P_0, R_S, T_0) # [kg/m^3]
v_e = stf.calc_v_e(T_0, R_S, gamma, P_e, P_0) # [m/s]

T_throat = stf.calc_T_throat(T_0, gamma) # [K]
P_throat = stf.calc_P_throat(P_0, gamma) # [Pa]
rho_throat = stf.calc_rho_throat(P_throat, R_S, T_throat) # [kg/m^3]
v_throat = stf.calc_v_throat(gamma, R_S, T_throat) # [m/s]
A_throat = stf.calc_A_throat(m_dot, rho_throat, v_throat) # [m^2]

M_e = stf.calc_M_e(P_e, P_0, gamma) # [unitless]
rho_e = stf.calc_rho_e(rho_0, M_e, gamma)
A_e = stf.calc_A_e(m_dot, rho_e, v_e) # [m^2]
T_e = stf.calc_T_e(T_0, gamma, M_e) # [K]

r_throat = stf.calc_r_throat(A_throat) # [m]
r_e = stf.calc_r_e(A_e) # [m]

A_throat_n = stf.calc_A_throat_n(A_throat, num_nozzle) # [m^2]
A_e_n = stf.calc_A_e_n(A_e, num_nozzle) # [m^2]
r_throat_n = stf.calc_r_throat_n(A_throat_n) # [m]
r_e_n = stf.calc_r_e_n(A_e_n) # [m]
dist_n = stf.calc_dist_n(r_throat_n, r_e_n) # [m]

P_nozzle_out_total = P_e * (1 + (gamma - 1) / 2 * M_e**2) ** (gamma / (gamma - 1))

## PLOT NOZZLE GEOMETRY

[X, Y, Z] = stf.plot_nozzle(A_e, A_throat, A_e, dist_n, dist_n, PLOT_AT_ONCE=PLOT_AT_ONCE, PLOT_WINDOW=pw)

## NOZZLE TABLE

nozzle_table = pd.DataFrame({
    "Parameter": [
        "Inlet Total Pressure (Pa)",
        "Exit Total Pressure (Pa)",
        "Ambient Pressure (Pa)",
        "Mass Flow Rate (kg/s)",
        "Inlet Total Temperature (K)",
        "Specific Heat Ratio (unitless)",
        "Gas Constant (J/(kg*K))",
        "Molar Mass (kg/mol)",
        "Number of Nozzles (unitless)",
        "Specific Gas Constant (J/(kg*K))",
        "Inlet Density (kg/m^3)",
        "Exit Velocity (m/s)",
        "Throat Temperature (K)",
        "Throat Pressure (Pa)",
        "Throat Density (kg/m^3)",
        "Throat Velocity (m/s)",
        "Throat Area (m^2)",
        "Exit Mach Number (unitless)",
        "Exit Density (kg/m^3)",
        "Exit Area (m^2)",
        "Exit Temperature (K)",
        "Throat Radius (m)",
        "Exit Radius (m)",
        "Nozzle Throat Area (m^2)",
        "Nozzle Exit Area (m^2)",
        "Nozzle Throat Radius (m)",
        "Nozzle Exit Radius (m)",
        "Nozzle Length (m)",
        "Nozzle Exit Total Pressure (Pa)"
    ],
    "Value": [
        P_0,
        P_e,
        P_A,
        m_dot,
        T_0,
        gamma,
        R,
        m_m,
        num_nozzle,
        R_S, 
        rho_0, 
        v_e, 
        T_throat, 
        P_throat, 
        rho_throat, 
        v_throat, 
        A_throat, 
        M_e, 
        rho_e, 
        A_e, 
        T_e, 
        r_throat, 
        r_e, 
        A_throat_n, 
        A_e_n, 
        r_throat_n, 
        r_e_n, 
        dist_n, 
        P_nozzle_out_total]
})

nozzle_table.to_csv("output_tables/nozzle_parameters.csv", index=False)

## TURBINE CALCULATIONS

rotor_radius = params['T_rotor_radius'] * params['in_to_m'] # [m]
blade_height = params['T_blade_height'] * params['mm_to_m'] # [m]
hub_radius = rotor_radius - blade_height # [m]
predicted_efficiency = 0.4 # [unitless]
turbine_horsepower = params['T_design_power_output'] / predicted_efficiency # [HP]
blade_solidity = eval(params['T_blade_solidity']) # [unitless]
turbine_rpm = params['T_design_rpm'] # [RPM]

mean_radius = (rotor_radius + hub_radius) / 2 # [m]
turbine_torque = turbine_horsepower * 5252 / turbine_rpm * params['lbft_to_nm'] # [N*m]
degree_of_reacion = 0 # [unitless]
num_blades = params['T_num_blades'] # [unitless]
chord = params['T_chord'] * params['cm_to_m'] # [m]
blade_spacing = 2 * np.pi * mean_radius / num_blades # [m]
beta_guess = 60 # [degrees] initial guess for relative angle, refined later
min_blade_thickness = 0.002 # [m] minimum blade thickness at the root
inlet_area = 0.0025 # [m^2] initial guess for inlet area
pitch = 2 * np.pi * mean_radius / num_blades # [m]

# velocity triangle calculations
v1, v2, w, u, alpha1, alpha2, beta = tf.rotor_back_calculate(
    turbine_rpm, 
    turbine_torque/num_blades, 
    m_dot/num_blades,
    np.deg2rad(beta_guess),
    mean_radius,
    v_e
)

alpha2 = np.mod(alpha2, 2 * np.pi) # ensure alpha2 is between 0 and 2pi

print(f"v1: {v1:.1f} m/s, v2: {v2:.1f} m/s, w: {w:.1f} m/s, u: {u:.1f} m/s, alpha1: {np.rad2deg(alpha1):.1f} deg, alpha2: {np.rad2deg(alpha2):.1f} deg, beta: {np.rad2deg(beta):.1f} deg")

print(type(blade_solidity))
tf.plot_velocity_triangles(v1, v2, u, w, w, chord, 0, beta, -beta, alpha1, -alpha2, PLOT_AT_ONCE=PLOT_AT_ONCE, PLOT_WINDOW=pw)

# Zweifel loading coefficient
zweifel_coeff = af.calc_zweifel_loading_coefficient(blade_solidity, np.deg2rad(beta), np.deg2rad(beta), w, w)
print(f"Zweifel loading coefficient: {zweifel_coeff:.3f} (target between 0.8 and 1.0)")

# blade geometry generation
kacker_okapuu_data = af.load_kacker_okapuu_data(path_str="/KackerOkapuu/kackerokapuu.pkl")
target_distance = 0.0028
rounds = 0.01785
x_lower, y_lower, x_upper, y_upper, distances, dist_plot_params = tbg.compute_geometry(
    chord, 
    np.degrees(beta), 
    target_distance, 
    blade_spacing, 
    rounds, 
    output_excel=True, 
    output_path=f"{output_table_path}design_table.xlsx"
)
tbg.plot_blade_geometry(x_lower, y_lower, x_upper, y_upper, PLOT_AT_ONCE=PLOT_AT_ONCE, PLOT_WINDOW=pw)

# plot the distances along the blade
print("Plotting distances along the blade...")
tbg.compute_distances_and_plot(*dist_plot_params, PLOT_AT_ONCE=PLOT_AT_ONCE, PLOT_WINDOW=pw)


# efficiency calculations
M_inlet_rel = w / np.sqrt(gamma * R_S * T_e)
M_inlet_abs = v1 / np.sqrt(gamma * R_S * T_e)
M_outlet_abs = v2 / np.sqrt(gamma * R_S * T_e)
C_p =  (10.8660+11.6259+12.9042)/3 * 0.238845896632776 # from textbook: 0.652 # Btu/(lb * deg F)
T_stator_inlet_imperial = T_0 * 9/5 # [R]

GAMMA = 1.1201
P_OUT_STATIC = 2.0684e+05 # [Pa]
TOTAL_PRESSURE_RATIO = 2.4132e+06 / P_0 # [unitless]
Q_FREESTREAM = 0.5 * rho_e * w**2
REYNOLDS_NUM = rho_e * w * chord / params['T_dynamic_viscosity'] # [unitless]
RH_RT = hub_radius / rotor_radius # [unitless]
PRESSURE_RATIO = 1 # [unitless] 
R_INLET = 733.98 # [J/(kg*K)] specific gas constant
P_INLET = 2.0684e+05 # [Pa] 

n_total, expected_hp = af.calculate_blade_design_efficiency(
    distances, 
    target_distance, 
    chord,
    blade_height,
    pitch,
    alpha1,
    alpha2,
    beta,
    C_p,
    turbine_horsepower,
    params['S_mass_flow_rate'],
    GAMMA,
    TOTAL_PRESSURE_RATIO,
    Q_FREESTREAM,
    REYNOLDS_NUM,
    RH_RT,
    PRESSURE_RATIO,
    R,
    T_e,
    T_stator_inlet_imperial,
    P_INLET,
    M_inlet_rel,
    M_inlet_rel, # outlet rel is the same
    kacker_okapuu_data)

print(f"Predicted turbine efficiency: {n_total:.3f} (expected horsepower: {expected_hp:.1f} HP)")
expected_torque  = expected_hp * 5252 / turbine_rpm * params['lbft_to_nm']
print(f"Expected torque: {expected_torque:.1f} Nm (design torque: {turbine_torque:.1f} Nm)")

# compute and plot mach and pressure distributions
areas = np.concatenate([distances, np.flip(distances)])
a_star = distances[0] / af.area_mach_relation(M_inlet_rel, gamma)
print(f"A* = {a_star:.6f}")
M_vec, P_vec = af.calculate_mach_pressure_distribution(areas, gamma, R, T_e, P_e, M_inlet_rel, M_inlet_rel, a_star)
af.plot_mach_pressure_distributions(M_vec, P_vec, PLOT_AT_ONCE=PLOT_AT_ONCE, PLOT_WINDOW=pw)

# 3D plots
if not USE_INPUT_GUI:
    plot_turbine = input("Plot 3D turbine model? (y/n): ")
    if plot_turbine.lower() == 'y':
        tf.plot_turbine(x_lower, y_lower, x_upper, y_upper, num_blades, hub_radius, chord, rotor_radius - hub_radius, min_blade_thickness)
else:
    plot_turbine = helpers.ask_yes_no("Plot 3D turbine model?", window_parent)
    if plot_turbine:
        tf.plot_turbine(x_lower, y_lower, x_upper, y_upper, num_blades, hub_radius, chord, rotor_radius - hub_radius, min_blade_thickness)

# write turbine results to csv
turbine_table = pd.DataFrame({
    "Parameter": [
        "Mean Radius (m)",
        "Blade Height (m)",
        "Hub Radius (m)",
        "Rotor Radius (m)",
        "Blade Solidity (unitless)",
        "Turbine RPM (RPM)",
        "Turbine Torque (Nm)",
        "Expected Torque (Nm)",
        "Degree of Reaction (unitless)",
        "Number of Blades (unitless)",
        "Chord (m)",
        "Blade Spacing (m)",
        "Inlet Area (m^2)",
        "Inlet Velocity (m/s)",
        "Outlet Velocity (m/s)",
        "Relative Velocity (m/s)",
        "Inlet Absolute Angle (deg)",
        "Outlet Absolute Angle (deg)",
        "Inlet Relative Angle (deg)",
        "Zweifel Loading Coefficient (unitless)",
        "Predicted Turbine Efficiency (unitless)",
        "Expected Horsepower (HP)"
    ],
    "Value": [
        mean_radius,
        blade_height,
        hub_radius,
        rotor_radius,
        blade_solidity,
        turbine_rpm,
        turbine_torque,
        expected_torque,
        degree_of_reacion,
        num_blades,
        chord,
        blade_spacing,
        inlet_area,
        v1,
        v2,
        w,
        np.rad2deg(alpha1),
        np.rad2deg(alpha2),
        np.rad2deg(beta),
        zweifel_coeff,
        n_total,
        expected_hp
    ]
})

turbine_table.to_csv("output_tables/turbine_parameters.csv", index=False)


# structural calculations
omega = turbine_rpm * 2 * np.pi / 60 # [rad/s]
height_blade = target_distance # [m]
rho_blade = 0.845; # [kg/m^3] density
length_blade = chord # [m]
height_bmin = min_blade_thickness # [m]
width_blade = rotor_radius - hub_radius # [m]
Z_blade = num_blades # [unitless]
mean_diameter = 2 * mean_radius # [m]

radius_turbine = sf.calc_radius_turbine(height_blade, hub_radius); # [m]
mass_blade = sf.calc_mass_blade(length_blade, width_blade, height_bmin, rho_blade); # [kg]
Force_centrifugal = sf.calc_Force_centrifugal(mass_blade, omega, radius_turbine); # [N]
stress_centrifugal = sf.calc_stress_centrifugal(rho_blade, height_blade, mean_diameter, omega); # [N/m^2]
Force_tangential = sf.calc_Force_tangential(m_dot, w, beta, w, -beta); # [N]
Force_axial = sf.calc_Force_axial(m_dot, v1, alpha1, w, -beta); # [N]
torque_blade = sf.calc_torque_blade(Force_tangential, height_blade); # [Nm]
torque_turbine = sf.calc_torque_turbine(Force_tangential, radius_turbine, Z_blade); # [Nm]
P = sf.calc_P(torque_turbine, omega); # [Nm/s]
Force_gas = sf.calc_Force_gas(Force_tangential, Force_axial); # [N]
Moment_Bending = sf.calc_Moment_Bending(height_blade, Z_blade, Force_gas); # [Nm]
I = sf.calc_I(length_blade, height_bmin); # [m^4]
stress_gas = sf.calc_stress_gas(height_bmin, Force_gas, width_blade, I); # [N/m^2]

# write results to csv
structure_table = pd.DataFrame({
    "Parameter": [
        "Force_centrifugal (N)",
        "stress_centrifugal (N/m^2)",
        "Force_tangential (N)",
        "Force_axial (N)",
        "torque_blade (Nm)",
        "torque_turbine (Nm)",
        "P (W)",
        "Force_gas (N)",
        "Moment_Bending (Nm)",
        "I (m^4)",
        "stress_gas (N/m^2)"
    ],
    "Value": [
        Force_centrifugal,
        stress_centrifugal,
        Force_tangential,
        Force_axial,
        torque_blade,
        torque_turbine,
        P,
        Force_gas,
        Moment_Bending,
        I,
        stress_gas
    ]
})

structure_table.to_csv("output_tables/turbine_results.csv", index=False)

"""
Off design calculations (vary mass flow rate)
"""

m_dot_offdesign = np.linspace(0.1 * m_dot, 1.5 * m_dot, 50) # [kg/s]
v_in_off_design = m_dot_offdesign / (rho_e * A_e) # [m/s]
off_design_params = [af.velocity_triangles_forward_calculate(v_off_design, alpha1, beta, mean_radius, m_dot_offdesign[i]) for i, v_off_design in enumerate(v_in_off_design)] # [(v2, alpha2, w1, U, torque)]

v2_offdesign = np.array([params[0] for params in off_design_params])
alpha2_offdesign = np.array([params[1] for params in off_design_params])
w_offdesign = np.array([params[2] for params in off_design_params])
u_offdesign = np.array([params[3] for params in off_design_params])
torque_offdesign = np.array([params[4] for params in off_design_params])

turbine_rpm = u_offdesign * 60 / (2 * np.pi * mean_radius)
# replace nan or singular values with 0
turbine_rpm = np.nan_to_num(turbine_rpm, nan=0.0, posinf=0.0, neginf=0.0)
turbine_horsepower = torque_offdesign * turbine_rpm / 5252 / params['lbft_to_nm']
m_dot_offdesign_lbm = m_dot_offdesign / params['lbm_to_kg']

# I printed to debug
# print("off_design_params", off_design_params)
# print("alpha2_offdesign:", alpha2_offdesign)
# print("w_offdesign:", w_offdesign)
# input("Press Enter to continue...")

# calculate efficiency for each off-design point
n_total_offdesign = []
expected_hp_offdesign = []
for i in range(len(m_dot_offdesign)):
    try:
        M_inlet_rel_offdesign = w_offdesign[i] / np.sqrt(gamma * R_S * T_e)

        n_total, expected_hp = af.calculate_blade_design_efficiency(
            distances, 
            target_distance, 
            chord,
            blade_height,
            pitch,
            alpha1,
            alpha2_offdesign[i],
            beta,
            C_p,
            turbine_horsepower[i],
            m_dot_offdesign[i],
            GAMMA,
            TOTAL_PRESSURE_RATIO,
            Q_FREESTREAM,
            REYNOLDS_NUM,
            RH_RT,
            PRESSURE_RATIO,
            R,
            T_e,
            T_stator_inlet_imperial,
            P_INLET,
            M_inlet_rel_offdesign,
            M_inlet_rel_offdesign, # outlet rel is the same
            kacker_okapuu_data,
            use_print=False)
    except Exception as e:
        print(f"Error calculating efficiency for mass flow rate {m_dot_offdesign[i]:.3f} kg/s: {e}")
        n_total = 0
        expected_hp = 0

    if np.isnan(n_total): # n_total < 0 or n_total > 1 or
        n_total = 0
        expected_hp = 0
    n_total_offdesign.append(n_total)
    expected_hp_offdesign.append(expected_hp)

# plot 3D plot of off-design turbine rpm and horsepower vs mass flow rate
fig = plt.figure(figsize=(10, 6))
ax = fig.add_subplot(111, projection='3d')
# plot as scatter
ax.scatter(m_dot_offdesign_lbm, turbine_rpm, expected_hp_offdesign, color='r', label='Data Points')

# plot best fit surface (quadratic)
# m_dot_mesh, rpm_mesh = np.meshgrid(m_dot_offdesign_lbm, turbine_rpm)
# coeffs = np.polyfit(m_dot_offdesign_lbm, turbine_rpm, 2)
# rpm_fit = np.polyval(coeffs, m_dot_mesh)
# coeffs_hp = np.polyfit(m_dot_offdesign_lbm, turbine_horsepower, 2)
# hp_fit = np.polyval(coeffs_hp, m_dot_mesh)
# ax.plot_surface(m_dot_mesh, rpm_mesh, hp_fit, alpha=0.5, color='b', label='Best Fit Surface')

# plot limit plane (HP = 0)
m_dot_plane = np.linspace(min(m_dot_offdesign_lbm), max(m_dot_offdesign_lbm), 10)
rpm_plane = np.linspace(min(turbine_rpm), max(turbine_rpm), 10)
m_dot_plane, rpm_plane = np.meshgrid(m_dot_plane, rpm_plane)
hp_plane = np.zeros_like(m_dot_plane)
ax.plot_surface(m_dot_plane, rpm_plane, hp_plane, alpha=0.3, color='gray', label='HP = 0 Plane')

# 100 HP plane
hp_100_plane = np.full_like(m_dot_plane, 100)
ax.plot_surface(m_dot_plane, rpm_plane, hp_100_plane, alpha=0.3, color='green', label='HP = 100 Plane')

ax.set_xlabel('Mass Flow Rate (lbm/s)')
ax.set_ylabel('Turbine RPM (RPM)')
ax.set_zlabel('Turbine Power (HP)')
ax.set_title('Off-Design Turbine Performance')
ax.legend()
if not PLOT_AT_ONCE:
    plt.show()
else:
    pw.addPlot("Off-Design Turbine Performance", fig)

# save fig to png
fig.savefig(f"{output_table_path}off_design_turbine_performance.png")

# create additional plot that hides all data points below the 100 HP plane
fig = plt.figure(figsize=(10, 6))
ax = fig.add_subplot(111, projection='3d')
# filter data points above 100 HP
mask = np.array(expected_hp_offdesign) >= 100
ax.scatter(np.array(m_dot_offdesign_lbm)[mask], np.array(turbine_rpm)[mask], np.array(expected_hp_offdesign)[mask], color='r', label='Data Points (HP >= 100)') 
ax.set_xlabel('Mass Flow Rate (lbm/s)')
ax.set_ylabel('Turbine RPM (RPM)')
ax.set_zlabel('Turbine Power (HP)')
ax.set_title('Off-Design Turbine Performance (HP >= 100)')
ax.legend()
if not PLOT_AT_ONCE:
    plt.show()
else:
    pw.addPlot("Off-Design Turbine Performance (HP >= 100)", fig)
# save fig to png
fig.savefig(f"{output_table_path}off_design_turbine_performance_hp100.png")
# save corresponding data to csv
off_design_data = pd.DataFrame({
    "Mass Flow Rate (lbm/s)": np.array(m_dot_offdesign_lbm)[mask],
    "Turbine RPM (RPM)": np.array(turbine_rpm)[mask],
    "Turbine Power (HP)": np.array(expected_hp_offdesign)[mask]
})
off_design_data.to_csv(f"{output_table_path}off_design_turbine_performance_hp100.csv", index=False)

# also plot mass flow rate vs efficiency
fig = plt.figure(figsize=(8, 5))
ax = fig.add_subplot(111)
ax.plot(m_dot_offdesign_lbm, n_total_offdesign, marker='o')
ax.set_xlabel('Mass Flow Rate (lbm/s)')
ax.set_ylabel('Turbine Efficiency (unitless)')
ax.set_title('Off-Design Turbine Efficiency')
ax.grid(True)
ax.set_ylim(0, 1)
if not PLOT_AT_ONCE:
    plt.show()
else:
    pw.addPlot("Off-Design Turbine Efficiency", fig)
# save fig to png
fig.savefig(f"{output_table_path}off_design_turbine_efficiency.png")

# plot masss flow rate vs torque
fig = plt.figure(figsize=(8, 5))
ax = fig.add_subplot(111)
ax.plot(m_dot_offdesign_lbm, torque_offdesign, marker='o', color='orange')
ax.set_xlabel('Mass Flow Rate (lbm/s)')
ax.set_ylabel('Turbine Torque (Nm)')
ax.set_title('Off-Design Turbine Torque')
ax.grid(True)
if not PLOT_AT_ONCE:
    plt.show()
else:
    pw.addPlot("Off-Design Turbine Torque", fig)
# save fig to png
fig.savefig(f"{output_table_path}off_design_turbine_torque.png")

"""
Off design calculations (vary mass flow rate and rpm)
"""

# m_dot_offdesign = np.linspace(0.1 * m_dot, 1.5 * m_dot, 50) # [kg/s]
# v_in_off_design = m_dot_offdesign / (rho_e * A_e) # [m/s]
# rpm_offdesign = np.linspace(0.5 * turbine_rpm, 1.5 * turbine_rpm, 50) # [RPM]
# u_offdesign = rpm_offdesign * 2 * np.pi * mean_radius / 60 # [m/s]

turbine_rpm = params['T_design_rpm'] # [RPM]
m_dot_range = np.linspace(0.1 * m_dot, 1.5 * m_dot, 50)
rpm_range = np.linspace(0.5 * turbine_rpm, 1.5 * turbine_rpm, 50)

# save mdot range and rpm range np arrays
np.savetxt(f"{output_table_path}off_design_mdot_range.csv", m_dot_range, delimiter=",")
np.savetxt(f"{output_table_path}off_design_rpm_range.csv", rpm_range, delimiter=",")

if not USE_INPUT_GUI:
    create_meshgrid = input("Create meshgrid for mass flow rate and rpm? (y/n): ").lower() == 'y'
else:
    create_meshgrid = helpers.ask_yes_no("Create meshgrid for mass flow rate and rpm?", window_parent)

if create_meshgrid:

    M_DOT, RPM = np.meshgrid(m_dot_range, rpm_range)

    V_IN = M_DOT / (rho_e * A_e)
    U = RPM * 2 * np.pi * mean_radius / 60

    torque_surface = np.zeros_like(M_DOT)
    hp_surface = np.zeros_like(M_DOT)
    alpha2_surface = np.zeros_like(M_DOT)
    w_surface = np.zeros_like(M_DOT)

    for i in range(M_DOT.shape[0]):
        for j in range(M_DOT.shape[1]):

            v_in = V_IN[i, j]
            rpm_val = RPM[i, j]
            u_val = U[i, j]
            mdot_val = M_DOT[i, j]

            v2, alpha2, w1, U_calc, torque = af.velocity_triangles_unsteady_forward_calculate(
                v_in,
                alpha1,
                beta,
                mean_radius,
                mdot_val,
                rpm_val,
                u_val
            )

            alpha2_surface[i, j] = alpha2
            w_surface[i, j] = w1
            torque_surface[i, j] = torque

            hp_surface[i, j] = torque * rpm_val / 5252 / params['lbft_to_nm']


    eff_surface = np.zeros_like(M_DOT)
    expected_hp_surface = np.zeros_like(M_DOT)

    for i in range(M_DOT.shape[0]):
        for j in range(M_DOT.shape[1]):

            try:
                M_inlet_rel = w_surface[i, j] / np.sqrt(gamma * R_S * T_e)

                n_total, expected_hp = af.calculate_blade_design_efficiency(
                    distances,
                    target_distance,
                    chord,
                    blade_height,
                    pitch,
                    alpha1,
                    alpha2_surface[i, j],
                    beta,
                    C_p,
                    hp_surface[i, j],
                    M_DOT[i, j],
                    GAMMA,
                    TOTAL_PRESSURE_RATIO,
                    Q_FREESTREAM,
                    REYNOLDS_NUM,
                    RH_RT,
                    PRESSURE_RATIO,
                    R,
                    T_e,
                    T_stator_inlet_imperial,
                    P_INLET,
                    M_inlet_rel,
                    M_inlet_rel,
                    kacker_okapuu_data,
                    use_print=False
                )

            except Exception:
                n_total = 0
                expected_hp = 0

            if np.isnan(n_total):
                n_total = 0
                expected_hp = 0

            eff_surface[i, j] = n_total
            expected_hp_surface[i, j] = expected_hp

    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')

    ax.plot_surface(RPM, M_DOT / params['lbm_to_kg'], hp_surface)

    ax.set_xlabel("RPM")
    ax.set_ylabel("Mass Flow (lbm/s)")
    ax.set_zlabel("Horsepower")
    ax.set_title("Off-Design Horsepower Surface")

    if not PLOT_AT_ONCE:
        plt.show()
    else:
        pw.addPlot("Off-Design Horsepower Surface", fig)

    # save data
    m_dot_surface_lbm = M_DOT / params['lbm_to_kg']

    rpm_flat = RPM.flatten()
    mdot_flat = m_dot_surface_lbm.flatten()
    hp_flat = hp_surface.flatten()
    expected_hp_flat = expected_hp_surface.flatten()
    eff_flat = eff_surface.flatten()
    alpha2_flat = alpha2_surface.flatten()
    w_flat = w_surface.flatten()
    torque_flat = torque_surface.flatten()

    mask = np.isfinite(hp_flat) & (hp_flat > 0)

    off_design_data = pd.DataFrame({
        "Mass Flow Rate (lbm/s)": mdot_flat[mask],
        "Turbine RPM (RPM)": rpm_flat[mask],
        "Torque (Nm)": torque_flat[mask],
        "Horsepower (HP)": hp_flat[mask],
        "Expected Horsepower (HP)": expected_hp_flat[mask],
        "Efficiency": eff_flat[mask],
        "Alpha2 (rad)": alpha2_flat[mask],
        "Relative Velocity W1 (m/s)": w_flat[mask]
    })
    off_design_data.to_csv(
        f"{output_table_path}off_design_turbine_performance_surface.csv",
        index=False
    )

    np.savetxt(f"{output_table_path}hp_surface.csv", hp_surface, delimiter=",")
    np.savetxt(f"{output_table_path}eff_surface.csv", eff_surface, delimiter=",")

    # create a 2D matrix of data and save to csv
    out = np.empty((len(rpm_range) + 1, len(m_dot_range) + 1))
    out[:] = np.nan

    # Fill headers
    out[0, 1:] = m_dot_range / params['lbm_to_kg'] # lbm/s
    out[1:, 0] = rpm_range

    # Fill z values
    out[1:, 1:] = hp_surface

    np.savetxt(f"{output_table_path}horsepower_surface.csv", out, delimiter=",")


## Nitrogen Cold Gas Calculations
def calc_cold_gas_performance():
    # put in func so that variables are local scope as I don't want to 
    # mess with the hot gas variables above
    gamma_N = 1.4
    R_N = 296.8 # J/(kg*K)
    rho_N = 1.2506 # kg/m^3

    Upstream_Temperature = 300 # K
    num_points = 100

    n2_standard_temp = 298; # K
    n2_standard_pressure = 101325; # Pa

    # estmate the inlet area from hot gas calculations
    # A_ratio = isentropic.calc_A_ratio(0.01, gamma)
    # A_inlet = A_throat / A_ratio

    # calculate mach numbers for cold gas
    # A_crit = A_inlet / isentropic.calc_A_ratio(0.01, gamma_N)
    # if A_crit <= A_throat:
    #     M_inlet = 1
    # else:
    #     M_inlet = iserntropic.calc_M_from_A_ratio(A_inlet / )
    # M_throat = isentropic.calc_M_from_A_ratio(A_throat / A_inlet, gamma_N, M1=0.01)
    # M_exit = isentropic.calc_M_from_A_ratio(A_e / A_throat, gamma_N, M1=M_throat)

  
    tescom_coeff_file = "./component_data/tescom26_1100_polyfit_coefficients.json"
    try:
        with open(tescom_coeff_file, 'r') as f:
            tescom_data = json.load(f)
            tescom_coeffs = tescom_data['Coefficients']
            tescom_pressure_range_psi = tescom_data['Static Pressure Range (Gauge)']
            tescom_flow_rate_range_scfm = tescom_data['Flow Rate Range (SCFM)']
            tescom_outlet_area_in2 = tescom_data['Outlet Area (in^2)']
            tescom_raw_data = tescom_data['Raw Data Points']
    except Exception as e:
        print(f"Error loading Tescom data: {e}\n Go to component_data and run the processing script")
        return

    # calculate Mach number at throat and exit for cold gas
    A_tescom = tescom_outlet_area_in2 * params['in_to_m']**2
    flow_rate_scfm = np.linspace(tescom_flow_rate_range_scfm[0], tescom_flow_rate_range_scfm[1], num_points)
    tescom_pressures_psig = np.polyval(tescom_coeffs, flow_rate_scfm)
    tescom_pressures_pa = tescom_pressures_psig * params['psi_to_pa'] + params['P_atm_pa']
    flow_rate_acfm = flow_rate_scfm / (tescom_pressures_pa / n2_standard_pressure) / (n2_standard_temp / Upstream_Temperature)
    flow_rate_m3_s = flow_rate_acfm * params['cfm_to_m3ps']
    
    # assuming upstream temperature, calculate density
    density_N2 = tescom_pressures_pa / (R_N * Upstream_Temperature) # ideal gas law

    mass_flow_rate = flow_rate_m3_s * density_N2 # kg/s
    mass_flow_rate_lbm_s = mass_flow_rate / params['lbm_to_kg'] # lbm/s
    velocity_tescom = mass_flow_rate / (density_N2 * A_tescom)
    # T_tescom = [tescom_pressures_pa[i] / (rho_N * R_N) for i in range(num_points)] # ideal gas law
    T_tescom = [Upstream_Temperature for i in range(num_points)] # assume constant temperature
    M_tescom = [velocity_tescom[i] / np.sqrt(gamma_N * R_N * T_tescom[i]) for i in range(num_points)]
    rho_0 = [density_N2[i] * isentropic.calc_rho0_rho_ratio(M_tescom[i], gamma_N) for i in range(num_points)]
    # T0 = [Upstream_Temperature * isentropic.calc_T0_T_ratio(M_tescom[i], gamma_N) for i in range(num_points)]
    T0 = [tescom_pressures_pa[i] / (rho_0[i] * R_N) for i in range(num_points)] 
    P_0 = [tescom_pressures_pa[i] * isentropic.calc_P0_P_ratio(M_tescom[i], gamma_N) for i in range(num_points)]

    # plot raw tescom data with polynomial fit
    # plt.figure(figsize=(8, 5))
    # plt.plot(flow_rate_scfm, tescom_pressures_psig, label='Tescom Static Pressure (Gauge)', color='orange')
    # plt.scatter(
    #     [point[0] for point in tescom_raw_data],
    #     [point[1] for point in tescom_raw_data],
    #     label='Tescom Raw Data Points',
    #     color='red',
    #     marker='x'
    # )    
    # plt.xlabel('Flow Rate (SCFM)')
    # plt.ylabel('Static Pressure (psig)')
    # plt.title('Tescom 26-1100 Static Pressure vs Flow Rate')
    # plt.grid()
    # plt.legend()
    # if not PLOT_AT_ONCE:
    #     plt.show()
    # else:
    #     pw.addPlot("Tescom 26-1100 Static Pressure vs Flow Rate", plt.gcf())

    # plot tescom data
    plt.figure(figsize=(8, 5))
    plt.plot(mass_flow_rate_lbm_s, M_tescom, label='Tescom Mach Number')
    plt.xlabel('Mass Flow Rate (lbm/s)')
    plt.ylabel('Mach Number')
    plt.legend()
    plt.title('Tescom 26-1100 Mach Number & Static Pressure vs Mass FLow Rate')
    plt.twinx()
    plt.plot(mass_flow_rate_lbm_s, tescom_pressures_pa, label='Tescom Static Pressure', color='orange')
    plt.ylabel('Static Pressure (Pa)')
    plt.grid()
    plt.legend()
    if not PLOT_AT_ONCE:
        plt.show()
    else:
        pw.addPlot("Tescom 26-1100 Mach Number vs Mass Flow Rate", plt.gcf())

    # calculate mach numbers
    # at throat
    M_throat = [isentropic.calc_M_from_A_ratio(A_tescom / A_throat, gamma_N, M1=M_tescom[i]) for i in range(num_points)]
    # at stator exit
    M_exit = [isentropic.calc_M_from_A_ratio(A_throat / A_e, gamma_N, M1=M_throat[i]) for i in range(num_points)]
    rho_exit = [rho_0[i] * isentropic.calc_rho_rho0_ratio(M_exit[i], gamma_N) for i in range(num_points)]
    T_e = [T0[i] * isentropic.calc_T_T0_ratio(M_exit[i], gamma_N) for i in range(num_points)]    
    P_e = [P_0[i] * isentropic.calc_P_P0_ratio(M_exit[i], gamma_N) for i in range(num_points)]
    V_exit = [M_exit[i] * np.sqrt(gamma_N * R_N * T_e[i]) for i in range(num_points)]

    # plot mach numbers
    plt.figure(figsize=(8, 5))
    plt.plot(mass_flow_rate_lbm_s, M_throat, label='Throat Mach Number')
    plt.plot(mass_flow_rate_lbm_s, M_exit, label='Exit Mach Number')
    plt.xlabel('Mass Flow Rate (lbm/s)')
    plt.ylabel('Mach Number')
    plt.title('Mach Number at Throat and Exit vs Mass Flow Rate')
    plt.grid()
    plt.legend()
    if not PLOT_AT_ONCE:
        plt.show()
    else:
        pw.addPlot("Mach Number at Throat and Exit vs Mass Flow Rate", plt.gcf())


    turbine_rpm = params['T_design_rpm'] # [RPM]
    rpm_range = np.linspace(0.5 * turbine_rpm, 1.5 * turbine_rpm, 50)

    # save mdot range and rpm range np arrays
    np.savetxt(f"{output_table_path}off_design_mdot_range_n2.csv", mass_flow_rate, delimiter=",")
    np.savetxt(f"{output_table_path}off_design_rpm_range_n2.csv", rpm_range, delimiter=",")

    if not USE_INPUT_GUI:
        create_meshgrid = input("Create meshgrid for mass flow rate and rpm? (y/n): ").lower() == 'y'
    else:
        create_meshgrid = helpers.ask_yes_no("Create meshgrid for mass flow rate and rpm?", window_parent)

    if create_meshgrid:

        M_DOT, RPM = np.meshgrid(mass_flow_rate, rpm_range)
        # print("M_DOT shape:", M_DOT.shape)
        # print("RPM shape:", RPM.shape)

        V_IN = np.meshgrid(V_exit, rpm_range)[0] # just repeat the exit velocity across the rpm range
        U = RPM * 2 * np.pi * mean_radius / 60
        T_e_mesh = np.meshgrid(T_e, rpm_range)[0] # repeat exit temperature across rpm range
        P_turbine_inlet_mesh = np.meshgrid(P_e, rpm_range)[0] # repeat exit pressure across rpm range

        torque_surface = np.zeros_like(M_DOT)
        hp_surface = np.zeros_like(M_DOT)
        alpha2_surface = np.zeros_like(M_DOT)
        w_surface = np.zeros_like(M_DOT)

        for i in range(M_DOT.shape[0]):
            for j in range(M_DOT.shape[1]):

                v_in = V_IN[i, j]
                rpm_val = RPM[i, j]
                u_val = U[i, j]
                mdot_val = M_DOT[i, j]

                v2, alpha2, w1, U_calc, torque = af.velocity_triangles_unsteady_forward_calculate(
                    v_in,
                    alpha1,
                    beta,
                    mean_radius,
                    mdot_val,
                    rpm_val,
                    u_val
                )

                alpha2_surface[i, j] = alpha2
                w_surface[i, j] = w1
                torque_surface[i, j] = torque

                hp_surface[i, j] = torque * rpm_val / 5252 / params['lbft_to_nm']


        eff_surface = np.zeros_like(M_DOT)
        expected_hp_surface = np.zeros_like(M_DOT)

        for i in range(M_DOT.shape[0]):
            for j in range(M_DOT.shape[1]):

                try:
                    M_inlet_rel = w_surface[i, j] / np.sqrt(gamma_N * R_N * T_e_mesh[i, j])
                    n2_total_pressure_ratio = isentropic.calc_total_pressure_ratio(M_inlet_rel, gamma_N)

                    n_total, expected_hp = af.calculate_blade_design_efficiency(
                        distances,
                        target_distance,
                        chord,
                        blade_height,
                        pitch,
                        alpha1,
                        alpha2_surface[i, j],
                        beta,
                        C_p,
                        hp_surface[i, j],
                        M_DOT[i, j],
                        gamma_N,
                        n2_total_pressure_ratio,
                        Q_FREESTREAM,
                        REYNOLDS_NUM,
                        RH_RT,
                        PRESSURE_RATIO,
                        R_N,
                        T_e_mesh[i, j],
                        T_stator_inlet_imperial,
                        P_turbine_inlet_mesh[i, j],
                        M_inlet_rel,
                        M_inlet_rel,
                        kacker_okapuu_data,
                        use_print=False
                    )

                except Exception:
                    n_total = 0
                    expected_hp = 0


                if np.isnan(n_total):
                    n_total = 0
                    expected_hp = 0

                eff_surface[i, j] = n_total
                expected_hp_surface[i, j] = expected_hp

        fig = plt.figure()
        ax = fig.add_subplot(111, projection='3d')

        ax.plot_surface(RPM, M_DOT / params['lbm_to_kg'], hp_surface)

        ax.set_xlabel("RPM")
        ax.set_ylabel("Mass Flow (lbm/s)")
        ax.set_zlabel("Horsepower")
        ax.set_title("Off-Design Horsepower Surface")

        if not PLOT_AT_ONCE:
            plt.show()
        else:
            pw.addPlot("Off-Design Horsepower Surface (N2)", fig)
            
        # hp contour
        plt.figure()
        cp = plt.contourf(RPM, M_DOT / params['lbm_to_kg'], hp_surface, levels=50, cmap='viridis')
        plt.colorbar(cp)
        plt.xlabel("RPM")
        plt.ylabel("Mass Flow (lbm/s)")
        plt.title("Off-Design Horsepower Contour (N2)")
        if not PLOT_AT_ONCE:
            plt.show()
        else:
            pw.addPlot("Off-Design Horsepower Contour (N2)", plt.gcf())

        # save data
        m_dot_surface_lbm = M_DOT / params['lbm_to_kg']

        rpm_flat = RPM.flatten()
        mdot_flat = m_dot_surface_lbm.flatten()
        hp_flat = hp_surface.flatten()
        expected_hp_flat = expected_hp_surface.flatten()
        eff_flat = eff_surface.flatten()
        alpha2_flat = alpha2_surface.flatten()
        w_flat = w_surface.flatten()
        torque_flat = torque_surface.flatten()

        mask = np.isfinite(hp_flat) & (hp_flat > 0)

        off_design_data = pd.DataFrame({
            "Mass Flow Rate (lbm/s)": mdot_flat[mask],
            "Turbine RPM (RPM)": rpm_flat[mask],
            "Torque (Nm)": torque_flat[mask],
            "Horsepower (HP)": hp_flat[mask],
            "Expected Horsepower (HP)": expected_hp_flat[mask],
            "Efficiency": eff_flat[mask],
            "Alpha2 (rad)": alpha2_flat[mask],
            "Relative Velocity W1 (m/s)": w_flat[mask]
        })
        off_design_data.to_csv(
            f"{output_table_path}off_design_turbine_performance_surface_n2.csv",
            index=False
        )

        np.savetxt(f"{output_table_path}hp_surface_n2.csv", hp_surface, delimiter=",")
        np.savetxt(f"{output_table_path}eff_surface_n2.csv", eff_surface, delimiter=",")

        # create a 2D matrix of data and save to csv
        out = np.empty((len(rpm_range) + 1, len(mass_flow_rate) + 1))
        out[:] = np.nan

        # Fill headers
        out[0, 1:] = mass_flow_rate / params['lbm_to_kg'] # lbm/s
        out[1:, 0] = rpm_range

        # Fill z values
        out[1:, 1:] = hp_surface

        np.savetxt(f"{output_table_path}horsepower_surface_n2.csv", out, delimiter=",")


if not USE_INPUT_GUI:
    run_n2_calcs = input("Run nitrogen cold gas performance calculations? (y/n): ").lower() == 'y'
else:
    run_n2_calcs = helpers.ask_yes_no("Run nitrogen cold gas performance calculations?", window_parent)

if run_n2_calcs:
    calc_cold_gas_performance()


if PLOT_AT_ONCE:
    # plt.show()
    pw.show()