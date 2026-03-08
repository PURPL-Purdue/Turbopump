import numpy as np
import yaml
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import pandas as pd
import os
import pickle


import turbine_design.structure_functions as sf
import turbine_design.stator_functions as stf
import turbine_design.analysis_functions as af
import turbine_design.turbine_blade_gen as tbg
import turbine_design.turbine_functions as tf

yaml_path = "../../params.yaml"
output_table_path = "output_tables/"

with open(yaml_path, 'r') as file:
    params = yaml.safe_load(file)

if os.path.exists(output_table_path) == False:
    os.makedirs(output_table_path)

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

[X, Y, Z] = stf.plot_nozzle(A_e, A_throat, A_e, dist_n, dist_n)

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
tf.plot_velocity_triangles(v1, v2, u, w, w, chord, 0, beta, -beta, alpha1, -alpha2)

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
tbg.plot_blade_geometry(x_lower, y_lower, x_upper, y_upper)

# plot the distances along the blade
print("Plotting distances along the blade...")
tbg.compute_distances_and_plot(*dist_plot_params)


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
af.plot_mach_pressure_distributions(M_vec, P_vec)

# 3D plots
plot_turbine = input("Plot 3D turbine model? (y/n): ")
if plot_turbine.lower() == 'y':
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
plt.show()
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
plt.show()
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
plt.show()
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
plt.show()
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

create_meshgrid = input("Create meshgrid for mass flow rate and rpm? (y/n): ")

if create_meshgrid.lower() == 'y':

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

    plt.show()

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
gamma_N = 1.4
R_N = 296.8 # J/(kg*K)
