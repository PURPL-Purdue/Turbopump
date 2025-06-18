import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import minimize, fsolve, fmin
from matplotlib.widgets import Slider
import random
import pandas as pd
import pickle
import os
import requests
import pathlib
from scipy.io import loadmat

"""
HELPERS
"""

def lerp(x1, y1, x2, y2, x):
    return (y2 - y1) / (x2 - x1) * (x-x1) + y1

num_points = 100

"""
CONSTANTS FROM HOT_FIRE_TURBINE.M (outputs)
"""

var_path = "../turbine_vars.mat"
var_data = loadmat(var_path)

INLET_MACH = float(var_data['M_inlet_rel'].item())
OUTLET_MACH = INLET_MACH
STATOR_EXIT_MACH = 2.2402
GAMMA = 1.1201 # -
T_INLET = 673.18 # Kelvin
R_INLET = 733.98 # specific gas constant J/(kgK)
P_INLET = 2.0684e+05 # Pa

DYNAMIC_VISCOSITY = 5.5522e-5 # Pa*s 0.55522
HORSE_POWER = float(var_data['horse_power'].item()) # hp
M_DOT = float(var_data['mass_flow'].item())  # kg/s
M_DOT_IMPERIAL = float(var_data['m_dot_imperial'].item()) # lb/s
C_p =  (10.8660+11.6259+12.9042)/3 * 0.238845896632776 # from textbook: 0.652 # Btu/(lb * deg F)
T_INLET_IMPERIAL = T_INLET * 9/5 # Rankine
T_INLET_STATOR_IMPERIAL = 876.05 * 9/5 # inlet temp before stator
PRESSURE_RATIO = 1 # Pin/Pout
P_OUT_STATIC = 2.0684e+05
P_OUT_TOTAL = P_OUT_STATIC * (1 + (GAMMA - 1) / 2 * STATOR_EXIT_MACH**2) ** (GAMMA / (GAMMA - 1))
TOTAL_PRESSURE_RATIO = 2.4132e+06 / P_OUT_STATIC # include stator
R_HUB = 0.0458 # hub radius m
R_TIP = 0.0508 # tip radius m
# P_CHAMBER = P_INLET / ((1 + (GAMMA - 1) / 2 * 2.2402**2)**(-GAMMA / (GAMMA-1))
P_CHAMBER = 2.41297e6
RH_RT = R_HUB / R_TIP
CHORD = 0.01 # m
NUM_BLADES = 35
ALPHA_IN = float(var_data['a_in_deg'].item()) # deg
ALPHA_OUT = float(var_data['a_out_deg'].item()) # deg 
BETA = float(var_data['b_deg'].item()) # deg
BLADE_HEIGHT = R_TIP - R_HUB # m
MEAN_RADIUS = (R_HUB + R_TIP) / 2 # m
PITCH = 2 * np.pi * MEAN_RADIUS / NUM_BLADES
W_REL_IN = float(var_data['w'].item()) # m/s
RHO_INLET = 0.418683 # kg/m^3
REYNOLDS_NUM = RHO_INLET * W_REL_IN * CHORD / DYNAMIC_VISCOSITY 
Q_FREESTREAM = 0.5 * RHO_INLET * W_REL_IN**2;

# print(INLET_MACH, HORSE_POWER, M_DOT, M_DOT_IMPERIAL, ALPHA_IN, ALPHA_OUT, BETA, W_REL_IN, sep='\n')

"""
thp   - Turbine horsepower (hp)
m_dot - Mass flow rate (lbm/s)
Cp    - Specific heat at constant pressure (Btu/lbm-R)
T0    - Inlet total temperature (R)
Rt    - Total pressure ratio (P_inlet / P_exit)
gamma - Specific heat ratio

Re - Reynolds number
MinH - relative Mach number at the hub at rotor inlet (can use M_in - ~5-10%)
Min - relative inlet Mach number 
Mout - relative outlet Mach number  
rH_rT - hub to tip radius ratio 
Pin_Pout - inlet to exit pressure ratio (1 b/c pure impulse)
c - chord
H - blade height
Lx - axial chord
alpha_in - inlet flow angle (degrees)
alpha_out - outlet flow angle (degrees)
alpha_m - mean flow angle (degrees)
beta - relative flow angle (degrees)
Kp, Kp_xi - empirical loss coefficients
Delta_Te - trailing edge temperature ratio
gamma - specific heat ratio
pitch_to_chord - ratio of pitch to chord
pitch - pitch of blades (2pir/N)
te_pitch - pitch of trailing edge (equal to pitch for our turbine)
throat - area of throat between blades (passage) (pitch - t_max)
"""

"""
EXPORT FUNCTIONS
"""

def export_parameters_to_excel(file_name, configurations, parameters):
    """
    Exports parameters to an Excel file in the required format.

    Parameters:
        file_name (str): The name of the Excel file to create.
        configurations (list of str): List of configuration names.
        parameters (dict): Dictionary where keys are parameter names and values are lists of parameter values.
    """
    # Prepare the data for the Excel sheet
    param_names = list(parameters.keys())  # Parameter names
    param_values = pd.DataFrame(parameters)  # Parameter values as a DataFrame

    # Insert configurations as the first column
    param_values.insert(0, 'Configurations', configurations)

    # Create an empty DataFrame for the top two blank rows
    blank_rows = pd.DataFrame([[""] * (len(param_names) + 1)] * 2)

    # Combine blank rows, parameter names (header), and parameter data
    header_row = [[""] + param_names]  # Parameter names start from B2
    full_data = pd.concat(
        [blank_rows, pd.DataFrame(header_row), param_values],
        ignore_index=True
    )

    # Save to an Excel file
    with pd.ExcelWriter(file_name, engine='openpyxl') as writer:
        full_data.to_excel(writer, index=False, header=False)

"""
ANALYSIS FUNCTIONS
"""

# Isentropic area-Mach relation function
def area_mach_relation(M, gamma):
    return (1 / M) * ((2 / (gamma + 1)) * (1 + (gamma - 1) / 2 * M**2))**((gamma + 1) / (2 * (gamma - 1)))

# Solve for Mach from area ratio using isentropic relations
def isentropic_mach_finder(area_ratio, gamma, M_inlet, M_outlet):
    func = lambda M: area_ratio - area_mach_relation(M, gamma)
    M_guess = (M_inlet + M_outlet) / 2
    Mach_solution, = fsolve(func, M_guess)
    return Mach_solution

# Calculate Mach and Pressure distributions
def calculate_mach_pressure_distribution(area_vec, gamma, R, T_inlet, P_inlet, M_inlet, M_outlet, a_star=None):
    Mach_vec = np.zeros_like(area_vec)
    P_vec = np.zeros_like(area_vec)

    Tt_inlet = T_inlet * (1 + ((gamma - 1) / 2) * M_inlet**2)
    Pt_inlet = P_inlet * (1 + ((gamma - 1) / 2) * M_inlet**2)**(gamma / (gamma - 1))

    if a_star is None:
        a_star = area_vec[0]

    for i, area in enumerate(area_vec):
        area_ratio = area / a_star
        Mach = isentropic_mach_finder(area_ratio, gamma, M_inlet, M_outlet)
        Mach_vec[i] = Mach
        T_static = Tt_inlet / (1 + ((gamma - 1) / 2) * Mach**2)
        P_vec[i] = Pt_inlet * (T_static / Tt_inlet)**(gamma / (gamma - 1))
    
    return Mach_vec, P_vec


# Plot Mach and Pressure distributions
def plot_mach_pressure_distributions(Mach_vec, P_vec):
    x = np.arange(1, len(Mach_vec) + 1)

    fig, axs = plt.subplots(2, 1, figsize=(10, 8))
    fig.suptitle("Analysis along Lower Surface")

    axs[0].plot(x, Mach_vec, 'b-')
    axs[0].grid(True)
    axs[0].set_xlabel("Lower Surface Coordinate")
    axs[0].set_ylabel("Mach Number")
    axs[0].set_title("Mach Number vs Lower (Pressure) Surface Coordinate")

    axs[1].plot(x, P_vec, 'g-')
    axs[1].grid(True)
    axs[1].set_xlabel("Lower Surface Coordinate")
    axs[1].set_ylabel("Static Pressure [Pa]")
    axs[1].set_title("Pressure Distribution along Pressure Surface")

    plt.tight_layout(rect=[0, 0, 1, 0.95])
    plt.show()

# Kantrowitz limit
def kantrowitz_limit(M0, gamma):
    term1 = M0 * ((gamma + 1) / (2 + (gamma - 1) * M0**2))**((gamma + 1) / (2 * (gamma - 1)))
    term2 = (((gamma + 1) * M0**2) / ((gamma - 1) * M0**2 + 2))**(-gamma / (gamma - 1))
    term3 = ((gamma + 1) / (2 * gamma * M0**2 - (gamma - 1)))**(-1 / (gamma - 1))
    return term1 * term2 * term3

def estimate_turbine_efficiency(thp, m_dot, Cp, T0, Rt, gamma):
    """
    TURBINE_EFFICIENCY Calculates the overall turbine efficiency

    Inputs:
        thp   - Turbine horsepower (hp)
        m_dot - Mass flow rate (lbm/s)
        Cp    - Specific heat at constant pressure (Btu/lbm-R)
        T0    - Inlet total temperature (R)
        Rt    - Total pressure ratio (P_inlet / P_exit)
        gamma - Specific heat ratio

    Output:
        eta_t - Turbine efficiency (unitless)
    """
    # Convert horsepower to ft-lbf/s: 1 hp = 550 ft-lbf/s
    thp_ftlbfs = thp * 550

    # Compute efficiency
    eta_t = thp_ftlbfs / (m_dot * Cp * T0 * (1 - (1 / Rt) ** ((gamma - 1) / gamma)) * 778)
    return eta_t

import math

def calc_relative_mach_number(V_abs, U, alpha, a):
    """
    V_abs: absolute velocity at inlet for M_rel_in or outlet for M_rel_out
    U: blade speed at rotor inlet radius for M_rel_in or outlet radius
    alpha: absolute flow angle for inlet (radians)
    a: speed of sound

    Returns:
        M_rel - relative Mach number
    """
    M_rel = np.sqrt(V_abs ** 2 + U ** 2 - 2 * V_abs * U * np.cos(alpha)) / a
    return M_rel

def kacker_okapuu_loss(Re, MinH, Min, Mout, rH_rT, Pin_Pout, Ptotal_rel,
                       c, H, Lx, alpha_in, alpha_out, alpha_m, beta, 
                       gamma, pitch_to_chord, pitch, te_pitch, throat, data):
    """
    Inputs:
        Re - Reynolds number
        MinH - relative Mach number at the hub at rotor inlet (can use M_in - ~5-10%)
        Min - relative inlet Mach number 
        Mout - relative outlet Mach number  
        rH_rT - hub to tip radius ratio 
        Pin_Pout - inlet to exit pressure ratio (1 b/c pure impulse)
        Ptotal_rel - total pressure in rotor relative frame
        c - chord
        H - blade height
        Lx - axial chord
        alpha_in - inlet flow angle (degrees)
        alpha_out - outlet flow angle (degrees)
        alpha_m - mean flow angle (degrees)
        beta - relative flow angle (degrees)
        Kp, Kp_xi - empirical loss coefficients
        Delta_Te - trailing edge temperature ratio
        gamma - specific heat ratio
        pitch_to_chord - ratio of pitch to chord
        pitch - pitch of blades (2pir/N)
        te_pitch - pitch of trailing edge (equal to pitch for our turbine)
        throat - area of throat between blades (passage) (pitch - t_max)

    Returns:
        Y_total - total loss coefficient
    """

    # Reynolds number correction factor
    if Re <= 2e5:
        xRe = (Re / 2e5) ** (-0.4)
    elif Re <= 1e6:
        xRe = 1.0
    else:
        xRe = (Re / 1e6) ** (-0.2)

    # Shock loss
    # term1 = (1 + (gamma - 1) / 2 * Min) ** (gamma / (gamma - 1))
    # term2 = (1 + (gamma - 1) / 2 * Mout) ** (gamma / (gamma - 1))
    # Y_shock = 0.75 * (MinH - 0.4) ** 1.75 * rH_rT * Pin_Pout * (term1 - term2)

    # alternative calcuation:
    dP_q1_hub = 0.75*(MinH-0.4)**1.75 # Eqn 4, this is at the hub
    dP_q1_shock = rH_rT * dP_q1_hub # Eqn 5
    Y_shock = dP_q1_shock * Pin_Pout * (1-(1+(gamma-1)/2*Min**2))/(1-(1+(gamma-1)/2*Mout**2)) # Eqn 6

    delta_p = (
            ((gamma + 1) * MinH**2) / ((gamma - 1) * MinH**2 + 2)
        ) ** (gamma / (gamma - 1)) * (
            (gamma + 1) / (2 * gamma * MinH**2 - (gamma - 1))
        ) ** (1 / (gamma - 1))
    dP_q1_shock = rH_rT * delta_p # Eqn 5
    if Min > 1.2:
        Y_shock = ((1 - dP_q1_shock) * Ptotal_rel) / Q_FREESTREAM
    
    print(f"Shock Loss: {Y_shock}")

    # Profile loss
    Yp_beta0 = data['Fig01_beta0'](pitch_to_chord, alpha_out)
    Yp_beta1_alpha2 = data['Fig02'](pitch_to_chord, alpha_out)
    t_max_c = data['Fig04'](np.abs(beta)+np.abs(alpha_out))
    
    Yp_amdc = (Yp_beta0 + np.abs(beta/alpha_out) *beta/alpha_out * (Yp_beta1_alpha2-Yp_beta0)) * ((t_max_c)/0.2)**(beta/alpha_out) # Eqn 2, AMDC = Ainley Mathieson Dunham Came
    
    YP_i0 = Yp_amdc
    
    estimate_with_mach = min(0.8, Mout)
    K1 = data['Fig08_K1'](estimate_with_mach)
    K2 = (Min/Mout)**2 
    Kp = 1-K2*(1-K1)
    Kp_xi = 1
    
    CFM = 1+60*(Mout-1)**2    # Eqn 9 
    
    YP = 0.914 * (2/3 * Kp * Kp_xi * YP_i0 + Y_shock)
    if Mout>1 and Mout < 1.4:
        YP = YP*CFM

    print(f"Profile Loss 0 Incedence Impulse Correlation: {YP_i0}")
    print(f"Profile Loss Corrected: {YP}")

    # Aspect ratio correction factor (assume χ_AR = 1)
    chi_AR = 1

    # Secondary loss
    alpha_in_rad = math.radians(alpha_in)
    alpha_out_rad = math.radians(alpha_out)
    alpha_m_rad = math.radians(alpha_m)
    term_angle = math.tan(alpha_in_rad) - math.tan(alpha_out_rad)
    beta_rad = math.radians(beta)

    # YS = (0.04 * (c / H) * chi_AR * (4 * term_angle ** 2) *
    #       (math.cos(alpha_out_rad) ** 2 / math.cos(alpha_m_rad)) *
    #       (math.cos(alpha_out_rad) / math.cos(alpha_in_rad)) *
    #       (1 - (Lx / H) ** 2 * (1 - Kp)))
    # alternative calculation:
    f_ar = (1-0.25*np.sqrt(2-H/c)) / (H/c) if H/c<=2 else 1/(H/c)
    alpham = np.arctan(0.5*(np.tan(alpha_in_rad) - np.tan(alpha_out_rad)))
    Cl_sc = 2*(np.tan(alpha_in_rad)+np.tan(alpha_out_rad))*np.cos(alpham)
    Ys_amdc = 0.0334 *f_ar *np.cos(alpha_out_rad)/np.cos(beta_rad) * (Cl_sc)**2 * np.cos(alpha_out_rad)**2 / np.cos(alpham)**3
    K3 = 1/(H/(Lx))**2     # Fig 13
    Ks = 1-K3*(1-Kp)        # Eqn 15
    YS = 1.2*Ys_amdc*Ks     # Eqn 16

    print(f"Secondary Flow Loss: {YS}")

    # Trailing edge loss (due to mixing of wakes and recovery pressure drive)
    # M2out = math.sqrt(Mout)
    # t1 = 1 + (gamma - 1) / 2 * M2out ** 2 * (1 - 1 / (1 - Delta_Te))
    # t2 = 1 + (gamma - 1) / 2 * M2out ** 2
    # YTe = (t1 ** (-gamma / (gamma - 1)) - 1) / (1 - t2 ** (-gamma / (gamma - 1)))
    # alternative calculation
    """
    non-dimensional trailing edge blockage parameter — how much of the throat is blocked by the trailing edge relative to the pitch. The higher this value:
        The thicker the trailing edge,
        The stronger the wake and mixing,
        The larger the trailing edge loss.
    """
    if np.abs(alpha_in-alpha_out)<5:
        delta_phi2 = data['Fig14_Impulse'](te_pitch*pitch / throat)
        '''
        pitch = 2pir/N, same radius so te_pitch = pitch, row.throat =pitch-t_max
        '''
    else:
        delta_phi2 = data['Fig14_Axial_Entry'](te_pitch*pitch / throat)
    
    YTe = (1-(gamma-1)/2 - Mout**2 * (1/(1-delta_phi2)-1))**(-gamma/(gamma-1))-1
    YTe = YTe / (1-(1+(gamma-1)/2*Mout**2)**(-gamma/(gamma-1)))


    if Mout > 1.2: 
        YTe = 0; # paper says ignore for supersonic
        print(f"Trailing Edge Losses Ignored (Accounted for in Profile Loss)")
    else:
        print(f"Trailing Edge Losses: {YTe}")
    

    # Total loss coefficient
    Y_total = xRe * YP + YS + YTe
    Y_total = xRe * YP + YTe

    print(f"Total Loss Coefficient: {Y_total}")

    return Y_total

"""
TURBINE BLADE GEOMETRY GENERATION
"""

# def solve_for_x1(c1, c2, r, y2, x2):
#     """
#     solves the equation (x+b)(r^2-(x+b)^2) = (y_2 - sqrt(r^2-(x+b)^2)) / (x_2 - x) for x.
#     """

#     print(c1, c2, r, y2, x2)

#     initial_guess = - c1 - r
#     def equation(x):
#         left = - (x + c1) / np.sqrt(r**2 - (x + c1)**2)
#         right = (y2 - (c2 + np.sqrt(r**2 - (x + c1)**2))) / (x2 - x)
#         return left - right

#     solution = fsolve(equation, initial_guess)
    
#     y1 = np.sqrt(r**2 - (solution[0] + c1)**2) + c2

#     return (solution[0], y1)

def solve_for_x1(c1, c2, r, y2, x2):
    """
    Solves the equation (x+b)(r^2-(x+b)^2) = (y_2 - sqrt(r^2-(x+b)^2)) / (x_2 - x) for x,
    using minimization for robustness.
    """
    # print(c1, c2, r, y2, x2)

    # Define the equation as the squared difference to minimize
    def equation(x):
        if np.abs(x + c1) >= r:
            return np.inf  # Ensure we don't evaluate outside the circle
        
        # Compute the terms
        left = - (x + c1) / np.sqrt(r**2 - (x + c1)**2)
        right = (y2 - (c2 + np.sqrt(r**2 - (x + c1)**2))) / (x2 - x)
        
        # Return squared difference
        return (left - right)**2

    # Set the initial guess and search range
    initial_guess = -c1 - r
    bounds = (-c1 - r, -c1)

    # Minimize the equation
    solution = fmin(equation, initial_guess, disp=False)

    # Calculate y1 based on the solution for x1
    x1 = solution[0]
    y1 = np.sqrt(r**2 - (x1 + c1)**2) + c2

    return (x1, y1)

def bezier_curve(control_points, t):
    n = len(control_points) - 1
    curve = np.zeros((len(t), 2))
    for i, (x, y) in enumerate(control_points):
        binomial_coeff = np.math.comb(n, i)
        curve[:, 0] += binomial_coeff * (1 - t) ** (n - i) * t ** i * x
        curve[:, 1] += binomial_coeff * (1 - t) ** (n - i) * t ** i * y
    return curve

def curvature(x, y):
    dx = np.gradient(x)
    dy = np.gradient(y)
    ddx = np.gradient(dx)
    ddy = np.gradient(dy)
    return (dx * ddy - dy * ddx) / (dx**2 + dy**2)**(3/2)

def compute_distances(x_lower, y_lower, upper_curve, target_distance, blade_spacing):
    y_lower = y_lower + blade_spacing;
    num_points = len(x_lower)
    # sampled_indices = random.sample(range(num_points // 2), min(20, num_points))
    sampled_indices = np.linspace(0, 2 * num_points//3, min(30, num_points), dtype=int)
    sse = 0

    dx = np.gradient(x_lower)
    dy = np.gradient(y_lower)

    for i in sampled_indices:
        x0, y0 = x_lower[i], y_lower[i]
        normal_dx = dy[i]
        normal_dy = -dx[i]
        normal_length = np.sqrt(normal_dx**2 + normal_dy**2)
        normal_dx /= normal_length
        normal_dy /= normal_length
        # normal_slope = normal_dy / normal_dx;
        normal_angle = np.angle(normal_dx + normal_dy*1j)

        min_distance = float('inf')
        min_slope_diff = float('inf')
        for bx, by in upper_curve:
            # slope = (by - y0) / (bx - x0)
            slope_angle = np.angle((bx - x0) + (by - y0) * 1j)
            # slope_diff = abs(slope_angle - normal_angle)
            slope_diff = min(abs(slope_angle - normal_angle), 2 * np.pi - abs(slope_angle - normal_angle))
            distance = abs((bx - x0) * normal_dx + (by - y0) * normal_dy)
            if slope_diff < min_slope_diff:
                min_distance = distance
                min_slope_diff = slope_diff

        sse += (min_distance - target_distance) ** 2

    return sse

def compute_arc(chord, beta, num_points):
    x_lower = np.linspace(0, chord, num_points)
    half_d = chord / 2
    lower_slope = np.tan(np.radians(beta))
    r = np.sqrt(half_d**2 + (half_d / lower_slope)**2)
    y_lower = np.sqrt(r**2 - (x_lower - half_d)**2) - np.sqrt(r**2 - half_d**2)
    y_prime = -(-half_d) / np.sqrt(r**2 - (-half_d)**2)
    return x_lower, y_lower, r, y_prime

def compute_upper_surface(params, target_distance, num_points, blade_spacing, radius, chord):
    x2, x3 = params
    c_d = chord / 2
    dist = c_d - x3

    r2 = radius - target_distance

    global upper_radius
    upper_radius = r2

    # circular section
    x_circ = np.linspace(x3, c_d + dist, num_points // 3)
    y_circ = np.sqrt(r2**2 - (x_circ - c_d)**2) - np.sqrt(radius**2 - c_d**2) + blade_spacing


    # bezier section & rounds
    # x_bez1 = np.linspace(x1, x3, num_points // 3)

    # x1 = 0
    # y1 = 0
    y3 = y_circ[0]
    slope = (-(c_d - dist) + c_d) / np.sqrt(-((c_d - dist)**2) + 2*(c_d - dist)*c_d + r2**2 - (c_d**2))
    y2 = slope * (x2 - c_d + dist) + y3

    global bez_y2, bez_x3
    bez_y2 = y2
    bez_x3 = x3

    ## rounds
    c1 = np.sqrt((y_prime_l * rounds_rad)**2 / (1 + y_prime_l ** 2))
    c2 = np.sqrt(rounds_rad**2 - c1**2)
    x1, y1 = solve_for_x1(c1, c2, rounds_rad, y2, x2)
    x_rounds1 = np.linspace(0, -c1 - rounds_rad, 50)
    x_rounds2 = np.linspace(-c1 - rounds_rad, x1, 50)
    y_rounds1 = c2 - np.sqrt(rounds_rad**2 - (x_rounds1 + c1)**2)
    y_rounds2 = c2 + np.sqrt(rounds_rad**2 - (x_rounds2 + c1)**2)

    control_points = [
        (x1, y1),
        (x2, y2),
        (x3, y3),
    ]
    t = np.linspace(0, 1, num_points // 3)
    bezier1 = bezier_curve(control_points, t)

    control_points2 = [
        (chord - x3, y3),
        (chord - x2, y2),
        (chord - x1, y1),
    ]
    bezier2 = bezier_curve(control_points2, t)

    x_upper = np.concatenate([x_rounds1, x_rounds2, bezier1[:, 0], x_circ, bezier2[:, 0]])
    y_upper = np.concatenate([y_rounds1, y_rounds2, bezier1[:, 1], y_circ, bezier2[:, 1]])
    
    upper = np.column_stack((x_upper, y_upper))
    return upper

def objective_fn(params, x_lower, y_lower, target_distance, num_points, blade_spacing, radius, chord):
    # x1, y1, x2, y2 = params  # Only optimize up to x = chord/2
    # x_mid, y_mid = x_lower[num_points // 2], y_lower[num_points // 2]
    # control_points = [
    #     (x_lower[0], y_lower[0]),
    #     (x1, y1),
    #     (x2, y2),
    #     (x_mid, y2),
    # ]
    # t = np.linspace(0, 1, num_points // 2)
    # bezier = bezier_curve(control_points, t)

    # bezier_curvature = curvature(bezier[:, 0], bezier[:, 1])
    # curvature_variance = np.var(bezier_curvature)

    upper_surface = compute_upper_surface(params, target_distance, num_points, blade_spacing, radius, chord)
    # len_upper_surf = upper_surface.shape[0]

    upper_surface_filtered = upper_surface[upper_surface[:, 0] > 0]
    len_upper_surf = upper_surface_filtered.shape[0]

    # Filter points where x > 0 for the lower surface
    x_lower_filtered = x_lower[:num_points // 2][x_lower[:num_points // 2] > 0]
    y_lower_filtered = y_lower[:num_points // 2][x_lower[:num_points // 2] > 0]

    # Compute the distance error
    distance_error = compute_distances(
        x_lower_filtered,
        y_lower_filtered,
        upper_surface_filtered[:len_upper_surf // 2, :],
        target_distance,
        blade_spacing
    )

    # distance_error = compute_distances(x_lower[:num_points // 2], y_lower[:num_points // 2], upper_surface[:len_upper_surf // 2, :], target_distance, blade_spacing)

    return distance_error # + curvature_variance

def optimize_upper_surface(x_lower, y_lower, target_distance, blade_spacing, radius, chord):
    num_points = len(x_lower)
    x_mid, y_mid = x_lower[num_points // 2], y_lower[num_points // 2]
    initial_guess = [1/3 * x_mid, 2/3 * x_mid]
    bounds = [
        (x_lower[0], x_mid), 
        (x_lower[0], x_mid),
    ]

    result = minimize(
        objective_fn,
        initial_guess,
        args=(x_lower, y_lower, target_distance, num_points, blade_spacing, radius, chord),
        bounds=bounds,
        method='SLSQP'
    )

    return result.x

def compute_distances_and_plot(x_lower, y_lower, upper_curve, target_distance, blade_spacing, ax):
    y_lower = y_lower + blade_spacing;
    num_points = len(x_lower)
    # sampled_indices = random.sample(range(num_points // 2), min(20, num_points))
    sampled_indices = np.linspace(0, 2 * num_points // 3, min(40, num_points), dtype=int)
    sse = 0

    dx = np.gradient(x_lower)
    dy = np.gradient(y_lower)

    distances = []

    for i in sampled_indices:
        x0, y0 = x_lower[i], y_lower[i]
        normal_dx = dy[i]
        normal_dy = -dx[i]
        normal_length = np.sqrt(normal_dx**2 + normal_dy**2)
        normal_dx /= normal_length
        normal_dy /= normal_length
        # normal_slope = normal_dy / normal_dx;
        normal_angle = np.angle(normal_dx + normal_dy*1j)

        min_distance = float('inf')
        min_slope_diff = float('inf')
        closest_point = None
        for bx, by in upper_curve:
            # slope = (by - y0) / (bx - x0)
            slope_angle = np.angle((bx - x0) + (by - y0) * 1j)
            # slope_diff = abs(slope_angle - normal_angle)
            slope_diff = min(abs(slope_angle - normal_angle), 2 * np.pi - abs(slope_angle - normal_angle))
            distance = abs((bx - x0) * normal_dx + (by - y0) * normal_dy)
            if slope_diff < min_slope_diff:
                min_distance = distance
                min_slope_diff = slope_diff
                closest_point = (bx, by)

        if closest_point:
            bx, by = closest_point
            ax.plot([x0, bx], [y0, by], 'g--', linewidth=0.5)

        distances.append(min_distance)
        sse += (min_distance - target_distance) ** 2


    A_STAR = distances[0] / area_mach_relation(INLET_MACH, GAMMA)
    print(f"A* = {A_STAR}")
    distances = np.concatenate([distances, np.flip(distances)])
    Mach_vec, P_vec = calculate_mach_pressure_distribution(distances, GAMMA, R_INLET, T_INLET, P_INLET, INLET_MACH, OUTLET_MACH, a_star=A_STAR)

    # plot_mach_pressure_distributions(Mach_vec, P_vec)

    # Check Kantrowitz limit
    kant_limit = kantrowitz_limit(INLET_MACH, GAMMA)
    max_area_ratio = np.max(distances) / A_STAR
    true_area_ratio = np.max(distances) / np.min(distances)
    kantrowitz_violation = max_area_ratio > kant_limit
    
    print(f"Kantrowitz Limit Violated: {kantrowitz_violation}, LIMIT = {kant_limit}, AREA_RATIO  = {max_area_ratio}, TRUE_AREA_RAATIO = {true_area_ratio}")

    # Efficiecy calculations
    blade_max_thickness = target_distance / 100 # convert to m
    
    n_th = estimate_turbine_efficiency(HORSE_POWER, M_DOT_IMPERIAL, C_p, T_INLET_STATOR_IMPERIAL, TOTAL_PRESSURE_RATIO, GAMMA)
    print(f"Efficiency from Isentropic Analysis: {n_th}")

    global data

    P_total_rotor_rel = P_vec[0] + Q_FREESTREAM

    y_loss = kacker_okapuu_loss(
        Re=REYNOLDS_NUM,
        MinH=(INLET_MACH * 0.95),
        Min=INLET_MACH,
        Mout=OUTLET_MACH,
        rH_rT=RH_RT,
        Pin_Pout=PRESSURE_RATIO,
        Ptotal_rel=P_total_rotor_rel,
        c=CHORD,
        H=BLADE_HEIGHT,
        Lx=CHORD,
        alpha_in=ALPHA_IN,
        alpha_out=ALPHA_OUT,
        alpha_m=((ALPHA_IN+ALPHA_OUT)/2),
        beta=BETA,
        gamma=GAMMA,
        pitch_to_chord=(PITCH/CHORD),
        pitch=PITCH,
        te_pitch=PITCH,
        throat=(PITCH - blade_max_thickness),
        data=data
    )
    print(f"Overall Losses: {y_loss}")

    n_total = n_th * (1 - y_loss)

    print(f"Overall Efficiency: {n_total}")
    print(f"Expected (After loss) Horsepower of Design: {HORSE_POWER * (1 - y_loss)} HP")
    print(f"Isentropic Efficiency: {(1 - y_loss)}")

    return sse


def update_plot(val):
    chord = chord_slider.val
    beta = beta_slider.val
    target_distance_percent = target_slider.val
    blade_spacing = blade_spacing_slider.val
    

    global y_prime_l, rounds_rad
    rounds_rad = round_slider.val

    x_lower, y_lower, radius, y_prime_l = compute_arc(chord, beta, num_points)

    params = optimize_upper_surface(x_lower, y_lower, target_distance_percent * blade_spacing, blade_spacing, radius, chord)

    full_upper = compute_upper_surface(params, target_distance_percent * blade_spacing, num_points, blade_spacing, radius, chord)

    lower_surface.set_data(x_lower, y_lower)
    upper_surface.set_data(full_upper[:, 0], full_upper[:, 1])
    lower_surface_guide.set_data(x_lower, y_lower + blade_spacing)

    for line in ax.lines[3:]:
        line.remove()

    upper_surface_filtered = full_upper[full_upper[:, 0] > 0]
    len_upper_surf = upper_surface_filtered.shape[0]

    # Filter points where x > 0 for the lower surface
    x_lower_filtered = x_lower[:num_points // 2][x_lower[:num_points // 2] > 0]
    y_lower_filtered = y_lower[:num_points // 2][x_lower[:num_points // 2] > 0]

    compute_distances_and_plot(x_lower[:num_points // 2], y_lower[:num_points // 2], upper_surface_filtered[:len_upper_surf // 2, :], target_distance_percent * blade_spacing, blade_spacing, ax)
    # compute_distances_and_plot(x_lower_filtered, y_lower_filtered, upper_surface_filtered[:len_upper_surf // 2, :], target_distance_percent * blade_spacing, blade_spacing, ax)

    fig.canvas.draw_idle()

    print(f"{radius}\n{rounds_rad}\n{params[0]}\n{bez_y2}\n{upper_radius}\n{blade_spacing}\n{beta}\n{chord}")

    configurations = ["TurbineConfig"]
    parameters = {
        "lower_rad": [radius],
        "fillet_rad": [rounds_rad],
        "bez_x": [params[0]],
        "bez_y": [bez_y2],
        "bez_x3": [bez_x3],
        "upper_rad": [upper_radius],
        "b_spacing": [blade_spacing],
        "l_angle": [beta],
        "chord": [chord],
    }
    # export_parameters_to_excel("design_table.xlsx", configurations, parameters)

def plot_geometry_with_sliders():
    global fig, ax, chord_slider, beta_slider, target_slider, blade_spacing_slider, round_slider, lower_surface, upper_surface, lower_surface_guide

    # initial_chord = 1.0
    # initial_beta = 30.0
    # initial_target_distance = 0.1
    # initial_blade_spacing = 1.0
    # initial_rounds = 0.005

    # initial_chord = 1.0
    # initial_beta = 31.8901619
    # initial_target_distance = 0.002776433 / 0.008
    # initial_blade_spacing = 0.8
    # initial_rounds = 0.0005

    """
    INITIAL PARAMTERS FROM VELOCITY TRIANGLES
    """

    # in cm
    initial_chord = 1.0
    initial_beta = BETA
    initial_target_distance = 0.0028 / 0.00822 # a percentage of blade spacing
    initial_blade_spacing = 0.8222
    initial_rounds = 0.01785 #0.0005

    global y_prime_l, rounds_rad
    rounds_rad = initial_rounds

    x_lower, y_lower, radius, y_prime_l = compute_arc(initial_chord, initial_beta, num_points)
    params = optimize_upper_surface(x_lower, y_lower, initial_target_distance * initial_blade_spacing, initial_blade_spacing, radius, initial_chord)

    full_upper = compute_upper_surface(params, initial_target_distance * initial_blade_spacing, num_points, initial_blade_spacing, radius, initial_chord)

    fig, ax = plt.subplots()
    plt.subplots_adjust(left=0.1, bottom=0.3)

    lower_surface, = ax.plot(x_lower, y_lower, 'b-', label='Lower Surface')
    upper_surface, = ax.plot(full_upper[:, 0], full_upper[:, 1], 'r-', label='Upper Surface')
    lower_surface_guide, = ax.plot(x_lower, y_lower + initial_blade_spacing, 'b-', label='Lower Surface Guide')

    upper_surface_filtered = full_upper[full_upper[:, 0] > 0]
    len_upper_surf = upper_surface_filtered.shape[0]

    compute_distances_and_plot(x_lower[:num_points // 2], y_lower[:num_points // 2], upper_surface_filtered[:len_upper_surf // 2, :], initial_target_distance * initial_blade_spacing, initial_blade_spacing, ax)


    ax.axis('equal')
    ax.grid(True)
    ax.set_xlabel('X')
    ax.set_ylabel('Y')
    ax.legend()

    ax_chord = plt.axes([0.1, 0.21, 0.8, 0.03])
    ax_beta = plt.axes([0.1, 0.16, 0.8, 0.03])
    ax_target = plt.axes([0.1, 0.11, 0.8, 0.03])
    ax_bspacing = plt.axes([0.1, 0.06, 0.8, 0.03])
    ax_rounds = plt.axes([0.1, 0.01, 0.8, 0.03])

    chord_slider = Slider(ax_chord, 'Chord', 0.5, 2.0, valinit=initial_chord)
    beta_slider = Slider(ax_beta, 'Beta', 0, 75, valinit=initial_beta)
    target_slider = Slider(ax_target, 'Target Area', 0.0, 1.0, valinit=initial_target_distance)
    blade_spacing_slider = Slider(ax_bspacing, 'Blade Spacing', 0.2, 3.0, valinit=initial_blade_spacing)
    round_slider = Slider(ax_rounds, 'Rounds', 0.0, 1.0, valinit=initial_rounds)

    chord_slider.on_changed(update_plot)
    beta_slider.on_changed(update_plot)
    target_slider.on_changed(update_plot)
    blade_spacing_slider.on_changed(update_plot)
    round_slider.on_changed(update_plot)

    plt.show()


def main():
    path = pathlib.Path(os.getcwd() + "/KackerOkapuu/kackerokapuu.pkl")
    global data
    data = None
    if not path.exists():
        print("Couldn't find data file. Downloading ... ")
        # url = "https://github.com/nasa/turbo-design/raw/main/references/Turbines/KackerOkapuu/kackerokapuu.pkl"
        # response = requests.get(url, stream=True)
        # with open(path.absolute(), mode="wb") as file:
        #     for chunk in response.iter_content(chunk_size=10 * 1024):
        #         file.write(chunk)
        
    with open(path.absolute(),'rb') as f:
        data = pickle.load(f) # type: ignore
        print(data)
    plot_geometry_with_sliders()

if __name__ == '__main__':
    main()