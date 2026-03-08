import numpy as np
from scipy.optimize import fsolve
import matplotlib.pyplot as plt
import numpy as np
import sympy as sp
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from scipy.optimize import fsolve, minimize
import math
import os
import pathlib
import requests
import pickle

def lerp(x1, y1, x2, y2, x):
    return (y2 - y1) / (x2 - x1) * (x-x1) + y1

def calc_blade_speed(radius, rpm):
    """Calculate blade speed given radius and RPM."""
    omega = rpm * 2 * np.pi / 60  # Convert RPM to angular velocity (rad/s)
    return np.array([radius * omega, 0])  # Blade speed (m/s)

def calc_zweifel_loading_coefficient(solidity, beta1, beta2, w1, w2):
    """Calculate Zweifel loading coefficient."""
    print(f"Solidity: {solidity}, beta1: {beta1}, beta2: {beta2}, w1: {w1}, w2: {w2}")
    phi_z = 2 * solidity * np.sin(beta2) * np.cos(beta2) * np.abs((w1/w2) - 1)
    return phi_z # target between 0.8 and 1.0

def prendtl_meyer(M, gamma):
    the_frac = np.sqrt((gamma + 1) / (gamma - 1))
    god_knows = np.atan(np.sqrt(gamma * (M**2 - 1) / (2 + gamma * (M**2 - 1))))
    literally_just_ref = np.atan(np.sqrt(M**2 - 1))
    nu = the_frac * god_knows - literally_just_ref
    return nu

def rotor_back_calculate(RPM, torque, mass_flow, beta, radius, V_in):
    """Back-calculate rotor parameters given operating conditions. Specifically solve for the velocity triangles."""
    U = calc_blade_speed(radius, RPM)[0]
    C = torque / (radius * mass_flow)

    v2, a1, a2, w, b = sp.symbols("v2 a1 a2 w b")

    """
    (i  ) v1 * sin(a1) == w * sin(b) + u
    (ii ) v1 * cos(a1) == w * cos(b)
    (iii) v2 * sin(a2) == w * sin(b) - u
    (iv ) v2 * cos(a2) == w * cos(b)
    (v  ) v1 * sin(a1) - v2 * sin(a2) == c
    (vi ) v1 * cos(a1) == v2 * cos(a2)
    """

    eq1 = sp.Eq(V_in * sp.cos(a1), w * sp.cos(b))
    eq2 = sp.Eq(V_in * sp.sin(a1), w * sp.sin(b) + U)
    eq3 = sp.Eq(v2 * sp.cos(a2), w * sp.cos(-b))
    eq4 = sp.Eq(v2 * sp.sin(a2), w * sp.sin(-b) + U)
    eq5 = sp.Eq(abs(V_in * sp.sin(a1) - v2 * sp.sin(a2)), C)

    solutions = sp.nsolve(
        [eq1, eq2, eq3, eq4, eq5],
        [v2, a1, a2, w, b],
        [V_in, beta + 10 * np.pi / 180, beta - 10 * np.pi / 180, V_in, beta],
    )

    v1 = V_in
    v2, a1, a2, w, b = map(float, solutions)
    a2 = -a2

    return v1, v2, w, U, a1, a2, b


def plot_velocity_triangles(
    v1, v2, u, w1, w2, chord_length, inclination_angle, beta1, beta2, alpha1, alpha2
):
    """Plot velocity triangles with angles."""
    plt.figure()
    plt.axis("equal")
    plt.grid(True)

    v1_vec = np.array([v1 * np.cos(alpha1), v1 * np.sin(alpha1)])
    v2_vec = np.array([v2 * np.cos(alpha2), v2 * np.sin(alpha2)])
    u1_vec = np.array([0, u])
    u2_vec = np.array([0, u])
    w1_vec = np.array([w1 * np.cos(beta1), w1 * np.sin(beta1)])
    w2_vec = np.array([w2 * np.cos(beta2), w2 * np.sin(beta2)])

    max_mag = max(
        np.linalg.norm(v1_vec),
        np.linalg.norm(v2_vec),
        np.linalg.norm(u1_vec),
        np.linalg.norm(u2_vec),
        np.linalg.norm(w1_vec),
        np.linalg.norm(w2_vec),
    )
    scale_factor = 2 * chord_length / max_mag

    v1_vec *= scale_factor
    v2_vec *= scale_factor
    u1_vec *= scale_factor
    u2_vec *= scale_factor
    w1_vec *= scale_factor
    w2_vec *= scale_factor

    plt.quiver(
        0,
        0,
        w1_vec[0],
        w1_vec[1],
        angles="xy",
        scale_units="xy",
        scale=1,
        color="r",
        label=f"w1 = {np.linalg.norm(w1_vec):.3f}",
    )
    plt.quiver(
        w1_vec[0],
        w1_vec[1],
        u1_vec[0],
        u1_vec[1],
        angles="xy",
        scale_units="xy",
        scale=1,
        color="b",
        label=f"u1 = {np.linalg.norm(u1_vec):.3f}",
    )
    plt.quiver(
        0,
        0,
        v1_vec[0],
        v1_vec[1],
        angles="xy",
        scale_units="xy",
        scale=1,
        color="g",
        label=f"v1 = {np.linalg.norm(v1_vec):.3f}",
    )

    chord_vec = chord_length * np.array(
        [np.cos(inclination_angle), np.sin(inclination_angle)]
    )
    chord_end = w1_vec + chord_vec

    plt.quiver(
        chord_end[0],
        chord_end[1],
        w2_vec[0],
        w2_vec[1],
        angles="xy",
        scale_units="xy",
        scale=1,
        color="r",
        label=f"w2 = {np.linalg.norm(w2_vec):.3f}",
    )
    plt.quiver(
        chord_end[0] + w2_vec[0],
        chord_end[1] + w2_vec[1],
        u2_vec[0],
        u2_vec[1],
        angles="xy",
        scale_units="xy",
        scale=1,
        color="b",
        label=f"u2 = {np.linalg.norm(u2_vec):.3f}",
    )
    plt.quiver(
        chord_end[0],
        chord_end[1],
        v2_vec[0],
        v2_vec[1],
        angles="xy",
        scale_units="xy",
        scale=1,
        color="g",
        label=f"v2 = {np.linalg.norm(v2_vec):.3f}",
    )

    plt.legend()
    plt.xlabel("X")
    plt.ylabel("Y")
    plt.title("Velocity Triangles with Chord Inclination")
    plt.show()

def velocity_triangles_forward_calculate(V_in, alpha1, beta, radius, mass_flow_rate):
    """Forward-calculate rotor parameters given inlet conditions."""
    # Calculate the outlet velocity and angles using the velocity triangle relations
    w1 = V_in * np.cos(alpha1) / np.cos(beta)
    U = V_in * np.sin(alpha1) - w1 * np.sin(beta)

    V_out = np.sqrt((w1 * np.cos(beta))**2 + (U - w1 * np.sin(beta))**2)
    alpha2 = np.arctan2(U - w1 * np.sin(beta), w1 * np.cos(beta))

    torque = (V_in * np.sin(alpha1) - V_out * np.sin(alpha2)) * mass_flow_rate * radius
    
    return V_out, alpha2, w1, U, torque

def velocity_triangles_unsteady_forward_calculate(V_in, alpha1, beta2, radius, mass_flow_rate, rpm, U):
    """Forward-calculate rotor parameters given inlet conditions including RPM (technically there is slip/incedence)"""
    V_theta = V_in * np.sin(alpha1)
    w_t = V_theta - U
    w_axial = V_in * np.cos(alpha1)
    w1 = np.sqrt(w_t**2 + w_axial**2)
    beta1 = np.arctan2(w_t, w_axial)
    V2_axial = w1 * np.cos(beta2)
    V2_t = w1 * np.sin(beta2) - U
    V_out = np.sqrt(V2_axial**2 + V2_t**2)
    alpha2 = np.arctan2(V2_t, V2_axial)
    torque = (V_in * np.sin(alpha1) - V_out * np.sin(alpha2)) * mass_flow_rate * radius
    return V_out, alpha2, w1, U, torque

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
                       gamma, pitch_to_chord, pitch, te_pitch, throat, Q_freestream, data, use_print=True):
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
        Y_shock = ((1 - dP_q1_shock) * Ptotal_rel) / Q_freestream
    
    if use_print:
        print(f"Shock Loss: {Y_shock}")

    # Profile loss
    try:
        Yp_beta0 = data['Fig01_beta0'](pitch_to_chord, alpha_out)
    except Exception as e:
        print("Error accessing Fig01_beta0 data:", e)
        Yp_beta0 = 0.02  # default fallback value
    try:
        Yp_beta1_alpha2 = data['Fig02'](pitch_to_chord, alpha_out)
    except Exception as e:
        print("Error accessing Fig02 data:", e)
        Yp_beta1_alpha2 = 0.02  # default fallback value
    try:
        figure_interp_val = np.abs(np.rad2deg(beta)) + np.abs(np.rad2deg(alpha_out))
        t_max_c = data['Fig04'](figure_interp_val)
    except Exception as e:
        print("Error accessing Fig04 data:", e)
        t_max_c = 0.02  # default fallback value

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

    if use_print:
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
    YS = abs(YS)  # ensure positive

    if use_print:
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
        if use_print:
            print(f"Trailing Edge Losses Ignored (Accounted for in Profile Loss)")
    else:
        if use_print:
            print(f"Trailing Edge Losses: {YTe}")
    

    # Total loss coefficient
    Y_total = xRe * YP + YS + YTe
    Y_total = xRe * YP + YTe

    if use_print:
        print(f"Total Loss Coefficient: {Y_total}")

    return Y_total

def calculate_blade_design_efficiency(distances, 
                                      blade_max_thickness, 
                                      CHORD,
                                      BLADE_HEIGHT,
                                      PITCH,
                                      ALPHA_IN,
                                      ALPHA_OUT,
                                      BETA,
                                      C_p,
                                      HORSE_POWER,
                                      M_DOT_IMPERIAL,
                                      GAMMA,
                                      TOTAL_PRESSURE_RATIO,
                                      Q_FREESTREAM,
                                      REYNOLDS_NUM,
                                      RH_RT,
                                      PRESSURE_RATIO,
                                      R_INLET,
                                      T_INLET,
                                      T_INLET_STATOR_IMPERIAL,
                                      P_INLET,
                                      INLET_MACH,
                                      OUTLET_MACH,
                                      data,
                                      use_print=True):
    A_STAR = distances[0] / area_mach_relation(INLET_MACH, GAMMA)
    if use_print:
        print(f"A* = {A_STAR}")
    distances = np.concatenate([distances, np.flip(distances)])
    Mach_vec, P_vec = calculate_mach_pressure_distribution(distances, GAMMA, R_INLET, T_INLET, P_INLET, INLET_MACH, OUTLET_MACH, a_star=A_STAR)

    # plot_mach_pressure_distributions(Mach_vec, P_vec)

    # Check Kantrowitz limit
    kant_limit = kantrowitz_limit(INLET_MACH, GAMMA)
    max_area_ratio = np.max(distances) / A_STAR
    true_area_ratio = np.max(distances) / np.min(distances)
    kantrowitz_violation = max_area_ratio < kant_limit
    
    if use_print:
        print(f"Kantrowitz Limit Violated: {kantrowitz_violation}, LIMIT = {kant_limit}, AREA_RATIO  = {max_area_ratio}, TRUE_AREA_RAATIO = {true_area_ratio}")

    # Efficiecy calculations
    blade_max_thickness = blade_max_thickness / 100 # convert to m
    
    n_th = estimate_turbine_efficiency(HORSE_POWER, M_DOT_IMPERIAL, C_p, T_INLET_STATOR_IMPERIAL, R_INLET, GAMMA)
    if use_print:
        print(f"Efficiency from Isentropic Analysis: {n_th}")

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
        Q_freestream=Q_FREESTREAM,
        data=data,
        use_print=use_print
    )
    if use_print:
        print(f"Overall Losses: {y_loss}")

    isentropic_efficiency = (1 - y_loss)
    n_total = n_th * isentropic_efficiency
    expected_hp = HORSE_POWER * isentropic_efficiency

    if use_print:
        print(f"Overall Efficiency: {n_total}")
        print(f"Expected (After loss) Horsepower of Design: {expected_hp} HP")
        print(f"Isentropic Efficiency: {isentropic_efficiency}")

    return n_total, expected_hp


def load_kacker_okapuu_data(path_str="/KackerOkapuu/kackerokapuu.pkl"):
    path = pathlib.Path(os.getcwd() + "/KackerOkapuu/kackerokapuu.pkl")
    # global data
    data = None
    if not path.exists():
        print("Couldn't find data file. Downloading ... ")
        url = "https://github.com/nasa/turbo-design/raw/main/references/Turbines/KackerOkapuu/kackerokapuu.pkl"
        response = requests.get(url, stream=True)
        with open(path.absolute(), mode="wb") as file:
            for chunk in response.iter_content(chunk_size=10 * 1024):
                file.write(chunk)
        
    with open(path.absolute(),'rb') as f:
        data = pickle.load(f) # type: ignore
        print(data)
    
    return data