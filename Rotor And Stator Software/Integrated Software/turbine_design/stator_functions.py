import math

def calc_R_S(R, m_m):  # [J/(kg*K)] (Specific Gas Constant)
    return (R / m_m) * 1000


def calc_rho_0(P_0, R_S, T_0):  # [kg/m^3]
    return P_0 / (R_S * T_0)


def calc_v_e(T_0, R_S, gamma, P_e, P_0):  # [m/s]
    return math.sqrt((T_0 * R_S) * ((2 * gamma) / (gamma - 1)) * (1 - (P_e / P_0) ** ((gamma - 1) / gamma)))


def calc_T_throat(T_0, gamma):  # [K]
    return T_0 / (1 + ((gamma - 1) / 2))


def calc_P_throat(P_0, gamma):  # [N/m^2]
    return P_0 / ((1 + ((gamma - 1) / 2)) ** ((gamma - 1) / gamma))


def calc_rho_throat(P_throat, R_S, T_throat):  # [kg/m^3]
    return P_throat / (R_S * T_throat)


def calc_v_throat(gamma, R_S, T_throat):  # [m/s]
    return math.sqrt(gamma * R_S * T_throat)


def calc_A_throat(m_dot, rho_throat, v_throat):  # [m^2]
    return m_dot / (rho_throat * v_throat)


def calc_M_e(P_e, P_0, gamma):
    P_ratio = P_e / P_0
    numerator = math.exp(math.log(P_ratio) / (-gamma / (gamma - 1))) - 1
    denominator = (gamma - 1) / 2
    return math.sqrt(numerator / denominator)


def calc_rho_e(rho_0, M_e, gamma):  # [kg/m^3]
    demon = 1 + ((gamma - 1) / 2) * M_e ** 2
    return rho_0 * demon ** (-1 / (gamma - 1))


def calc_A_e(m_dot, rho_e, v_e):  # [m^2]
    return m_dot / (rho_e * v_e)


def calc_T_e(T_0, gamma, M_e):  # [K]
    demon = 1 + ((gamma - 1) / 2) * M_e ** 2
    return T_0 / demon


def calc_r_throat(A_throat):  # [m]
    return math.sqrt(A_throat / math.pi)


def calc_r_e(A_e):  # [m]
    return math.sqrt(A_e / math.pi)


def calc_dist(r_throat, r_e):  # [m]
    return (r_e - r_throat) / math.tan(math.radians(15))


def calc_A_throat_n(A_throat, n):  # [m^2]
    return A_throat / n


def calc_A_e_n(A_e, n):  # [m^2]
    return A_e / n


def calc_r_throat_n(A_throat_n):  # [m]
    return math.sqrt(A_throat_n / math.pi)


def calc_r_e_n(A_e_n):  # [m]
    return math.sqrt(A_e_n / math.pi)


def calc_dist_n(r_throat_n, r_e_n):  # [m]
    return (r_e_n - r_throat_n) / math.tan(math.radians(15))


def calc_F_thrust(m_dot, v_e):  # [N]
    return m_dot * v_e


## NOZZLE PLOTTING FUNCTIONS

import numpy as np
import matplotlib.pyplot as plt

def exp_scale(x):
    return (np.exp(x) - 1) / (np.exp(1) - 1)

def plot_nozzle(A_inlet, A_throat, A_exit, inlet_len, outlet_len):
    num_points = 100
    r_inlet = np.sqrt(A_inlet / np.pi)
    r_throat = np.sqrt(A_throat / np.pi)
    r_exit = np.sqrt(A_exit / np.pi)
    throat_len = 0.1 * inlet_len
    x_total = inlet_len + throat_len + outlet_len
    x = np.linspace(0, x_total, num_points + 1)
    thetas = np.linspace(0, 2 * np.pi, 1000)
    X, Thetas = np.meshgrid(x, thetas)

    R = np.zeros_like(X)
    # Inlet region
    inlet_mask = X <= inlet_len
    R[inlet_mask] = r_inlet + (r_throat - r_inlet) * exp_scale(X[inlet_mask] / inlet_len)
    # Throat region
    throat_mask = (X > inlet_len) & (X <= (inlet_len + throat_len))
    R[throat_mask] = r_throat
    # Outlet region
    outlet_mask = (X > (inlet_len + throat_len)) & (X <= x_total)
    R[outlet_mask] = r_throat + ((r_exit - r_throat) / outlet_len) * (X[outlet_mask] - (inlet_len + throat_len))

    Y = R * np.cos(Thetas)
    Z = R * np.sin(Thetas)
    C = X

    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')
    surf = ax.plot_surface(X, Y, Z, facecolors=plt.cm.viridis(C / x_total), shade=True)
    ax.set_xlabel("length [m]")
    ax.set_ylabel("width [m]")
    ax.set_zlabel("height [m]")
    ax.set_title("Plot of nozzle geometry")
    plt.colorbar(surf, ax=ax, shrink=0.5, aspect=5)
    plt.show()

    minZ = np.min(Z)
    maxZ = np.max(Z)
    plt.figure()
    plt.contour(X, Y, Z, levels=np.linspace(minZ, maxZ / 2, 10))
    plt.xlabel("length [m]")
    plt.ylabel("width [m]")
    plt.title("Contour Plot for nozzle")
    plt.show()

    return X, Y, Z

def save_to_stl(X, Y, Z, filename):
    # Optional: Requires numpy-stl and mesh generation from surface
    # This is a placeholder; actual STL export from surface mesh is nontrivial.
    print(f"Surface export to {filename} not implemented. Use a mesh generator like trimesh or numpy-stl.")
