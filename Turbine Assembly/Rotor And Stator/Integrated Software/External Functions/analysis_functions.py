import numpy as np
from scipy.optimize import fsolve
import matplotlib.pyplot as plt

def calculate_mach_pressure_distribution(area_vec, gamma, R, T_inlet, P_inlet, M_inlet, M_outlet):
    # Initialize Mach number and pressure arrays
    Mach_vec = np.zeros_like(area_vec)
    P_vec = np.zeros_like(area_vec)

    # Total temperature and pressure at the inlet
    Tt_inlet = T_inlet * (1 + ((gamma - 1) / 2) * M_inlet**2)
    Pt_inlet = P_inlet * (1 + ((gamma - 1) / 2) * M_inlet**2)**(gamma / (gamma - 1))

    # Calculate Mach number along the blade
    for i in range(len(area_vec)):
        area_ratio = area_vec[i] / area_vec[0]

        # Assume isentropic and numerically solve for Mach number from area-Mach equation
        Mach_vec[i] = isentropic_mach_finder(area_ratio, gamma, M_inlet, M_outlet)

        # Calculate static temperature and pressure using isentropic relations
        T_static = Tt_inlet / (1 + ((gamma - 1) / 2) * Mach_vec[i]**2)
        P_vec[i] = Pt_inlet * (T_static / Tt_inlet)**(gamma / (gamma - 1))

    return Mach_vec, P_vec

def isentropic_mach_finder(area_ratio, gamma, M_inlet, M_outlet):
    # Initial guess for Mach number as the average of inlet and outlet Mach numbers
    Mach_guess = (M_inlet + M_outlet) / 2

    # Function for the isentropic area-Mach relation
    def func(M):
        return area_ratio - (1 / M) * ((2 / (gamma + 1)) * (1 + ((gamma - 1) / 2) * M**2))**((gamma + 1) / (2 * (gamma - 1)))

    # Solve for Mach number using fsolve
    Mach = fsolve(func, Mach_guess)

    return Mach[0]  # fsolve returns a list, return the first element

def plot_mach_pressure_distributions(Mach_vec, P_vec):
    x = np.arange(1, len(Mach_vec) + 1)  # lower coordinate

    plt.figure()
    plt.suptitle("Analysis along Lower Surface")
    
    # Mach Number Plot
    plt.subplot(2, 1, 1)
    plt.plot(x, Mach_vec, 'b-')
    plt.grid(True)
    plt.xlabel("Lower Surface Coordinate")
    plt.ylabel("Mach Number")
    plt.title("Mach Number vs Lower (Pressure) Surface Coordinate")

    # Pressure Distribution Plot
    plt.subplot(2, 1, 2)
    plt.plot(x, P_vec, 'g-')
    plt.grid(True)
    plt.xlabel("Lower Surface Coordinate")
    plt.ylabel("Static Pressure [Pa]")
    plt.title("Pressure Distribution along Pressure Surface")

    plt.tight_layout()
    plt.show()

def kantrowitz_limit(M0, gamma):
    term1 = M0 * ((gamma + 1) / (2 + (gamma - 1) * M0**2))**((gamma + 1) / (2 * (gamma - 1)))
    term2 = (((gamma + 1) * M0**2) / (((gamma - 1) * M0**2) + 2))**(-gamma / (gamma - 1))
    term3 = ((gamma + 1) / (2 * gamma * M0**2 - (gamma - 1)))**(-1 / (gamma - 1))
    
    A_ratio = term1 * term2 * term3
    return A_ratio

def calc_relative_mach(U, V, theta, gamma, R, T):
    a = np.sqrt(gamma * R * T)
    M_rel = np.sqrt((V * np.cos(theta) - U) ** 2 + (V * np.sin(theta)) ** 2) / a
    return M_rel

def calculate_inlet_area(U, V, theta, gamma, R, T, A_exit, M_exit):
    M_rel = calc_relative_mach(U, V, theta, gamma, R, T)

    A_star_inlet_ratio = area_mach_relation(M_rel, gamma)
    A_star_exit_ratio = area_mach_relation(M_exit, gamma)

    A_inlet = A_exit * (A_star_inlet_ratio / A_star_exit_ratio)
    return A_inlet

def area_mach_relation(M, gamma):
    A_ratio = (1 / M) * ((2 / (gamma + 1)) * (1 + (gamma - 1) / 2 * M**2))**((gamma + 1) / (2 * (gamma - 1)))
    return A_ratio