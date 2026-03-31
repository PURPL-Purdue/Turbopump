"""
Functions for common isentropic relations for compressible flow
Author: Amanjyoti Mridha
"""

import numpy as np
from scipy.optimize import fsolve, least_squares, minimize

def calc_A_ratio(M, gamma):
    """
    Calculates A/A* for given M and gamma for istentropic flow
    """
    return ((gamma + 1) / 2) ** (-(gamma + 1) / (2 * (gamma - 1))) * (1 + (gamma - 1) / 2 * M**2)**((gamma + 1) / (2 * (gamma - 1))) / M

def calc_M_from_A_star_ratio(A_ratio, gamma, subsonic=True):
    """
    Calculates M for given A/A* and gamma for istentropic flow
    If solution not found, returns 1
    """
    
    def equation(M):
        return calc_A_ratio(M, gamma) - A_ratio
    
    if subsonic:
        M_initial_guess = 0.5
        bounds = (0, 1)
    else:
        M_initial_guess = 2.0
        bounds = (1, np.inf)
    
    # M_solution, = fsolve(equation, M_initial_guess)
    result = least_squares(equation, M_initial_guess, bounds=bounds)
    if not result.success:
        print("Warning: Could not find a solution for M from A/A* ratio. Check inputs.")
        return 1
    return result.x[0]

def calc_M_from_A_ratio(A_ratio, gamma, M1):
    """
    Calculates M2 for given A1/A2, will return 1 if flow is choked
    """
    if M1 == 1:
        if A_ratio == 1:
            return 1;
        elif A_ratio > 1: # converging (A1>A2) -> subsonic (not very realistic case)
            return calc_M_from_A_star_ratio(A_ratio, gamma, subsonic=True)
        else: # A1 < A2 -> diverging -> supersonic
            # going from choked up to supersonic
            return calc_M_from_A_star_ratio(A_ratio, gamma, subsonic=False)
    
    A1_A_star = calc_A_ratio(M1, gamma)
    if M1 < 1 and A1_A_star <= A_ratio: # A* > A2 => A1/A* < A1/A2
        return 1.0  # Choked flow
    elif M1 > 1 and A1_A_star <= A_ratio: # A* > A2 => A1/A* < A1/A2
        return 1.0  # Choked flow but unrealistic to approach from supersonic
    else:
        A2_A_star = A1_A_star / A_ratio
        return calc_M_from_A_star_ratio(A2_A_star, gamma, subsonic=(M1 < 1))

def calc_P0_P_ratio(M, gamma):
    """
    Calculates P0/P for given M and gamma for istentropic flow
    """
    return (1 + (gamma - 1) / 2 * M**2) ** (gamma / (gamma - 1))

def calc_T0_T_ratio(M, gamma):
    """
    Calculates T0/T for given M and gamma for istentropic flow
    """
    return 1 + (gamma - 1) / 2 * M**2

def calc_rho0_rho_ratio(M, gamma):
    """
    Calculates rho0/rho for given M and gamma for istentropic flow
    """
    return (1 + (gamma - 1) / 2 * M**2) ** (1 / (gamma - 1))

def calc_P_P0_ratio(M, gamma):
    """
    Calculates P/P0 for given M and gamma for istentropic flow
    """
    return 1 / calc_P0_P_ratio(M, gamma)

def calc_T_T0_ratio(M, gamma):
    """
    Calculates T/T0 for given M and gamma for istentropic flow
    """
    return 1 / calc_T0_T_ratio(M, gamma)

def calc_rho_rho0_ratio(M, gamma):
    """
    Calculates rho/rho0 for given M and gamma for istentropic flow
    """
    return 1 / calc_rho0_rho_ratio(M, gamma)

def calc_G_func(M, gamma):
    """
    Calculates theoretical G function for M and gamma
    """
    return (1 + gamma * M**2) / (1 + (gamma - 1) / 2 * M**2)**(gamma / (gamma - 1))

def calculate_D_func(M, gamma):
    """
    Calculates theoretical D function for M and gamma
    """
    return M / (1 + (gamma - 1) / 2 * M**2)**((gamma + 1) / (2 * (gamma - 1)))

def calculate_N_func(M, gamma):
    """
    Calculates theoretical N function for M and gamma
    """
    return calculate_D_func(M, gamma) / calc_G_func(M, gamma)

