import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import minimize, fsolve, fmin
from matplotlib.widgets import Slider
import random
import pandas as pd

num_points = 100

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