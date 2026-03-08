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
import math

from .helpers import export_parameters_to_excel


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
        binomial_coeff = math.comb(n, i)
        curve[:, 0] += binomial_coeff * (1 - t) ** (n - i) * t ** i * x
        curve[:, 1] += binomial_coeff * (1 - t) ** (n - i) * t ** i * y
    return curve

def curvature(x, y):
    dx = np.gradient(x)
    dy = np.gradient(y)
    ddx = np.gradient(dx)
    ddy = np.gradient(dy)
    return (dx * ddy - dy * ddx) / (dx**2 + dy**2)**(3/2)

def compute_distance_error(x_lower, y_lower, upper_curve, target_distance, blade_spacing):
    y_lower = y_lower + blade_spacing;
    num_points = len(x_lower)
    # sampled_indices = random.sample(range(num_points // 2), min(20, num_points))
    sampled_indices = np.linspace(0, 2 * num_points//3, min(30, num_points), dtype=int)
    sse = 0

    dx = np.gradient(x_lower)
    dy = np.gradient(y_lower)

    distances = []

    # find index of upper curve where x is closest to x_lower[-1]
    upper_end_index = np.argmin(np.abs(upper_curve[:, 0] - x_lower[-1]))

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
        for bx, by in upper_curve[: upper_end_index]:
            # slope = (by - y0) / (bx - x0)
            slope_angle = np.angle((bx - x0) + (by - y0) * 1j)
            # slope_diff = abs(slope_angle - normal_angle)
            slope_diff = min(abs(slope_angle - normal_angle), 2 * np.pi - abs(slope_angle - normal_angle))
            distance = abs((bx - x0) * normal_dx + (by - y0) * normal_dy)
            if slope_diff < min_slope_diff:
                min_distance = distance
                min_slope_diff = slope_diff

        sse += (min_distance - target_distance) ** 2
        distances.append(min_distance)

    return sse, distances

def compute_arc(chord, beta, num_points):
    x_lower = np.linspace(0, chord, num_points)
    half_d = chord / 2
    lower_slope = np.tan(np.radians(beta))
    r = np.sqrt(half_d**2 + (half_d / lower_slope)**2)
    y_lower = np.sqrt(r**2 - (x_lower - half_d)**2) - np.sqrt(r**2 - half_d**2)
    y_prime = -(-half_d) / np.sqrt(r**2 - (-half_d)**2)
    return x_lower, y_lower, r, y_prime

def compute_upper_surface(params, target_distance, num_points, blade_spacing, radius, chord, y_prime_l, rounds_rad):
    x2, x3 = params
    c_d = chord / 2
    dist = c_d - x3

    r2 = radius - target_distance

    computed_parameters = {}

    computed_parameters["upper_radius"] = r2

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

    computed_parameters["bez_y2"] = y2
    computed_parameters["bez_x3"] = x3

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

    x_mid = 0.5 * (np.min(x_circ) + np.max(x_circ))
    x_rounds2_reflected = 2 * x_mid - x_rounds2[::-1]
    x_rounds1_reflected = 2 * x_mid - x_rounds1[::-1]

    y_rounds2_reversed = y_rounds2[::-1]
    y_rounds1_reversed = y_rounds1[::-1]

    x_upper = np.concatenate([x_rounds1, x_rounds2, bezier1[:, 0], x_circ, bezier2[:, 0], x_rounds2_reflected, x_rounds1_reflected])
    y_upper = np.concatenate([y_rounds1, y_rounds2, bezier1[:, 1], y_circ, bezier2[:, 1], y_rounds2_reversed, y_rounds1_reversed])
    
    # x_upper = np.concatenate([x_rounds1, x_rounds2, bezier1[:, 0], x_circ, bezier2[:, 0]])
    # y_upper = np.concatenate([y_rounds1, y_rounds2, bezier1[:, 1], y_circ, bezier2[:, 1]])

    upper = np.column_stack((x_upper, y_upper))
    return upper, computed_parameters

def objective_fn(params, x_lower, y_lower, target_distance, num_points, blade_spacing, radius, chord, y_prime_l, rounds_rad):
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

    upper_surface, param_dict = compute_upper_surface(params, target_distance, num_points, blade_spacing, radius, chord, y_prime_l, rounds_rad)
    # len_upper_surf = upper_surface.shape[0]

    upper_surface_filtered = upper_surface[upper_surface[:, 0] > 0]
    len_upper_surf = upper_surface_filtered.shape[0]

    # Filter points where x > 0 for the lower surface
    x_lower_filtered = x_lower[:num_points // 2][x_lower[:num_points // 2] > 0]
    y_lower_filtered = y_lower[:num_points // 2][x_lower[:num_points // 2] > 0]

    # Compute the distance error
    distance_error, distances = compute_distance_error(
        x_lower_filtered,
        y_lower_filtered,
        upper_surface_filtered[:len_upper_surf // 2, :],
        target_distance,
        blade_spacing
    )

    # distance_error = compute_distances(x_lower[:num_points // 2], y_lower[:num_points // 2], upper_surface[:len_upper_surf // 2, :], target_distance, blade_spacing)

    return distance_error # + curvature_variance

def optimize_upper_surface(x_lower, y_lower, target_distance, blade_spacing, radius, chord, y_prime_l, rounds_rad):
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
        args=(x_lower, y_lower, target_distance, num_points, blade_spacing, radius, chord, y_prime_l, rounds_rad),
        bounds=bounds,
        method='SLSQP'
    )

    return result.x

def compute_distances_and_plot(x_lower, y_lower, upper_curve, target_distance, blade_spacing) -> tuple[float, list[float]]:
    
    # ---- Create figure with 2 rows of subplots ----
    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(8, 10), height_ratios=[2, 1])
    plt.subplots_adjust(hspace=0.35)

    # ---- Plot upper & lower surfaces on ax1 ----
    ax1.plot(x_lower, y_lower, 'b-', label='Lower Surface')
    ax1.plot(upper_curve[:, 0], upper_curve[:, 1], 'r-', label='Upper Surface')
    ax1.plot(x_lower, y_lower + blade_spacing, 'b--', label='Lower Surface Guide')
    ax1.axis('equal')
    ax1.grid(True)
    ax1.set_xlabel('X/C')
    ax1.set_ylabel('Y/C')
    ax1.legend()

    # ---- compute distances ----
    y_lower = y_lower + blade_spacing
    num_points = len(x_lower)
    sampled_indices = np.linspace(0, num_points - 1, min(40, num_points), dtype=int)
    sse = 0

    dx = np.gradient(x_lower)
    dy = np.gradient(y_lower)

    distances = []

    upper_end_index = np.argmin(np.abs(upper_curve[:, 0] - x_lower[-1]))

    for i in sampled_indices:
        x0, y0 = x_lower[i], y_lower[i]
        normal_dx = dy[i]
        normal_dy = -dx[i]
        normal_length = np.sqrt(normal_dx**2 + normal_dy**2)
        normal_dx /= normal_length
        normal_dy /= normal_length
        normal_angle = np.angle(normal_dx + normal_dy*1j)

        min_distance = float('inf')
        min_slope_diff = float('inf')
        closest_point = None

        for bx, by in upper_curve[0: upper_end_index]:
            slope_angle = np.angle((bx - x0) + (by - y0) * 1j)
            slope_diff = min(abs(slope_angle - normal_angle),
                             2 * np.pi - abs(slope_angle - normal_angle))
            distance = abs((bx - x0) * normal_dx + (by - y0) * normal_dy)
            if slope_diff < min_slope_diff:
                min_distance = distance
                min_slope_diff = slope_diff
                closest_point = (bx, by)

        if closest_point:
            bx, by = closest_point
            ax1.plot([x0, bx], [y0, by], 'g--', linewidth=0.5)

        distances.append(min_distance)
        sse += (min_distance - target_distance) ** 2

    # ---- Plot distances on ax2 ----
    ax2.plot(x_lower[sampled_indices], distances, 'ko-', label='Computed Distances')
    ax2.axhline(target_distance, color='r', linestyle='--', label='Target Distance')
    ax2.set_xlabel('X/C (Lower Surface)')
    ax2.set_ylabel('Distance (% of Chord)')
    ax2.legend()
    ax2.grid(True)

    plt.show()

    return sse, distances



def plot_geometry_with_sliders(num_points = 1000, beta = np.deg2rad(60)):
    """
    INITIAL PARAMTERS FROM VELOCITY TRIANGLES
    """

    # in cm
    initial_chord = 1.0
    initial_beta = beta
    initial_target_distance = 0.0028 / 0.00822 # a percentage of blade spacing
    initial_blade_spacing = 0.8222
    initial_rounds = 0.01785 #0.0005

    x_lower, y_lower, radius, y_prime_l = compute_arc(initial_chord, initial_beta, num_points)
    params = optimize_upper_surface(x_lower, y_lower, initial_target_distance * initial_blade_spacing, initial_blade_spacing, radius, initial_chord, y_prime_l, initial_rounds)

    full_upper, param_dict = compute_upper_surface(params, initial_target_distance * initial_blade_spacing, num_points, initial_blade_spacing, radius, initial_chord, y_prime_l, initial_rounds)

    fig, ax = plt.subplots()
    plt.subplots_adjust(left=0.1, bottom=0.3)

    lower_surface, = ax.plot(x_lower, y_lower, 'b-', label='Lower Surface')
    upper_surface, = ax.plot(full_upper[:, 0], full_upper[:, 1], 'r-', label='Upper Surface')
    lower_surface_guide, = ax.plot(x_lower, y_lower + initial_blade_spacing, 'b-', label='Lower Surface Guide')

    upper_surface_filtered = full_upper[full_upper[:, 0] > 0]
    len_upper_surf = upper_surface_filtered.shape[0]

    sse, distances = compute_distances_and_plot(x_lower[:num_points // 2], y_lower[:num_points // 2], upper_surface_filtered[:len_upper_surf // 2, :], initial_target_distance * initial_blade_spacing, initial_blade_spacing, ax)


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

    def update_plot(val):
        chord = chord_slider.val
        beta = beta_slider.val
        target_distance_percent = target_slider.val
        blade_spacing = blade_spacing_slider.val
        rounds_rad = round_slider.val

        global x_lower, y_lower
        x_lower, y_lower, radius, y_prime_l = compute_arc(chord, beta, num_points)

        params = optimize_upper_surface(x_lower, y_lower, target_distance_percent * blade_spacing, blade_spacing, radius, chord, y_prime_l, rounds_rad)

        full_upper, param_dict = compute_upper_surface(params, target_distance_percent * blade_spacing, num_points, blade_spacing, radius, chord, y_prime_l, initial_rounds)

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

        global sse, distances
        sse, distances = compute_distances_and_plot(x_lower[:num_points // 2], y_lower[:num_points // 2], upper_surface_filtered[:len_upper_surf // 2, :], target_distance_percent * blade_spacing, blade_spacing, ax)
        # compute_distances_and_plot(x_lower_filtered, y_lower_filtered, upper_surface_filtered[:len_upper_surf // 2, :], target_distance_percent * blade_spacing, blade_spacing, ax)

        fig.canvas.draw_idle()

        print(f"{radius}\n{rounds_rad}\n{params[0]}\n{param_dict["bez_y2"]}\n{param_dict["upper_radius"]}\n{blade_spacing}\n{beta}\n{chord}")

        configurations = ["TurbineConfig"]
        parameters = {
            "lower_rad": [radius],
            "fillet_rad": [rounds_rad],
            "bez_x": [params[0]],
            "bez_y": param_dict["bez_y2"],
            "bez_x3": param_dict["bez_x3"],
            "upper_rad": param_dict["upper_radius"],
            "b_spacing": [blade_spacing],
            "l_angle": [beta],
            "chord": [chord],
        }
        export_parameters_to_excel("design_table.xlsx", configurations, parameters)

    chord_slider.on_changed(update_plot)
    beta_slider.on_changed(update_plot)
    target_slider.on_changed(update_plot)
    blade_spacing_slider.on_changed(update_plot)
    round_slider.on_changed(update_plot)

    plt.show()

    x_upper, y_upper = full_upper[:, 0], full_upper[:, 1]

    return x_lower, y_lower, x_upper, y_upper, distances

def compute_geometry(chord, beta, target_distance, blade_spacing, rounds, num_points = 1000, output_excel=False, output_path="design_table.xlsx") -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray, list[float], list]:
    """
    INITIAL PARAMTERS FROM VELOCITY TRIANGLES
    """

    # normalize by chord
    target_distance = target_distance / chord
    blade_spacing = blade_spacing / chord
    rounds = rounds
    chord_norm = 1.0

    x_lower, y_lower, radius, y_prime_l = compute_arc(chord_norm, beta, num_points)
    params = optimize_upper_surface(x_lower, y_lower, target_distance, blade_spacing, radius, chord_norm, y_prime_l, rounds)

    full_upper, param_dict = compute_upper_surface(params, target_distance, num_points, blade_spacing, radius, chord_norm, y_prime_l, rounds)

    upper_surface_filtered = full_upper[full_upper[:, 0] > 0]
    len_upper_surf = upper_surface_filtered.shape[0]

    dist_plot_params = [x_lower[:num_points // 2], y_lower[:num_points // 2], upper_surface_filtered[:len_upper_surf // 2, :], target_distance, blade_spacing]
    sse, distances = compute_distance_error(*dist_plot_params)

    x_upper, y_upper = full_upper[:, 0], full_upper[:, 1]

    configurations = ["TurbineConfig"]
    parameters = {
        "lower_rad": [radius],
        "fillet_rad": [rounds],
        "bez_x": [params[0]],
        "bez_y": param_dict["bez_y2"],
        "bez_x3": param_dict["bez_x3"],
        "upper_rad": param_dict["upper_radius"],
        "b_spacing": [blade_spacing],
        "l_angle": [beta],
        "chord": [chord],
    }
    if output_excel:
        export_parameters_to_excel(output_path, configurations, parameters)

    return x_lower, y_lower, x_upper, y_upper, distances, dist_plot_params


def plot_blade_geometry(x_lower, y_lower, x_upper, y_upper):
    fig, ax = plt.subplots()
    ax.plot(x_lower, y_lower, 'b-', label='Lower Surface')
    ax.plot(x_upper, y_upper, 'r-', label='Upper Surface')
    ax.axis('equal')
    ax.grid(True)
    ax.set_xlabel('X/C')
    ax.set_ylabel('Y/C')
    ax.legend()
    plt.show()

def plot_distance_distribution(x_lower, distances):
    fig, ax = plt.subplots()
    ax.plot(x_lower[:len(distances)], distances, 'g-', label='Distance Distribution')
    ax.axhline(y=np.mean(distances), color='r', linestyle='--', label='Mean Distance')
    ax.axis('equal')
    ax.grid(True)
    ax.set_xlabel('X')
    ax.set_ylabel('Distance to Upper Surface')
    ax.legend()
    plt.show()
