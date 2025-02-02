import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import minimize
from matplotlib.widgets import Slider
import random

num_points = 100

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
    sampled_indices = random.sample(range(num_points // 2), min(20, num_points))
    sse = 0

    for i in sampled_indices:
        x0, y0 = x_lower[i], y_lower[i]
        dx = np.gradient(x_lower)
        dy = np.gradient(y_lower)
        normal_dx = -dy[i]
        normal_dy = dx[i]
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
            slope_diff = abs(slope_angle - normal_angle)
            distance = abs((bx - x0) * normal_dx + (by - y0) * normal_dy)
            if slope_diff < min_slope_diff:
                min_distance = distance

        sse += (min_distance - target_distance) ** 2

    return sse

def compute_arc(chord, beta, num_points):
    x_lower = np.linspace(0, chord, num_points)
    half_d = chord / 2
    lower_slope = np.tan(np.radians(beta))
    r = np.sqrt(half_d**2 + (half_d / lower_slope)**2)
    y_lower = np.sqrt(r**2 - (x_lower - half_d)**2) - np.sqrt(r**2 - half_d**2)
    return x_lower, y_lower, r

def compute_upper_surface(params, target_distance, num_points, blade_spacing, radius, chord):
    x2, x3 = params
    c_d = chord / 2
    dist = c_d - x3

    r2 = radius - target_distance

    # circular section
    x_circ = np.linspace(x3, c_d + dist, num_points // 3)
    y_circ = np.sqrt(r2**2 - (x_circ - c_d)**2) - np.sqrt(radius**2 - c_d**2) + blade_spacing

    # bezier section
    # x_bez1 = np.linspace(x1, x3, num_points // 3)
    x1 = 0
    y1 = 0
    y3 = y_circ[0]
    slope = (-(c_d - dist) + c_d) / np.sqrt(-((c_d - dist)**2) + 2*(c_d - dist)*c_d + r2**2 - (c_d**2))
    y2 = slope * (x2 - c_d + dist) + y3

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

    x_upper = np.concatenate([bezier1[:, 0], x_circ, bezier2[:, 0]])
    y_upper = np.concatenate([bezier1[:, 1], y_circ, bezier2[:, 1]])
    
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
    len_upper_surf = upper_surface.shape[0]

    distance_error = compute_distances(x_lower[:num_points // 2], y_lower[:num_points // 2], upper_surface[:len_upper_surf // 2, :], target_distance, blade_spacing)

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

def update_plot(val):
    chord = chord_slider.val
    beta = beta_slider.val
    target_distance_percent = target_slider.val
    blade_spacing = blade_spacing_slider.val

    x_lower, y_lower, radius = compute_arc(chord, beta, num_points)
    params = optimize_upper_surface(x_lower, y_lower, target_distance_percent * blade_spacing, blade_spacing, radius, chord)

    # half_upper = compute_upper_surface(params, target_distance_percent * blade_spacing, num_points, blade_spacing, radius, chord)
    # full_upper = np.vstack([
    #     half_upper,
    #     np.flipud(np.column_stack((chord - half_upper[:, 0], half_upper[:, 1])))
    # ])
    full_upper = compute_upper_surface(params, target_distance_percent * blade_spacing, num_points, blade_spacing, radius, chord)

    lower_surface.set_data(x_lower, y_lower)
    upper_surface.set_data(full_upper[:, 0], full_upper[:, 1])
    lower_surface_guide.set_data(x_lower, y_lower + blade_spacing)

    fig.canvas.draw_idle()


def plot_geometry_with_sliders():
    global fig, ax, chord_slider, beta_slider, target_slider, blade_spacing_slider, lower_surface, upper_surface, lower_surface_guide

    initial_chord = 1.0
    initial_beta = 30.0
    initial_target_distance = 0.1
    initial_blade_spacing = 1.0

    x_lower, y_lower, radius = compute_arc(initial_chord, initial_beta, num_points)
    params = optimize_upper_surface(x_lower, y_lower, initial_target_distance * initial_blade_spacing, initial_blade_spacing, radius, initial_chord)

    # t = np.linspace(0, 1, 100)
    # bezier = bezier_curve(bezier_control_points, t)

    # half_upper = compute_upper_surface(params, initial_target_distance * initial_blade_spacing, initial_blade_spacing, radius, initial_chord)
    # full_upper = np.vstack([
    #     half_upper,
    #     np.flipud(np.column_stack((initial_chord - half_upper[:, 0], half_upper[:, 1])))
    # ])
    full_upper = compute_upper_surface(params, initial_target_distance * initial_blade_spacing, num_points, initial_blade_spacing, radius, initial_chord)

    fig, ax = plt.subplots()
    plt.subplots_adjust(left=0.1, bottom=0.3)

    lower_surface, = ax.plot(x_lower, y_lower, 'b-', label='Lower Surface')
    upper_surface, = ax.plot(full_upper[:, 0], full_upper[:, 1], 'r-', label='Upper Surface')
    lower_surface_guide, = ax.plot(x_lower, y_lower + initial_blade_spacing, 'b-', label='Lower Surface Guide')

    ax.axis('equal')
    ax.grid(True)
    ax.set_xlabel('X')
    ax.set_ylabel('Y')
    ax.legend()

    ax_chord = plt.axes([0.1, 0.2, 0.8, 0.03])
    ax_beta = plt.axes([0.1, 0.15, 0.8, 0.03])
    ax_target = plt.axes([0.1, 0.1, 0.8, 0.03])
    ax_bspacing = plt.axes([0.1, 0.05, 0.8, 0.03])

    chord_slider = Slider(ax_chord, 'Chord', 0.5, 2.0, valinit=initial_chord)
    beta_slider = Slider(ax_beta, 'Beta', 0, 75, valinit=initial_beta)
    target_slider = Slider(ax_target, 'Target Area', 0.0, 1.0, valinit=initial_target_distance)
    blade_spacing_slider = Slider(ax_bspacing, 'Blade Spacing', 0.2, 3.0, valinit=initial_target_distance)

    chord_slider.on_changed(update_plot)
    beta_slider.on_changed(update_plot)
    target_slider.on_changed(update_plot)
    blade_spacing_slider.on_changed(update_plot)

    plt.show()

plot_geometry_with_sliders()
