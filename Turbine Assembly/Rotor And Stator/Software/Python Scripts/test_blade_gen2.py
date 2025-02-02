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

def compute_distances(x_lower, y_lower, bezier, target_distance, blade_spacing):
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
        for bx, by in bezier:
            # slope = (by - y0) / (bx - x0)
            slope = np.angle((bx - x0) + (by - y0) * 1j)
            slope_diff = abs(slope - normal_angle)
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
    return x_lower, y_lower

def symmetric_objective_function(params, x_lower, y_lower, target_distance, num_points, blade_spacing):
    x1, y1, x2, y2 = params  # Only optimize up to x = chord/2
    x_mid, y_mid = x_lower[num_points // 2], y_lower[num_points // 2]
    control_points = [
        (x_lower[0], y_lower[0]),
        (x1, y1),
        (x2, y2),
        (x_mid, y2),
    ]
    t = np.linspace(0, 1, num_points // 2)
    bezier = bezier_curve(control_points, t)

    bezier_curvature = curvature(bezier[:, 0], bezier[:, 1])
    curvature_variance = np.var(bezier_curvature)

    distance_error = compute_distances(x_lower[:num_points // 2], y_lower[:num_points // 2], bezier, target_distance, blade_spacing)

    return distance_error + curvature_variance

def optimize_upper_surface(x_lower, y_lower, target_distance, blade_spacing):
    num_points = len(x_lower)
    x_mid, y_mid = x_lower[num_points // 2], y_lower[num_points // 2]
    initial_guess = [x_lower[num_points // 4], y_lower[num_points // 4], x_lower[num_points // 2 - 1], y_lower[num_points // 2 - 1] + blade_spacing / 2]
    bounds = [
        (x_lower[0], x_mid), 
        (min(y_lower), y_mid + blade_spacing), 
        (x_lower[num_points // 4], x_mid),  
        (y_mid, y_mid + blade_spacing)
    ]

    result = minimize(
        symmetric_objective_function,
        initial_guess,
        args=(x_lower, y_lower, target_distance, num_points, blade_spacing),
        bounds=bounds,
        method='SLSQP'
    )

    x1, y1, x2, y2 = result.x
    return [
        (x_lower[0], y_lower[0]),
        (x1, y1),
        (x2, y2),
        (x_mid, y2),
    ]

def update_plot(val):
    chord = chord_slider.val
    beta = beta_slider.val
    target_distance = target_slider.val
    blade_spacing = blade_spacing_slider.val

    x_lower, y_lower = compute_arc(chord, beta, num_points)
    bezier_control_points = optimize_upper_surface(x_lower, y_lower, target_distance * blade_spacing, blade_spacing)

    t = np.linspace(0, 1, num_points // 2)
    bezier_half = bezier_curve(bezier_control_points, t)
    bezier_full = np.vstack([
        bezier_half,
        np.flipud(np.column_stack((chord - bezier_half[:, 0], bezier_half[:, 1])))
    ])

    lower_surface.set_data(x_lower, y_lower)
    upper_surface.set_data(bezier_full[:, 0], bezier_full[:, 1])
    lower_surface_guide.set_data(x_lower, y_lower + blade_spacing)

    fig.canvas.draw_idle()


def plot_geometry_with_sliders():
    global fig, ax, chord_slider, beta_slider, target_slider, blade_spacing_slider, lower_surface, upper_surface, lower_surface_guide

    initial_chord = 1.0
    initial_beta = 30.0
    initial_target_distance = 0.1
    initial_blade_spacing = 1.0

    x_lower, y_lower = compute_arc(initial_chord, initial_beta, num_points)
    bezier_control_points = optimize_upper_surface(x_lower, y_lower, initial_target_distance * initial_blade_spacing, initial_blade_spacing)

    t = np.linspace(0, 1, 100)
    bezier = bezier_curve(bezier_control_points, t)

    fig, ax = plt.subplots()
    plt.subplots_adjust(left=0.1, bottom=0.3)

    lower_surface, = ax.plot(x_lower, y_lower, 'b-', label='Lower Surface')
    upper_surface, = ax.plot(bezier[:, 0], bezier[:, 1], 'r-', label='Upper Surface')
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
