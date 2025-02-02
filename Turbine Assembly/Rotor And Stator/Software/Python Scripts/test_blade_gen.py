import numpy as np
import pygame

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.widgets import Slider

# Function to generate blade geometry
def generate_blade_geometry(chord, beta, A_inlet, B_spacing, blade_thickness):
    theta_l = np.radians(beta)
    num_points = 100

    # Lower surface guide
    x_lower = np.linspace(0, chord, num_points)
    y_guide_lower = np.tan(theta_l) * x_lower
    y_guide_lower[x_lower > chord / 2] = (
        np.tan(-theta_l) * (x_lower[x_lower > chord / 2] - chord / 2) +
        np.tan(theta_l) * (chord / 2)
    )

    # Lower surface 
    half_d = chord/2
    lower_slope = np.tan(theta_l)
    r = np.sqrt((half_d) ^ 2 + (half_d / lower_slope)^2)
    y_lower = np.sqrt(r^2 - (x_lower - half_d))^2 - np.sqrt(r^2 - half_d^2)


    # Upper surface guide curve (simplified for visualization)
    x_upper = np.linspace(0, chord, num_points)
    y_guide_upper = y_lower + B_spacing

    # Bézier curve for upper surface
    x1 = A_inlet * np.cos(theta_l)
    y1 = -A_inlet * np.sin(theta_l)
    x2 = chord / 2
    y2 = max(y_upper) + blade_thickness

    t = np.linspace(0, 1, num_points)
    B_x = (1 - t) ** 2 * x1 + 2 * (1 - t) * t * chord / 2 + t ** 2 * x2
    B_y = (1 - t) ** 2 * y1 + 2 * (1 - t) * t * y2 + t ** 2 * y2

    # Concatenate upper surface points
    x_upper = np.concatenate([x_upper, B_x])
    y_upper = np.concatenate([y_upper, B_y])

    return x_lower, y_lower, x_upper, y_upper

def calc_sq_err_min_dist(x_lower, y_lower, x0, y0, x1, slope, x2, y2, A_inlet, do_plot):
    num_points = 100
    y1 = y0 + slope * (x1 - x0)
    line_x = np.linspace(x0, x1, num_points)
    line_y = np.linspace(y0, y1, num_points)

    x_m = (y2 - y1) / slope + x1
    y_m = y2

    t = np.linspace(0, 1, num_points)
    bezier_x = (1 - t)**2 * x1 + 2 * (1 - t) * t * x_m + t**2 * x2
    bezier_y = (1 - t)**2 * y1 + 2 * (1 - t) * t * y_m + t**2 * y2

    x_upper = np.concatenate([line_x, bezier_x])
    y_upper = np.concatenate([line_y, bezier_y])

    # Numerical derivatives using finite differences
    dx_dt = np.gradient(bezier_x, t)
    dy_dt = np.gradient(bezier_y, t)
    d2x_dt2 = np.gradient(dx_dt, t)
    d2y_dt2 = np.gradient(dy_dt, t)

    # Curvature calculation
    curvature = (dx_dt * d2y_dt2 - dy_dt * d2x_dt2) / (dx_dt**2 + dy_dt**2)**(3/2)
    curvature[np.isnan(curvature) | np.isinf(curvature)] = 0
    penalty = np.sqrt(np.var(curvature))

    if do_plot:
        plt.figure()
        plt.plot(t, curvature, 'b-')
        plt.xlabel('t')
        plt.ylabel('Curvature κ(t)')
        plt.title('Curvature along the Bezier Curve')
        print(f'Penalty (variance of curvature): {penalty:.6f}')

    sample_size = round(0.3 * len(x_lower))
    sample_indices = np.sort(np.random.choice(len(x_lower), sample_size, replace=False))
    x_lower_sample = x_lower[sample_indices]
    y_lower_sample = y_lower[sample_indices]

    min_distances = np.zeros(len(x_lower_sample))
    min_dist_x = np.zeros(len(x_lower_sample))

    sse = 0

    for i, (x0, y0) in enumerate(zip(x_lower_sample, y_lower_sample)):
        if i < len(x_lower_sample) - 1:
            dx = x_lower_sample[i + 1] - x0
            dy = y_lower_sample[i + 1] - y0
        else:
            dx = x0 - x_lower_sample[i - 1]
            dy = y0 - y_lower_sample[i - 1]

        # Calculate normal (inwards/right)
        normal_dx = -dy
        normal_dy = dx

        # Normalize
        normal_length = np.sqrt(normal_dx**2 + normal_dy**2)
        normal_dx /= normal_length
        normal_dy /= normal_length

        projection_length = 0
        old_min = float('inf')
        min_distance = 99
        min_dist_idx = 0

        # Project the normal until it "intersects" the upper surface
        while min_distance < old_min:
            old_min = min_distance
            projection_length += 1e-5

            x_normal = x0 + normal_dx * projection_length
            y_normal = y0 + normal_dy * projection_length

            distances = np.sqrt((x_upper - x_normal)**2 + (y_upper - y_normal)**2)
            min_distance = np.min(distances)
            min_dist_idx = np.argmin(distances)

        min_distances[i] = old_min
        min_dist_x[i] = x_upper[min_dist_idx]

        sse += (old_min - A_inlet)**2

    err = sse

    if do_plot:
        print(f"Min SSE: {sse:.4f}")

        plt.figure()
        plt.plot(x_lower, y_lower, 'r-', label="Lower Surface")
        plt.plot(x_upper, y_upper, 'b-', label="Upper Surface")
        plt.title("Blade Channel Optimization")
        plt.axis('equal')
        plt.grid(True)
        plt.xlabel("Blade Coordinate [m]")
        plt.ylabel("Vertical Component [m]")
        plt.legend()
        plt.show()

    return err


# Function to calculate and display the live plot
def update_plot(val):
    chord = chord_slider.val
    beta = beta_slider.val
    B_spacing = spacing_slider.val
    blade_thickness = thickness_slider.val

    guide, surface = generate_blade_geometry(
        chord, beta, A_inlet, B_spacing, blade_thickness
    )

    x_lower, y_lower, x_upper, y_upper = surface
    x_guide_lower, y_guide_lower, x_guide_upper, y_guide_upper = guide
    

    # Update plot data
    lower_surface.set_data(x_lower, y_lower)
    upper_surface.set_data(x_upper, y_upper)

    # Recalculate normals for "area calculation"
    for normal in normals:
        normal.remove()
    normals.clear()
    for i in range(0, len(x_lower), 5):
        normal_dx = -(y_upper[i] - y_lower[i])
        normal_dy = (x_upper[i] - x_lower[i])
        length = np.sqrt(normal_dx ** 2 + normal_dy ** 2)
        normal_dx /= length
        normal_dy /= length
        normal = ax.arrow(
            x_lower[i], y_lower[i], 0.1 * normal_dx, 0.1 * normal_dy,
            head_width=0.01, head_length=0.02, color='green', alpha=0.7
        )
        normals.append(normal)

    fig.canvas.draw_idle()

# Parameters
A_inlet = 0.1  # Example parameter, adjust as needed
initial_chord = 1.0
initial_beta = 10.0
initial_B_spacing = 0.1
initial_blade_thickness = 0.05

# Create the figure and axis
fig, ax = plt.subplots()
plt.subplots_adjust(left=0.1, bottom=0.3)
ax.axis('equal')
ax.grid(True)
ax.set_xlabel("x [m]")
ax.set_ylabel("y [m]")
ax.set_title("Blade Geometry with Adjustable Parameters")

# Initial plot
x_lower, y_lower, x_upper, y_upper = generate_blade_geometry(
    initial_chord, initial_beta, A_inlet, initial_B_spacing, initial_blade_thickness
)
lower_surface, = ax.plot(x_lower, y_lower, 'b-', label="Lower Surface")
upper_surface, = ax.plot(x_upper, y_upper, 'r-', label="Upper Surface")
normals = []

# Add sliders
ax_chord = plt.axes([0.1, 0.2, 0.8, 0.03])
ax_beta = plt.axes([0.1, 0.15, 0.8, 0.03])
ax_spacing = plt.axes([0.1, 0.1, 0.8, 0.03])
ax_thickness = plt.axes([0.1, 0.05, 0.8, 0.03])

chord_slider = Slider(ax_chord, 'Chord', 0.5, 2.0, valinit=initial_chord)
beta_slider = Slider(ax_beta, 'Beta (°)', 0, 45, valinit=initial_beta)
spacing_slider = Slider(ax_spacing, 'Blade Spacing', 0.05, 0.2, valinit=initial_B_spacing)
thickness_slider = Slider(ax_thickness, 'Thickness', 0.02, 0.1, valinit=initial_blade_thickness)

# Update plot when sliders are adjusted
chord_slider.on_changed(update_plot)
beta_slider.on_changed(update_plot)
spacing_slider.on_changed(update_plot)
thickness_slider.on_changed(update_plot)

plt.legend()
plt.show()