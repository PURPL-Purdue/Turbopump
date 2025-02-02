import numpy as np
import sympy as sp
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from scipy.optimize import fsolve, minimize


def calc_blade_speed(radius, rpm):
    """Calculate blade speed given radius and RPM."""
    omega = rpm * 2 * np.pi / 60  # Convert RPM to angular velocity (rad/s)
    return np.array([radius * omega, 0])  # Blade speed (m/s)


def rotor_back_calculate(RPM, torque, mass_flow, beta, radius, V_in):
    """Back-calculate rotor parameters given operating conditions."""
    U = calc_blade_speed(radius, RPM)[0]
    C = torque / (radius * mass_flow)

    v2, a1, a2, w, b = sp.symbols("v2 a1 a2 w b")

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
        label="w1",
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
        label="u1",
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
        label="v1",
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
        label="w2",
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
        label="u2",
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
        label="v2",
    )

    plt.legend()
    plt.xlabel("X")
    plt.ylabel("Y")
    plt.title("Velocity Triangles with Chord Inclination")
    plt.show()


def calculate_blade_efficiency(mass_flow_rate, C1, C2, W1, W2, beta1, beta2, alpha1):
    kb = W2 / W1  # single-stage, single-rotor impulse turbine

    # power transferred to rotor blades
    Pb = (
        (1 / 4)
        * mass_flow_rate
        * C1**2
        * (1 + kb * np.cos(beta2) / np.cos(beta1))
        * np.cos(alpha1) ** 2
    )

    # input power
    input_power = (1 / 2) * mass_flow_rate * C1**2

    # blade efficiency
    efficiency = Pb / input_power

    return efficiency


def curvature_function(y2, m, x1, x2, y1):
    x_m = (y2 - y1) / m + x1
    return abs((2 * x2 - 2 * x_m) * (2 * y1 - 2 * y2)) / (
        8 * ((x2 - x_m) ** 2) ** (3 / 2)
    )


def generate_blade_geom(c, beta, A_inlet, B_spacing, blade_thickness):
    theta_l = np.radians(beta)
    num_points = 100

    # Lower Surface
    x1 = np.linspace(0, c / 2, num_points)
    bottom_curve1 = np.tan(theta_l) * x1
    x2 = np.linspace(c / 2, c, num_points)
    bottom_curve2 = np.tan(-theta_l) * (x2 - c / 2) + np.tan(theta_l) * (c / 2)

    r = c / 2 * np.sin(theta_l)
    C1 = c / 2 - np.sqrt((np.tan(theta_l) ** 2 * r**2) / (1 + np.tan(theta_l) ** 2))
    C2 = c / 2 + np.sqrt((np.tan(theta_l) ** 2 * r**2) / (1 + np.tan(theta_l) ** 2))
    x1 = np.linspace(0, C1, num_points)
    x2 = np.linspace(C1, C2, num_points)
    x3 = np.linspace(C2, c, num_points)
    curve_DC_1 = np.tan(theta_l) * x1
    curve_CC = np.sqrt(r**2 - (x2 - c / 2) ** 2)
    curve_DC_2 = np.tan(-theta_l) * (x3 - c / 2) + np.tan(theta_l) * (c / 2)

    x_lower = np.concatenate([x1, x2, x3])
    y_lower = np.concatenate([curve_DC_1, curve_CC, curve_DC_2])

    # Upper Surface Guide Curves
    x_upper1 = np.linspace(0, A_inlet * np.cos(theta_l), num_points)
    upper_curve1 = -1 / np.tan(theta_l) * x_upper1
    x_upper2 = np.linspace(c - A_inlet * np.cos(theta_l), c, num_points)
    upper_curve2 = 1 / np.tan(theta_l) * (x_upper2 - c)

    # Bézier Curve
    x1 = A_inlet * np.cos(theta_l)
    y1 = -1 / np.tan(theta_l) * x1
    x2 = c / 2
    m = (upper_curve1[-1] - (-B_spacing)) / (A_inlet * np.cos(theta_l))

    y2_initial_guess = r - B_spacing
    y2 = fsolve(
        lambda y2: curvature_function(y2, m, x1, x2, y1) - (1 / (r - A_inlet)),
        y2_initial_guess,
    )[0]
    y2 = r - B_spacing + blade_thickness
    max_thickness = r

    x_m = (y2 - y1) / m + x1
    t = np.linspace(0, 1, num_points)
    B1 = (1 - t) ** 2 * x1 + 2 * (1 - t) * t * x_m + t**2 * x2
    B2 = (1 - t) ** 2 * y1 + 2 * (1 - t) * t * y2 + t**2 * y2

    x_upper = np.concatenate([x_upper1, B1, np.flip(c - B1), x_upper2])
    y_upper = np.concatenate([upper_curve1, B2, np.flip(B2), upper_curve2]) + B_spacing

    # Plot
    plt.figure()
    plt.plot(x_lower, y_lower, "k-", linewidth=2)
    plt.plot(x_upper, y_upper, "k-", linewidth=2)
    plt.xlabel("x [m]")
    plt.ylabel("y [m]")
    plt.title("Blade Geometry Parametric Curves")
    plt.grid(True)
    plt.axis("equal")
    plt.show()

    return max_thickness, x_lower, y_lower, x_upper, y_upper


def generate_blade_geom_constant_area_v2(c, beta, A_inlet, B_spacing, blade_thickness):
    theta_l = np.radians(beta)  # Convert to radians

    num_points = 100

    # Lower Surface Guide Curves
    x1 = np.linspace(0, c / 2, num_points)
    bottom_curve1 = np.tan(theta_l) * x1
    x2 = np.linspace(c / 2, c, num_points)
    bottom_curve2 = np.tan(-theta_l) * (x2 - c / 2) + np.tan(theta_l) * (c / 2)

    # Lower Surface
    half_d = c / 2
    lower_slope = np.tan(theta_l)
    r = np.sqrt(half_d**2 + (half_d / lower_slope) ** 2)
    x_lower = np.linspace(0, c, num_points)
    y_lower = np.sqrt(r**2 - (x_lower - half_d) ** 2) - np.sqrt(r**2 - half_d**2)

    # Upper Surface Guide Curves
    x_upper1 = np.linspace(0, A_inlet * np.cos(theta_l), num_points)
    upper_curve1 = -1 / np.tan(theta_l) * x_upper1
    x_upper2 = np.linspace(c - A_inlet * np.cos(theta_l), c, num_points)
    upper_curve2 = 1 / np.tan(theta_l) * (x_upper2 - c)

    # Solve for y2 such that curvature k = 1
    y2 = r - B_spacing + blade_thickness
    max_thickness = r

    x_m = (y2 - upper_curve1[-1]) / (
        (upper_curve1[-1] + B_spacing) / A_inlet * np.cos(theta_l)
    ) + x1[-1]

    def calc_sq_err_min_dist(vars):
        x1_opt, y2_opt = vars
        return np.abs((2 * x2[-1] - 2 * x_m) * (2 * upper_curve1[-1] - 2 * y2_opt)) / (
            8 * ((x2[-1] - x_m) ** 2) ** (3 / 2)
        )

    initial_guess = [x1[-1], y2]
    result = minimize(calc_sq_err_min_dist, initial_guess, method="Nelder-Mead")
    optimal_x1, optimal_y2 = result.x

    # Parametric Bézier Curve
    t = np.linspace(0, 1, num_points)
    B1 = (1 - t) ** 2 * optimal_x1 + 2 * (1 - t) * t * x_m + t**2 * x2[-1]
    B2 = (
        (1 - t) ** 2 * upper_curve1[-1]
        + 2 * (1 - t) * t * optimal_y2
        + t**2 * optimal_y2
    )

    # Upper Surface
    x_upper = np.concatenate([x_upper1, B1, np.flip(c - B1), x_upper2])
    y_upper = np.concatenate(
        [upper_curve1 - B_spacing, B2, np.flip(B2), upper_curve2 - B_spacing]
    )

    # Plotting
    plt.figure()
    plt.plot(x1, bottom_curve1, "b--", linewidth=2)
    plt.plot(x2, bottom_curve2, "b--", linewidth=2)
    plt.plot(x_lower, y_lower, "k-", linewidth=2)
    plt.plot(x_upper, y_upper, "k-", linewidth=2)
    plt.xlabel("x [m]")
    plt.ylabel("y [m]")
    plt.title("Blade Geometry Parametric Curves")
    plt.grid()
    plt.axis("equal")
    plt.show()

    # Calculate cross-sectional area
    cross_sectional_area = np.trapz(y_upper - y_lower, x_upper)

    return max_thickness, cross_sectional_area, x_lower, y_lower, x_upper, y_upper


def calc_sq_err_min_dist(x_lower, y_lower, x0, y0, x1, slope, x2, y2, A_inlet, do_plot):
    num_points = 100
    y1 = y0 + slope * (x1 - x0)
    line_x = np.linspace(x0, x1, num_points)
    line_y = np.linspace(y0, y1, num_points)

    x_m = (y2 - y1) / slope + x1
    y_m = y2

    t = np.linspace(0, 1, num_points)
    bezier_x = (1 - t) ** 2 * x1 + 2 * (1 - t) * t * x_m + t**2 * x2
    bezier_y = (1 - t) ** 2 * y1 + 2 * (1 - t) * t * y_m + t**2 * y2

    x_upper = np.concatenate([line_x, bezier_x])
    y_upper = np.concatenate([line_y, bezier_y])

    # Numerical derivatives using finite differences
    dx_dt = np.gradient(bezier_x, t)  # First derivative of x
    dy_dt = np.gradient(bezier_y, t)  # First derivative of y
    d2x_dt2 = np.gradient(dx_dt, t)  # Second derivative of x
    d2y_dt2 = np.gradient(dy_dt, t)  # Second derivative of y

    # Curvature calculation
    curvature = (dx_dt * d2y_dt2 - dy_dt * d2x_dt2) / (dx_dt**2 + dy_dt**2) ** (3 / 2)
    curvature[np.isnan(curvature) | np.isinf(curvature)] = 0
    penalty = np.sqrt(np.var(curvature))

    if do_plot:
        plt.figure()
        plt.plot(t, curvature, "b-")
        plt.xlabel("t")
        plt.ylabel("Curvature κ(t)")
        plt.title("Curvature along the Bezier Curve")
        print(f"Penalty (variance of curvature): {penalty:.6f}")

    sample_size = round(0.3 * len(x_lower))
    sample_indices = np.sort(np.random.choice(len(x_lower), sample_size, replace=False))
    x_lower = x_lower[sample_indices]
    y_lower = y_lower[sample_indices]

    min_distances = np.zeros(len(x_lower))
    min_dist_x = np.zeros(len(x_lower))

    sse = 0

    for i in range(len(x_lower)):
        x0 = x_lower[i]
        y0 = y_lower[i]

        # Approx tangent with finite differences
        if i < len(x_lower) - 1:
            dx = x_lower[i + 1] - x0
            dy = y_lower[i + 1] - y0
        else:
            dx = x0 - x_lower[i - 1]
            dy = y0 - y_lower[i - 1]

        # Calculate normal (inwards/right)
        normal_dx = -dy
        normal_dy = dx

        # Normalize
        normal_length = np.sqrt(normal_dx**2 + normal_dy**2)
        normal_dx /= normal_length
        normal_dy /= normal_length

        # Initialize projection length and distance
        projection_length = 0
        old_min = np.inf
        min_distance = 99
        min_dist_idx = 0

        # Project the normal until it "intersects" the upper surface
        while min_distance < old_min:
            old_min = min_distance
            projection_length += 1e-5

            x_normal = x0 + normal_dx * projection_length
            y_normal = y0 + normal_dy * projection_length

            distances = np.sqrt((x_upper - x_normal) ** 2 + (y_upper - y_normal) ** 2)
            min_distance = np.min(distances)
            min_dist_idx = np.argmin(distances)

        min_distances[i] = old_min
        min_dist_x[i] = x_upper[min_dist_idx]

        sse += (old_min - A_inlet) ** 2

    err = sse  # + penalty (if you want to include the penalty in the error)

    if do_plot:
        print(f"Min SSE: {sse:.4f}")

        plt.figure()
        plt.plot(x_lower, y_lower, "r-")
        plt.plot(x_upper, y_upper, "b-")
        plt.title("Blade Channel Optimization")
        plt.axis("equal")
        plt.grid(True)
        plt.xlabel("Blade Coordinate [m]")
        plt.ylabel("Vertical Component [m]")
        plt.show()

    return err


def calc_and_plot_min_dist(x_lower, y_lower, x_upper, y_upper, A_inlet):
    sample_size = round(0.3 * len(x_lower))
    sample_indices = np.sort(np.random.choice(len(x_lower), sample_size, replace=False))
    x_lower = x_lower[sample_indices]
    y_lower = y_lower[sample_indices]
    print(len(x_lower))

    plt.figure()

    # Plot the lower and upper surfaces
    plt.subplot(2, 1, 1)
    plt.plot(x_lower, y_lower, 'k-', linewidth=2, label='Lower Surface')
    plt.plot(x_upper, y_upper, 'b-', linewidth=2, label='Upper Surface')

    upper_distances = np.sqrt(np.diff(x_upper)**2 + np.diff(y_upper)**2)
    intersection_threshold = np.min(upper_distances)

    min_distances = np.zeros(len(x_lower))
    min_dist_x = np.zeros(len(x_lower))

    for i in range(len(x_lower)):
        x0 = x_lower[i]
        y0 = y_lower[i]

        # Approx tangent with finite differences
        if i < len(x_lower) - 1:
            dx = x_lower[i + 1] - x0
            dy = y_lower[i + 1] - y0
        else:
            dx = x0 - x_lower[i - 1]
            dy = y0 - y_lower[i - 1]

        # Calculate normal (inwards/right)
        normal_dx = -dy
        normal_dy = dx

        # Normalize the normal vector
        normal_length = np.sqrt(normal_dx**2 + normal_dy**2)
        normal_dx /= normal_length
        normal_dy /= normal_length

        # Initialize projection length and distance
        projection_length = 0
        old_min = np.inf
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

        plt.plot([x0, x_upper[min_dist_idx]], [y0, y_upper[min_dist_idx]], 'r--')
        min_distances[i] = np.sqrt((x0 - x_upper[min_dist_idx])**2 + (y0 - y_upper[min_dist_idx])**2)
        min_dist_x[i] = x_upper[min_dist_idx]

    plt.xlabel('X Coordinate [m]')
    plt.ylabel('Y Coordinate [m]')
    plt.title('Lower and Upper Surfaces with Normals')
    plt.legend(loc='best')
    plt.grid(True)
    plt.axis('equal')

    # Plot minimum distances
    plt.subplot(2, 1, 2)
    plt.plot(min_dist_x, min_distances, 'g-o', linewidth=1.5)
    plt.xlabel('Upper Curve X [m]')
    plt.ylabel('Minimum Distance to Lower Surface [m]')
    plt.title('Minimum Distance from Upper Surface to Lower Surface')
    plt.grid(True)

    plt.show()

    return min_distances


def calc_cross_sectional_area(x_lower, y_lower, x_upper, y_upper):
    plt.figure()
    plt.plot(x_lower, y_lower, 'k-', linewidth=2, label='Lower Surface')
    plt.plot(x_upper, y_upper, 'k-', linewidth=2, label='Upper Surface')

    upper_surf_area = 0
    lower_surf_area = 0

    # Calculate area for lower surface using trapezoidal rule
    for i in range(len(x_lower) - 1):
        x0 = x_lower[i]
        y0 = y_lower[i]
        x1 = x_lower[i + 1]
        y1 = y_lower[i + 1]

        dx = x1 - x0
        dA = dx * (y1 + y0) / 2
        lower_surf_area += dA

    # Calculate area for upper surface using trapezoidal rule
    for i in range(len(x_upper) - 1):
        x0 = x_upper[i]
        y0 = y_upper[i]
        x1 = x_upper[i + 1]
        y1 = y_upper[i + 1]

        dx = x1 - x0
        dA = dx * (y1 + y0) / 2
        upper_surf_area += dA

    # Cross sectional area (in square meters)
    cross_sectional_area = upper_surf_area - lower_surf_area
    cross_sectional_area *= 10000  # Convert to cm^2

    # Fill the cross-sectional area
    plt.fill_betweenx(np.concatenate([y_upper, y_lower[::-1]]), 
                      np.concatenate([x_upper, x_lower[::-1]]), 
                      color='yellow', alpha=0.3, label=f'Cross Sectional Area = {cross_sectional_area:.2f} cm²')

    plt.xlabel('X Coordinate [m]')
    plt.ylabel('Y Coordinate [m]')
    plt.title('Cross Sectional Area')
    plt.legend(loc='best')
    plt.grid(True)
    plt.axis('equal')
    plt.show()

    return cross_sectional_area


def plot_turbine(x_lower, y_lower, x_upper, y_upper, num_blades, hub_radius, chord, blade_radius, offset):
    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')
    ax.set_box_aspect([1, 1, 1])  # Equal scaling for all axes

    r1 = hub_radius
    r2 = blade_radius

    # Generate a mesh grid for the cylinder part (hub)
    theta, Z = np.meshgrid(np.linspace(0, 2 * np.pi, 100), np.linspace(0, chord, 100))
    X_cyl = r1 * np.cos(theta)
    Y_cyl = r1 * np.sin(theta)

    # Plot the cylindrical hub
    ax.plot_surface(X_cyl, Y_cyl, Z, alpha=0.3, rstride=1, cstride=1, edgecolor='none')

    # Generate the cap of the hub
    theta_cap = np.linspace(0, 2 * np.pi, 100)
    X_cap = r1 * np.cos(theta_cap)
    Y_cap = r1 * np.sin(theta_cap)

    # Plot the top and bottom caps of the hub
    ax.fill(X_cap, Y_cap, 0, 'b', edgecolor='none')
    ax.fill(X_cap, Y_cap, chord, 'b', edgecolor='none')

    # Generate the blade geometry
    x = np.concatenate([x_lower, np.flip(x_upper)])
    y = np.concatenate([y_lower, np.flip(y_upper)]) - offset
    z = np.linspace(r1 - r2 / 4, r1 + r2, 100)

    X, Z = np.meshgrid(x, z)
    Y, Z = np.meshgrid(y, z)

    # Plot the blades
    for k in range(num_blades):
        theta_blade = 2 * np.pi * k / num_blades

        X_rotated = Y * np.cos(theta_blade) - Z * np.sin(theta_blade)
        Y_rotated = Y * np.sin(theta_blade) + Z * np.cos(theta_blade)

        # Plot the rotated blades
        ax.plot_surface(X_rotated, Y_rotated, X, facecolors=plt.cm.viridis(np.random.rand()), edgecolor='none')

    # Set plot labels and title
    ax.set_xlabel('X [m]')
    ax.set_ylabel('Y [m]')
    ax.set_zlabel('Z [m]')
    ax.set_title('Plot of Turbine')

    plt.show()


def prendtl_meyer(M, gamma):
    the_frac = np.sqrt((gamma + 1) / (gamma - 1))
    god_knows = np.atan(np.sqrt(gamma * (M^2 - 1) / (2 + gamma * (M^2 - 1))))
    literally_just_ref = np.atan(np.sqrt(M^2 - 1))
    nu = the_frac * god_knows - literally_just_ref
    return nu
