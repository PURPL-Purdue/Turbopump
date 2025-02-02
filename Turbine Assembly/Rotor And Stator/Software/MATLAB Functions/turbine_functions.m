function funs = turbine_functions
    % function that returns a struct of function handles
    funs.calc_blade_speed = @calc_blade_speed;
    funs.rotor_back_calculate = @rotor_back_calculate;
    funs.plot_velocity_triangles_angles = @plot_velocity_triangles_angles;
    funs.calculate_blade_efficiency = @calculate_blade_efficiency;
    funs.generate_blade_geom = @generate_blade_geom;
    funs.generate_blade_geom_constant_area = @generate_blade_geom_constant_area;
    funs.calc_and_plot_min_dist = @calc_and_plot_min_dist;
    funs.calc_cross_sectional_area = @calc_cross_sectional_area;
    funs.plot_turbine = @plot_turbine;
    funs.curvature_function = @curvature_function;
    funs.prendtl_meyer = @prendtl_meyer;
end

% function definitions
function U = calc_blade_speed(radius, rpm)
    % rpm to angular velocity (rad/s)
    omega = rpm * 2 * pi / 60;

    U = [radius * omega, 0];  % blade speed (m/s)
end

function [v1, v2, w, U, a1, a2, b] = rotor_back_calculate(RPM, torque, mass_flow, beta, radius, V_in)
    % Convert RPM to angular velocity (rad/s)
    U = calc_blade_speed(radius, RPM);
    U = U(1);
    
    C = torque / radius / mass_flow;
    % C = 4 * U;
    %fprintf("C param: %.3f\n", C);

    syms v2 a1 a2 w u c b v1

    % eq1 = v1 * sin(a1) == w * sin(b) + u
    % eq2 = v1 * cos(a1) == w * cos(b)
    % eq3 = v2 * sin(a2) == w * sin(b) - u
    % eq4 = v2 * cos(a2) == w * cos(b)
    % eq5 = v1 * sin(a1) - v2 * sin(a2) == c
    % eq6 = v1 * cos(a1) == v2 * cos(a2)


    eq6 = [v1 * cos(a1) ; v1 * sin(a1)] == [w * cos(b) ; w * sin(b) + u];
    eq7 = [v2 * cos(a2) ; v2 * sin(a2)] == [w * cos(-b) ; w * sin(-b) + u];
    eq8 = abs(v1 * sin(a1) - v2 * sin(a2)) == c;

    % equations = subs([eq1, eq2, eq3, eq4, eq5], [v1, b, u, c], [V_in, deg2rad(beta), U, C])
    equations = subs([eq6(1), eq6(2), eq7(1), eq7(2), eq8], [v1, u, c], [V_in, U, C]);
    % equations = subs([eq1, eq2, eq3, eq4, eq5], [v1, u, c], [V_in, U, C])

    [v2, a1, a2, w, b] = vpasolve(equations, [v2, a1, a2, w, b], [V_in; beta + 10 * pi/180; beta - 10 * pi/180; V_in; beta]);
    % [v2, a1, a2, w] = solve(equations, [v2, a1, a2, w], [v1, b + 10 * pi / 180, b + 10 * pi / 180, v1])
    
    v1 = V_in;
    v2 = double(v2);
    a1 = double(a1);
    a2 = double(a2);
    w = double(w);
    b = double(b);
    % [a1, a2, v2, w] = solve([eq6(1), eq6(2), eq7(1), eq7(2), eq8], [a1, a2, v2, w])

    % go back to abs angle measures
    a2 = -a2;
    % beta_out = -beta_out;
    
    % disp results
    %fprintf('Inlet Absolute Velocity (V_in): %.2f m/s\n', v1);
    %fprintf('Outlet Absolute Velocity (V_out): %.2f m/s\n', v2);
    %fprintf('Inlet Relative Velocity (W_in): %.2f m/s\n', w);
    %fprintf('Outlet Relative Velocity (W_out): %.2f m/s\n', w);
    %fprintf('Inlet Abs Angle (a_in): %.2f deg\n', rad2deg(a1));
    %fprintf('Outlet Abs Angle (a_out): %.2f deg\n', rad2deg(a2));
    %fprintf("Beta: %.3f\n", rad2deg(b));
    %fprintf('Turbine velocity (u): %.2f m/s\n', U);
end

function plot_velocity_triangles_angles(v1, v2, u, w1, w2, chord_length, inclination_angle, beta1, beta2, alpha1, alpha2)
    % Plot w1, u1 with the tail of u1 on the head of w1, and v1 connecting tail of w1 to head of u1
    figure;
    hold on;
    axis equal;
    grid on;

    v1 = [v1 * cos(alpha1), v1 * sin(alpha1)];
    v2 = [v2 * cos(alpha2), v2 * sin(alpha2)];
    u1 = [0, u];
    u2 = [0, u];
    w1 = [w1 * cos(beta1), w1 * sin(beta1)];
    w2 = [w2 * cos(beta2), w2 * sin(beta2)];

    mag_v1 = norm(v1);
    mag_v2 = norm(v2);
    mag_u1 = norm(u1);
    mag_u2 = norm(u2);
    mag_w1 = norm(w1);
    mag_w2 = norm(w2);
    
    % i rescale here so that the diagram is readable
    max_magnitude = max([mag_v1, mag_v2, mag_u1, mag_u2, mag_w1, mag_w2]);
    scale_factor = 2 * chord_length / max_magnitude;
    
    % Scale all vectors
    v1 = v1 * scale_factor;
    v2 = v2 * scale_factor;
    u1 = u1 * scale_factor;
    u2 = u2 * scale_factor;
    w1 = w1 * scale_factor;
    w2 = w2 * scale_factor;
    
    quiver(0, 0, w1(1), w1(2), 0, 'r', 'LineWidth', 2);
    quiver(w1(1), w1(2), u1(1), u1(2), 0, 'b', 'LineWidth', 2);
    quiver(0, 0, v1(1), v1(2), 0, 'g', 'LineWidth', 2);
    
    theta = inclination_angle;
    chord_vec = chord_length * [cos(theta), sin(theta)];
    quiver(w1(1), w1(2), chord_vec(1), chord_vec(2), 0, 'k', 'LineWidth', 2);
    
    chord_end = w1 + chord_vec;
    
    quiver(chord_end(1), chord_end(2), w2(1), w2(2), 0, 'r', 'LineWidth', 2);
    quiver(chord_end(1) + w2(1), chord_end(2) + w2(2), u2(1), u2(2), 0, 'b', 'LineWidth', 2);
    quiver(chord_end(1), chord_end(2), v2(1), v2(2), 0, 'g', 'LineWidth', 2);
    
    text(w1(1)/2, w1(2)/2, sprintf('w1 = %.2f m/s', mag_w1), 'FontSize', 12);
    text(w1(1) + u1(1)/2, w1(2) + u1(2)/2, sprintf('u1 = %.2f m/s', mag_u1), 'FontSize', 12);
    text(v1(1)/2, v1(2)/2, sprintf('v1 = %.2f m/s', mag_v1), 'FontSize', 12);
    text(chord_end(1) + w2(1)/2, chord_end(2) + w2(2)/2, sprintf('w2 = %.2f m/s', mag_w2), 'FontSize', 12);
    text(chord_end(1) + w2(1) + u2(1)/2, chord_end(2) + w2(2) + u2(2)/2, sprintf('u2 = %.2f m/s', mag_u2), 'FontSize', 12);
    text(chord_end(1) + v2(1)/2, chord_end(2) + v2(2)/2, sprintf('v2 = %.2f m/s', mag_v2), 'FontSize', 12);
    
    xlabel('X');
    ylabel('Y');
    title('Velocity Triangles with Chord Inclination');
    hold off;
end

function efficiency = calculate_blade_efficiency(mass_flow_rate, C1, C2, W1, W2, beta1, beta2, alpha1)    
    kb = W2/W1;  % single-stage, single-rotor impulse turbine
    
    % power transferred to rotor blades
    Pb = (1/4) * mass_flow_rate * C1^2 * (1 + kb * cos(beta2)/cos(beta1)) * cos(alpha1)^2;
    
    % inpute power
    input_power = (1/2) * mass_flow_rate * C1^2;
    
    % blade efficiency
    efficiency = Pb / input_power;
end

function [max_thickness, cross_sectional_area, x_lower, y_lower, x_upper, y_upper, areas] = generate_blade_geom(c, beta, A_inlet, B_spacing, blade_thickness)
    theta_l = beta * pi / 180;

    num_points = 100;
    
    figure;
    hold on;
    
    %%% Lower Surface

    % Lower Surface Guide Curves
    x1 = linspace(0, c/2, num_points);
    bottom_curve1 = tan(theta_l) * x1;
    x2 = linspace(c/2, c, num_points);
    bottom_curve2 = tan(-theta_l) * (x2 - c/2) + tan(theta_l) * (c/2);
    plot(x1, bottom_curve1, 'b--', 'LineWidth', 2);
    plot(x2, bottom_curve2, 'b--', 'LineWidth', 2);

    % Lower Surface
    r = c/2 * sin(theta_l);
    C1 = c/2 - sqrt((tan(theta_l)^2 * r^2) / (1 + tan(theta_l) ^ 2));
    C2 = c/2 + sqrt((tan(theta_l)^2 * r^2) / (1 + tan(theta_l) ^ 2));
    x1 = linspace(0, C1, num_points);
    x2 = linspace(C1, C2, num_points);
    x3 = linspace(C2, c, num_points);
    curve_DC_1 = tan(theta_l) * x1;
    curve_CC = sqrt(r^2 - (x2 - c/2) .^ 2);
    curve_DC_2 = tan(-theta_l) * (x3 - c/2) + tan(theta_l) * (c/2);
    h1 = plot(x1, curve_DC_1, 'k-', 'LineWidth', 2);
    h2 = plot(x2, curve_CC, 'k-', 'LineWidth', 2);
    h3 = plot(x3, curve_DC_2, 'k-', 'LineWidth', 2);
    x_lower = [x1, x2, x3];
    y_lower = [curve_DC_1, curve_CC, curve_DC_2];

    % Upper Surface Guide Curves
    % Area Normals
    x_upper1 = linspace(0, A_inlet * cos(theta_l), num_points);
    upper_curve1 = -1 / tan(theta_l) * x_upper1;
    x_upper2 = linspace(c - A_inlet * cos(theta_l), c, num_points);
    upper_curve2 = 1 / tan(theta_l) * (x_upper2 - c);
    plot(x_upper1, upper_curve1, 'r--', 'LineWidth', 2);
    plot(x_upper2, upper_curve2, 'r--', 'LineWidth', 2);
    
    % Quadratic Bezier Curve
    x1 = A_inlet * cos(theta_l);
    y1 = -1 / tan(theta_l) * x1;
    x2 = c / 2;
    m = (upper_curve1(end) - (-B_spacing)) / (A_inlet * cos(theta_l));
    
    % Guide Tangents
    u_x_1 = linspace(0, c/2, num_points);
    u_a_1 = m * u_x_1 - B_spacing;
    u_x_2 = linspace(c/2, c, num_points);
    u_a_2 = m * (c - u_x_2) - B_spacing;
    plot(u_x_1, u_a_1, 'r--', 'LineWidth', 2);
    plot(u_x_2, u_a_2, 'r--', 'LineWidth', 2);
    
    % Solve for y2 such that k = 1
    y2_initial_guess = r - B_spacing;
    % y2 = fzero(@(y2) curvature_function(y2, m, x1, x2, y1) - 1, y2_initial_guess);
    % c_gap = 
    % k = NaN;
    % while isnan(k)
    %     y2 = fzero(@(y2) curvature_function(y2, m, x1, x2, y1) - (1/(r-A_inlet)), y2_initial_guess);
    %     k = abs((2*x2 - 2*x_m)*(2*y1 - 2*y2))/(8*((x2 - x_m)^2)^(3/2));
    % end
    y2 = fzero(@(y2) curvature_function(y2, m, x1, x2, y1) - (1/(r-A_inlet)), y2_initial_guess); 

    y2 = r - B_spacing + blade_thickness;
    max_thickness = r;

    x_m = (y2 - y1) / m + x1;
    k = abs((2*x2 - 2*x_m)*(2*y1 - 2*y2))/(8*((x2 - x_m)^2)^(3/2));
    % fprintf("k = %.4f\n", k)
    
    % Parametric Bézier Guide curve
    t = linspace(0, 1, num_points);
    B1 = (1 - t).^2 * x1 + 2 * (1 - t) .* t * x_m + t.^2 * x2;
    B2 = (1 - t).^2 * y1 + 2 * (1 - t) .* t * y2 + t.^2 * y2;
    plot(B1, B2, 'r--', 'LineWidth', 2);
    plot(c - B1, B2, 'r--', 'LineWidth', 2);
    
    % Upper surface
    % Parametric Bézier curve
    h4 = plot(B1, B2 + B_spacing, 'k-', 'LineWidth', 2);
    h5 = plot(c - B1, B2 + B_spacing, 'k-', 'LineWidth', 2);
    
    % Upper Surface Tails
    u_x_1 = linspace(0, x1, num_points);
    u_a_1 = m * u_x_1;
    u_x_2 = linspace(c - x1, c, num_points);
    u_a_2 = m * (c - u_x_2);
    h6 = plot(u_x_1, u_a_1, 'k-', 'LineWidth', 2);
    h7 = plot(u_x_2, u_a_2, 'k-', 'LineWidth', 2);

    B1_post = flip(c - B1);
    B2_post = flip(B2);
    x_upper = [u_x_1, B1, B1_post, u_x_2];
    y_upper = [u_a_1 - B_spacing, B2, B2_post, u_a_2 - B_spacing];
    
    handles = [h1, h2, h3, h4, h5, h6, h7];
    
    % shifted_handles = gobjects(size(handles));
    for i = 1:length(handles)
        x_data = get(handles(i), 'XData');
        y_data = get(handles(i), 'YData');
        % y_shifted = y_data - shift_amount;
        plot(x_data, y_data - B_spacing, 'g--', 'LineWidth', 2);
        plot(x_data, y_data + B_spacing, 'g--', 'LineWidth', 2);
        % shifted_handles(i) = plot(x_data, y_shifted, 'g--', 'LineWidth', 2);
        % set(shifted_handles(i), 'Color', get(handles(i), 'Color'));
    end
    
    % Label and axes
    xlabel('x [m]');
    ylabel('y [m]');
    title('Blade Geomtetry Parametric Curves');
    grid on;
    axis equal;
    hold off
    
    % Save Blade to SVG
    new_fig = figure;
    new_ax = axes(new_fig);
    axis(new_ax, 'equal');
    copyobj([h1, h2, h3, h4, h5, h6, h7], new_ax);
    set(new_ax, 'Visible', 'off'); 
    set(new_ax, 'XTick', []); 
    set(new_ax, 'YTick', []);
    set(new_ax, 'Box', 'off');
    % exportgraphics(new_ax, 'blade.svg', 'ContentType', 'vector');
    saveas(new_fig, 'cgrotor.svg');

    % calculate minimum spacing
    areas = calc_and_plot_min_dist(x_lower, y_lower, x_upper, y_upper, A_inlet);

    % calculate cross sectional area
    y_upper = y_upper + B_spacing;
    cross_sectional_area = calc_cross_sectional_area(x_lower, y_lower, x_upper, y_upper);


    %% EXPORT GEOMETRY KEYPOINTS
    % Initialize table to store key points
    % key_points_table = table([], [], [], [], 'VariableNames', {'Type', 'X', 'Y', 'Additional'});

    % %%% Lower Surface

    % % Line segment DC_1
    % key_points_table = [key_points_table; table({'Line'}, [x1(1); x1(end)], [curve_DC_1(1); curve_DC_1(end)], {''; ''})];

    % % Circular arc CC
    % arc_center = [c/2, 0];
    % arc_radius = r;
    % arc_start_angle = atan2(curve_CC(1) - arc_center(2), x2(1) - arc_center(1));
    % arc_end_angle = atan2(curve_CC(end) - arc_center(2), x2(end) - arc_center(1));
    % key_points_table = [key_points_table; table({'Arc'}, arc_center(1), arc_center(2), {arc_radius, arc_start_angle, arc_end_angle})];

    % % Line segment DC_2
    % key_points_table = [key_points_table; table({'Line'}, [x3(1); x3(end)], [curve_DC_2(1); curve_DC_2(end)], {''; ''})];

    % %%% Upper Surface

    % % Bézier curve
    % key_points_table = [key_points_table; table({'Bezier'}, [x1; x_m; x2], [y1; y2; y2], {''; ''})];

    % % Upper tail segments
    % key_points_table = [key_points_table; table({'Line'}, [u_x_1(1); u_x_1(end)], [u_a_1(1); u_a_1(end)], {''; ''})];
    % key_points_table = [key_points_table; table({'Line'}, [u_x_2(1); u_x_2(end)], [u_a_2(1); u_a_2(end)], {''; ''})];

    % % Save key points to a CSV file
    % writetable(key_points_table, 'key_points.csv');
end

function [max_thickness, cross_sectional_area, x_lower, y_lower, x_upper, y_upper, areas] = generate_blade_geom_constant_area(c, beta, A_inlet, B_spacing, blade_thickness)
    theta_l = beta * pi / 180;

    num_points = 100;
    
    figure;
    hold on;
    
    %%% Lower Surface

    % Lower Surface Guide Curves
    x1 = linspace(0, c/2, num_points);
    bottom_curve1 = tan(theta_l) * x1;
    x2 = linspace(c/2, c, num_points);
    bottom_curve2 = tan(-theta_l) * (x2 - c/2) + tan(theta_l) * (c/2);
    plot(x1, bottom_curve1, 'b--', 'LineWidth', 2);
    plot(x2, bottom_curve2, 'b--', 'LineWidth', 2);

    % Lower Surface
    r = c/2 * sin(theta_l);
    C1 = c/2 - sqrt((tan(theta_l)^2 * r^2) / (1 + tan(theta_l) ^ 2));
    C2 = c/2 + sqrt((tan(theta_l)^2 * r^2) / (1 + tan(theta_l) ^ 2));
    x1 = linspace(0, C1, num_points);
    x2 = linspace(C1, C2, num_points);
    x3 = linspace(C2, c, num_points);
    curve_DC_1 = tan(theta_l) * x1;
    curve_CC = sqrt(r^2 - (x2 - c/2) .^ 2);
    curve_DC_2 = tan(-theta_l) * (x3 - c/2) + tan(theta_l) * (c/2);
    h1 = plot(x1, curve_DC_1, 'k-', 'LineWidth', 2);
    h2 = plot(x2, curve_CC, 'k-', 'LineWidth', 2);
    h3 = plot(x3, curve_DC_2, 'k-', 'LineWidth', 2);
    x_lower = [x1, x2, x3];
    y_lower = [curve_DC_1, curve_CC, curve_DC_2];

    % Upper Surface Guide Curves
    % Area Normals
    x_upper1 = linspace(0, A_inlet * cos(theta_l), num_points);
    upper_curve1 = -1 / tan(theta_l) * x_upper1;
    x_upper2 = linspace(c - A_inlet * cos(theta_l), c, num_points);
    upper_curve2 = 1 / tan(theta_l) * (x_upper2 - c);
    plot(x_upper1, upper_curve1, 'r--', 'LineWidth', 2);
    plot(x_upper2, upper_curve2, 'r--', 'LineWidth', 2);
    
    % Quadratic Bezier Curve
    x1 = A_inlet * cos(theta_l);
    y1 = -1 / tan(theta_l) * x1;
    x2 = c / 2;
    m = (upper_curve1(end) - (-B_spacing)) / (A_inlet * cos(theta_l));
    
    % Guide Tangents
    u_x_1 = linspace(0, c/2, num_points);
    u_a_1 = m * u_x_1 - B_spacing;
    u_x_2 = linspace(c/2, c, num_points);
    u_a_2 = m * (c - u_x_2) - B_spacing;
    plot(u_x_1, u_a_1, 'r--', 'LineWidth', 2);
    plot(u_x_2, u_a_2, 'r--', 'LineWidth', 2);
    
    % Solve for y2 such that k = 1
    y2 = r - B_spacing + blade_thickness;
    max_thickness = r;

    x_m = (y2 - y1) / m + x1;
    k = abs((2*x2 - 2*x_m)*(2*y1 - 2*y2))/(8*((x2 - x_m)^2)^(3/2));
    fprintf("k = %.4f\n", k)

    half_len = floor(length(x_lower) / 2);
    x_lower_first_half = x_lower(1:half_len);
    y_lower_first_half = y_lower(1:half_len)
    slope = m;
    objectiveFunction = @(vars) calc_sq_err_min_dist( ...
        x_lower_first_half, y_lower_first_half, 0, -B_spacing, vars(1), slope, x2, vars(2), A_inlet, false);
    initial_guess = [x1, y2];
    [optimal_vars, fval] = fminsearch(objectiveFunction, initial_guess);
    
    optimal_x1 = optimal_vars(1);
    % optimal_x_m = optimal_vars(2);
    % optimal_y_m = optimal_vars(3);
    optimal_y2 = optimal_vars(2);

    fprintf('Optimal x1: %.4f\n', optimal_x1);
    % fprintf('Optimal x_m: %.4f\n', optimal_x_m);
    % fprintf('Optimal y_m: %.4f\n', optimal_y_m);
    fprintf('Optimal y2: %.4f\n', optimal_y2);
    fprintf('Minimum penalty: %.4f\n', fval);
    
    % Parametric Bézier Guide curve
    t = linspace(0, 1, num_points);
    B1 = (1 - t).^2 * x1 + 2 * (1 - t) .* t * x_m + t.^2 * x2;
    B2 = (1 - t).^2 * y1 + 2 * (1 - t) .* t * y2 + t.^2 * y2;
    plot(B1, B2, 'r--', 'LineWidth', 2);
    plot(c - B1, B2, 'r--', 'LineWidth', 2);
    
    % Upper surface
    % Parametric Bézier curve
    h4 = plot(B1, B2 + B_spacing, 'k-', 'LineWidth', 2);
    h5 = plot(c - B1, B2 + B_spacing, 'k-', 'LineWidth', 2);
    
    % Upper Surface Tails
    u_x_1 = linspace(0, x1, num_points);
    u_a_1 = m * u_x_1;
    u_x_2 = linspace(c - x1, c, num_points);
    u_a_2 = m * (c - u_x_2);
    h6 = plot(u_x_1, u_a_1, 'k-', 'LineWidth', 2);
    h7 = plot(u_x_2, u_a_2, 'k-', 'LineWidth', 2);

    B1_post = flip(c - B1);
    B2_post = flip(B2);
    x_upper = [u_x_1, B1, B1_post, u_x_2];
    y_upper = [u_a_1 - B_spacing, B2, B2_post, u_a_2 - B_spacing];
    
    handles = [h1, h2, h3, h4, h5, h6, h7];
    
    % shifted_handles = gobjects(size(handles));
    for i = 1:length(handles)
        x_data = get(handles(i), 'XData');
        y_data = get(handles(i), 'YData');
        % y_shifted = y_data - shift_amount;
        plot(x_data, y_data - B_spacing, 'g--', 'LineWidth', 2);
        plot(x_data, y_data + B_spacing, 'g--', 'LineWidth', 2);
        % shifted_handles(i) = plot(x_data, y_shifted, 'g--', 'LineWidth', 2);
        % set(shifted_handles(i), 'Color', get(handles(i), 'Color'));
    end
    
    % Label and axes
    xlabel('x [m]');
    ylabel('y [m]');
    title('Blade Geomtetry Parametric Curves');
    grid on;
    axis equal;
    hold off
    
    % Save Blade to SVG
    new_fig = figure;
    new_ax = axes(new_fig);
    axis(new_ax, 'equal');
    copyobj([h1, h2, h3, h4, h5, h6, h7], new_ax);
    set(new_ax, 'Visible', 'off'); 
    set(new_ax, 'XTick', []); 
    set(new_ax, 'YTick', []);
    set(new_ax, 'Box', 'off');
    % exportgraphics(new_ax, 'blade.svg', 'ContentType', 'vector');
    saveas(new_fig, 'cgrotor.svg');

    % calculate minimum spacing
    areas = calc_and_plot_min_dist(x_lower, y_lower, x_upper, y_upper, A_inlet);

    % calculate cross sectional area
    y_upper = y_upper + B_spacing;
    cross_sectional_area = calc_cross_sectional_area(x_lower, y_lower, x_upper, y_upper);

    calc_sq_err_min_dist( ...
        x_lower_first_half, y_lower_first_half, 0, -B_spacing, optimal_x1, slope, x2, y2, A_inlet, true);


    %% EXPORT GEOMETRY KEYPOINTS
    % Initialize table to store key points
    % key_points_table = table([], [], [], [], 'VariableNames', {'Type', 'X', 'Y', 'Additional'});

    % %%% Lower Surface

    % % Line segment DC_1
    % key_points_table = [key_points_table; table({'Line'}, [x1(1); x1(end)], [curve_DC_1(1); curve_DC_1(end)], {''; ''})];

    % % Circular arc CC
    % arc_center = [c/2, 0];
    % arc_radius = r;
    % arc_start_angle = atan2(curve_CC(1) - arc_center(2), x2(1) - arc_center(1));
    % arc_end_angle = atan2(curve_CC(end) - arc_center(2), x2(end) - arc_center(1));
    % key_points_table = [key_points_table; table({'Arc'}, arc_center(1), arc_center(2), {arc_radius, arc_start_angle, arc_end_angle})];

    % % Line segment DC_2
    % key_points_table = [key_points_table; table({'Line'}, [x3(1); x3(end)], [curve_DC_2(1); curve_DC_2(end)], {''; ''})];

    % %%% Upper Surface

    % % Bézier curve
    % key_points_table = [key_points_table; table({'Bezier'}, [x1; x_m; x2], [y1; y2; y2], {''; ''})];

    % % Upper tail segments
    % key_points_table = [key_points_table; table({'Line'}, [u_x_1(1); u_x_1(end)], [u_a_1(1); u_a_1(end)], {''; ''})];
    % key_points_table = [key_points_table; table({'Line'}, [u_x_2(1); u_x_2(end)], [u_a_2(1); u_a_2(end)], {''; ''})];

    % % Save key points to a CSV file
    % writetable(key_points_table, 'key_points.csv');
end

function [max_thickness, cross_sectional_area, x_lower, y_lower, x_upper, y_upper, areas] = generate_blade_geom_constant_area_v2(c, beta, A_inlet, B_spacing, blade_thickness)
    theta_l = beta * pi / 180;

    num_points = 100;
    
    figure;
    hold on;
    
    %%% Lower Surface

    % Lower Surface Guide Curves
    x1 = linspace(0, c/2, num_points);
    bottom_curve1 = tan(theta_l) * x1;
    x2 = linspace(c/2, c, num_points);
    bottom_curve2 = tan(-theta_l) * (x2 - c/2) + tan(theta_l) * (c/2);
    plot(x1, bottom_curve1, 'b--', 'LineWidth', 2);
    plot(x2, bottom_curve2, 'b--', 'LineWidth', 2);

    % Lower Surface
    half_d = c / 2;
    lower_slope = tan(theta_l);
    r = sqrt((half_d) ^ 2 + (half_d / lower_slope)^2);
    x_lower = linspace(0, c, num_points);
    y_lower = sqrt(r^2 - (x - half_d))^2 - sqrt(S^2 - half_d^2);
    h1 = plot(x_lower, y_lower, 'k-', 'LineWidth', 2);

    % Upper Surface Guide Curves
    % Area Normals
    x_upper1 = linspace(0, A_inlet * cos(theta_l), num_points);
    upper_curve1 = -1 / tan(theta_l) * x_upper1;
    x_upper2 = linspace(c - A_inlet * cos(theta_l), c, num_points);
    upper_curve2 = 1 / tan(theta_l) * (x_upper2 - c);
    plot(x_upper1, upper_curve1, 'r--', 'LineWidth', 2);
    plot(x_upper2, upper_curve2, 'r--', 'LineWidth', 2);
    
    % Quadratic Bezier Curve
    x1 = A_inlet * cos(theta_l);
    y1 = -1 / tan(theta_l) * x1;
    x2 = c / 2;
    m = (upper_curve1(end) - (-B_spacing)) / (A_inlet * cos(theta_l));
    
    % Guide Tangents
    u_x_1 = linspace(0, c/2, num_points);
    u_a_1 = m * u_x_1 - B_spacing;
    u_x_2 = linspace(c/2, c, num_points);
    u_a_2 = m * (c - u_x_2) - B_spacing;
    plot(u_x_1, u_a_1, 'r--', 'LineWidth', 2);
    plot(u_x_2, u_a_2, 'r--', 'LineWidth', 2);
    
    % Solve for y2 such that k = 1
    y2 = r - B_spacing + blade_thickness;
    max_thickness = r;

    x_m = (y2 - y1) / m + x1;
    k = abs((2*x2 - 2*x_m)*(2*y1 - 2*y2))/(8*((x2 - x_m)^2)^(3/2));
    fprintf("k = %.4f\n", k)

    half_len = floor(length(x_lower) / 2);
    x_lower_first_half = x_lower(1:half_len);
    y_lower_first_half = y_lower(1:half_len)
    slope = m;
    objectiveFunction = @(vars) calc_sq_err_min_dist( ...
        x_lower_first_half, y_lower_first_half, 0, -B_spacing, vars(1), slope, x2, vars(2), A_inlet, false);
    initial_guess = [x1, y2];
    [optimal_vars, fval] = fminsearch(objectiveFunction, initial_guess);
    
    optimal_x1 = optimal_vars(1);
    % optimal_x_m = optimal_vars(2);
    % optimal_y_m = optimal_vars(3);
    optimal_y2 = optimal_vars(2);

    fprintf('Optimal x1: %.4f\n', optimal_x1);
    % fprintf('Optimal x_m: %.4f\n', optimal_x_m);
    % fprintf('Optimal y_m: %.4f\n', optimal_y_m);
    fprintf('Optimal y2: %.4f\n', optimal_y2);
    fprintf('Minimum penalty: %.4f\n', fval);
    
    % Parametric Bézier Guide curve
    t = linspace(0, 1, num_points);
    B1 = (1 - t).^2 * x1 + 2 * (1 - t) .* t * x_m + t.^2 * x2;
    B2 = (1 - t).^2 * y1 + 2 * (1 - t) .* t * y2 + t.^2 * y2;
    plot(B1, B2, 'r--', 'LineWidth', 2);
    plot(c - B1, B2, 'r--', 'LineWidth', 2);
    
    % Upper surface
    % Parametric Bézier curve
    h4 = plot(B1, B2 + B_spacing, 'k-', 'LineWidth', 2);
    h5 = plot(c - B1, B2 + B_spacing, 'k-', 'LineWidth', 2);
    
    % Upper Surface Tails
    u_x_1 = linspace(0, x1, num_points);
    u_a_1 = m * u_x_1;
    u_x_2 = linspace(c - x1, c, num_points);
    u_a_2 = m * (c - u_x_2);
    h6 = plot(u_x_1, u_a_1, 'k-', 'LineWidth', 2);
    h7 = plot(u_x_2, u_a_2, 'k-', 'LineWidth', 2);

    B1_post = flip(c - B1);
    B2_post = flip(B2);
    x_upper = [u_x_1, B1, B1_post, u_x_2];
    y_upper = [u_a_1 - B_spacing, B2, B2_post, u_a_2 - B_spacing];
    
    handles = [h1, h2, h3, h4, h5, h6, h7];
    
    % shifted_handles = gobjects(size(handles));
    for i = 1:length(handles)
        x_data = get(handles(i), 'XData');
        y_data = get(handles(i), 'YData');
        % y_shifted = y_data - shift_amount;
        plot(x_data, y_data - B_spacing, 'g--', 'LineWidth', 2);
        plot(x_data, y_data + B_spacing, 'g--', 'LineWidth', 2);
        % shifted_handles(i) = plot(x_data, y_shifted, 'g--', 'LineWidth', 2);
        % set(shifted_handles(i), 'Color', get(handles(i), 'Color'));
    end
    
    % Label and axes
    xlabel('x [m]');
    ylabel('y [m]');
    title('Blade Geomtetry Parametric Curves');
    grid on;
    axis equal;
    hold off
    
    % Save Blade to SVG
    new_fig = figure;
    new_ax = axes(new_fig);
    axis(new_ax, 'equal');
    copyobj([h1, h2, h3, h4, h5, h6, h7], new_ax);
    set(new_ax, 'Visible', 'off'); 
    set(new_ax, 'XTick', []); 
    set(new_ax, 'YTick', []);
    set(new_ax, 'Box', 'off');
    % exportgraphics(new_ax, 'blade.svg', 'ContentType', 'vector');
    saveas(new_fig, 'cgrotor.svg');

    % calculate minimum spacing
    areas = calc_and_plot_min_dist(x_lower, y_lower, x_upper, y_upper, A_inlet);

    % calculate cross sectional area
    y_upper = y_upper + B_spacing;
    cross_sectional_area = calc_cross_sectional_area(x_lower, y_lower, x_upper, y_upper);

    calc_sq_err_min_dist( ...
        x_lower_first_half, y_lower_first_half, 0, -B_spacing, optimal_x1, slope, x2, y2, A_inlet, true);


    %% EXPORT GEOMETRY KEYPOINTS
    % Initialize table to store key points
    % key_points_table = table([], [], [], [], 'VariableNames', {'Type', 'X', 'Y', 'Additional'});

    % %%% Lower Surface

    % % Line segment DC_1
    % key_points_table = [key_points_table; table({'Line'}, [x1(1); x1(end)], [curve_DC_1(1); curve_DC_1(end)], {''; ''})];

    % % Circular arc CC
    % arc_center = [c/2, 0];
    % arc_radius = r;
    % arc_start_angle = atan2(curve_CC(1) - arc_center(2), x2(1) - arc_center(1));
    % arc_end_angle = atan2(curve_CC(end) - arc_center(2), x2(end) - arc_center(1));
    % key_points_table = [key_points_table; table({'Arc'}, arc_center(1), arc_center(2), {arc_radius, arc_start_angle, arc_end_angle})];

    % % Line segment DC_2
    % key_points_table = [key_points_table; table({'Line'}, [x3(1); x3(end)], [curve_DC_2(1); curve_DC_2(end)], {''; ''})];

    % %%% Upper Surface

    % % Bézier curve
    % key_points_table = [key_points_table; table({'Bezier'}, [x1; x_m; x2], [y1; y2; y2], {''; ''})];

    % % Upper tail segments
    % key_points_table = [key_points_table; table({'Line'}, [u_x_1(1); u_x_1(end)], [u_a_1(1); u_a_1(end)], {''; ''})];
    % key_points_table = [key_points_table; table({'Line'}, [u_x_2(1); u_x_2(end)], [u_a_2(1); u_a_2(end)], {''; ''})];

    % % Save key points to a CSV file
    % writetable(key_points_table, 'key_points.csv');
end


function err = calc_sq_err_min_dist(x_lower, y_lower, x0, y0, x1, slope, x2, y2, A_inlet, do_plot)
    num_points = 100;
    y1 = y0 + slope * (x1 - x0);
    line_x = linspace(x0, x1, num_points);
    line_y = linspace(y0, y1, num_points);
    % slope = (y1 - y0) / (x1 - x0);
    x_m = (y2 - y1) / slope + x1;
    y_m = y2;

    t = linspace(0, 1, num_points);
    bezier_x = (1 - t).^2 * x1 + 2 * (1 - t) .* t * x_m + t.^2 * x2;
    bezier_y = (1 - t).^2 * y1 + 2 * (1 - t) .* t * y_m + t.^2 * y2;

    x_upper = [line_x, bezier_x];    y_upper = [line_y, bezier_y];

    % Numerical derivatives using finite differences
    dx_dt = gradient(bezier_x, t);  % First derivative of x
    dy_dt = gradient(bezier_y, t);  % First derivative of y
    d2x_dt2 = gradient(dx_dt, t);   % Second derivative of x
    d2y_dt2 = gradient(dy_dt, t);   % Second derivative of y

    % Curvature calculation
    curvature = (dx_dt .* d2y_dt2 - dy_dt .* d2x_dt2) ./ (dx_dt.^2 + dy_dt.^2).^(3/2);
    curvature(isnan(curvature) | isinf(curvature)) = 0;
    penalty = sqrt(var(curvature));

    if do_plot
        figure;
        plot(t, curvature, 'b-');
        xlabel('t');
        ylabel('Curvature \kappa(t)');
        title('Curvature along the Bezier Curve');

        fprintf('Penalty (variance of curvature): %.6f\n', penalty);
    end

    sample_size = round(0.3 * length(x_lower));
    sample_indices = sort(randperm(length(x_lower), sample_size));
    x_lower = x_lower(sample_indices);
    y_lower = y_lower(sample_indices);

    min_distances = zeros(1, length(x_lower));
    min_dist_x = zeros(1, length(x_lower));

    sse = 0;

    for i = 1:length(x_lower)
        x0 = x_lower(i);
        y0 = y_lower(i);

        % approx tangent w/ finite differences
        if i < length(x_lower)
            dx = x_lower(i + 1) - x0;
            dy = y_lower(i + 1) - y0;
        else
            dx = x0 - x_lower(i - 1);
            dy = y0 - y_lower(i - 1);
        end

        % calc normal (inwards/right)
        normal_dx = -dy;
        normal_dy = dx;

        % normalize
        normal_length = sqrt(normal_dx^2 + normal_dy^2);
        normal_dx = normal_dx / normal_length;
        normal_dy = normal_dy / normal_length;

        % Initialize projection length and distance
        projection_length = 0;
        old_min = inf;
        min_distance = 99;
        min_dist_idx = 1;

        % project the normal until it "intersects" the upper surface
        while min_distance < old_min
            old_min = min_distance;
            projection_length = projection_length + 1e-5;

            x_normal = x0 + normal_dx * projection_length;
            y_normal = y0 + normal_dy * projection_length;

            distances = sqrt((x_upper - x_normal).^2 + (y_upper - y_normal).^2);
            
            [min_distance, min_dist_idx] = min(distances);
        end

        min_distances(i) = old_min; %sqrt((x0 - x_upper(min_dist_idx))^2 + (y0 - y_upper(min_dist_idx))^2);
        min_dist_x(i) = x_upper(min_dist_idx);

        sse = sse + (old_min - A_inlet) ^ 2;
    end

    err = sse; % + penalty;

    if do_plot
        fprintf("Min SSE: %.4f", sse);

        figure;
        hold on;
        plot(x_lower, y_lower, 'r-');
        plot(x_upper, y_upper, 'b-');
        title("Blade Channel Optimization")
        axis equal;
        grid on;
        xlabel("Blade Coordinate [m]");
        ylabel("Vertical Component [m]");
    end
end

function [min_distances] = calc_and_plot_min_dist(x_lower, y_lower, x_upper, y_upper, A_inlet) 
    sample_size = round(0.3 * length(x_lower));
    sample_indices = sort(randperm(length(x_lower), sample_size));
    x_lower = x_lower(sample_indices);
    y_lower = y_lower(sample_indices);
    disp(length(x_lower))

    figure;
    subplot(2, 1, 1);
    hold on;

    plot(x_lower, y_lower, 'k-', 'LineWidth', 2, 'DisplayName', 'Lower Surface');
    plot(x_upper, y_upper, 'b-', 'LineWidth', 2, 'DisplayName', 'Upper Surface');

    upper_distances = sqrt(diff(x_upper).^2 + diff(y_upper).^2);
    intersection_threshold = min(upper_distances);

    min_distances = zeros(1, length(x_lower));
    min_dist_x = zeros(1, length(x_lower));

    for i = 1:length(x_lower)
        x0 = x_lower(i);
        y0 = y_lower(i);

        % approx tangent w/ finite differences
        if i < length(x_lower)
            dx = x_lower(i + 1) - x0;
            dy = y_lower(i + 1) - y0;
        else
            dx = x0 - x_lower(i - 1);
            dy = y0 - y_lower(i - 1);
        end

        % calc normal (inwards/right)
        normal_dx = -dy;
        normal_dy = dx;

        % normalize
        normal_length = sqrt(normal_dx^2 + normal_dy^2);
        normal_dx = normal_dx / normal_length;
        normal_dy = normal_dy / normal_length;

        % % project
        % projection_length = A_inlet; 
        % x_normal = x0 + normal_dx * projection_length;
        % y_normal = y0 + normal_dy * projection_length;

        % % approx intersection with min dist from curr point to curve
        % distances = sqrt((x_upper - x_normal).^2 + (y_upper - y_normal).^2);
        % [min_distances(i), idx] = min(distances);
        % min_dist_x(i) = x_upper(idx);

        % Initialize projection length and distance
        projection_length = 0;
        old_min = inf;
        min_distance = 99;
        min_dist_idx = 1;

        % project the normal until it "intersects" the upper surface
        while min_distance < old_min
            old_min = min_distance;
            projection_length = projection_length + 1e-5;

            x_normal = x0 + normal_dx * projection_length;
            y_normal = y0 + normal_dy * projection_length;

            distances = sqrt((x_upper - x_normal).^2 + (y_upper - y_normal).^2);
            
            [min_distance, min_dist_idx] = min(distances);
        end

        plot([x0, x_upper(min_dist_idx)], [y0, y_upper(min_dist_idx)], 'r--');
        min_distances(i) = sqrt((x0 - x_upper(min_dist_idx))^2 + (y0 - y_upper(min_dist_idx))^2);
        min_dist_x(i) = x_upper(min_dist_idx);
        %fprintf('Min at iteration %d is %.3f at x=%.3f\n', i, min_distances(i), min_dist_x(i));
    end

    xlabel('X Coordinate [m]');
    ylabel('Y Coordinate [m]');
    title('Lower and Upper Surfaces with Normals');
    legend('Lower Surface', 'Upper Surface', 'Location', 'best');
    grid on;
    axis equal;

    subplot(2, 1, 2);
    plot(min_dist_x, min_distances, 'g-o', 'LineWidth', 1.5);
    xlabel('Upper Curve X [m]');
    ylabel('Minimum Distance to Lower Surface [m]');
    title('Minimum Distance from Upper Surface to Lower Surface');
    grid on;

    %fprintf('Overall minimum distance is %.5f\n', min(min_distances));
end

function [cross_sectional_area] = calc_cross_sectional_area(x_lower, y_lower, x_upper, y_upper) 
    figure;
    hold on;

    plot(x_lower, y_lower, 'k-', 'LineWidth', 2, 'DisplayName', 'Lower Surface');
    plot(x_upper, y_upper, 'k-', 'LineWidth', 2, 'DisplayName', 'Upper Surface');

    upper_surf_area = 0;
    lower_surf_area = 0;

    for i = 1:(length(x_lower) - 1)
        x0 = x_lower(i);
        y0 = y_lower(i);
        x1 = x_lower(i+1);
        y1 = y_lower(i+1);

        dx = x1 - x0;
        dA = dx * (y1 + y0) / 2;
        lower_surf_area = lower_surf_area + dA;
    end

    for i = 1:(length(x_upper) - 1)
        x0 = x_upper(i);
        y0 = y_upper(i);
        x1 = x_upper(i+1);
        y1 = y_upper(i+1);

        dx = x1 - x0;
        dA = dx * (y1 + y0) / 2;
        upper_surf_area = upper_surf_area + dA;
    end

    cross_sectional_area = upper_surf_area - lower_surf_area;
    cross_sectional_area = cross_sectional_area * 10000;

    fill([x_upper, fliplr(x_lower)], [y_upper, fliplr(y_lower)], 'y', 'FaceAlpha', 0.3); % plot cross sectional area

    xlabel('X Coordinate [m]');
    ylabel('Y Coordinate [m]');
    title('Cross Sectional Area');
    legend('Lower Surface', 'Upper Surface', "Cross Sectional Area = " + cross_sectional_area + " cm^2", 'Location', 'best');
    grid on;
    axis equal;

    %fprintf('Cross Sectional Area %.5f cm^2\n', cross_sectional_area);
end

function plot_turbine(x_lower, y_lower, x_upper, y_upper, num_blades, hub_radius, chord, blade_radius, offset)
    figure;
    hold on;

    r1 = hub_radius; 
    r2 = blade_radius; 

    [theta, Z] = meshgrid(linspace(0, 2*pi, 100), linspace(0, chord, 100)); % Adjust height range for cylinder
    X_cyl = r1 * cos(theta);
    Y_cyl = r1 * sin(theta);

    surf(X_cyl, Y_cyl, Z, 'FaceAlpha', 0.3, 'EdgeColor', 'none');

    theta_cap = linspace(0, 2*pi, 100);
    X_cap = r1 * cos(theta_cap);
    Y_cap = r1 * sin(theta_cap);

    fill3(X_cap, Y_cap, zeros(size(X_cap)), 'b', 'EdgeColor', 'none');
    % 'FaceAlpha', 0.5
    fill3(X_cap, Y_cap, chord * ones(size(X_cap)), 'b', 'EdgeColor', 'none');

    x = [x_lower, fliplr(x_upper)];
    y = [y_lower, fliplr(y_upper)] - offset;
    z = linspace(r1 - r2/4, r1 + r2, 100);

    [X, Z] = meshgrid(x, z);
    [Y, Z] = meshgrid(y, z);

    % surf(Y, Z, X, 'FaceColor', 'interp', 'EdgeColor', 'none');

    for k = 1:num_blades
        theta = 2 * pi * k / num_blades;

        X_rotated = Y * cos(theta) - Z * sin(theta);
        Y_rotated = Y * sin(theta) + Z * cos(theta);
        
        surf(X_rotated, Y_rotated, X, 'FaceColor', 'interp', 'EdgeColor', 'none');

        % x_rot = y * cos(theta) - z * sin(theta);
        % y_rot = y * sin(theta) + z * cos(theta);
        % fill3(x_rot, y_rot, x)
        % surf(X_rotated, Y_rotated, X - r2 / 2, 'FaceColor', 'interp', 'EdgeColor', 'none'); 
    end

    hold off;
    axis equal;
    xlabel('X [m]'); 
    ylabel('Y [m]'); 
    zlabel('Z [m]');
    title('Plot of turbine');

end

function k = curvature_function(y2, m, x1, x2, y1)
    x_m = (y2 - y1) / m + x1;

    k = abs((2*x2 - 2*x_m)*(2*y1 - 2*y2))/(8*((x2 - x_m)^2)^(3/2));
end

function nu=prendtl_meyer(M, gamma) 
    the_frac = sqrt((gamma + 1) / (gamma - 1));
    god_knows = atan(sqrt(gamma * (M^2 - 1) / (2 + gamma * (M^2 - 1))));
    literally_just_ref = atan(sqrt(M^2 - 1));
    nu = the_frac * god_knows - literally_just_ref;
end


