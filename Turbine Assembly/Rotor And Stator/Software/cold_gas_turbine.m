clear;
clc;

% main
% load('../../params.mat')

% Y down is progression of turbine stages
% Turbines rotate X ->

%% Cold Gas Testing Constants:
rotor_radius = 0.04; % [m]
hub_radius = 0.02; % [m]
mass_flow_n2 = 0.05; % [kg/s]
shaft_power = 10; % [kW]
gamma_n2 = 1.4; % specific heat ratio

turbine_rpm = 10000; % [rpm]
input_velocity = NaN; % TODO: From Andrew
nozzle_radius = NaN; % TODO: From Andrew
% TODO: Andrew's constants for nozzle:

%% Conversion Factors
kw_to_hp = 1.34102209; % [hp/kw]

% Input parameters
radius = (rotor_radius + hub_radius) / 2; % meters
horse_power = shaft_power * kw_to_hp % HP
torque = 746 * horse_power/(2 * pi * turbine_rpm/60) % N*m
degree_of_reaction = 0;
num_blades = 10;
chord = 0.02; % [m]
blade_spacing = 2 * pi * hub_radius / num_blades; % [m]
V_in = 564; % [m/2]
beta = 60; % deg
stagger_angle = 0; % [deg]

[v1, v2, w, u, a1, a2, b] = rotorBackCalculate(turbine_rpm, torque / num_blades, mass_flow_n2, deg2rad(beta), radius, V_in);
plot_velocity_triangles_angles(v1, v2, u, w, w, chord, 0, b, -b, a1, -a2);

efficiency = calculate_blade_efficiency(mass_flow_n2, v1, v2, w, w, b, b, a1);
fprintf("isentropic efficiency: %.4f\n", efficiency)  

% v_i = prendtl_meyer()

generate_blade_geom(chord, rad2deg(b), 0.005, blade_spacing)
% generate_blade_geom(0.04, 60, 0.01, 0.0126)

function U = calc_blade_speed(radius, rpm)
    % rpm to angular velocity (rad/s)
    omega = rpm * 2 * pi / 60;

    U = [radius * omega, 0];  % blade speed (m/s)
end

function [v1, v2, w, U, a1, a2, b] = rotorBackCalculate(RPM, torque, mass_flow, beta, radius, V_in)
    % Convert RPM to angular velocity (rad/s)
    U = calc_blade_speed(radius, RPM);
    U = U(1);
    
    C = torque / radius / mass_flow;
    % C = 4 * U;
    fprintf("C param: %.3f\n", C);

    syms v2 a1 a2 w u c b v1

    eq1 = v1 * sin(a1) == w * sin(b) + u
    eq2 = v1 * cos(a1) == w * cos(b)
    eq3 = v2 * sin(a2) == w * sin(b) - u
    eq4 = v2 * cos(a2) == w * cos(b)
    eq5 = v1 * sin(a1) - v2 * sin(a2) == c
    % eq6 = v1 * cos(a1) == v2 * cos(a2)


    eq6 = [v1 * cos(a1) ; v1 * sin(a1)] == [w * cos(b) ; w * sin(b) + u]
    eq7 = [v2 * cos(a2) ; v2 * sin(a2)] == [w * cos(-b) ; w * sin(-b) + u]
    eq8 = abs(v1 * sin(a1) - v2 * sin(a2)) == c

    % equations = subs([eq1, eq2, eq3, eq4, eq5], [v1, b, u, c], [V_in, deg2rad(beta), U, C])
    equations = subs([eq6(1), eq6(2), eq7(1), eq7(2), eq8], [v1, u, c], [V_in, U, C])
    % equations = subs([eq1, eq2, eq3, eq4, eq5], [v1, u, c], [V_in, U, C])

    [v2, a1, a2, w, b] = vpasolve(equations, [v2, a1, a2, w, b], [V_in; beta + 10 * pi/180; beta - 10 * pi/180; V_in; beta])
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
    fprintf('Inlet Absolute Velocity (V_in): %.2f m/s\n', v1);
    fprintf('Outlet Absolute Velocity (V_out): %.2f m/s\n', v2);
    fprintf('Inlet Relative Velocity (W_in): %.2f m/s\n', w);
    fprintf('Outlet Relative Velocity (W_out): %.2f m/s\n', w);
    fprintf('Inlet Abs Angle (a_in): %.2f deg\n', rad2deg(a1));
    fprintf('Outlet Abs Angle (a_out): %.2f deg\n', rad2deg(a2));
    fprintf("Beta: %.3f\n", rad2deg(b));
    fprintf('Turbine velocity (u): %.2f m/s\n', U);
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
    
    text(w1(1)/2, w1(2)/2, 'w1', 'FontSize', 12);
    text(w1(1) + u1(1)/2, w1(2) + u1(2)/2, 'u1', 'FontSize', 12);
    text(v1(1)/2, v1(2)/2, 'v1', 'FontSize', 12);
    text(chord_end(1) + w2(1)/2, chord_end(2) + w2(2)/2, 'w2', 'FontSize', 12);
    text(chord_end(1) + w2(1) + u2(1)/2, chord_end(2) + w2(2) + u2(2)/2, 'u2', 'FontSize', 12);
    text(chord_end(1) + v2(1)/2, chord_end(2) + v2(2)/2, 'v2', 'FontSize', 12);
    
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


function generate_blade_geom(c, beta, A_inlet, B_spacing)
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
    y2 = r - B_spacing + 0.005

    x_m = (y2 - y1) / m + x1;
    k = abs((2*x2 - 2*x_m)*(2*y1 - 2*y2))/(8*((x2 - x_m)^2)^(3/2));
    fprintf("k = %.4f\n", k)
    
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
    xlabel('x');
    ylabel('y');
    title('Blade Geomtetry Parametric Curves');
    grid on;
    axis equal;
    hold off
    
    % Save Blade to SVG
    new_fig = figure;
    new_ax = axes(new_fig);
    copyobj([h1, h2, h3, h4, h5, h6, h7], new_ax);
    set(new_ax, 'Visible', 'off'); 
    set(new_ax, 'XTick', []); 
    set(new_ax, 'YTick', []);
    set(new_ax, 'Box', 'off');
    % exportgraphics(new_ax, 'blade.svg', 'ContentType', 'vector');
    saveas(new_fig, 'cgrotor.svg');
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