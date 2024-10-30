% TODO: add mass flow optimization
% TODO: prendtl-meyer based US  
% TODO: spans 

%% COLD GAS CONSTANTS

P_0 = 3447378.6466; % [N/m^2, 500psi]
P_e = 101352.93221; % [N/m^2, 14.7psi]
P_A = P_e;
m_dot = 0.05; % [kg/s]
T_0 = 250; % [K]
gamma = 1.4; % (Specific heat ratio of Air)
R = 8.3145; % [J/(mol*K)] (Universal Gas Constant)
m_m = 28.96; % [g/mol] (Molar Mass of Air)
n = 2; % (number of nozzles)

%% TURBINE CONSTANTS

rotor_radius = 0.04; % [m]
hub_radius = 0.02; % [m]
mass_flow_n2 = m_dot; % [kg/s]
shaft_power = 10; % [kW]
gamma_n2 = 1.4; % specific heat ratio

turbine_rpm = 10000; % [rpm]

%% Conversion Factors
kw_to_hp = 1.34102209; % [hp/kw]

% Input parameters
radius = (rotor_radius + hub_radius) / 2; % [m]
horse_power = shaft_power * kw_to_hp; % [HP]
torque = 746 * horse_power/(2 * pi * turbine_rpm/60); % N*m
degree_of_reaction = 0;
num_blades = 10;
chord = 0.02; % [m]
blade_spacing = 2 * pi * hub_radius / num_blades; % [m]
beta = 60; % initial estimate [deg] - refined in vtriangle
% V_in = 564; % [m/2]

%% COLD GAS NOZZLE CALCULATIONS

R_S = calc_R_S(R, m_m) % The actual value is 296.1, idk why its giving a different number
rho_0 = calc_rho_0(P_0, R_S, T_0) % [kg/m^3]
v_e = calc_v_e(T_0, R_S, gamma, P_e, P_0) % [m/s]


T_throat = calc_T_throat(T_0, gamma) % [K]
P_throat = calc_P_throat(P_0, gamma) % [N/m^2]
rho_throat = calc_rho_throat(P_throat, R_S, T_throat) % [kg/m^3]
v_throat = calc_v_throat(gamma, R_S, T_throat) % [m/s]
A_throat = calc_A_throat(m_dot, rho_throat, v_throat) % [m^2]


M_e = calc_M_e(P_e, P_0, gamma)
rho_e = calc_rho_e(rho_0, M_e, gamma) % [kg/m^3]
A_e = calc_A_e(m_dot, rho_e, v_e) % [m^2]
T_e = calc_T_e(T_0, gamma, M_e) % [K]


r_throat = calc_r_throat(A_throat) % [m]
r_e = calc_r_e(A_e) % [m]


dist = calc_dist(r_throat, r_e) % [m]
F_thrust = calc_F_thrust(m_dot, v_e) % [N]


A_throat_n = calc_A_throat_n(A_throat, n) % [m^2]
A_e_n = calc_A_e_n(A_e, n) % [m^2]
r_throat_n = calc_r_throat_n(A_throat_n) % [m]
r_e_n = calc_r_e_n(A_e_n) % [m]
dist_n = calc_dist_n(r_throat_n, r_e_n) % [m]

%% PLOT NOZZLE GEOMETRY

[X,Y,Z] = plot_nozzle(A_e, A_throat, A_e, dist_n, dist_n);


%% TURBINE CALCULATONS

min_blade_thickness = 0.005; % [m]
inlet_area = 0.005; % [m]

[v1, v2, w, u, a1, a2, b] = rotorBackCalculate(turbine_rpm, torque / num_blades, mass_flow_n2, deg2rad(beta), radius, v_e);
plot_velocity_triangles_angles(v1, v2, u, w, w, chord, 0, b, -b, a1, -a2);

efficiency = calculate_blade_efficiency(mass_flow_n2, v1, v2, w, w, b, b, a1);
fprintf("isentropic efficiency: %.4f\n", efficiency)

[max_blade_thickness, cross_sectional_area] = generate_blade_geom(chord, rad2deg(b), inlet_area, blade_spacing, min_blade_thickness)

%% STRUCTURE CALCULATIONS

% These values are random rough values %
w = turbine_rpm / 60 * 2 * pi % [rad/s] (Conversion of 10,000 rpm to rads)
height_blade = max_blade_thickness % [m]
rho_blade = 0.845 % [kg/m^3]
Length_blade = chord % [m]
height_bmin = min_blade_thickness % [m]
width_blade = rotor_radius - hub_radius % [m]
Z_blade = num_blades % [Number of blades]
mean_diameter = radius * 2; % [m] mean blade diam

radius_turbine = calc_radius_turbine(height_blade, hub_radius) % [m]
mass_blade = calc_mass_blade(Length_blade, width_blade, height_bmin, rho_blade) % [kg]
Force_centrifugal = calc_Force_centrifugal(mass_blade, w, radius_turbine) % [N]
stress_centrifugal = calc_stress_centrifugal(rho_blade, height_blade, mean_diameter, w) % [N/m^2]
Force_tangential = calc_Force_tangential(mass_flow_n2, w, b, w, -b) % [N]
Force_axial = calc_Force_axial(mass_flow_n2, v1, a1, w, -b) % [N]
torque_blade = calc_torque_blade(Force_tangential, height_blade) % [Nm]
torque_turbine = calc_torque_turbine(Force_tangential, radius_turbine, Z_blade) % [Nm]
P = calc_P(torque_turbine, w) % [Nm/s]
Force_gas = calc_Force_gas(Force_tangential, Force_axial) % [N]
Moment_Bending = calc_Moment_Bending(height_blade, Z_blade, Force_gas) % [Nm]
I = calc_I(Length_blade, height_bmin) % [m^4]
stress_gas = calc_stress_gas(height_bmin, Force_gas, width_blade, I) % [N/m^2]

%% STRUCTURE FUNCTIONS

function radius_turbine = calc_radius_turbine(height_blade, radius_hub) % [m]
    radius_turbine = 0.5*height_blade+radius_hub;
end
function mass_blade = calc_mass_blade(Length_blade, width_blade, height_bmin, rho_blade) % [kg]
    mass_blade = Length_blade*width_blade*height_bmin*rho_blade;
end
function Force_centrifugal = calc_Force_centrifugal(mass_blade, w, radius_turbine) % [N]
    Force_centrifugal = mass_blade*(w^2)*radius_turbine;
end
function stress_centrifugal = calc_stress_centrifugal(radius_turbine, rho_blade, height_blade, w) % [N/m^2]
    stress_centrifugal = radius_turbine*rho_blade*height_blade*(w^2);
end
function Force_tangential = calc_Force_tangential(m_dot, V_1, beta_1, V_2, beta_2) % [N] [Wants angles in radians]
    Force_tangential = m_dot*(V_1*cos(beta_1)+V_2*cos(beta_2));
end
function Force_axial = calc_Force_axial(m_dot, C_1, alpha_1, V_2, beta_2) % [N]
    Force_axial = m_dot*(C_1*sin(alpha_1)+V_2*sin(beta_2));
end
function torque_blade = calc_torque_blade(Force_tangential, height_blade) % [Nm]
    torque_blade = 0.5*Force_tangential*height_blade;
end
function torque_turbine = calc_torque_turbine(Force_tangential, radius_turbine, Z_blade) % [Nm]
    torque_turbine = Force_tangential*radius_turbine*Z_blade;
end
function P = calc_P(torque_turbine, w) % [Nm/s]
    P = torque_turbine*w;
end
function Force_gas = calc_Force_gas(Force_tangential, Force_axial) % [N]
    Force_gas = sqrt((Force_tangential^2)+(Force_axial^2));
end
function Moment_Bending = calc_Moment_Bending(height_blade, Z_blade, Force_gas) % [Nm]
    Moment_Bending = (height_blade/(2*Z_blade))*Force_gas;
end
function I = calc_I(Length_blade, height_bmin) % [m^4]
    I = (1/12)*(Length_blade^3)*(height_bmin);
end
function stress_gas = calc_stress_gas(height_bmin, Force_gas, width_b, I) % [N/m^2]
    stress_gas = (0.5*height_bmin*Force_gas*width_b)/I;
end

%% NOZZLE FUNCTIONS

function R_S = calc_R_S(R, m_m) % [J/(kg*K)] (Specific Gas Constant)
    R_S = (R/m_m)*1000;
end
function rho_0 = calc_rho_0(P_0, R_S, T_0) % [kg/m^3]
    rho_0 = P_0/(R_S*T_0);
end
function v_e = calc_v_e(T_0, R_S, gamma, P_e, P_0) % [m/s]
    v_e = sqrt((T_0*R_S)*((2*gamma)/(gamma-1))*(1-(P_e/P_0)^((gamma-1)/gamma)));
end
function T_throat = calc_T_throat(T_0, gamma) % [K]
    T_throat = T_0/(1+((gamma-1)/2));
end
function P_throat = calc_P_throat(P_0, gamma) % [N/m^2]
    P_throat = P_0/((1+((gamma-1)/2)^((gamma-1)/gamma)));
end
function rho_throat = calc_rho_throat(P_throat, R_S, T_throat) % [kg/m^3]
    rho_throat = P_throat/(R_S*T_throat);
end
function v_throat = calc_v_throat(gamma, R_S, T_throat) % [m/s]
    v_throat = sqrt((gamma*R_S*T_throat));
end
function A_throat = calc_A_throat(m_dot, rho_throat, v_throat) % [m^2]
    A_throat = m_dot/(rho_throat*v_throat);
end

function M_e = calc_M_e(P_e, P_0, gamma)
    P_ratio = P_e/P_0;
    numerator = exp(log(P_ratio)/(-gamma/(gamma-1)))-1;
    denominator = (gamma-1)/2;
    M_e = sqrt(numerator/denominator);
end
function rho_e = calc_rho_e(rho_0, M_e, gamma) % [kg/m^3]
    demon = (1+((gamma-1)/2)*M_e^2);
    rho_e = rho_0*demon^(-1/(gamma-1));
end
function A_e = calc_A_e(m_dot, rho_e, v_e) % [m^2]
    A_e = m_dot/(rho_e*v_e);
end
function T_e = calc_T_e(T_0, gamma, M_e) % [K]
    demon = (1+((gamma-1)/2)*M_e^2);
    T_e = T_0/demon;
end
function r_throat = calc_r_throat(A_throat) % [m]
    r_throat = sqrt(A_throat/pi);
end
function r_e = calc_r_e(A_e) % [m]
    r_e = sqrt(A_e/pi);
end
function dist = calc_dist(r_throat, r_e) % [m]
    dist = (r_e-r_throat)/(tand(15));
end
function A_throat_n = calc_A_throat_n(A_throat, n) % [m^2]
    A_throat_n = A_throat/n;
end
function A_e_n = calc_A_e_n(A_e, n) % [m^2]
    A_e_n = A_e/n;
end
function r_throat_n = calc_r_throat_n(A_throat_n) % [m]
    r_throat_n = sqrt(A_throat_n/pi);
end
function r_e_n = calc_r_e_n(A_e_n) % [m]
    r_e_n = sqrt(A_e_n/pi);
end
function dist_n = calc_dist_n(r_throat_n, r_e_n) % [m]
    dist_n = (r_e_n-r_throat_n)/(tand(15));
end
function F_thrust = calc_F_thrust(m_dot, v_e) % [N]
    F_thrust = m_dot*v_e;
end


%% NOZZLE PLOTTING FUNCTIONS

function [X,Y,Z] = plot_nozzle(A_inlet, A_throat, A_exit, inlet_len, outlet_len)
    figure;
    num_points = 100;
    r_inlet = sqrt(A_inlet/pi);
    r_throat = sqrt(A_throat/pi);
    r_exit = sqrt(A_exit/pi);
    % outlet_len = (r_exit - r_throat) * tand(angle);
    throat_len = 0.1 * inlet_len;
    x_total = inlet_len + throat_len + outlet_len;
    x_step = x_total / num_points;
    [X, Thetas] = meshgrid(0:x_step:x_total, 0:(2*pi/1000):2*pi);
    R = (r_inlet + (r_throat - r_inlet) .* exp_scale(X ./ inlet_len)) .* (X <= inlet_len);
    R = R + r_throat .* (inlet_len < X & X <= (inlet_len + throat_len));
    R = R + (r_throat + ((r_exit - r_throat)/outlet_len .* (X - (inlet_len + throat_len)) )) .* ((inlet_len + throat_len) < X & X <= x_total);
    Y = R .* cos(Thetas);
    Z = R .* sin(Thetas);
    C = X;
    surf(X,Y,Z,C)
    shading interp
    axis square
    colorbar
end

function percentage=exp_scale(x)
    percentage = (exp(x) - 1) ./ (exp(1) - 1);
end

function save_to_stl(X, Y, Z, filename)
    [faces, vertices] = surf2patch(X, Y, Z, 'triangles');
    tr = triangulation(faces, vertices);
    stlwrite(tr, filename);    
    disp(['Surface exported to ', filename]);
end

%% TURBINE GENERATION FUNCTIONS

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


function [max_thickness, cross_sectional_area] = generate_blade_geom(c, beta, A_inlet, B_spacing, blade_thickness)
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
    y2 = r - B_spacing + blade_thickness
    max_thickness = r;

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

    % calculate minimum spacing
    calc_and_plot_min_dist(x_lower, y_lower, x_upper, y_upper, A_inlet);

    % calculate cross sectional area
    cross_sectional_area = calc_cross_sectional_area(x_lower, y_lower, x_upper, y_upper + B_spacing);
end

function calc_and_plot_min_dist(x_lower, y_lower, x_upper, y_upper, A_inlet) 
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
        fprintf('Min at iteration %d is %.3f at x=%.3f\n', i, min_distances(i), min_dist_x(i));
    end

    xlabel('X Coordinate');
    ylabel('Y Coordinate');
    title('Lower and Upper Surfaces with Normals');
    legend('Lower Surface', 'Upper Surface', 'Location', 'best');
    grid on;
    axis equal;

    subplot(2, 1, 2);
    plot(min_dist_x, min_distances, 'g-o', 'LineWidth', 1.5);
    xlabel('Upper Curve X');
    ylabel('Minimum Distance to Lower Surface');
    title('Minimum Distance from Upper Surface to Lower Surface');
    grid on;

    fprintf('Overall minimum distance is %.5f\n', min(min_distances));
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

    fill([x_upper, fliplr(x_lower)], [y_upper, fliplr(y_lower)], 'y', 'FaceAlpha', 0.3); % plot cross sectional area

    xlabel('X Coordinate');
    ylabel('Y Coordinate');
    title('Lower and Upper Surfaces Area');
    legend('Lower Surface', 'Upper Surface', "Cross Sectional Area = " + cross_sectional_area, 'Location', 'best');
    grid on;
    axis equal;

    fprintf('Cross Sectional Area %.5f\n', cross_sectional_area);
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