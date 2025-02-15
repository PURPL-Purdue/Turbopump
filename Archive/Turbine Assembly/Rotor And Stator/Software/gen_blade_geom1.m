clear;
clc;

% main
% load('../../params.mat')

c = 5.9;
beta = 52;
A_inlet = 1.4;
B_spacing = 2.6;

generate_blade_geom(c, beta, A_inlet, B_spacing)

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
    y2_initial_guess = 1.6;
    y2 = fzero(@(y2) curvature_function(y2, m, x1, x2, y1) - 1, y2_initial_guess);

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
    saveas(new_fig, 'rotor.svg');
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