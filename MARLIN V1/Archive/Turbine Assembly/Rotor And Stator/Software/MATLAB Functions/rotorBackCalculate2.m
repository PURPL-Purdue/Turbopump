function [v1, v2, w1, w2, a1, a2, U] = rotorBackCalculate2(RPM, torque, mass_flow, beta_in, beta_out, alpha_in, radius)
    % Convert RPM to angular velocity (rad/s)
    U = calc_blade_speed(radius, RPM);
    U = U(1);
    
    C = torque / radius / mass_flow;
    % C = 4 * U;
    fprintf("C param: %.3f\n", C);
    fprintf("Betas: %.3f, %.3f\n", beta_in, beta_out);

    % calculate w1, a1, and v1
    w1 = U / (cos(beta_in) * tan(alpha_in) - sin(beta_in));
    a1 = alpha_in;
    v1 = w1 * cos(beta_in) / cos(alpha_in);

    % calculate w2, a2, and v2
    w2 = abs((C + v1 * sin(alpha_in) - U) / sin(beta_out));
    a2 = atan2((w2 * sin(beta_out) + U) , (w2 * cos(beta_out)));
    v2 = w2 * cos(beta_out) / cos(a2);
    
    % disp results
    fprintf('Inlet Absolute Velocity (V_in): %.2f m/s\n', v1);
    fprintf('Outlet Absolute Velocity (V_out): %.2f m/s\n', v2);
    fprintf('Inlet Relative Velocity (W_in): %.2f m/s\n', w1);
    fprintf('Outlet Relative Velocity (W_out): %.2f m/s\n', w2);
    fprintf('Inlet Relative Angle (a_in): %.2f deg\n', rad2deg(alpha_in));
    fprintf('Outlet Relative Angle (a_out): %.2f deg\n', rad2deg(a2));
    fprintf('Turbine velocity (u): %.2f m/s\n', U);
end