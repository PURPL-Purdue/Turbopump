function [v1, v2, w1, w2, a1, a2, U] = rotorBackCalculate(RPM, torque, mass_flow, beta_in, beta_out, radius)
    % Convert RPM to angular velocity (rad/s)
    U = calc_blade_speed(radius, RPM);
    U = U(1);
    
    C = torque / radius / mass_flow;
    % C = 4 * U;
    fprintf("C param: %.3f\n", C);
    fprintf("Betas: %.3f, %.3f\n", beta_in, beta_out);

    % beta_out = -beta_out;

    % w1 = (C) / (sin(beta_in) + tan(beta_out) * cos(beta_in)) + 3 * U;
    w1 = (C) / (tan(beta_out) * cos(beta_in) - sin(beta_in)) + 3 * U;
    a1 = atan2((w1 * sin(beta_in) + U) , (w1 * cos(beta_in)));
    v1 = w1 * cos(beta_in) / cos(a1);

    % w2 = -(v1 * sin(a1) - C + U) / sin(beta_out);
    w2 = w1 * cos(beta_in) / cos(beta_out);
    a2 = atan2((w2 * sin(beta_out) - U) , (w2 * cos(beta_out)));
    v2 = w2 * cos(beta_out) / cos(a2);

    % go back to abs angle measures
    % a2 = -a2;
    % beta_out = -beta_out;
    
    % disp results
    fprintf('Inlet Absolute Velocity (V_in): %.2f m/s\n', v1);
    fprintf('Outlet Absolute Velocity (V_out): %.2f m/s\n', v2);
    fprintf('Inlet Relative Velocity (W_in): %.2f m/s\n', w1);
    fprintf('Outlet Relative Velocity (W_out): %.2f m/s\n', w2);
    fprintf('Inlet Relative Angle (a_in): %.2f deg\n', rad2deg(a1));
    fprintf('Outlet Relative Angle (a_out): %.2f deg\n', rad2deg(a2));
    fprintf('Turbine velocity (u): %.2f m/s\n', U);
end