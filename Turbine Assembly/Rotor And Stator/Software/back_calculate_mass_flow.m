horse_power = 100;
rpm = 50000;
radius = 0.0375;
torque = horse_power * 5252 / rpm * 1.355817; % N*m
a_in = 35.63709471;
a_out = 23.12163596;
Beta = 29.76722393;
u = 196.3495408;

data = readtable("HotGas-nozzle_values.csv");


M_e = calc_mach_from_

[v1, v2, w, U, a1, a2, b] = rotorBackCalculate(rpm, torque, radius, Beta, u, v_guess, a_in, a_out);

function [v1, v2, w, U, a1, a2, b] = rotorBackCalculate(RPM, torque, radius, beta, U, v_guess, a_in, a_out)
    syms v2 a1 a2 w u c b v1 m_dot

    % eq1 = [v1 * cos(a1) ; v1 * sin(a1)] == [w * cos(b) ; w * sin(b) + u]
    % eq2 = [v2 * cos(a2) ; v2 * sin(a2)] == [w * cos(-b) ; w * sin(-b) + u]
    % eq3 = abs(v1 * sin(a1) - v2 * sin(a2)) == torque / radius / m_dot

    % % equations = subs([eq1, eq2, eq3, eq4, eq5], [v1, b, u, c], [V_in, deg2rad(beta), U, C])
    % equations = subs([eq1(1), eq1(2), eq2(1), eq2(2), eq3], [a1, a2, u, b], [a_in, a_out, U, beta])
    % % equations = subs([eq1, eq2, eq3, eq4, eq5], [v1, u, c], [V_in, U, C])

    % [v1, v2, w, m_dot] = vpasolve(equations, [v1, v2, w, m_dot], [v_guess; v_guess; v_guess; 4])
    % % [v2, a1, a2, w] = solve(equations, [v2, a1, a2, w], [v1, b + 10 * pi / 180, b + 10 * pi / 180, v1])

    eq4 = v1 * cos(a1) / cos(b) == w;
    eq5 = (v1 * sin(a1) - u) / sin(b) == w;

    eqns = subs([eq4, eq5], [a1, b, u], [a_in, beta, U]);
    [v1, w] = vpasolve(eqns, [v1, w], [v_guess, v_guess]);

    v1 = double(v1)
    w = double(w)
    m_dot = torque / (radius * abs(2 * w * sin(beta)))
    
    % v1 = double(v1)
    % v2 = double(v2)
    % a1 = double(a1)
    % a2 = double(a2)
    % w = double(w)
    % m_dot = double(m_dot)

    % % [a1, a2, v2, w] = solve([eq6(1), eq6(2), eq7(1), eq7(2), eq8], [a1, a2, v2, w])

    % % go back to abs angle measures
    % a2 = -a2;
    % % beta_out = -beta_out;
    
    % % disp results
    % fprintf('Inlet Absolute Velocity (V_in): %.2f m/s\n', v1);
    % fprintf('Outlet Absolute Velocity (V_out): %.2f m/s\n', v2);
    % fprintf('Inlet Relative Velocity (W_in): %.2f m/s\n', w);
    % fprintf('Outlet Relative Velocity (W_out): %.2f m/s\n', w);
    % fprintf('Inlet Abs Angle (a_in): %.2f deg\n', rad2deg(a1));
    % fprintf('Outlet Abs Angle (a_out): %.2f deg\n', rad2deg(a2));
    % fprintf("Beta: %.3f\n", rad2deg(b));
    % fprintf('Turbine velocity (u): %.2f m/s\n', U);
end