P_0 = 3447378.6466; % [N/m^2, 500psi]
P_e = 101352.93221; % [N/m^2, 14.7psi]
P_A = P_e;

T_0 = 293; % [K]
gamma = 1.4; % (Specific heat ratio of Air)
R = 8.3145; % [J/(mol*K)] (Universal Gas Constant)
m_m = 28.02; % [g/mol] (Molar Mass of Air)

rpm = 50000;
radius = 0.0375;
a_in = 35.63709471;
a_out = 23.12163596;
Beta = 29.76722393;
u = 196.3495408;

nozzle_data = readtable("ColdGas-nozzle_values.csv");
nozzle_dict = struct();
for i = 1:height(nozzle_data)
    nozzle_dict.(nozzle_data.Variable{i}) = nozzle_data.Value(i);
end

A_throat = nozzle_dict.A_throat;
A_e = nozzle_dict.A_e;
area_ratio = A_e / A_throat;
R_air = nozzle_dict.R_S;
T_e = nozzle_dict.T_e;
rho_e = nozzle_dict.rho_e;

M_e = isentropicMachFinder(area_ratio, gamma)
v = mps_per_mach(M_e, gamma, R_air, T_e)

horse_power_test = 100;
torque_test = horse_power_test * 5252 / rpm * 1.355817;
m_dot = solveMassFlow(rpm, torque_test, radius, deg2rad(Beta), u, v, deg2rad(a_in), deg2rad(a_out), rho_e, false);

hp_values = 50:200;
m_dot_values = zeros(size(hp_values));
m_dot_values_imp = zeros(size(hp_values));
v_dot_values = zeros(size(hp_values));

for i = 1:length(hp_values)
    horse_power = hp_values(i);
    torque = horse_power * 5252 / rpm * 1.355817;
    [m_dot, m_dot_imperial, v_dot] = solveMassFlow(rpm, torque, radius, deg2rad(Beta), u, v, deg2rad(a_in), deg2rad(a_out), rho_e, true);
    m_dot_values(i) = m_dot;
    m_dot_values_imp(i) = m_dot_imperial;
    v_dot_values(i) = v_dot;
end

figure;
plot(hp_values, m_dot_values, 'LineWidth', 1.5)
xlabel('Horsepower')
ylabel('Mass Flow Rate (kg/s)')
title('Mass Flow Rate vs. Horsepower')
grid on


figure;
plot(hp_values, m_dot_values_imp, 'LineWidth', 1.5)
xlabel('Horsepower')
ylabel('Mass Flow Rate (lb/s)')
title('Mass Flow Rate vs. Horsepower')
grid on

figure;
plot(hp_values, v_dot_values, 'LineWidth', 1.5)
xlabel('Horsepower')
ylabel('Volumetric Flow Rate (m^3/min)')
title('Volumetric Flow Rate vs. Horsepower')
grid on

function [m_dot, m_dot_imperial, v_dot] = solveMassFlow(RPM, torque, radius, beta, U, v, a_in, a_out_guess, rho, hide)
    kgps_to_lbps = 2.2046244201838;
    
    a1 = a_in;

    % syms v2 a2 w m_dot

    % eq6 = [v * cos(a1); v * sin(a1)] == [w * cos(beta); w * sin(beta) + U];
    % eq7 = [v2 * cos(a2); v2 * sin(a2)] == [w * cos(-beta); w * sin(-beta) + U];
    % eq8 = abs(v * sin(a1) - v2 * sin(a2)) == torque / radius / m_dot;

    % equations = subs([eq6(1), eq7(1), eq7(2), eq8]);
    % solution = vpasolve(equations, [v2, a2, w, m_dot], [v; a_out_guess; v; 4]);
    % m_dot = double(solution.m_dot)

    w = v * cos(a_in) / cos(beta);
    v2sina2 = w * sin(-beta) + U;
    m_dot = torque / radius / abs(v * sin(a1) - v2sina2); 
    m_dot_imperial = m_dot * kgps_to_lbps;
    v_dot = m_dot / rho * 60;

    if (~hide)
        fprintf("Mass Flow:\n");
        fprintf("\t%.5f kg/s\n", m_dot);
        fprintf("\t%.5f lb/s\n", m_dot_imperial);
        fprintf("\t%.5f m^3/min\n", v_dot);
    end
end


function Mach = isentropicMachFinder(area_ratio, gamma)
    options = optimset('Display', 'off');
    Mach_guess = 1.2;    
    % ((gamma + 1) / 2) ^ (-(gamma + 1) / (2 * (gamma - 1))) * (1 + (gamma - 1) / 2 * M^2) ^ ((gamma + 1) / (2 * (gamma - 1))) / M
    % (1/M) * ((2/(gamma+1)) * (1 + ((gamma-1)/2) * M^2))^((gamma+1)/(2*(gamma-1)))
    func = @(M) area_ratio - ((gamma + 1) / 2) ^ (-(gamma + 1) / (2 * (gamma - 1))) * (1 + (gamma - 1) / 2 * M^2) ^ ((gamma + 1) / (2 * (gamma - 1))) / M;
    Mach = fsolve(func, Mach_guess, options);
end

function v = mps_per_mach(Mach, gamma, R, T)
    a = sqrt(gamma * R * T);
    v = Mach * a;
end