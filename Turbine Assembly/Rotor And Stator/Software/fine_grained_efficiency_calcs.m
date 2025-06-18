% Example Inputs
Re = 3e5;
MinH = 1.2;
M2in = 1.5;
M2out = 2.0;
rH_rT = 0.9;
Pin_Pout = 1.05;
L = 0.04;       % chord [m]
H = 0.06;       % blade height [m]
Lx = 0.02;      % axial chord [m]
alpha_in = 70;  % degrees
alpha_out = 30; % degrees
alpha_m = 50;   % degrees
Kp = 1.0;       % profile correction factor
Kp_xi = 1.0;    % assumed constant
YP_i0 = 0.02;   % base profile loss
Delta_Te = 0.05;
gamma = 1.33;

Y = kacker_okapuu_loss(Re, MinH, M2in, M2out, rH_rT, Pin_Pout, ...
                       L, H, Lx, alpha_in, alpha_out, alpha_m, ...
                       Kp, Kp_xi, YP_i0, Delta_Te, gamma);
fprintf('Total loss coefficient Y = %.5f\n', Y);


function eta_t = estimate_turbine_efficiency(thp, m_dot, Cp, T0, Rt, gamma)
    %TURBINE_EFFICIENCY Calculates the overall turbine efficiency
    %
    % Inputs:
    %   thp   - Turbine horsepower (hp)
    %   m_dot - Mass flow rate (lbm/s)
    %   Cp    - Specific heat at constant pressure (Btu/lbm-R)
    %   T0    - Inlet total temperature (R)
    %   Rt    - Total pressure ratio (P_inlet / P_exit)
    %   gamma - Specific heat ratio
    %
    % Output:
    %   eta_t - Turbine efficiency (unitless)

        % Convert horsepower to ft-lbf/s: 1 hp = 550 ft-lbf/s
        thp_ftlbfs = thp * 550;

        % Compute efficiency
        eta_t = thp_ftlbfs / (m_dot * Cp * T0 * (1 - (1/Rt)^((gamma - 1)/gamma)));
end


function M_rel = calc_relative_mach_number(V_abs, U, alpha, a)
    % V_abs: absolute velocity at inlet for M_rel_in or outlet for M_rel_out
    % U: blade speed at rotor inlet radius for M_rel_in or outlet radius
    % alpha: absolute flow angle for inlet 
    M_rel = sqrt(V_abs ^ 2 + U^2 - 2 * V_abs * U * cos(alpha)) / a
end

function Y_total = kacker_okapuu_loss(Re, MinH, Min, Mout, rH_rT, Pin_Pout, ...
                                       c, H, Lx, alpha_in, alpha_out, alpha_m, ...
                                       Kp, Kp_xi, Delta_Te, gamma)
    % Re - Reynolds number
    % MinH - relative Mach number at the hub at rotor inlet (can use M_in - ~5-10%)
    % Min - relative inlet Mach number squared
    % Mout - relative outlet Mach number squared 
    % rH_rT - hub to tip radius ratio
    % Pin_Pout - inlet to exit pressure ratio
    % c - chord
    % H - blade height
    % Lx - axial chord
    % alpha_in - inlet 
    
    

    %% Reynolds number correction factor
    if Re <= 2e5
        xRe = (Re / (2e5))^(-0.4);
    elseif Re <= 1e6
        xRe = 1.0;
    else
        xRe = (Re / (1e6))^(-0.2);
    end

    %% Baseline loss for blade with zero incidence (momentum boundary layer estimate)


    %% Shock loss
    term1 = (1 + (gamma - 1)/2 * Min^2)^(gamma / (gamma - 1));
    term2 = (1 + (gamma - 1)/2 * Mout^2)^(gamma / (gamma - 1));
    Y_shock = 0.75 * (MinH - 0.4)^1.75 * (rH_rT) * (Pin_Pout) * (term1 - term2);

    %% Profile loss
    YP = 0.914 * (2/3 * Kp * Kp_xi * YP_i0 + Y_shock);

    %% Aspect ratio correction factor (assume χ_AR = 1 unless given)
    chi_AR = 1;  

    %% Secondary loss
    alpha_in_rad = deg2rad(alpha_in);
    alpha_out_rad = deg2rad(alpha_out);
    alpha_m_rad = deg2rad(alpha_m);

    term_angle = tan(alpha_in_rad) - tan(alpha_out_rad);
    YS = 0.04 * (L / H) * chi_AR * (4 * term_angle^2) * ...
        (cos(alpha_out_rad)^2 / cos(alpha_m_rad)) * ...
        (cos(alpha_out_rad) / cos(alpha_in_rad)) * ...
        (1 - (Lx / H)^2 * (1 - Kp));

    %% Trailing edge loss
    t1 = 1 + (gamma - 1)/2 * M2out^2 * (1 - 1 / (1 - Delta_Te));
    t2 = 1 + (gamma - 1)/2 * M2out^2;
    YTe = (t1^(-gamma / (gamma - 1)) - 1) / (1 - t2^(-gamma / (gamma - 1)));

    %% Total loss coefficient
    Y_total = xRe * YP + YS + YTe;
end


function Y_s = shock_loss(M1, gamma)
    % Estimate total pressure loss coefficient from a normal shock
    % M1: Mach number before shock
    % gamma: specific heat ratio (e.g., 1.4 for air)

    % Total pressure ratio across normal shock
    P0_ratio = ((1 + ((gamma - 1)/2) * M1^2) / ...
            (gamma * M1^2 - (gamma - 1)/2))^(gamma / (gamma - 1)) * ...
            ((gamma + 1) * M1^2 / ((gamma - 1) * M1^2 + 2))^(1 / (gamma - 1));

    Y_s = 1 - P0_ratio;
end

function Y_v = viscous_loss(V_rel, rho, chord, mu)
    % Estimate viscous profile loss coefficient using flat plate analogy
    % V_rel: relative velocity [m/s]
    % rho: density [kg/m^3]
    % chord: blade chord length [m]
    % mu: dynamic viscosity [Pa·s]

    Re = (rho * V_rel * chord) / mu;
    f = 0.316 / Re^0.25; % Blasius correlation
    Y_v = f; % approximate as loss coefficient
end

function P_ideal = ideal_power(mdot, V_exit)
    % Ideal power from stator kinetic energy (impulse stage)
    % V_exit: absolute velocity from stator exit (m/s)

    P_ideal = mdot * (V_exit^2) / 2;
end

function eta_stage = stage_efficiency(P_actual, P_ideal)
    % Overall stage efficiency
    eta_stage = P_actual / P_ideal;
end

function eta_rotor = total_loss_efficiency(Y_s, Y_v)
    % Efficiency considering shock and viscous losses
    % Assume multiplicative loss instead of additive
    eta_rotor = (1 - Y_s) * (1 - Y_v);
end

function Y_v = estimate_profile_loss(delta_theta_deg, t_over_c, Re_rel, M_rel)
    % ESTIMATE_PROFILE_LOSS Estimate profile (viscous) loss coefficient for turbine rotor
    % 
    % Inputs:
    %   delta_theta_deg : Flow turning angle across blade (degrees)
    %   t_over_c        : Blade thickness-to-chord ratio (t/c)
    %   Re_rel          : Reynolds number based on chord and relative velocity (optional)
    %   M_rel           : Relative inlet Mach number to rotor (optional)
    %
    % Output:
    %   Y_v             : Estimated profile (viscous) loss coefficient

    % if nargin < 3 || isempty(Re_rel)
    %     Re_rel = 1e6; % assume default if unknown
    % end
    % if nargin < 4 || isempty(M_rel)
    %     M_rel = 1; % assume supersonic if unknown
    % end

    % Convert turning angle to radians
    delta_theta = deg2rad(delta_theta_deg);

    % Empirical constants from cascade data
    K = 0.014;  % Loss due to turning
    C = 0.02;   % Loss due to thickness

    % Basic profile loss model
    Y_v_basic = K * delta_theta^2 + C * t_over_c;

    % Reynolds number correction (optional)
    Re_correction = (1e6 / Re_rel)^0.2;  % Mild increase in loss with low Re
    Y_v = Y_v_basic * Re_correction;

    % Supersonic adjustment (shock-boundary layer interaction)
    if M_rel > 1.0
        supersonic_factor = 1 + 0.5 * (M_rel - 1);
        Y_v = Y_v * supersonic_factor;
        warning('Supersonic flow detected (M = %.2f). Shock-boundary interaction may increase profile losses.', M_rel);
    end
end

function zeta_sf = secondary_flow_loss_model(V1, V2, c, nu, h, s, alpha1, alpha2, beta)
    % Calculate Secondary Flow loss at tips and root
    % V1 and V2 are relative velocity for rotor stage
    % c is chord len
    % nu is dynamic viscosity
    % h is height
    % s is spacing
    % 
    K1 = 4.65;
    K2 = 0.675; 
    Re_c = V2 * c / nu;
    c_fp = 0.074 * Re_c ^ (-1/5);
    vel_ratio = V1 / V2;
    f = 0.2 * sum([1, vel_ratio, vel_ratio^2, vel_ratio^3, vel_ratio^4]);
    c_f = c_fp * (f ^ 0.8);
    zeta_sf = c / h * 2 * c_f * (K1 + K2 * (C_L * c / s)^2 * sin(alpha1) / (sin(alpha2)^2));
end

function n_sf = static_efficiency(tau, omega, W, C_p, T_0, P4, P0, gamma)
    % tau is torque
    % omega is angular velocity
    % W is work
    % p4 / p0 is static to total pressure ratio accross the turbine

    n_sf = (tau * omega / W) / (C_p * T_0 * (1 - (P4 / P0) ^ ((gamma - 1) / gamma)));
end


%%%%0.38 OF ratio
%  stiffness of 5 psi
% chatter 