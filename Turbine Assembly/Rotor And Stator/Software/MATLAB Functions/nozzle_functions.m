function funs = nozzle_functions
    % function that returns a struct of function handles
    funs.calc_R_S = @calc_R_S;
    funs.calc_rho_0 = @calc_rho_0;
    funs.calc_v_e = @calc_v_e;
    funs.calc_T_throat = @calc_T_throat;
    funs.calc_P_throat = @calc_P_throat;
    funs.calc_rho_throat = @calc_rho_throat;
    funs.calc_v_throat = @calc_v_throat;
    funs.calc_A_throat = @calc_A_throat;
    funs.calc_M_e = @calc_M_e;
    funs.calc_rho_e = @calc_rho_e;
    funs.calc_A_e = @calc_A_e;
    funs.calc_T_e = @calc_T_e;
    funs.calc_r_throat = @calc_r_throat;
    funs.calc_r_e = @calc_r_e;
    funs.calc_dist = @calc_dist;
    funs.calc_A_throat_n = @calc_A_throat_n;
    funs.calc_A_e_n = @calc_A_e_n;
    funs.calc_r_throat_n = @calc_r_throat_n;
    funs.calc_r_e_n = @calc_r_e_n;
    funs.calc_dist_n = @calc_dist_n;
    funs.calc_F_thrust = @calc_F_thrust;
end

% function definitions:
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
