function funs = analysis_functions
    % Function that returns a struct of function handles
    funs.calculateMachPressureDistribution = @calculateMachPressureDistribution;
    funs.isentropicMachFinder = @isentropicMachFinder;
    funs.plotMachPressureDistributions = @plotMachPressureDistributions;
    funs.kantrowitz_limit = @kantrowitz_limit;
    funs.calc_relative_mach = @calc_relative_mach;
end

% function definitions
function [Mach_vec, P_vec] = calculateMachPressureDistribution(area_vec, gamma, R, T_inlet, P_inlet, M_inlet, M_outlet)
    Mach_vec = zeros(size(area_vec));
    P_vec = zeros(size(area_vec));
    
    % total temperature and pressure at inlet
    Tt_inlet = T_inlet * (1 + ((gamma - 1) / 2) * M_inlet^2);
    Pt_inlet = P_inlet * (1 + ((gamma - 1) / 2) * M_inlet^2)^(gamma / (gamma - 1));

    % calc mach number along the blade 
    for i = 1:length(area_vec)
        area_ratio = area_vec(i) / area_vec(1);
        
        % assume isentropic and numerically solve for mach number from area - mach equation (Andrew uses for nozzle)
        Mach_vec(i) = isentropicMachFinder(area_ratio, gamma, M_inlet, M_outlet);
        
        % calc static temperature and pressure using isentropic relations
        T_static = Tt_inlet / (1 + ((gamma - 1) / 2) * Mach_vec(i)^2);
        P_vec(i) = Pt_inlet * (T_static / Tt_inlet)^(gamma / (gamma - 1));
    end
end

function Mach = isentropicMachFinder(area_ratio, gamma, M_inlet, M_outlet)
    options = optimset('Display', 'off');
    Mach_guess = (M_inlet + M_outlet) / 2; % guess as avg mach
    
    % function for the isentropic area-Mach relation
    func = @(M) area_ratio - (1/M) * ((2/(gamma+1)) * (1 + ((gamma-1)/2) * M^2))^((gamma+1)/(2*(gamma-1)));
    
    % sovle for mach
    Mach = fsolve(func, Mach_guess, options);
end

function plotMachPressureDistributions(Mach_vec, P_vec)
    x = [1:1:length(Mach_vec)]; % lower coordinate  

    figure;
    sgtitle("Analysis along Lower Surface")
    hold on
    subplot(2, 1, 1);
    plot(x, Mach_vec, 'b-')
    grid on
    xlabel("Lower Surface Coordinate")
    ylabel("Mach Number")
    title("Mach Number vs Lower (Pressure) Surface Coordinate")

    subplot(2, 1, 2)
    plot(x, P_vec, 'g-')
    grid on
    xlabel("Lower Surface Coordinate")
    ylabel("Static Pressure [Pa]")
    title("Pressure Distribution along Pressure Surface ")
end

function A_ratio = kantrowitz_limit(M0, gamma)
    term1 = M0 * ( (gamma + 1) / (2 + (gamma - 1) * M0^2) )^((gamma + 1) / (2 * (gamma - 1)));
    term2 = ( ( (gamma + 1) * M0^2 ) / ( (gamma - 1) * M0^2 + 2 ) )^(-gamma / (gamma - 1));
    term3 = ( (gamma + 1) / (2 * gamma * M0^2 - (gamma - 1)) )^(-1 / (gamma - 1));
    
    A_ratio = term1 * term2 * term3;
end

function M_rel = calc_relative_mach(U, V, theta, gamma, R, T)
    a = sqrt(gamma * R * T);
    M_rel = sqrt((V * cos(theta) - U) ^ 2 + (V * sin(theta)) ^ 2) / a;
end

function A_inlet = calculate_inlet_area(U, V, theta, gamma, R, T, A_exit, M_exit)
    M_rel = calc_relative_mach(U, V, theta, gamma, R, T);

    A_star_inlet_ratio = area_mach_relation(M_rel, gamma);
    A_star_exit_ratio = area_mach_relation(M_exit, gamma);

    A_inlet = A_exit * (A_star_inlet_ratio / A_star_exit_ratio);
end

function A_ratio = area_mach_relation(M, gamma)
    A_ratio = (1 / M) * ((2 / (gamma + 1)) * (1 + (gamma - 1) / 2 * M^2))^((gamma + 1) / (2 * (gamma - 1)));
end
