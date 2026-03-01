%% ========================================================================
%  1-D TRANSIENT CYLINDRICAL HEAT TRANSFER SOLVER (ROCKET WALL)
%
%  Author : Rafael Macia Titos
%
%  Physics:
%  - 1-D radial transient conduction in cylindrical coordinates
%  - Gas-side convection via Bartz correlation
%  - Outer-wall natural convection to ambient
%  - Fully implicit finite-volume formulation (TDMA)
%
%  Governing equation:
%    rho*cp*dT/dt = (1/r)*d/dr ( k*r*dT/dr )
%% ========================================================================

clear; clc;

%% -------------------- INPUT DATA ----------------------------------------
param   = readtable("turbopump_properties.csv");
contour = importdata("contour.csv");

param.Station = categorical(param.Station);

section = input( ...
    'Section to Analyze: Injector, Combustor, Throat, Exit \n','s');

%% -------------------- GAS PROPERTIES -------------------------------------
% T_g   : Gas static temperature [K]
% Pr    : Prandtl number [-]
% mu    : Dynamic viscosity [Pa·s]
% cp_g  : Gas specific heat at constant pressure [J/(kg·K)]
% rhog  : Gas density [kg/m^3]
% gam   : Ratio of specific heats [-]
% Ma    : Local Mach number [-]

T_g  = param.Temperature_K(param.Station == section);       % [K]
Pr   = param.Prandtl(param.Station == section);             % [-]
mu   = param.Viscosity_Pa_s(param.Station == section);      % [Pa·s]
cp_g = param.Cp_J_kgK(param.Station == section);            % [J/(kg·K)]
rhog = param.Density_kg_m3(param.Station == section);       % [kg/m^3]
gam  = param.Gamma(param.Station == section);               % [-]
Ma   = param.Mach(param.Station == section);                % [-]

% Chamber pressure (taken at injector)
Pc = param.Pressure_Pa(param.Station == "Injector");        % [Pa]

%% -------------------- GEOMETRY -------------------------------------------
% Nozzle contour radius converted from mm to m
r_contour = contour(:,2) / 1000;                             % [m]

% Internal diameter per station
param.Diameter = 2 * [ ...
    r_contour(1);        % Injector
    r_contour(1);        % Combustor
    min(r_contour);      % Throat
    r_contour(end) ];    % Exit

D_in   = param.Diameter(param.Station == section);           % [m]
A_flow = pi * (D_in/2)^2;                                    % [m^2]

% Throat reference quantities
Throat_D = 2 * min(r_contour);                               % [m]
Throat_A = pi * (Throat_D/2)^2;                              % [m^2]
Cstar    = 1796.7;                                           % [m/s]

%% -------------------- WALL MATERIAL --------------------------------------
% Material properties:
% [ k (W/m·K), rho (kg/m^3), cp (J/kg·K) ]

if section ~= "Throat"
    MaterialProp = [46.6, 7850, 477];     % AISI 8630 steel
else
    MaterialProp = [380, 7900, 385];      % Copper
end

k_w = MaterialProp(1);                    % Thermal conductivity [W/m·K]
rho = MaterialProp(2);                    % Density [kg/m^3]
cp  = MaterialProp(3);                    % Specific heat [J/kg·K]

% Wall thickness per station
param.Thickness = [0.023012; 0.023012; 0.038265; 0.021856];
L = param.Thickness(param.Station == section);               % [m]

%% -------------------- RADIAL DISCRETIZATION ------------------------------
% 1-D finite-volume discretization in cylindrical coordinates

n  = 1000;            % Number of radial control volumes [-]
dr = L / n;           % Radial step size [m]

r_inner = D_in / 2;   % Inner wall radius [m]
r_outer = r_inner + L;% Outer wall radius [m]

% Cell faces and centers
r_f     = linspace(r_inner, r_outer, n+1)';                   % [m]
r_nodes = 0.5*(r_f(1:end-1) + r_f(2:end));                    % [m]

% Geometric quantities (per unit axial length)
Af  = 2*pi*r_f;                                               % Face area [m]
Vol = pi*(r_f(2:end).^2 - r_f(1:end-1).^2);                   % Volume [m^2]

%% -------------------- BOUNDARY CONDITIONS --------------------------------
h_nat = 5;              % Natural convection coefficient [W/m^2·K]
T_amb = 300;            % Ambient temperature [K]

tf = 2.0;               % Final time [s]
dt = 5e-4;              % Time step [s]

phi = 300 * ones(n,1);  % Initial wall temperature [K]

% Material-dependent melting temperature
param.melt = [1733; 1733; 1358; 1733];   % [K]
T_melt = param.melt(param.Station == section);

%% -------------------- GAS-SIDE HEAT TRANSFER ------------------------------
% Bartz correlation for turbulent rocket nozzle flow
% h_g in [W/m^2·K]

hG_base = (0.026 / Throat_D^0.2) * ...
          (mu^0.2 * cp_g / Pr^0.6) * ...
          (Pc / Cstar)^0.8 * ...
          (Throat_A / A_flow)^0.9;

% Adiabatic wall (recovery) temperature [K]
Taw = T_g * (1 + (gam-1)/2 * Ma^2 * Pr^(1/3));

%% -------------------- NUMERICAL PRECOMPUTATION ----------------------------
% Conductive conductances at faces [W/K]
aw = k_w * Af(1:end-1) / dr;     % West faces
ae = k_w * Af(2:end)   / dr;     % East faces

% Transient storage coefficient [W/K]
ap0 = rho * cp * Vol / dt;

%% -------------------- VISUALIZATION SETUP --------------------------------
writerObj = VideoWriter('Cylindrical_Heat_Transfer.mp4','MPEG-4');
writerObj.FrameRate = 30;
open(writerObj);

fig = figure('Color','w');
hP = plot(r_nodes - r_inner, phi, 'r','LineWidth',2);
grid on;
axis([0 L 280 3200]);
yline(T_melt,'--k','Melting Point','LabelVerticalAlignment','bottom');

xlabel('Wall depth (m)');
ylabel('Temperature (K)');

%% -------------------- TRANSIENT SOLVER -----------------------------------
% Fully implicit finite-volume formulation:
%
%   ap0*T_i^(n+1) =
%     aw*T_{i-1}^(n+1) + ae*T_{i+1}^(n+1)
%     + ap0*T_i^n + source terms
%
% Boundary conditions are applied via equivalent source terms.

Melt = false;

for t = 0:dt:tf

    %% Gas-side heat transfer update (Bartz sigma correction)
    Tw = phi(1);                                   % Inner wall temp [K]
    term_Ma = 1 + (gam-1)/2 * Ma^2;

    sigma = 1 / ( ...
        (0.5*(Tw/T_g)*term_Ma + 0.5)^0.68 * ...
         term_Ma^0.12 );

    hG = hG_base * sigma;                          % [W/m^2·K]

    %% Assemble coefficient matrix
    d = ap0 + aw + ae;

    % Inner wall: gas convection
    d(1) = ap0(1) + ae(1) + hG*Af(1);

    % Outer wall: ambient convection
    d(end) = ap0(end) + aw(end) + h_nat*Af(end);

    % Off-diagonals
    a_diag = -ae(1:end-1);
    b_diag = -aw(2:end);

    %% Right-hand side
    rhs = ap0 .* phi;
    rhs(1)   = rhs(1)   + hG    * Af(1)   * Taw;
    rhs(end) = rhs(end) + h_nat * Af(end) * T_amb;

    %% Thomas Algorithm (TDMA)
    dp = d; cp_vec = rhs;
    for i = 2:n
        m = b_diag(i-1) / dp(i-1);
        dp(i)     = dp(i)     - m * a_diag(i-1);
        cp_vec(i) = cp_vec(i) - m * cp_vec(i-1);
    end

    phi(n) = cp_vec(n) / dp(n);
    for i = n-1:-1:1
        phi(i) = (cp_vec(i) - a_diag(i)*phi(i+1)) / dp(i);
    end

    %% Visualization and melting check
    if mod(round(t/dt),10) == 0
        set(hP,'YData',phi);
        title(sprintf('t = %.2f s | h_g = %.0f W/m²K | T_w = %.0f K', ...
              t, hG, phi(1)));
        drawnow;
        writeVideo(writerObj,getframe(fig));

        if phi(1) > T_melt && ~Melt
            fprintf('WARNING: Wall melting at t = %.3f s\n', t);
            Melt = true;
        end
    end
end

close(writerObj);
disp('Simulation complete.');