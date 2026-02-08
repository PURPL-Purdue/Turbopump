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
runtime = input('Runtime? \n');
MaterialChange = input('Is there any change in material?(0=no / 1=yes) \n');
if MaterialChange == 0
    L = input('Thickness? (mm)\n')/1000;
    mtype1 = input('Material? (1=copper/2=steelAISI8630/3=graphite) \n');
    mtype2 = mtype1; % Set material 2 same as 1
    L1 = L; L2 = 0;  % Entire length is material 1
else
    mtype1 = input('Material 1 (Inner)? (1=copper/2=steelAISI8630/3=graphite) \n');
    L1 = input('Thickness 1? (mm)\n')/1000;
    mtype2 = input('Material 2 (Outer)? (1=copper/2=steelAISI8630/3=graphite) \n');
    L2 = input('Thickness 2? (mm)\n')/1000;
    L = L1 + L2;
end
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
% [ kw (W/m·K), rho (kg/m^3), cp (J/kg·K), Polynomial Max Temp ]
materials = cell(3,1);
materials{1} = [507.735, -0.540163, 0.000706823, -3.13884*10^-7; 9075, -0.512415, 0 , 0; 188.548, 1.00219, -0.00132148, 5.75175*10^-7;1200,0 ,0, 0];   % Copper
materials{2} = [46.6,0 , 0, 0; 7850, 0, 0, 0; 475, 0, 0, 0;5000, 0, 0, 0]; % AISI 8630
materials{3} = [168, 0, 0 ,0; 2500, 0 ,0, 0; 717, 0, 0, 0;5000 0 0 0];  % Graphite
Mat1 = materials{mtype1};
Mat2 = materials{mtype2};

% Wall thickness per station [Injector, Combustor, Throat, Exit]
%param.Thickness = [0.0228854; 0.0228854; 0.03; 0.021856];
%L = param.Thickness(param.Station == section);               % [m]

%% -------------------- RADIAL DISCRETIZATION ------------------------------
% 1-D finite-volume discretization in cylindrical coordinates

n  = 1000;            % Number of radial control volumes [-]
dr = L / n;           % Radial step size [m]

r_inner = D_in / 2;   % Inner wall radius [m]
r_outer = r_inner + L;% Outer wall radius [m]

% Cell faces and centers
r_f     = linspace(r_inner, r_outer, n+1)';                   % [m]
r_nodes = 0.5*(r_f(1:end-1) + r_f(2:end));                    % [m]
n_match = round(L1 / dr);

% Geometric quantities (per unit axial length)
Af  = 2*pi*r_f;                                               % Face area [m]
Vol = pi*(r_f(2:end).^2 - r_f(1:end-1).^2);                   % Volume [m^2]

%% -------------------- BOUNDARY CONDITIONS --------------------------------
h_nat = 5;              % Natural convection coefficient [W/m^2·K]
T_amb = 300;            % Ambient temperature [K]

tf = 60;               % Final time [s]
dt = 1e-3;              % Time step [s]

phi = 300 * ones(n,1);  % Initial wall temperature [K]

% Material-dependent melting temperature
param.melt = [1733; 1733; 3500; 1733];   % [K]
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
aw_pre = Af(1:end-1) / dr;     % West faces
ae_pre = Af(2:end)   / dr;     % East faces

%% -------------------- VISUALIZATION SETUP --------------------------------
writerObj = VideoWriter('Cylindrical_Heat_Transfer.mp4','MPEG-4');
writerObj.FrameRate = 30;
open(writerObj);

fig = figure('Color','w');
hP = plot(r_nodes - r_inner, phi, 'r','LineWidth',2);
grid on;
axis([0 L 280 4000]);
yline(T_melt,'--k','Melting Point','LabelVerticalAlignment','bottom');

xlabel('Wall depth (m)');
ylabel('Temperature (K)');

Melt = false;

for t = 0:dt:tf

    % 1. Gas-side Heat Transfer Logic (Boundary Condition Switch)
    Tw = phi(1);
    if t <= runtime
        % FORCED CONVECTION (Bartz Correlation)
        term_Ma = 1 + (gam-1)/2 * Ma^2;
        sigma = 1 / ((0.5*(Tw/T_g)*term_Ma + 0.5)^0.68 * term_Ma^0.12);
        h_inner = hG_base * sigma; 
        T_ref_inner = Taw; % Gas drives toward recovery temperature
    else
        % NATURAL CONVECTION (After 2 seconds)
        h_inner = h_nat; % Example value for stagnant hot gas [W/m²K]
        T_ref_inner = T_amb; % Assuming gas cools toward ambient or a set T
    end

    %% 2. Evaluate Properties with Material-Specific Clamping
    % Initialize vectors for the whole domain
    k_nodes = zeros(n,1); 
    rho_v   = zeros(n,1); 
    cp_v    = zeros(n,1);

    % --- SEGMENT 1 (Inner Material) ---
    idx1 = 1:n_match;
    T_limit1 = Mat1(4,1);
    T_prop1 = min(phi(idx1), T_limit1); % Apply material 1 limit

    k_nodes(idx1) = Mat1(1,1) + Mat1(1,2)*T_prop1 + Mat1(1,3)*T_prop1.^2 + Mat1(1,4)*T_prop1.^3;
    rho_v(idx1)   = Mat1(2,1) + Mat1(2,2)*T_prop1 + Mat1(2,3)*T_prop1.^2 + Mat1(2,4)*T_prop1.^3;
    cp_v(idx1)    = Mat1(3,1) + Mat1(3,2)*T_prop1 + Mat1(3,3)*T_prop1.^2 + Mat1(3,4)*T_prop1.^3;

    % --- SEGMENT 2 (Outer Material) ---
    if n_match < n
        idx2 = n_match+1:n;
        T_limit2 = Mat2(4,1);
        T_prop2 = min(phi(idx2), T_limit2); % Apply material 2 limit

        k_nodes(idx2) = Mat2(1,1) + Mat2(1,2)*T_prop2 + Mat2(1,3)*T_prop2.^2 + Mat2(1,4)*T_prop2.^3;
        rho_v(idx2)   = Mat2(2,1) + Mat2(2,2)*T_prop2 + Mat2(2,3)*T_prop2.^2 + Mat2(2,4)*T_prop2.^3;
        cp_v(idx2)    = Mat2(3,1) + Mat2(3,2)*T_prop2 + Mat2(3,3)*T_prop2.^2 + Mat2(3,4)*T_prop2.^3;
    end

    % 3. Calculate Face Conductance (Harmonic Mean)
    k_face = (2 * k_nodes(1:end-1) .* k_nodes(2:end)) ./ (k_nodes(1:end-1) + k_nodes(2:end));
    cond_int = (Af(2:n) .* k_face) / dr;
    % 4. Assemble TDMA Coefficients

    % Update storage coefficient
    ap0 = (rho_v .* cp_v .* Vol) / dt;
    ae = zeros(n,1); aw = zeros(n,1);
    ae(1:n-1) = cond_int;
    aw(2:n)   = cond_int;
    d = ap0 + aw + ae;

    % --- INNER BOUNDARY (Switchable Convection) ---
    U_inner = 1 / (1/h_inner + (dr/2)/k_nodes(1));
    d(1) = ap0(1) + ae(1) + U_inner * Af(1);

    % --- OUTER BOUNDARY (Adiabatic) ---
    d(end) = ap0(end) + aw(end);
    %U_outer = 1 / (1/h_nat + (dr/2)/k_nodes(end));
    %d(end) = ap0(end) + aw(end) + U_outer * Af(end); 


    a_diag = -ae(1:end-1);
    b_diag = -aw(2:end);

    % 5. Right-Hand Side (RHS)
    rhs = ap0 .* phi;
    rhs(1) = rhs(1) + U_inner * Af(1) * T_ref_inner;
    % rhs(end) has no addition because it is adiabatic
    %rhs(end) = rhs(end) + U_outer * Af(end) * T_amb; %if Natural Convection

    % 6. Thomas Algorithm (TDMA) Solver
    dp = d; cp_vec = rhs;
    for i = 2:n
        m = b_diag(i-1) / dp(i-1);
        dp(i) = dp(i) - m * a_diag(i-1);
        cp_vec(i) = cp_vec(i) - m * cp_vec(i-1);
    end
    phi(n) = cp_vec(n) / dp(n);
    for i = n-1:-1:1
        phi(i) = (cp_vec(i) - a_diag(i)*phi(i+1)) / dp(i);
    end

     % 7. Visualization and Check
    if mod(round(t/dt),10) == 0
        set(hP,'YData',phi);
        title(sprintf('t = %.2f s | T_w = %.0f K | T_{Oring} %.0f K', t, phi(1), phi(end)));
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

%% ========================================================================
% Version History
%V1: Initial State
%V2: Added dynamic material properties (Temperature Dependant)[Only Copper]
%V3: Added Runtime Input to Check Thermal Soaking 
%V4: Multiple Material & Thickness and Material Selection [Steel and Graphite added (constant)]