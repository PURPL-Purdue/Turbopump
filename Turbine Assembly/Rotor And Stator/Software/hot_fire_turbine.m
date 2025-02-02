% TODO: add mass flow optimization
% TODO: prendtl-meyer based US  
% TODO: spans 
% TODO: calc area ratio relative to throat for relative mach -> inlet area then optimize to minimize change in inlet area

%% INCLUDE EXTERNAL FUNCTIONS

addpath(genpath('MATLAB Functions'))

%% CONVERSION FACTORS 
lb_to_kg = 0.45359237;
bar_to_pa = 100000;
psi_to_pa = 6894.75729;
in_to_cm = 2.54; %1 in = 2.54 cm

%% HOT GAS CONSTANTS

P_0 = 350 * psi_to_pa; % 24.132 * bar_to_pa; % [N/m^2, 500psi]
P_e = 30 * psi_to_pa; % [N/m^2, 14.7psi]
P_A = P_e;
m_dot_imperial = 0.81; % 1.087; % pounds / second 
m_dot = m_dot_imperial * lb_to_kg; % [kg/s]
T_0 = 876.05; % [K]
gamma = 1.1201; % (Specific heat ratio )
R = 8.3145; % [J/(mol*K)] (Universal Gas Constant)
m_m = 11.328; % [g/mol] (Molar Mass )
n = 8; % (number of nozzles)

%% TURBINE CONSTANTS

rotor_radius = 2 * in_to_cm / 100; % 0.04; % [m]
hub_radius = rotor_radius - 0.005; % 0.035; % [m] % 0.02
mass_flow = m_dot; % [kg/s]
shaft_power = 200; % 150; % [kW]
blade_solidity = 1 / 0.8; % c / s
turbine_rpm = 50000; % [rpm]

%% Conversion Factors
kw_to_hp = 1.34102209; % [hp/kw]

% Input parameters
radius = (rotor_radius + hub_radius) / 2; % [m]
horse_power = 200; % shaft_power * kw_to_hp; % [HP]
torque = horse_power * 5252 / turbine_rpm * 1.355817; % 746 * horse_power/(2 * pi * turbine_rpm/60); % N*m
degree_of_reaction = 0;
num_blades = 35; % 10
chord = 0.01; % [m]
blade_spacing = 2 * pi * hub_radius / num_blades; % [m]
beta = 60; % initial estimate [deg] - refined in vtriangle
% V_in = 564; % [m/2]
min_blade_thickness = 0.002; % [m] % 0.01
inlet_area = 0.0025; % [m]

%% GAS CONSTANTS TABLE

% gas parameters
Variable = {'P_0'; 'P_e'; 'P_A'; 'm_dot'; 'T_0'; 'gamma'; 'R'; 'm_m'};
Value = [P_0; P_e; P_A; m_dot; T_0; gamma; R; m_m];
Units = {'N/m^2'; 'N/m^2'; 'N/m^2'; 'kg/s'; 'K'; '-'; 'J/(mol*K)'; 'g/mol'};
Description = {
    'Total pressure (500 psi)';
    'Exit pressure (14.7 psi)';
    'Ambient pressure';
    'Mass flow rate';
    'Total temperature';
    'Specific heat ratio';
    'Universal gas constant';
    'Molar mass'
};

T = table(Variable, Value, Units, Description);
writetable(T, 'OutputTables/FinalHotGas-gas_values.csv');
disp(T);

%% HOT GAS NOZZLE CALCULATIONS
nozzle_funs = nozzle_functions();

R_S = nozzle_funs.calc_R_S(R, m_m); % The actual value is 296.1, idk why its giving a different number
rho_0 = nozzle_funs.calc_rho_0(P_0, R_S, T_0); % [kg/m^3]
v_e = nozzle_funs.calc_v_e(T_0, R_S, gamma, P_e, P_0); % [m/s]


T_throat = nozzle_funs.calc_T_throat(T_0, gamma); % [K]
P_throat = nozzle_funs.calc_P_throat(P_0, gamma); % [N/m^2]
rho_throat = nozzle_funs.calc_rho_throat(P_throat, R_S, T_throat); % [kg/m^3]
v_throat = nozzle_funs.calc_v_throat(gamma, R_S, T_throat); % [m/s]
A_throat = nozzle_funs.calc_A_throat(m_dot, rho_throat, v_throat); % [m^2]


M_e = nozzle_funs.calc_M_e(P_e, P_0, gamma);
rho_e = nozzle_funs.calc_rho_e(rho_0, M_e, gamma); % [kg/m^3]
A_e = nozzle_funs.calc_A_e(m_dot, rho_e, v_e); % [m^2]
T_e = nozzle_funs.calc_T_e(T_0, gamma, M_e); % [K]


r_throat = nozzle_funs.calc_r_throat(A_throat); % [m]
r_e = nozzle_funs.calc_r_e(A_e); % [m]


dist = nozzle_funs.calc_dist(r_throat, r_e); % [m]
F_thrust = nozzle_funs.calc_F_thrust(m_dot, v_e); % [N]


A_throat_n = nozzle_funs.calc_A_throat_n(A_throat, n); % [m^2]
A_e_n = nozzle_funs.calc_A_e_n(A_e, n); % [m^2]
r_throat_n = nozzle_funs.calc_r_throat_n(A_throat_n); % [m]
r_e_n = nozzle_funs.calc_r_e_n(A_e_n); % [m]
dist_n = nozzle_funs.calc_dist_n(r_throat_n, r_e_n); % [m]

%% PLOT NOZZLE GEOMETRY

nozzle_plot_funs = nozzle_plot_functions();

[X,Y,Z] = nozzle_plot_funs.plot_nozzle(A_e, A_throat, A_e, dist_n, dist_n);

%% NOZZLE TABLE


% Nozzle Parameters
Variable = {
    'R_S'; 'rho_0'; 'v_e'; 'T_throat'; 'P_throat'; 'rho_throat'; 'v_throat';
    'A_throat'; 'M_e'; 'rho_e'; 'A_e'; 'T_e'; 'r_throat'; 'r_e'; 'dist'; 
    'F_thrust'; 'A_throat_n'; 'A_e_n'; 'r_throat_n'; 'r_e_n'; 'dist_n'
};
Value = [
    R_S; rho_0; v_e; T_throat; P_throat; rho_throat; v_throat;
    A_throat; M_e; rho_e; A_e; T_e; r_throat; r_e; dist;
    F_thrust; A_throat_n; A_e_n; r_throat_n; r_e_n; dist_n
];
Units = {
    'J/(kg·K)'; 'kg/m^3'; 'm/s'; 'K'; 'N/m^2'; 'kg/m^3'; 'm/s';
    'm^2'; '-'; 'kg/m^3'; 'm^2'; 'K'; 'm'; 'm'; 'm'; 
    'N'; 'm^2'; 'm^2'; 'm'; 'm'; 'm'
};
Description = {
    'Specific gas constant'; 'Stagnation density'; 'Exit velocity'; 
    'Throat temperature'; 'Throat pressure'; 'Throat density'; 'Throat velocity';
    'Throat area'; 'Exit Mach number'; 'Exit density'; 'Exit area'; 
    'Exit temperature'; 'Throat radius'; 'Exit radius'; 'Nozzle length'; 
    'Thrust force'; 'Throat area per nozzle'; 'Exit area per nozzle'; 
    'Throat radius per nozzle'; 'Exit radius per nozzle'; 'Nozzle length per nozzle'
};

T = table(Variable, Value, Units, Description);
writetable(T, 'OutputTables/FinalHotGas-nozzle_values.csv');
disp(T);

%% TURBINE CALCULATONS

turbine_funs = turbine_functions();
analysis_funs = analysis_functions();


v_e = 1666.5847; % m/s
[v1, v2, w, u, a1, a2, b] = turbine_funs.rotor_back_calculate(turbine_rpm, torque / num_blades, mass_flow / num_blades, deg2rad(60), radius, v_e);
turbine_funs.plot_velocity_triangles_angles(v1, v2, u, w, w, chord, 0, b, -b, a1, -a2);

% a2 = mod(a2, 360);
v2_comp =  [v2 * cos(a2), v2 * sin(a2)];
v2 = sqrt(v2_comp(1) ^ 2 + v2_comp(2) ^ 2);
a2 = atan2d(v2_comp(1), v2_comp(2));

c_blade_sp = 0.8 * chord % 0.6 - 1.2
c_num_blades = 2 * pi * radius / c_blade_sp
% c_area = mass_flow / c_num_blades / r_e / (v1 * sin(a1))
c_area = A_e / c_num_blades
c_arc_theta = c_area / (0.5 * (rotor_radius^2 - hub_radius^2))
c_arc_len = radius * c_arc_theta


efficiency = turbine_funs.calculate_blade_efficiency(mass_flow, v1, v2, w, w, b, b, a1);
fprintf("isentropic efficiency: %.4f\n", efficiency)

[max_blade_thickness, cross_sectional_area, x_lower, y_lower, x_upper, y_upper, areas] = turbine_funs.generate_blade_geom_constant_area(chord, rad2deg(b), inlet_area, blade_spacing, min_blade_thickness);
turbine_funs.plot_turbine(x_lower, y_lower, x_upper, y_upper, num_blades, hub_radius, chord, rotor_radius - hub_radius, max_blade_thickness / 2)
areas = areas * (rotor_radius - hub_radius);

[Mach_vec, P_vec] = analysis_funs.calculateMachPressureDistribution(areas, gamma, R, T_e, P_e, M_e, M_e);
analysis_funs.plotMachPressureDistributions(Mach_vec, P_vec); 

%% TURBINE SPAN CALCS

% num_spans = 10;
% span_positions = linspace(hub_radius, rotor_radius, num_spans);
% cross_sectional_areas = zeros(1, num_spans);
% efficiencies = zeros(1, num_spans);

% for i = 1:num_spans
%     radius = span_positions(i);

%     [v1, v2, w, u, a1, a2, b] = turbine_funs.rotor_back_calculate(turbine_rpm, torque / num_blades / num_spans, mass_flow / num_blades / num_spans, deg2rad(60), radius, v_e);
    
%     turbine_funs.plot_velocity_triangles_angles(v1, v2, u, w, w, chord, 0, b, -b, a1, -a2);

%     efficiencies(i) = turbine_funs.calculate_blade_efficiency(mass_flow, v1, v2, w, w, b, b, a1);

%     [max_blade_thickness, cross_sectional_area, x_lower, y_lower, x_upper, y_upper, area] = ...
%         turbine_funs.generate_blade_geom_constant_area(chord, rad2deg(b), inlet_area, blade_spacing, min_blade_thickness);


%     cross_sectional_areas(i) = cross_sectional_area;
%     areas = area * (radius - hub_radius);
    
%     % turbine_funs.plot_turbine(x_lower, y_lower, x_upper, y_upper, num_blades, radius, chord, rotor_radius - hub_radius, max_blade_thickness / 2);

%     [Mach_vec, P_vec] = analysis_funs.calculateMachPressureDistribution(areas, gamma, R, T_e, P_e, M_e, M_e);
%     analysis_funs.plotMachPressureDistributions(Mach_vec, P_vec);
% end

% % Print overall efficiency (average over all spans)
% average_efficiency = mean(efficiencies);
% fprintf("Average isentropic efficiency: %.4f\n", average_efficiency)


%% TURBINE TABLE

% Turbine values
Variable = {
    'rotor_radius'; 'hub_radius'; 'mass_flow'; 'shaft_power'; 'turbine_rpm';
    'radius'; 'horse_power'; 'torque'; 'degree_of_reaction'; 'num_blades';
    'chord'; 'blade_spacing'; 'min_blade_thickness'; 'inlet_area'; 
    'V_in'; 'V_out'; 'W_in'; 'W_out'; 'a_in'; 'a_out'; 'Beta'; 'u'
};
Value = [
    rotor_radius; hub_radius; mass_flow; shaft_power; turbine_rpm;
    radius; horse_power; torque; degree_of_reaction; num_blades;
    chord; blade_spacing; min_blade_thickness; inlet_area;
    v1; v2; w; w; rad2deg(a1); rad2deg(a2); rad2deg(b); u
];
Units = {
    'm'; 'm'; 'kg/s'; 'kW'; 'rpm';
    'm'; 'HP'; 'N*m'; '-'; '-';
    'm'; 'm'; 'm'; 'm²';
    'm/s'; 'm/s'; 'm/s'; 'm/s'; '°'; '°'; '°'; 'm/s'
};
Description = {
    'Rotor radius'; 'Hub radius'; 'Mass flow rate'; 'Shaft power'; 'Turbine RPM';
    'Average radius'; 'Horsepower'; 'Torque'; 'Degree of reaction'; 'Number of blades';
    'Blade chord length'; 'Blade spacing'; 'Minimum blade thickness'; 'Inlet area';
    'Inlet absolute velocity'; 'Outlet absolute velocity'; 'Inlet relative velocity'; 'Outlet relative velocity';
    'Inlet absolute angle'; 'Outlet absolute angle'; 'Blade angle (Beta)'; 'Turbine velocity'
};

T = table(Variable, Value, Units, Description);
writetable(T, 'OutputTables/FinalHotGas-turbine_values.csv');
disp(T);

%% Main Input Output Table (Turbine Geometry Specific Values)
Variable = {
    'Inputs';
    'mass flow';
    'turbine diameter';
    'turbine power';
    'annulus width';
    'solidity';
    '';
    'Outputs';
    'blade spacing';
    'num blades';
    'chord';
    'inlet area';
    'absolute inlet velocity';
    'relative inlet angle';
    'absolute inlet angle';
};

% TODO FINISH THIS AND MAKE SURE ANGLES ARE POS
Value = {
    '';
    mass_flow;
    2 * rotor_radius;
    horse_power;
    rotor_radius - hub_radius;
    blade_solidity;
    '';
    '';
    c_blade_sp;
    ceil(c_num_blades);
    chord;
    c_arc_len; % area at mean rad
    v1;
    rad2deg(b);
    rad2deg(a1)
};
Units = {
    ''; 
    'kg/s'; 'm'; 'HP'; 'm'; '-';
    ''; '';
    'm'; '-'; 'm'; 'm'; 'm/s'; 'deg (°)'; 'deg (°)' 
};
Description = {
    '';
    'mass flow'; 'outer diameter'; 'shaft power'; 'annulus width'; 'blade solidity';
    '';'';
    'blade spacing'; 'num bldades'; 'chord'; 'mean rad arc length'; 'inlet absolute velocity'; 'inlet relative angle'; 'inlet absolute angle'
};

T = table(Variable, Value, Units, Description);
writetable(T, 'OutputTables/FinalHotGas-IO.csv');
disp(T);

%% STRUCTURE CALCULATIONS

structure_funs = structure_functions();

% These values are random rough values %
omega = turbine_rpm / 60 * 2 * pi; % [rad/s] (Conversion of 10,000 rpm to rads)
height_blade = max_blade_thickness; % [m]
rho_blade = 0.845; % [kg/m^3]
Length_blade = chord; % [m]
height_bmin = min_blade_thickness; % [m]
width_blade = rotor_radius - hub_radius; % [m]
Z_blade = num_blades; % [Number of blades]
mean_diameter = radius * 2; % [m] mean blade diam

radius_turbine = structure_funs.calc_radius_turbine(height_blade, hub_radius); % [m]
mass_blade = structure_funs.calc_mass_blade(Length_blade, width_blade, height_bmin, rho_blade); % [kg]
Force_centrifugal = structure_funs.calc_Force_centrifugal(mass_blade, omega, radius_turbine); % [N]
stress_centrifugal = structure_funs.calc_stress_centrifugal(rho_blade, height_blade, mean_diameter, omega); % [N/m^2]
Force_tangential = structure_funs.calc_Force_tangential(mass_flow, w, b, w, -b); % [N]
Force_axial = structure_funs.calc_Force_axial(mass_flow, v1, a1, w, -b); % [N]
torque_blade = structure_funs.calc_torque_blade(Force_tangential, height_blade); % [Nm]
torque_turbine = structure_funs.calc_torque_turbine(Force_tangential, radius_turbine, Z_blade); % [Nm]
P = structure_funs.calc_P(torque_turbine, omega); % [Nm/s]
Force_gas = structure_funs.calc_Force_gas(Force_tangential, Force_axial); % [N]
Moment_Bending = structure_funs.calc_Moment_Bending(height_blade, Z_blade, Force_gas); % [Nm]
I = structure_funs.calc_I(Length_blade, height_bmin); % [m^4]
stress_gas = structure_funs.calc_stress_gas(height_bmin, Force_gas, width_blade, I); % [N/m^2]

%% STRUCTURE TABLES

% structure table
Variable = {
    'radius_turbine'; 'mass_blade'; 'Force_centrifugal'; 'stress_centrifugal';
    'Force_tangential'; 'Force_axial'; 'torque_blade'; 'torque_turbine';
    'P'; 'Force_gas'; 'Moment_Bending'; 'I'; 'stress_gas'
};
Value = [
    radius_turbine; mass_blade; Force_centrifugal; stress_centrifugal;
    Force_tangential; Force_axial; torque_blade; torque_turbine;
    P; Force_gas; Moment_Bending; I; stress_gas
];
Units = {
    'm'; 'kg'; 'N'; 'N/m^2';
    'N'; 'N'; 'N*m'; 'N*m';
    'W'; 'N'; 'N*m'; 'm^4'; 'N/m^2'
};
Description = {
    'Turbine radius'; 'Blade mass'; 'Centrifugal force on blade'; 'Centrifugal stress';
    'Tangential force'; 'Axial force'; 'Blade torque'; 'Turbine torque';
    'Power'; 'Gas force'; 'Bending moment'; 'Second moment of area (I)'; 'Gas stress'
};

T = table(Variable, Value, Units, Description);
writetable(T, 'OutputTables/FinalHotGas-structure_values.csv');
disp(T);