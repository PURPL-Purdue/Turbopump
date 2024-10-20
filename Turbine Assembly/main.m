%% THIS SCRIPT LOADS ALL DESIGN CONSTANTS INTO A FILE CALLED "params.mat"
params_file = 'params.mat';

%% Turbine Constants
% TODO: Add final turbopump constants

%% Cold Gas Testing Constants:
rotor_radius = 0.04; % [m]
hub_radius = 0.02; % [m]
mass_flow_n2 = 0.17; % [kg/s]
shaft_power = 10; % [kW]
gamma_n2 = 1.4; % specific heat ratio

turbine_rpm = 10000; % [rpm]
input_velocity = NaN; % TODO: From Andrew
nozzle_radius = NaN; % TODO: From Andrew
% TODO: Andrew's constants for nozzle:

%% Conversion Factors
kw_to_hp = 1.34102209; % [hp/kw]

%% Save Parameters
% save(params_file, ...
%     "rotor_radius", "hub_radius", "mass_flow_n2", "shaft_power", "gamma_n2", "nozzle_radius", "input_velocity", "turbine_rpm", "shaft_power", ...
%     "kw_to_hp");
save(params_file)