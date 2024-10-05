%% THIS SCRIPT LOADS ALL DESIGN CONSTANTS INTO A FILE CALLED "params.mat"
params_file = 'params.mat';

%% Turbine Constants
% TODO: Add final turbopump constants

%% Cold Gas Testing Constants:
rotor_radius = 0.04; % [m]
mass_flow_n2 = 0.17; % [kg/s]
shaft_power = 10; % [kW]
turbine_rpm = 10000; % [rpm]
input_velocity = NaN; % TODO: From Andrew
nozzle_radius = NaN; % TODO: From Andrew
% TODO: Andrew's constants for nozzle:


%% Save Parameters
save(params_file, ...
    "rotor_radius", "input_velocity", "turbine_rpm", "shaft_power", "mass_flow_n2");