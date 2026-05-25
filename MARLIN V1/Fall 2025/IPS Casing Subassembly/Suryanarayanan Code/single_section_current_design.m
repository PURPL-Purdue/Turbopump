%% input
format long;
seal_radius = .75; %inches
sealNum = 2; %number of seals in the seal stack up

clearance = .004; %inches (radial)
tooth_width = .005; %inches
cavity_width = .031; %inches

teeth = 6;
loxPin = 1000; %psi
loxRecycleP = 500; %psi
loxRho = 1000; % density kg/m^3
loxMu = 1.0016e-3; %Pa s



%% change units
loxPin = 6894.76*loxPin; %Pa
loxRecycleP = 6894.76*loxRecycleP; %Pa

seal_radius = seal_radius * .0245; %m
clearance = clearance * .0245; %m
tooth_width = tooth_width * .0245; %m
cavity_width = cavity_width * .0245; %m
tooth_pitch = tooth_width + cavity_width; %m

%% Calculation
[mdotloxRecyle,ReloxRecycle] = laby_leakage(clearance, tooth_pitch, tooth_width, seal_radius, teeth, loxPin, loxRecycleP, loxRho, loxMu); %kg/s
QloxRecycle = mdotloxRecyle / loxRho;
mdotloxRecyle_imperial = mdotloxRecyle*2.2046; %lbm/s
QloxRecycle_imperial = QloxRecycle*264.171*60; %gal/m;


