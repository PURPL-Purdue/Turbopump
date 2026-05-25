%% input
format long;
seal_radii_min = .315; %inches
step = .011; %inches
sealNum = 5; %number of seals in the seal stack up

clearance = .002; %inches
tooth_width = .01; %inches
cavity_width = .031; %inches

teeth = 5;
loxPin = 1000; %psi
loxRecycleP = 50; %psi
loxRho = 1140; % density kg/m^3
loxMu = 1.95e-4; %Pa s

atmP = 14.7; %psi
N2inPressure = 1000; %psi

fuelPin = 1000; %psi
fuelRecycleP = 60; %psi
fuelRho = 800; %kg/m^3
fuelMu = 1.92e-3; %Pa s

%% change units
loxPin = 6894.76*loxPin; %Pa
loxRecycleP = 6894.76*loxRecycleP; %Pa
N2inPressure = 6894.76*N2inPressure; %Pa
fuelPin = 6894.76*fuelPin; %Pa
fuelRecycleP = 6894.76*fuelRecycleP; %Pa
atmP = atmP * 6894.76; %Pa

step = step * .0245; %m
seal_radii_min = seal_radii_min * .0245; %m
clearance = clearance * .0245; %m
tooth_width = tooth_width * .0245; %m
cavity_width = cavity_width * .0245; %m
tooth_pitch = tooth_width + cavity_width; %m

seal_radii = [];
for i=1:sealNum
    seal_radii(i) = seal_radii_min + (i-1)*step;
end

[mdotloxRecyle,ReloxRecycle] = laby_leakage(clearance, tooth_pitch, tooth_width, seal_radii(1), teeth, loxPin, loxRecycleP, loxRho, loxMu); %kg/s
QloxRecycle = mdotloxRecyle / loxRho;
mdotloxRecyle = mdotloxRecyle*2.2046; %lbm/s
QloxRecycle = QloxRecycle*264.171; %gal/s;

[mdotLoxDrain, ReLoxDrain] = laby_leakage(clearance, tooth_pitch, tooth_width, seal_radii(2), teeth, loxRecycleP, atmP, loxRho, loxMu);
QloxDrain = mdotLoxDrain / loxRho;
mdotLoxDrain = mdotLoxDrain*2.2046; %lbm/s
QloxDrain = QloxDrain*264.1721; %gal/s

[mdotfuelRecycle,ReFuelRecycle] = laby_leakage(clearance, tooth_pitch, tooth_width, seal_radii(6), teeth, fuelPin, fuelRecycleP, fuelRho, fuelMu); %kg/s
QfuelRecycle = mdotfuelRecycle / fuelRho *264.172; %gal/s
mdotfuelRecycle = mdotfuelRecycle * 2.2046; %lbm/s

[mdotfuelDrain, ReFuelDrain] = laby_leakage(clearance, tooth_pitch, tooth_width, seal_radii(5), teeth, fuelRecycleP, atmP, fuelRho, fuelMu); %kg/s
QfuelDrain = mdotfuelDrain /fuelRho * 264.172;   %gal/s
mdotfuelDrain = mdotfuelDrain * 2.2046;   %lbm/s

fprintf("\n\nReynolds Numbers:\nLOx Recycle = %f\nLox Drain = %f\nFuel Recycle = %f\nFuel Drain = %f\n",ReloxRecycle, ReLoxDrain, ReFuelRecycle, ReFuelDrain)
fprintf("\nLeakage Rates (lbm/s):\n LOx Recycle = %f\nLox Drain = %f\nFuel Recycle = %f\nFuel Drain = %f\n",mdotloxRecyle, mdotLoxDrain, mdotfuelRecycle, mdotfuelDrain)
fprintf("\nLeakage Rates (gal/s):\n LOx Recycle = %f\nLox Drain = %f\nFuel Recycle = %f\nFuel Drain = %f\n",QloxRecycle, QloxDrain, QfuelRecycle, QfuelDrain)

