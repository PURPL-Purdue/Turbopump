%input
format long;
seal_radius = .3346; %inches


clearance = .002; %inches
tooth_width = .01; %inches
cavity_width = .031; %inches
teeth = 10;

atmP = 14.7; %psi

rp1Pin = 60; %psi
rp1Pdrain = 30; %psi
rp1Rho = 800; %kg/m^3
rp1Mu = 1.92e-3; %Pa s




%change units

rp1Pin = 6894.76*rp1Pin;
rp1Pdrain = 6894.76*rp1Pdrain;
atmP = atmP * 6894.76;

seal_radius = seal_radius*.0245;
clearance = clearance * .0245;
tooth_width = tooth_width * .0245;
cavity_width = cavity_width * .0245;
tooth_pitch = tooth_width + cavity_width;


mdotrp1 = laby_leakage(clearance, tooth_pitch, tooth_width, seal_radius, teeth, rp1Pin, rp1Pdrain, rp1Rho, rp1Mu); %kg/s
Qrp1 = mdotrp1 / rp1Rho;
mdotrp1 = mdotrp1*2.2046 %lbm/s
Qrp1 = Qrp1*264.171 %gal/s;

