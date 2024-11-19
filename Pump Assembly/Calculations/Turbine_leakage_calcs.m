%input
format long;
seal_radii_min = .257; %inches


clearance = .002; %inches
tooth_width = .01; %inches
cavity_width = .031; %inches
teeth = 10;
loxPin = 1000; %psi
loxPdrain = 50; %psi
loxPvent = 1000; %psi
loxRho = 1140; %kg/m^3
loxMu = 1.95e-4; %Pa s

atmP = 14.7; %psi

rp1Pin = 1000; %psi
rp1Pdrain = 60; %psi
rp1Pvent = 1000; %psi
rp1Rho = 800; %kg/m^3
rp1Mu = 1.92e-3; %Pa s




%change units
loxPin = 6894.76*loxPin;
loxPdrain = 6894.76*loxPdrain;
loxPvent = 6894.76*loxPvent;
rp1Pin = 6894.76*rp1Pin;
rp1Pdrain = 6894.76*rp1Pdrain;
rp1Pvent = 6894.76*rp1Pvent;
atmP = atmP * 6894.76;

step = step * .0245;
seal_radii_min = seal_radii_min * .0245;
clearance = clearance * .0245;
tooth_width = tooth_width * .0245;
cavity_width = cavity_width * .0245;
tooth_pitch = tooth_width + cavity_width;

seal_radii = [];
for i=1:sealNum
    seal_radii(i) = seal_radii_min + (i-1)*step;
end

mdotlox1 = laby_leakage(clearance, tooth_pitch, tooth_width, seal_radii(1), teeth, loxPin, loxPdrain, loxRho, loxMu); %kg/s
Qlox1 = mdotlox1 / loxRho;
mdotlox1 = mdotlox1*2.2046 %lbm/s
Qlox1 = Qlox1*264.171; %gal/s;

mdotlox2 = laby_leakage(clearance, tooth_pitch, tooth_width, seal_radii(2), teeth, loxPdrain, atmP, loxRho, loxMu);
Qlox2 = mdotlox2 / loxRho;
mdotlox2 = mdotlox2*2.2046 %lbm/s
Qlox2 = Qlox2*264.1721; %gal/s

mdotRp11 = laby_leakage(clearance, tooth_pitch, tooth_width, seal_radii(6), teeth, rp1Pin, rp1Pdrain, rp1Rho, rp1Mu); %kg/s
Qrp11 = mdotRp11 / rp1Rho *264.172;%gal/s
mdotRp11 = mdotRp11 * 2.2046 %lbm/s



mdotRp12 = laby_leakage(clearance, tooth_pitch, tooth_width, seal_radii(5), teeth, rp1Pdrain, atmP, rp1Rho, rp1Mu); %kg/s
Qrp11 = mdotRp12 /rp1Rho * 264.172;%gal/s
mdotRp12 = mdotRp12 * 2.2046   %lbm/s
