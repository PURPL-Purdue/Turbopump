import numpy as np
import matplotlib.pyplot as mpl
from CoolProp.CoolProp import PropsSI
import CoolProp
import yaml

## Conversion factors ##
psi2Pa = 6894.76
kgpm3tolbpft3 = 0.062428
ft3togal = 7.48052
R2K = 0.555556
m2in = 39.3700787
gpm2m3s = 0.00006309019640344
N2lbf = 0.224809

#Define the head and Q for both props
H_kero = 1422.51
H_LOx = 1494.31

Q_kero = 66.49 #
Q_lox = 86.92 #

#Define temps for each surrogate
water_T = 527.67 #(R) just calculated 20 C in Rankine (room temp)
LN2_T = PropsSI('T','P', 101325,'Q', 0, 'Nitrogen')/R2K #(R) Saturation temp at atmospheric pressure (this is what prop will be conditioned to)

#Then get rho for each (surrogates! not the actual props) using coolprop
water_rho = PropsSI('D','T', water_T*R2K,'P', 150*psi2Pa, 'Water') #(kg/m^3) Keeping this SI to use in orifice equation
LN2_rho = PropsSI('D','T', LN2_T*R2K,'P', 200*psi2Pa, 'Nitrogen') #(kg/m^3) Keeping this SI to use in orifice equation

#Calculate orifices
water_P = H_kero * water_rho * kgpm3tolbpft3 / 144 + 150 #(psi) Water outlet pressure
LN2_P = H_LOx * LN2_rho * kgpm3tolbpft3 / 144 + 200 #(psi) LN2 outlet pressure
print(water_P)
print(LN2_P)

water_CdA = Q_kero * gpm2m3s / ((2*(water_P - 14.7)*psi2Pa/water_rho)**0.5) #(m^2)
LN2_CdA = Q_kero * gpm2m3s / ((2*(LN2_P - 14.7)*psi2Pa/LN2_rho)**0.5) #(m^2)

water_diam = 2*np.sqrt(water_CdA/(0.9*np.pi))*m2in #(in) Using a C_d of 0.9
LN2_diam = 2*np.sqrt(LN2_CdA/(0.9*np.pi))*m2in #(in) Using a C_d of 0.9

#Calculate reaction force
water_mdot = Q_kero * gpm2m3s * water_rho #(kg/s)
LN2_mdot = Q_kero * gpm2m3s * water_rho #(kg/s)

water_vexit = Q_kero * gpm2m3s / (water_CdA / 0.9) #(m/s)
LN2_vexit = Q_lox * gpm2m3s / (LN2_CdA / 0.9) #(m/s)

water_F = water_mdot * water_vexit * N2lbf #(lbf)
LN2_F = LN2_mdot * LN2_vexit * N2lbf #(lbf)

#Outputs
print(f'Water orifice diameter: {water_diam:.3f} in, reaction force: {water_F:.1f} lbf')
print(f'LN2 orifice diameter: {LN2_diam:.3f} in, reaction force: {LN2_F:.1f} lbf')