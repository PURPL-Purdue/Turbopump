import numpy as npy



mToFeet = 3.28084  
lsToGalMin = 15.85032

volumetric_flow = 26.81081
NPSH = 10 * mToFeet
NPSH_r = 150
Q = 1.69 * lsToGalMin # volumetric flow converted from l/s, gal/min 
n = 40000 # roations speed, rpm 
rho = 0.3 # hub to tip ratio 
H = 140 * mToFeet # head rise meters to feet 
g = 32.2






N_ss = (n * Q)**(1/2)/ (NPSH**(3/4)) # Suction specific speed (Huzel and Huang 192 
N_ssP = N_ss / (1 - rho**2)**(1/2) # Corrected Suction specific speed 
phi = (3574/N_ssP) / ((1 + (1 + 6 * (3574 / N_ssP) ** 2) ** (1/2)) / 2)  
D_t = 0.37843 * (Q/ ((1-rho**2) * n * phi)) ** (1/3)
D_h = D_t * rho 
## TODO Velocity triangle
delta_H = H 
u = (npy.pi * D_t * n) / 60
eta_bl = .85
C_u = (g * delta_H) / (eta_bl * u)
C_m = u * phi
gamma = npy.atan(C_m/(u-C_u))
beta = gamma/0.575
## TODO delta = 2 * npy.pi * r * npy.tan(gamma)  needs velocity triangle

Ns  = (n * (Q) ** 0.5) / (delta_H) ** .75
N = 3 
S = (npy.pi* D_h) / 3 
sigma = 2.5 
C = sigma * S 
# TODO Beta =  outlet blade angle-inlet blade angle


# USE PINT!!!! OR ELSE IM GOING TO KILL SOMEONE