#this code is based on equation presented in this website:
#https://www.mathworks.com/help/hydro/ref/annularorifice.html

import math

c = .002 #clearance inches
r = .334 #inner radius inches
R = c + r #outer radius inches
pa = 800 #psi in
pb = 50 #psi out
nu = 8.5035e-5 #kinematic viscosity ft^2/s
l = .3 #shaft length inches
rho = .029 #density in lbm/in^3

nu = nu/144 #convert unit to in^2/s for final equation

volflow = (math.pi * ((R-r)**3) * (R) * (pa - pb)) / (6 * nu * rho * l)
massflow = rho * volflow
print(massflow)