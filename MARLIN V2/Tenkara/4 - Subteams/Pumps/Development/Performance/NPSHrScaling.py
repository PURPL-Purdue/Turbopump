from pint import Quantity as Q_
import numpy as np
from matplotlib import pyplot as plt
import CoolProp.CoolProp as cp

waterVaporPressure = Q_(cp.PropsSI('P','T',293,'Q', 0, 'Water'), 'Pa')
waterDensity = Q_(cp.PropsSI('D','T', 293, 'P', 101325, 'Water'), 'kg/m^3')
g = Q_(9.81, 'm/s^2')

vaporHead = waterVaporPressure / (waterDensity * g)

N_1 = Q_(35000, 'rpm')
N_2_arr = Q_(np.linspace(0, 35000, 100), 'rpm')
NPSHi_1 = Q_(114, 'm')
NPSHi_2_arr = []
NPSHi_pressure_arr = []

closestN = Q_(0, 'rpm')
best_dist_from_curve = Q_(1000, 'm')

def findNPSHr(N1,N2,NPSH1):
    NPSH2 = ((N2 / N1).to('dimensionless') ** 2) * NPSH1
    return NPSH2

def findPress(density, head, g):
    p = density * head * g
    return p

for N2 in N_2_arr:
    NPSH_loop = findNPSHr(N_1,N2,NPSHi_1)
    press_loop = findPress(waterDensity, NPSH_loop, g)

    NPSHi_2_arr.append(NPSH_loop)
    NPSHi_pressure_arr.append(press_loop.to('psi'))

    if(np.abs(NPSH_loop - vaporHead) < best_dist_from_curve):
        closestN = N2
        best_dist_from_curve = np.abs(findNPSHr(N_1,closestN,NPSHi_1) - vaporHead)



NPSH_vals = [j.magnitude for j in NPSHi_2_arr]
press_vals = [p.magnitude for p in NPSHi_pressure_arr]

plt.plot(N_2_arr, NPSH_vals, label='NPSHi')
plt.plot(N_2_arr, press_vals, label='NPSHi Pressure')
plt.plot([closestN.magnitude], [findNPSHr(N_1,closestN,NPSHi_1).magnitude], 'o', label='Water vapor head')
plt.xlabel('Speed (rpm)')
plt.ylabel('NSPHi (m)')
plt.legend()
plt.grid()
plt.show()