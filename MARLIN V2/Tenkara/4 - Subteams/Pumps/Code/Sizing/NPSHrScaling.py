from pint import Quantity as Q_
import numpy as np
from matplotlib import pyplot as plt
import CoolProp.CoolProp as cp

waterVaporPressure = Q_(cp.PropsSI('P','T',293,'Q', 0, 'Water'), 'Pa')
waterDensity = Q_(cp.PropsSI('D','T', 293, 'P', 101325, 'Water'), 'kg/m^3')
g = Q_(9.81, 'm/s^2')

vaporHead = waterVaporPressure / (waterDensity * g)
atmHead = Q_(101325,'Pa') / (waterDensity * g)

N_1 = Q_(35000, 'rpm')
N_2_arr = Q_(np.linspace(0, 35000, 100), 'rpm')
NPSHi_1 = Q_(114, 'm')
P_1 = Q_(13,'kW')
Q_1 = Q_(1.78, 'liter / second')
P_2_arr = []
NPSHi_2_arr = []
NPSHi_pressure_arr = []
Q_2_arr = []
mdot_2_arr = []

vaporN = Q_(0, 'rpm')
atmN = Q_(0, 'rpm')
best_dist_from_curve_vap = Q_(1000, 'm')
best_dist_from_curve_atm = Q_(1000, 'm')

def findNPSHr(N1,N2,NPSH1):
    NPSH2 = ((N2 / N1).to('dimensionless') ** 2) * NPSH1
    return NPSH2

def findPress(density, head, g):
    p = density * head * g
    return p

def findPower(N1,N2, P1):
    P2 = P1 / ((N1/N2).to('dimensionless') ** 3)
    return P2

def findQ(N1,N2, Q1):
    Q2 = ((N2/N1).to('dimensionless')) * Q1
    return Q2

for N2 in N_2_arr:
    # Calc new values for current index in rpm array
    NPSH_loop = findNPSHr(N_1,N2,NPSHi_1)
    press_loop = findPress(waterDensity, NPSH_loop, g)
    P_loop = findPower(N_1, N2, P_1)
    Q_loop = findQ(N_1, N2, Q_1)
    mdot_loop = (Q_loop * waterDensity).to('kg/second')

    # Add new values to output arrays
    NPSHi_2_arr.append(NPSH_loop)
    NPSHi_pressure_arr.append(press_loop.to('psi'))
    P_2_arr.append(P_loop)
    Q_2_arr.append(Q_loop)
    mdot_2_arr.append(mdot_loop)

    if(np.abs(NPSH_loop - vaporHead) < best_dist_from_curve_vap):
        vaporN = N2
        best_dist_from_curve_vap = np.abs(findNPSHr(N_1,vaporN,NPSHi_1) - vaporHead)

    if(np.abs(NPSH_loop - atmHead) < best_dist_from_curve_atm):
        atmN = N2
        best_dist_from_curve_atm = np.abs(findNPSHr(N_1,atmN,NPSHi_1) - atmHead)



NPSH_vals = [p.magnitude for p in NPSHi_2_arr]
press_vals = [p.magnitude for p in NPSHi_pressure_arr]
power_vals = [P.magnitude for P in P_2_arr]
Q_2_arr = [p.magnitude for p in Q_2_arr]
mdot_2_arr = [p.magnitude for p in mdot_2_arr]

plt.close('all')

plt.figure()
plt.plot(N_2_arr, NPSH_vals, label='NPSHi (m)')
plt.plot(N_2_arr, press_vals, label='NPSHi Pressure (PSI)')
plt.plot([vaporN.magnitude], [findNPSHr(N_1,vaporN,NPSHi_1).magnitude], 'o', label='Water vapor head')
plt.plot([atmN.magnitude], [findNPSHr(N_1,atmN,NPSHi_1).magnitude], 'o', label='Atmospheric pressure head')
plt.xlabel('Speed (rpm)')
plt.ylabel('NSPHi (m) / PSI')
plt.legend()
plt.grid()

plt.figure()
plt.axvline(x=atmN.magnitude, color='Gray', linestyle='--', linewidth=2)
plt.plot(N_2_arr, power_vals, label='Required Power (kW)', color='Blue')
plt.xlabel('Speed (rpm)')
plt.ylabel('Required Power (kW)')
plt.twinx()
plt.plot(N_2_arr, mdot_2_arr, label='Required Mass Flow (kg/s)', color='Orange')
plt.ylabel('Required Mass Flow (kg/s)')
plt.legend()
plt.grid()

plt.show(block=False)