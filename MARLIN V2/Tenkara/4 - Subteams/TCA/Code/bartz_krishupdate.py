"""
Reads Outputs/propreties.csv and calculates the film coefficient in function of axial position.
"""

import pandas as pd
from pint import UnitRegistry
import numpy as np
import yaml
import matplotlib.pyplot as plt
from scipy.interpolate import CubicSpline

#######################################################################################
# Read properties.csv and read dataframe
# =====================================================================================

df = pd.read_csv('Outputs/properties.csv')
ureg = UnitRegistry()

with open('Inputs/TCA_params.yaml') as f:
    p = yaml.safe_load(f)
    pc = p['chamber_pressure'] * ureg.psi

x = df['x_m'].values * ureg.m
y = df['y_m'].values * ureg.m
gamma = df['gamma'].values
Cp = df['Cp [(KJ/kg-K)]'].values * ureg.kilojoule / (ureg.kilogram * ureg.kelvin)
k = df['k [W/m-K]'].values * ureg.watt / (ureg.meter * ureg.kelvin)
MW = df['MW [kg/kmol]'].values * ureg.kilogram / ureg.kmole
T = df['T_local [K]'].values * ureg.kelvin
P = df['P_local [bar]'].values * ureg.bar
mu = df['Viscosity [Pa*s]'].values * ureg.pascal * ureg.second
Pr = df['Prandtl Number'].values
M = df['Mach'].values
Tw = df['Tw [K]'].values * ureg.kelvin   

rt = np.interp(0.0, x, y)
Dt = 2 * rt
At = np.pi * rt ** 2

A = np.pi * y ** 2

cstar = p['cstar'] * ureg.meter / ureg.second

########################################################################################
# Calculate the gas-side heat transfer coefficient for one poistion as a function of the local properties
# ====================================================================================

def calc_hg_point(Dt, mu, Cp, Pr, pc, cstar, rc, Tw, T, gamma, M):
    """
    Calculates gas-side heat transfer coefficient hg
    at one point in the thrust chamber.

    Inputs:
        Dt     = throat diameter [m]
        mu     = dynamic viscosity [Pa*s]
        Cp     = specific heat [J/kg-K]
        Pr     = Prandtl number
        pc     = chamber pressure [Pa]
        cstar  = characteristic velocity [m/s]
        rc     = throat radius of curvature [m]
        Tw     = local wall temperature [K]
        T      = local static gas temperature [K]
        gamma  = specific heat ratio
        M      = local Mach number

    Returns:
        hg
    """

    Dt = Dt.to(ureg.meter)
    mu = mu.to(ureg.pascal * ureg.second)
    Cp = Cp.to(ureg.joule / (ureg.kilogram * ureg.kelvin))
    pc = pc.to(ureg.pascal)
    cstar = cstar.to(ureg.meter / ureg.second)
    rc = rc.to(ureg.meter)
    Tw = Tw.to(ureg.kelvin)
    T = T.to(ureg.kelvin)

    # isentropic flow area ratio
    areaRatio = (1/M) * ((2/(gamma+1)) * (1 + ((gamma-1)/2)*(M**2)))**((gamma+1)/(2*(gamma-1)))
    
    isentropic_ff = 1 + ((gamma - 1) / 2) * M**2

    sigma = 1 / ((0.5 * ((Tw / T) + 0.5)**0.68)* isentropic_ff**0.12)

    hg = ((0.026 / Dt**0.2) * (mu**0.2 * Cp / Pr**0.6) * (pc / cstar)**0.8 * (Dt / rc)**0.1 * (areaRatio)**0.9 * sigma)

    return hg

########################################################################################
# Calculate the slope of the viscosity-temperature curve
# ====================================================================================

T_vals = T.to(ureg.kelvin).magnitude
mu_vals = mu.to(ureg.pascal * ureg.second).magnitude

w = np.gradient(np.log(mu_vals), np.log(T_vals))


########################################################################################
# Calculate the radius of curvature of the nozzle wall
# ====================================================================================

spline = CubicSpline(x, y)

dy_dx = spline(0.0, 1)
d2y_dx2 = spline(0.0, 2)

rc = (1 + dy_dx**2)**1.5 / abs(d2y_dx2)

#######################################################################################
# Calculate heat transfer coefficient and recovery temperature
# =====================================================================================



isentropic_ff = 1 + (gamma - 1) / 2 * M ** 2
sigma = 1/((0.5*((Tw/T) + 0.5)**(.68))*isentropic_ff**(.12))
hg = (0.026/Dt**0.2)*(mu**0.2*Cp.to(ureg.joule/(ureg.kilogram*ureg.kelvin))/Pr**0.6)*(pc.to(ureg.pascal)/cstar)**0.8*(Dt/rc)**0.1*(At/A)**0.9*sigma

Tr = T * (1 + ((gamma - 1) / 2 * Pr**(1/3) * M ** 2))

#print(Tr.to(ureg.kelvin))
#print(hg.to(ureg.watt/(ureg.meter**2*ureg.kelvin)))


#plt.plot(x.magnitude, hg.magnitude)

#plt.xlabel("Axial position, x [m]")
#plt.ylabel("Gas-side heat transfer coefficient, $h_g$ [W/m²K]")
#plt.grid(True)
#plt.show()

#plt.plot(x.magnitude, Tr.magnitude)

#plt.xlabel("Axial position, x [m]")
#plt.ylabel("Recovery Temperature, $T_r$ [K]")
#plt.grid(True)
#plt.show()

#######################################################################################
# Fit heat transfer coefficient and recovery temperature to a spline
# =====================================================================================


#Getting correct units and set them to magnitudes
x_data = x.to(ureg.meter).magnitude
hg_data = hg.to(ureg.watt / (ureg.meter**2 * ureg.kelvin)).magnitude
Tr_data = Tr.to(ureg.kelvin).magnitude

#Deleting the repeated indices
x_unique, indices = np.unique(x_data, return_index=True)
hg_unique = hg_data[indices]
Tr_unique = Tr_data[indices]

#Generating the splines
splinehg = CubicSpline(x_unique, hg_unique)
hg_fit = splinehg(x_data)

splineTr = CubicSpline(x_unique, Tr_unique)
Tr_fit = splineTr(x_data)

#######################################################################################
# Write heat transfer coefficient and recovery temperature to csv file
# =====================================================================================

df2 = pd.DataFrame({
    'x_m': x.magnitude,
    'hg_W_m2K': hg.to(ureg.watt/(ureg.meter**2*ureg.kelvin)).magnitude,
    'Tr_K': Tr.to(ureg.kelvin).magnitude})

#df2.to_csv('Outputs/heat_transfer_parameters.csv', index=False)

#######################################################################################
# Fitting the spline to get a larger database of points
# =====================================================================================

x_ansys = np.linspace(x_data.min(), x_data.max(), 500)

hg_ansys = splinehg(x_ansys)
Tr_ansys = splineTr(x_ansys)

np.savetxt(
    "Outputs/ansys_input2.csv",
    np.column_stack((x_ansys, hg_ansys,Tr_ansys)),
    delimiter=","
)