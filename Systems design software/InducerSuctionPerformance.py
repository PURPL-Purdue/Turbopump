import numpy as np
import matplotlib.pyplot as mpl
from CoolProp.CoolProp import PropsSI
import CoolProp
import yaml


## Conversion factors ##
psi2Pa = 6894.76
Ntolbf = 0.224809
lbtokg = 0.453592
kgpm3tolbpft3 = 0.062428
ft3togal = 7.48052
R2K = 0.555556
fttom = 0.3048
ft2tom2 = 0.092903
ft3tom3 = 0.0283168
gEnglish = 32.174 #(ft/s^2)

## YAML import ##
with open('pumps_params.yaml') as file:
    pumps_params = yaml.safe_load(file)

with open('rotordynamics_params.yaml') as file:
    rotordynamics_params = yaml.safe_load(file)

## Rotordynamics var initialize ##
N = rotordynamics_params['shaft_speed'] #(RPM)

## Systems variables initialize ##
lox_temp = pumps_params['lox_inlet_design_temperature'] #(R)
kero_temp = pumps_params['kero_inlet_design_temperature'] #(R)

lox_inlet_press = pumps_params['lox_inlet_design_pressure'] #(psi)
kero_inlet_press = pumps_params['kero_inlet_design_pressure'] #(psi)

lox_outlet_press = pumps_params['lox_outlet_design_pressure'] #(psi)
kero_outlet_press = pumps_params['kero_outlet_design_pressure'] #(psi)

lox_mdot = pumps_params['lox_design_mdot'] * 1 #(lbm/s)
kero_mdot = pumps_params['kero_design_mdot'] #(lbm/s)

lox_vapor_P = PropsSI('P','T', lox_temp*R2K,'Q',0,'oxygen') / psi2Pa #(psi)
kero_vapor_P = 0.031 #(psi)

lox_rho = kgpm3tolbpft3 * PropsSI('D','T', lox_temp*R2K,'P', 14.7*psi2Pa, 'oxygen') #(lbm/ft^3)
kero_rho = 50.3 #(lbm/ft^3)

kero_vdot = kero_mdot / (kero_rho) #(ft^3/s)
lox_vdot = lox_mdot / (lox_rho) #(ft^3/s)

lox_NPSHa = 144 * (lox_inlet_press - lox_vapor_P) / lox_rho

kero_NPSHa = 144 * (kero_inlet_press - kero_vapor_P) / kero_rho

## Pump hardware parameters ##
lox_inducer_D1 = 1.13 / 12 #(ft)
lox_inducer_DH1 = 0.41 / 12 #(ft)
lox_inducer_DH2 = 0.63 / 12 #(ft) 
lox_inducer_E_fraction = 0.1 #(dH_inducer/dH_total)

lox_A1 = np.pi * (((lox_inducer_D1 / 2)**2) - ((lox_inducer_DH1 / 2)**2)) #(ft^2)
lox_A2 = np.pi * (((lox_inducer_D1 / 2)**2) - ((lox_inducer_DH2 / 2)**2)) #(ft^2)

kero_D2 = 1.16 / 12 #(ft)
kero_DH2 = 0.63 / 12 #(ft)

kero_A2 = np.pi * (((kero_D2 / 2)**2) - ((kero_DH2 / 2)**2)) #(ft^2)

## Nss limit prediction ##
# LOX inducer
lox_inducer_headrise = 144 * (100) / lox_rho #(ft) As per the pumps CDR, 100psi
lox_impeller_NPSHa = lox_NPSHa + lox_inducer_headrise

lox_u1 = lox_inducer_D1 * np.pi * N / 60 #(ft/s)
lox_tau1 = 2 * (gEnglish) * (lox_NPSHa) / (lox_u1 ** 2)

lox_c_m1 = (lox_vdot) / (lox_A1) #(ft/s)
lox_phi1 = lox_c_m1 / lox_u1
lox_inducer_Nss_performance = 8147 * (lox_phi1 ** 0.5) * (lox_tau1 ** (-3/4)) * ((1 - (lox_inducer_DH1 / lox_inducer_D1)) ** 0.5)

lox_inducer_Nss = N * (((lox_mdot / lox_rho) * ft3togal * 60)** 0.5) / ((lox_NPSHa)**0.75)

# LOX impeller
lox_u2 = lox_inducer_D1 * np.pi * N / 60 #(ft/s)
lox_tau2 = 2 * (gEnglish) * (lox_impeller_NPSHa) / (lox_u2 ** 2)

lox_c_m2 = (lox_vdot) / (lox_A2) #(ft/s)
lox_phi2 = lox_c_m2 / lox_u2
lox_impeller_Nss_performance = 8147 * (lox_phi2 ** 0.5) * (lox_tau2 ** (-3/4)) * ((1 - (lox_inducer_DH2 / lox_inducer_D1)) ** 0.5)

lox_impeller_Nss = N * (((lox_mdot / lox_rho) * ft3togal * 60)** 0.5) / ((lox_NPSHa)**0.75)

# Kero
kero_u2 = kero_D2 * np.pi * N / 60 #(ft/s)
kero_tau2 = 2 * (gEnglish) * (kero_NPSHa) / (kero_u2 ** 2)

kero_Nss = N * (((kero_mdot / kero_rho) * ft3togal * 60)** 0.5) / ((kero_NPSHa)**0.75)

kero_c_m2 = (kero_vdot) / (kero_A2) #(ft/s)
kero_phi2 = kero_c_m2 / kero_u2
kero_impeller_Nss_performance = 8147 * (kero_phi2 ** 0.5) * (kero_tau2 ** (-3/4)) * ((1 - (kero_DH2 / kero_D2)) ** 0.5)


## Outuputs ##
print(f"Tip speed: {lox_u1:.2f}")
print(f"Meridional speed: {lox_c_m1:.2f}")
print(f"Inlet energy ratio: {lox_NPSHa / ((lox_c_m1** 2) / (2 * gEnglish)):.2f}\n")

print(f"Predicted LOX inducer Nss limit: {lox_inducer_Nss_performance:.2f}")
print(f"LOX inducer Nss: {lox_inducer_Nss:.2f}")
print(f"LOX NPSHa: {lox_NPSHa:.2f}")
print(f"Phi: {lox_phi1:.3f}\n")

print(f"Predicted kero impeller Nss limit: {kero_impeller_Nss_performance:.2f}")
print(f"Kero Nss: {kero_Nss:.2f}")
print(f"Kero NPSHa: {kero_NPSHa:.2f}")
print(f"Phi: {kero_phi2:.3f}\n")