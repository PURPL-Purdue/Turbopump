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


## YAML import ##
with open('pumps_params.yaml') as file:
    pumps_params = yaml.safe_load(file)

with open('rotordynamics_params.yaml') as file:
    rotordynamics_params = yaml.safe_load(file)

with open('TCA/TCA_params.yaml') as file:
    TCA_params = yaml.safe_load(file)


## Var initialize ##
N = rotordynamics_params['shaft_speed'] #(RPM)

lox_temp = pumps_params['lox_inlet_design_temperature'] #(R)
kero_temp = pumps_params['kero_inlet_design_temperature'] #(R)

lox_mdot = pumps_params['lox_design_mdot'] #(lbm/s)
kero_mdot = pumps_params['kero_design_mdot'] #(lbm/s)

lox_vapor_P = PropsSI('P','T', lox_temp*R2K,'Q',0,'oxygen')/psi2Pa #(psi)
lox_rho = kgpm3tolbpft3 * PropsSI('D','T', lox_temp*R2K,'P', 14.7*psi2Pa, 'oxygen') #(lbm/ft^3)

kero_vapor_P = 0.031 #(psi)
kero_rho = 50.3 #(lbm/ft^3)

#Pulling from the TCA YAML
TCA_mdot = TCA_params['turbopump_mdot'] #(lbm/s)
TCA_of = TCA_params['oxidizer_fuel_ratio'] 
lox_pressure = TCA_params["Injector inlet pressure"]["LOX"] #(psi)
kero_pressure = TCA_params["Injector inlet pressure"]["RP1"] #(psi)
lox_inlet_P = pumps_params['lox_inlet_design_pressure'] #(psi)
kero_inlet_P = pumps_params['kero_inlet_design_pressure'] #(psi)

## OK, now set the params for the test campaign ##
loxPump_water_T = 50 + 459.67 #(R)
loxPump_LN2_T = PropsSI('T','P', 14.7*psi2Pa,'Q',0,'water')/R2K #(R)

keroPump_water_T = 50 + 459.67 #(R)

kero_vdot = 60 * ft3togal * ((TCA_mdot / (TCA_of + 1)) + 1) / kero_rho #(ft^3/s) 
lox_vdot = 60 * ft3togal * ((TCA_mdot - TCA_mdot / (TCA_of + 1)) + 0.25) / lox_rho #(ft^3/s)

lox_headrise = (lox_pressure - lox_inlet_P) * 144 / lox_rho
kero_headrise = (kero_pressure - kero_inlet_P) * 144 / kero_rho

print(f'LOx headrise: {lox_headrise:.2f} ft, Kero headrise: {kero_headrise:.2f} ft')
print(f'LOx volumetric flowrate: {lox_vdot:.2f} gpm, Kero volumetric flowrate: {kero_vdot:.2f} gpm')