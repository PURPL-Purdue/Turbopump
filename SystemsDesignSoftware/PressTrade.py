## SUMMARY ##
# I put together this trade script to figure out if our kerosene pumps will need an inducer
# by running an Nss sweep based on inlet pressures.This is also the first script to make use 
# of the "source of truth" YAML that we are implementing into this project and acts as a 
# reference of how to use it. If you have any questions about this script, reach out to 
# Alejandro Diaz Contreras or current turbopump team leadership.

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

## Calcs ##
lox_tank_P = np.linspace(14.7, 50, 100)
lox_NPSHa = 144 * (lox_tank_P - lox_vapor_P) / lox_rho
lox_Nss = N * (((lox_mdot / lox_rho) * ft3togal * 60)** 0.5) / ((lox_NPSHa)**0.75)

kero_tank_P = np.linspace(14.7, 50, 100)
kero_NPSHa = 144* (kero_tank_P - kero_vapor_P) / kero_rho
kero_Nss = N * (((kero_mdot / kero_rho) * ft3togal * 60)** 0.5) / ((kero_NPSHa)**0.75)

## Display ##
mpl.plot(lox_tank_P, lox_Nss, label="LOX Nss", color="blue")
mpl.plot(kero_tank_P, kero_Nss, label="Kero Nss", color="red")
mpl.title("Nss for given inlet pressures")
mpl.xlabel("Inlet pressure (psi)")
mpl.ylabel("Nss")
mpl.legend()
mpl.xlim(14.7, 50)
mpl.ylim(0, 80000)
mpl.axhline(y=10000, color='black', linestyle="-", label="10,000 Nss")
mpl.axhline(y=55000, color='black', linestyle="-", label="55,000 Nss")
mpl.show()