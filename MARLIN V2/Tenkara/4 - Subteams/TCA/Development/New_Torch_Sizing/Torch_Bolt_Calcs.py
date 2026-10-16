# TCA Igniter Sizing Code
# 11/28/2025 : Created - Louis DeSano

# units standard: every variable is stored as SI, convert when inputting to external libraries if needed
# rocketCEA: uses English Engineering
# pyfluids: uses SI with Celsius

import numpy as np
import matplotlib.pyplot as plt
from rocketcea.cea_obj import CEA_Obj
import cea as cea
from pyfluids import Fluid, FluidsList, Input
import pandas as pd
import yaml

## Unit Conversions & Constants ##
n2lbf = 4.44822     # [N/lbf]

psi2Pa = 6894.76    # [psi/Pa]

lbm2kg = 0.453592   #[lbm/kg]

ft2m = 0.3048       # [ft/m]
m2in = 39.3701      # [m/in]

bar2psi = 14.503773773    # [psi/bar]
bar2pa = 100000           # [pa/bar]

R2K = 0.555556      # [R/K] (5/9)

# Constants
G0 = 9.81           # [m/s^2]
P_ATM = 101325      # [Pa]
T_AMB_CELSIUS = 20 # [deg C]
T_AMB_KELVIN = 293.15 # [deg K]

###############

with open("C:/Users/igoto/Downloads/GH1/Turbopump/MARLIN V2/Tenkara/4 - Subteams/TCA/Inputs/TCA_params.yaml") as f:
    tca_yaml = yaml.safe_load(f)

#pressure contact area defined until inner oring channel

max_radius = 0.45     # [in]
max_radius = max_radius / m2in     # [m]

max_area = np.pi * (max_radius)**2

Safety_factor = 2     #factor of safety for bolts

tca_chamber_pressure = tca_yaml['chamber_pressure']    # [psi]

tca_chamber_pressure = tca_chamber_pressure * psi2Pa   # [Pa]

pres_force = tca_chamber_pressure * max_area       # [Pa/m^2 = N]

print(f'The force acting is {pres_force} N')

######High strength steel bolts max yield strength
max_yield = 150000   # [psi]
max_yield = max_yield * psi2Pa
number_bolts = 3    #number of bolts

min_dia = np.sqrt((number_bolts / np.pi) * (Safety_factor * pres_force) / (max_yield))

print(f"The minimum minor diameter of the bolts should be {min_dia * m2in} inches")