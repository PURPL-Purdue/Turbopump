# Torch_Analysis
import numpy as np
import matplotlib.pyplot as plt
from rocketcea.cea_obj import CEA_Obj
import cea as cea
from pyfluids import Fluid, FluidsList, Input
import pandas as pd

## Define torch parameters
At = 1     # Throat area
S_f = 1     # Fuel injector size
S_o = 1     # Oxidizer injector size
p_c = 1     # Torch Chamber Pressure
T_c = 1     # Torch Chamber Temperature

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

