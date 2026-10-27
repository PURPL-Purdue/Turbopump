import numpy as np
import math
import matplotlib.pyplot as plt
from rocketcea.cea_obj import CEA_Obj
import cea as cea
from pyfluids import Fluid, FluidsList, Input
import pandas as pd
import yaml


###############################
#DEFINE CONSTANTS
###############################

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

R_OX = 0.25984      #[kJ / kg - K]
R_F = .51828        #[kJ / kg - K]

R_OXJ = 259.84      #[J / kg - K]
R_FJ = 518.28        #[J / kg - K]

###############################
#DEFINE GLOBAL VARIABLES
###############################

FEED_TEMP = 293.15    #Degrees Kelvin
OX_DIA = 0.063          #inches
FUEL_DIA = 0.059        #inches
CD = 0.6
K_F = 1.31  # ch4 specific heat ratio
K_OX = 1.40 # gox specific heat ratio  








###############################
#DEFINE SUBFUNCTIONS
###############################






###############################
#MAIN FUNCTION
###############################

def main():

    G_ox = np.sqrt(K_OX / R_OXJ) * ((K_OX + 1) / 2) ** (-(K_OX + 1) / (2 * (K_OX - 1)))
    G_f = np.sqrt(K_F / R_FJ) * ((K_F + 1) / 2) ** (-(K_F + 1) / (2 * (K_F - 1)))

    min_of = 0.5
    max_of = 2.5
    of_step = 0.05

    min_mdot = 15      #g/s
    max_mdot = 30      #g/s
    mdot_step = 2.5

    min_pf = 150      #psia
    max_pf = 600      #psia
    pf_step = 5

    pc_target = 200.0     #psia
    
    of_sweep = np.arange(min_of, (max_of + of_step), of_step)
    mdot_sweep = np.arange(min_mdot, (max_mdot + mdot_step), mdot_step)
    pf_sweep = np.arange(min_pf, (max_pf + pf_step), pf_step) * psi2Pa

    for of in of_sweep:
        for pf in pf_sweep:
            


    df = pd.DataFrame({'O/F': [], 'MDOT (g/s)': [], 'Ox MDOT (g/s)': [], 'Fuel MDOT (g/s)': [], 'Choked Ox feed pressure (psia)': [], 'Choked Fuel feed pressure (psia)': [], 'Min ox stiffness': [], 'Min fuel stiffness': []})


    for mdot in mdot_sweep:
        for of in of_sweep:
            mdot_ox = (of / (of + 1)) * mdot
            mdot_f = mdot - mdot_ox

            stiffness_ox = 200        #%
            stiffness_f = 200         #%

            c_ratio_ox = (2 / (K_OX+1)) ** (K_OX / (K_OX - 1)) 
            c_ratio_f = (2 / (K_F+1)) ** (K_F / (K_F - 1)) 

            A_ox = np.pi * ((OX_DIA / (m2in * 2)))**2
            A_f = np.pi * ((FUEL_DIA / (m2in * 2)))**2

            p_ox = mdot_ox * np.sqrt(FEED_TEMP * R_OXJ / K_OX) / (CD * A_ox * ((K_OX + 1) / 2) ** (-(K_OX + 1) / (2 * (K_OX - 1))))
            p_ox = mdot_ox * np.sqrt(FEED_TEMP * R_OXJ / K_OX) / (CD * A_ox * ((K_OX + 1) / 2) ** (-(K_OX + 1) / (2 * (K_OX - 1))))
            p_f = mdot_f * np.sqrt(FEED_TEMP * R_FJ / K_F) / (CD * A_f * ((K_F + 1) / 2) ** (-(K_F + 1) / (2 * (K_F - 1))))

            p_line_ox = (1 + stiffness_ox) * c_ratio_ox**-1 * pc_target
            p_line_f = (1 + stiffness_f) * c_ratio_f**-1 * pc_target

    return