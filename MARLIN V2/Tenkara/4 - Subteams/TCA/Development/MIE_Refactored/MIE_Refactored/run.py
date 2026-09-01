######################################################################################
# Author: Louis DeSano
# Latest Revision: April 7, 2026
#
# This set of functions determines:
#   - Minimum ignition energy
#   - Chemical Power Addition
#   - Igniter mass flow rate
#
# Outputs: (Found in Results Folder Created at Location of setpoints.yaml)
# vs_OF.png/csv: 
#   - Torch Energy Release, Combustion Temperature, and MIE mdot vs. Torch OF
# reaction.png: 
#   - Well Stirred Reactor Ignition vs. Time
#
# Usage:
# Call run.py in command line with the setpoints yaml as its only argument:
# example: python3 run.py setpoints.yaml
#
# Dependencies:
# Libraries
#   - Cantera, Coolprop, Numpy, Matplotlib, RocketCEA
# Files:
# - MIE.py, WSR.py, TCA_params.yaml
######################################################################################

# Import Dependent Libraries
import cantera as ct
import CoolProp as cp
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import time
from pyfluids import Fluid, FluidsList, Input
from rocketcea.cea_obj import CEA_Obj
import yaml
import sys
import warnings
import os
warnings.filterwarnings("ignore")

IN2M = 0.0254
FT2IN = 12

# Import Combustion Codes
from WSR import Auto_Igntition_Temp
from MIE import chain_reaction_power, plot_OF_sweep

#YOU WERE GETTING ENTHALPY FOR IPA AT EACH STATE FROM THE WEBSITE
"""
import Torch_OF"""

#Yaml Safe Open Function
def load_config(config_path):
    with open(config_path, 'r') as f:
        # Use safe_load for security when dealing with untrusted input
        config_data = yaml.safe_load(f)
    return config_data

# Parse Setpoints From YAML
TCA_config = load_config('TCA_params.yaml')

OF = TCA_config['oxidizer_fuel_ratio']
mdot_main = TCA_config['turbopump_mdot']
Pc_main = TCA_config['tca_chamber_pressure']
lstar = TCA_config['characteristic_length']
Dt = TCA_config['tca_throat_diameter']
expansion_ratio = TCA_config['tca_expansion_ratio']
eta_cstar = TCA_config['c_star_efficiency']

# read in simulation configs
sim_config_path = sys.argv[1]
#sim_config_path = "Refactored_MIE/IPA/setpoints.yaml"
sim_config = load_config(sim_config_path)
config_dir = os.path.join(os.path.dirname(sim_config_path), 'Results')
os.makedirs(config_dir, exist_ok=True)

MCA_inlet = sim_config["MCA_inlet"]
mech_MCA = sim_config["mech_MCA"]
cea_fuel = sim_config["CEA_fuel"]
cea_ox = sim_config["CEA_ox"]

mech_torch = sim_config["mech_torch"]
fuel_torch = sim_config["torch_fuel"]
ox_torch = sim_config["torch_ox"]
OF_torch = sim_config["OF_torch"]
pc_torch = sim_config["pc_torch"]
bisect_count = sim_config["WSR_bisection"]

# determine stay time for WSR model
# get density of chamber 
cea = CEA_Obj(oxName=cea_ox, fuelName=cea_fuel)
rho_main = cea.get_Chamber_Density(Pc=Pc_main, MR=OF, eps=1) / (12**3)
cstar = cea.get_Cstar(Pc=Pc_main, MR=OF) * eta_cstar * FT2IN
At = np.pi * (Dt/2)**2
Vc = lstar * At
ts = Vc * rho_main / mdot_main
print(f"Stay Time [s]: {ts:0.5f}")


# Well Stirred Reactor Model to Find AIT of Propellants for Chamber conditions:
# Stay Time, Chamber Pressure, OF Ratio
T_ign = Auto_Igntition_Temp(bisect_count=bisect_count, mdot=mdot_main, OF=OF, residence_time=ts, Pc=Pc_main, mech=mech_MCA, prop_temps=MCA_inlet, save_path=config_dir)

if sim_config['OF_sweep']:
    OF_min = sim_config['OF_min']
    OF_max = sim_config['OF_max']
    OF_step = sim_config['OF_step']
    OF_set = np.arange(OF_min,OF_max,OF_step)
    energy = []
    temp = []
    mdot = []

    for OF_torch in OF_set:
        mdot_torch, P_ig, Torch_E_release, Tc_torch, _ = chain_reaction_power(config=sim_config, T_ign=T_ign, 
                    mech_torch=mech_torch, fuel_torch=fuel_torch, ox_torch=ox_torch, OF_torch=OF_torch, pc_torch=pc_torch, 
                    mech_MCA=mech_MCA, MCA_INLET=MCA_inlet, mdot_MCA=mdot_main, OF_MCA=OF, pc_MCA=Pc_main)
        energy.append(Torch_E_release * 1e-6)
        temp.append(Tc_torch)
        mdot.append(mdot_torch)

    plot_OF_sweep(OF_set=OF_set, E_ign=energy, T_comb=temp, mdot=mdot, file_path=config_dir)

    # save data set to sheet
    df = pd.DataFrame({
        "OF_torch": OF_set,
        "Energy (MJ)": energy,
        "Temperature (K)": temp,
        "mdot (kg/s)": mdot
    })
    sheet_path = os.path.join(config_dir, "vs_OF.csv")
    df.to_csv(sheet_path, index=False)

else:
    mdot_torch, P_ig, _, _, _ = chain_reaction_power(config=sim_config, T_ign=T_ign, 
            mech_torch=mech_torch, fuel_torch=fuel_torch, ox_torch=ox_torch, OF_torch=OF_torch, pc_torch=pc_torch, 
            mech_MCA=mech_MCA, MCA_INLET=MCA_inlet, mdot_MCA=mdot_main, OF_MCA=OF, pc_MCA=Pc_main)
