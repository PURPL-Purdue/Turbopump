## GG V2.0 Film Cooling Sizing Script
# Authors: Miles Cohen & Lily Brideou
# Updated on: 9/23/2026
# Description: Script for sizing film cooling in the Gas Generator V2.0
# Method: Given input mass flow rates and temperatures, this script
# uses CEA to model the combustion excluding any fild cooling, then 
# assumes the film cooling and combustion products are fully gassous and well
# mixed and reaches thermal equilibrium before exiting the combustion chamber.

import numpy as np
import pandas as pd
import yaml 
from rocketcea.cea_obj import CEA_Obj
import CoolProp as CP
import os

def_path = os.path.join(r'\MARLIN V2\Tenkara\4 - Subteams\Gas Generator\Inputs\GG_hardware_definitions.yaml')
config = yaml.safe_load(open(def_path, 'r'))

def mdot(cda, pup, pdows, rho):
    return cda * np.sqrt(2 * rho * (pup - pdows))

def getCEA_h(Pc, OF):
    cea = CEA_Obj(oxname = config['GG_ox_type'], fuelname = config['GG_fuel_type'])
    H_list = cea.get_Enthalpies(Pc, OF,1, 0, 0)
    return cea

def IPA_Cp_gas(T):
    return config['cp_fuel_g_m'] * T + config['cp_fuel_g_b']

