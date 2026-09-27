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
import math
import matplotlib.pyplot as plt

#######################################
#Reading in Hardware Definition
#######################################
#def_path = os.path.join(r'\MARLIN V2\Tenkara\4 - Subteams\Gas Generator\Inputs\GG_hardware_definitions.yaml')
#def_path = 'C:\Users\miles\OneDrive\Documents\GitHub\Turbopump\MARLIN V2\Tenkara\4 - Subteams\Gas Generator\Inputs\GG_hardware_definitions.yaml'
script_dir = os.path.dirname(os.path.abspath(__file__))
parent_dir = os.path.dirname(script_dir)  # Go up one level
def_path = os.path.join(parent_dir, 'Inputs', 'GG_hardware_definition.yaml')
print(def_path)
config = yaml.safe_load(open(def_path, 'r'))
#######################################
# Function Definitions
#######################################
cea = CEA_Obj(config['GG_fuel_type'], config['GG_ox_type'], config['GG_fuel_type']) 


def mdot(cda, pup, Pc, rho):
    #print(f'cda: {cda} pup: {pup} Pc: {Pc} rho: {rho}')
    Pc = np.nan_to_num(Pc, nan=0.0)
    Pc = 0 if Pc is None else Pc
    mdot_out = cda * np.sqrt(2 * rho * max((pup - Pc),0))
    return float(mdot_out)

def UpdateCEA_Cp(Pc, mdot_ox, mdot_IPA):
    OF = mdot_ox/max(mdot_IPA, 0.0000000001)
    Cp = 4.1868 * cea.get_Chamber_Cp(Pc/psi2pa, OF, 1, 0)
    return Cp

def IPA_Cp_gas(T):
    #Returns the Cp of IPA in kJ/kg K
    #input is Kelvin
    Cp = config['cp_fuel_g_m'] * T + config['cp_fuel_g_b']
    return 1000 * Cp / IPA_MM
#Feed Pressures

def CEA_init_state(mdot_ox, mdot_ipa, Pc):
    OF = mdot_ox/max(mdot_ipa, 0.0000000001)
    Cp = 4.1868 * cea.get_Chamber_Cp(Pc/psi2pa, OF, 1,0)
    #print(f'Initial state - OF: {OF}, Cp: {Cp}')
    h = -2.326 * cea.get_Enthalpies(Pc/psi2pa, OF, 1, 0, 0)[0]
    T = 0.555555 * cea.get_Tcomb(Pc/psi2pa, OF)
    mass = mdot_ipa + mdot_ox
    
    H = mass * h
    return [mass, H, T, Cp]

def equilibrium(state1, state2):
    M1, M2 = state1[0], state2[0]
    H1, H2 = state1[1], state2[1]
    T1, T2 = state1[2], state2[2]
    Cp1, Cp2 = state1[3], state2[3]
    Mf = M1+M2
    Hf = H1+H2
    Cp = (M1*Cp1 + M2*Cp2)/(M1+M2)
    Tf = Hf/(Mf*Cp)
    return [Mf, Hf, Tf, Cp]

def Pc(LOx_feed, IPA_feed, Pc_in, pcA):
    #Recursivly solves for chamber pressure
    mdot_ox = mdot(config['LOx_Cda'],  LOx_feed, Pc_in, config['LOx_rho'])
    mdot_ipa = mdot(config['fuel_Cda'], IPA_feed, Pc_in, config['IPA_rho'])
    OF = mdot_ox/max(mdot_ipa, 0.000000000001)
    #print(f'mdots {mdot_ox} and {mdot_ipa}')
    cstar = 0.3048 * cea.get_Cstar(Pc_in/psi2pa, OF)
    Pc_f = cstar * (mdot_ox+mdot_ipa) / config['stator_throat_area']
    pcA.append(Pc_f)
    #print(len(pcA))
    #if len(pcA) == 500:
       #plt.plot(pcA) 
       #plt.show()
       #plt.pause(1000)
    if len(pcA) >= 1000:
        raise RuntimeError("Pressure solver did not converge")
    # d(error) / dt threshhold 
    if abs(Pc_in-Pc_f) < 3000:
        #plt.plot(pcA)
        #plt.show()
        return Pc_f
    
    else:
        #print(f'CHamber pressure {Pc_f}')
        error = Pc_f - Pc_in
        new_Pc = Pc_in + 0.05 * error
        return Pc(LOx_feed, IPA_feed, new_Pc, pcA)


def Heat_addition(state1, E_in, flag = 0):
    M1 = state1[0]
    H1 = state1[1]
    T1 = state1[2]
    Cp1 = state1[3]
    Tf =  T1 + (E_in)/(M1*Cp1)
    
    
    return[M1, H1+E_in, Tf, Cp1]

def IPA_satTemp(Pc):
    Pc = Pc / 100000
    A = 4.57795
    B = 1221.423
    C = -87.474
    # Inverted Antoine equation solving for T in K
    return (B / (A - math.log10(Pc))) - C

#####################################
# Variable Definition
#####################################

step_size = 0.1
psi2pa = 6894.757
# Feed pressures
LOx_feed_pressure = 3430000     #pa
IPA_feed_pressure = 3430000     #pa
# IPA inlet temp
IPA_temp = 298                   # K
#IPA_H_mol = 318.2                # Enthalpy at 298K kJ/mol
IPA_MM = 60.1                    # g/mol
#IPA_H = 1000 * IPA_H_mol / IPA_MM# kJ/kg
# Chamber Pressure
Pc_guess = 2400000              #pa  
Pc_array = []
film_temp = []
gas_temp = []

##########################################
# Define input states
##########################################
# Solve for chamber pressure
Pc_in = Pc(LOx_feed_pressure, IPA_feed_pressure, Pc_guess, Pc_array)
# Solve for mass flows into the chamber
mdot_ox = mdot(config['LOx_Cda'], LOx_feed_pressure, Pc_in, config['LOx_rho'])
mdot_IPA = mdot(config['fuel_Cda'], IPA_feed_pressure, Pc_in, config['IPA_rho'])
mdot_film = mdot(config['film_Cda'], IPA_feed_pressure, Pc_in, config['IPA_rho'])
print(f'IPA core mdot: {mdot_IPA:.3f}    Ox mdot: {mdot_ox:.3f}     film mdot: {mdot_film:.3f}')
print(f'core OF Ratio: {(mdot_ox/mdot_IPA):.3f}')
print(f'Total OF Ratio: {(mdot_ox/(mdot_IPA + mdot_film)):.3f}')
print(f'film cooling %: {100*(mdot_film / (mdot_film + mdot_ox+mdot_IPA)):.3f}')
# Get combustion gas properties
combustion_gas = CEA_init_state(mdot_ox, mdot_IPA, Pc_in)
# Define film cooling initial properties
Cp_IPA_l = 1000 * config['cp_IPA_l'] / IPA_MM             # 1000 * kJ/mol K * mol/ g = kJ/kg K
film = [mdot_film, 0, IPA_temp, Cp_IPA_l]
delta_T_init = abs(combustion_gas[2] - film[2])
film_temp.append(film[2])
gas_temp.append(combustion_gas[2])
#####################################
# Begin Simulation
#####################################
# subcooled liquid heating of the IPA
T_dif = IPA_satTemp(Pc_in) - film[2]
##print(f'Temperature difference: {T_dif:.3f} K')
E_dif = T_dif * film[3] * film[0]
#print(f'Energy difference: {E_dif:.3f} kJ')
film = Heat_addition(film, E_dif)
combustion_gas = Heat_addition(combustion_gas, -E_dif)
film_temp.append(film[2])
gas_temp.append(combustion_gas[2])
# vaporization of the IPA
E_vap = film[0] * (1000 * config['h_vap_IPA'] / IPA_MM)
film[1] = film[1] + E_vap
combustion_gas = Heat_addition(combustion_gas, -E_vap)
film_temp.append(film[2])
gas_temp.append(combustion_gas[2])
film[3] = IPA_Cp_gas(film[2])
step_num = 0
while abs(film[2] - combustion_gas[2]) > 1 and step_num < 2000:
    E_trans = 0.5 * step_size * (combustion_gas[2] - film[2])
    #print(f'E_trans {E_trans:.3f}')
    film = Heat_addition(film, E_trans)
    combustion_gas = Heat_addition(combustion_gas, -E_trans)
    film_temp.append(film[2])
    gas_temp.append(combustion_gas[2])
    film[3] = IPA_Cp_gas(film[2])
    step_num += 1
#print('flag 7')
OutPut_gas = [film[0] + combustion_gas[0], film[1] + combustion_gas[1], film[2],
               (film[0]*film[3]+combustion_gas[0]*combustion_gas[3])/(film[0]+combustion_gas[0])]

plt.plot(film_temp)
plt.plot(gas_temp)
plt.legend(['film','combustion gas'])
plt.grid()
plt.show()

print(f'Output mdot:  {OutPut_gas[0]:.3f} kg/s')
print(f'Output Temp:  {OutPut_gas[2]:.3f} K')
print(f'Chamber Pressure:  {Pc_in/psi2pa:.1f} psi')
    
    


