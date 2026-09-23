# ## Adaptation of Hunstville Ignition Energy Paper for Turbopump TCA
import numpy as np
import cantera as ct
import CoolProp as cp
import pandas as pd
import yaml
import matplotlib.pyplot as plt
import os

PSI2PA = 6894.76 # psi to pa
R2K = 0.555556 # Rankine to Kelvin
BTU2J = 1055.06 # BTU to J
LBM2KG = 0.453592 # lbm to kg

T_amb = 298 # K
# function to get masss averaged chemical potential of mixture 
def mass_avg_chem_pot(mixture):
    mu = mixture.chemical_potentials      # J/kmol
    W  = mixture.molecular_weights        # kg/kmol
    Y  = mixture.Y                        # mass fractions
    return np.dot(Y, mu / W)              # J/kg mixture

def load_config(config_path):
    with open(config_path, 'r') as f:
        # Use safe_load for security when dealing with untrusted input
        config_data = yaml.safe_load(f)
    return config_data

def chain_reaction_power(config, T_ign, mech_torch, fuel_torch, ox_torch, OF_torch, pc_torch, mech_MCA, MCA_INLET, mdot_MCA, OF_MCA, pc_MCA):

    # Torch
    pc_torch *= PSI2PA
    # MCA
    mdot_MCA *= LBM2KG
    pc_MCA *= PSI2PA
    
    # ### Find Change in enthalpy required to reach Autoigniton Temperature for MCA Propellants
    
    # parse inlet state temperatures [K]
    MCA_fuel,MCA_ox = MCA_INLET.keys()
    T1_f = MCA_INLET[MCA_fuel]
    T1_ox_MCA = MCA_INLET[MCA_ox]


    if config['cp_flag_fuel']: 
        cp_fuel_MCA = config['cp_fuel']
        h1_f_MCA = cp.CoolProp.PropsSI("H", "T", T1_f, "P", pc_MCA, cp_fuel_MCA)
        h2_f_MCA = cp.CoolProp.PropsSI("H", "T", T_ign, "P", pc_MCA, cp_fuel_MCA)
    else:
        h1_f_MCA = config['h1_f']
        h2_f_MCA = config['h2_f']


    if config['cp_flag_ox']: 
        cp_ox_MCA = config['cp_ox']
        h1_ox_MCA = cp.CoolProp.PropsSI("H", "T", T1_ox_MCA, "P", pc_MCA, cp_ox_MCA)
        h2_ox_MCA = cp.CoolProp.PropsSI("H", "T", T_ign, "P", pc_MCA, cp_ox_MCA)
    else:
        h1_ox_MCA = config['h1_ox']
        h2_ox_MCA = config['h2_ox']
       

    dh_f = h2_f_MCA - h1_f_MCA
    dh_ox_MCA = h2_ox_MCA - h1_ox_MCA

    df = pd.DataFrame([dh_f*1e-6, dh_ox_MCA*1e-6], [MCA_fuel, MCA_ox], columns=["Delta h [MJ/kg]"])

    dh_mix = (1/(OF_MCA+1)) * dh_f + (OF_MCA/(OF_MCA+1)) * dh_ox_MCA # J/kg
    print("MCA Propellants: \n State 1: Injection \n State 2: Autoignition Temperature\n")
    print(df)
    print(f"\nMass Averaged dh: {dh_mix * 10**-6:0.3f} [MJ/kg]")


    # find difference in chemical potentials in torch reaction to find energy released from reaction
    torch_rxn = ct.Solution(mech_torch)
    torch_rxn.TPY = T_amb, pc_torch, {ox_torch: OF_torch/(OF_torch+1), fuel_torch:(1/(OF_torch+1)) } 

    # get starting chemical potentials, equillibriate reaction and find difference
    pre_rxn = mass_avg_chem_pot(torch_rxn)
    torch_rxn.equilibrate("HP")
    post_rxn = mass_avg_chem_pot(torch_rxn)

    combustion_energy_torch = np.abs(post_rxn - pre_rxn)


    print("Torch Energy Release:\n")
    print(f"Flame Temp: {torch_rxn.T:0.1f} [K]")
    print(f"Combustion Energy: {combustion_energy_torch* 10**-6:0.3f} [MJ/kg]")


    # same process for MCA Combustion energy
    # find difference in chemical potentials in torch reaction to find energy released from reaction
    MCA_rxn = ct.Solution(mech_MCA)
    MCA_rxn.TPY = T_ign, pc_MCA, {MCA_ox: OF_MCA/(OF_MCA+1), MCA_fuel:(1/(OF_MCA+1))} 

    # get starting chemical potentials, equillibriate reaction and find difference
    pre_rxn = mass_avg_chem_pot(MCA_rxn)
    MCA_rxn.equilibrate("HP")
    post_rxn = mass_avg_chem_pot(MCA_rxn)

    combustion_energy_MCA = np.abs(post_rxn - pre_rxn)

    print("MCA Energy Release:\n")
    print(f"Flame Temp: {MCA_rxn.T:0.1f} [K]")
    print(f"Combustion Energy: {combustion_energy_MCA* 10**-6:0.3f} [MJ/kg]")


    # Compute percent of MCA that needs to ignite for chain reaction
    x = dh_mix/combustion_energy_MCA 
    P_ig = mdot_MCA * x * dh_mix

    print(f"{x*100:0.3f}% of mdot_MCA to cause chain reaction")
    print(f"Torch Power Required: {P_ig*10**-3:0.3f} [kW]")


    nu = 0.2 # Heat Transfer Efficiency from torch exhaust to MCA Propellants
    mdot_torch = P_ig / (combustion_energy_torch * nu)
    mdot_torch_perfect = P_ig / (combustion_energy_torch)


    print(f"Assuming Perfect Heat Transfer:")
    print(f" mdot Torch: {mdot_torch_perfect:0.5f} [kg/s]\n")

    print(f"Assuming Heat Transfer Efficiency of {nu*100:0.2f}%")
    print(f" mdot Torch: {mdot_torch:0.5f} [kg/s]")
    print(f" mdot Torch: {mdot_torch/LBM2KG:0.5f} [lbm/s]")

    return mdot_torch, P_ig, combustion_energy_torch, torch_rxn.T, x

def plot_OF_sweep(OF_set, E_ign, T_comb, mdot, file_path):
    fig, axs = plt.subplots(3, 1, figsize=(10, 10), sharex=True)

    axs[0].plot(OF_set, E_ign)
    axs[0].set_ylabel("Energy [MJ/kg]")
    axs[0].set_title("Energy vs OF")
    axs[0].grid(True)

    axs[1].plot(OF_set, T_comb)
    axs[1].set_ylabel("Temp [K]")
    axs[1].set_title("Temp vs OF")
    axs[1].grid(True)

    axs[2].plot(OF_set, mdot)
    axs[2].set_ylabel("mdot [kg/s]")
    axs[2].set_title("mdot vs OF")
    axs[2].set_xlabel("OF")
    axs[2].grid(True)

    fig_path = os.path.join(file_path, "vs_OF.png")
    plt.savefig(fig_path)
    plt.close()

# ### Tabled for now: Injecting pilot flame into a wsr
# #### -> can be either (likely both) torch combustion products -> torch ignition or torch combustion products -> MCA 
# #### -> may be able to incorporate flame kernel radius into this

# 


