import numpy as np
import math
import matplotlib.pyplot as plt
import cea as cea
from pyfluids import Fluid, FluidsList, Input
from rocketcea.cea_obj import CEA_Obj
import pandas as pd
import yaml

cea = CEA_Obj(fuelName='CH4', oxName='GOX')

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

psi2Pa = 6894.76     # [psi/Pa]
lbm2kg = 0.453592    # [lbm/kg]
m2in = 39.3701       # [m/in]
bar2pa = 100000.0    # [Pa/bar]

r2f = -458.67

FEED_TEMP = 293.15   # [K] (20 deg C)
CD = 0.8
CSTAR_EFF = 0.90

K_F = 1.30           # CH4 specific heat ratio
K_OX = 1.40          # GOX specific heat ratio
R_OXJ = 259.84       # [J / kg-K]
R_FJ = 518.28        # [J / kg-K]

OX_DIA = 0.063       # [in]
FUEL_DIA = 0.063     # [in]

D_TORCH = 0.268      # [in]
P_OX = 570
P_F = 550

A_torch = np.pi * ((D_TORCH / (m2in * 2)))**2
A_ox = np.pi * ((OX_DIA / (m2in * 2)))**2
A_f = np.pi * ((FUEL_DIA / (m2in * 2)))**2

def choked_pressure(mdot, area, k, T, R, Cd=1.0):
    """Calculates upstream stagnation pressure required to choke mdot through an orifice."""
    factor = (2.0 / (k + 1.0)) ** ((k + 1.0) / (2.0 * (k - 1.0)))
    return mdot / (Cd * area * np.sqrt((k / (R * T)) * factor))

def choked_mdot(pressure, area, k, T, R, Cd=1.0):
    """Calculates choked mass flow given upstream stagnation pressure."""
    factor = (2.0 / (k + 1.0)) ** ((k + 1.0) / (2.0 * (k - 1.0)))
    return pressure * Cd * area * np.sqrt((k / (R * T)) * factor)

#fuel_PYF = Fluid(FluidsList.Methane).with_state(Input.pressure(p_line), Input.temperature(FEED_TEMP))
#ox_PYF = Fluid(FluidsList.Oxygen).with_state(Input.pressure(p_line), Input.temperature(FEED_TEMP))

#rho_f = fuel_PYF.density            # [kg/m^3]
#rho_ox = ox_PYF.density             # [kg/m^3]
#sonic_f = fuel_PYF.sound_speed      # [m/s]
#sonic_ox = ox_PYF.sound_speed       # [m/s]

mdot_set = 0.0783 * lbm2kg     #lbm/s to kg/s
OF_set = 1.5
md_ox_set = OF_set * mdot_set / (OF_set + 1)
md_f_set = mdot_set - md_ox_set

#p_ox_set = md_ox_set * np.sqrt(FEED_TEMP * R_OXJ / K_OX) / (CD * A_ox * ((K_OX + 1) / 2) ** (-(K_OX + 1) / (2 * (K_OX - 1))))
#p_f_set = md_f_set * np.sqrt(FEED_TEMP * R_FJ / K_F) / (CD * A_f * ((K_F + 1) / 2) ** (-(K_F + 1) / (2 * (K_F - 1))))

#ox_PYF = Fluid(FluidsList.Oxygen).with_state(Input.pressure(p_ox_set), Input.temperature(FEED_TEMP))

#rho_ox = ox_PYF.density             # [kg/m^3]
#p_ox_set = (md_ox_set)**2 / ((CD * A_ox)**2 * (K_OX * rho_ox * ( 2 / (K_OX + 1))**((K_OX + 1)/(K_OX - 1))))

#print(f"The setpoint feed pressures are {(p_ox_set / psi2Pa) : .3f} for Ox and {(p_f_set / psi2Pa) : .3f} for fuel")

#P_OX = p_ox_set / psi2Pa
#P_F= p_f_set / psi2Pa

pox_arr = [P_OX, (P_OX * 1.15), P_OX, (P_OX * 0.85), P_OX, (P_OX * 1.10), (P_OX * 0.90)]
pf_arr = [P_F, P_F, (P_F * 0.85), P_F, (P_F * 1.15), (P_F * 0.90), (P_F * 1.10)]

pox_arr = np.array(pox_arr) * psi2Pa
pf_arr = np.array(pf_arr) * psi2Pa

mdot_oxs = pox_arr / (np.sqrt(FEED_TEMP * R_OXJ / K_OX) / (CD * A_ox * ((K_OX + 1) / 2) ** (-(K_OX + 1) / (2 * (K_OX - 1)))))
mdot_fs = pf_arr / (np.sqrt(FEED_TEMP * R_FJ / K_F) / (CD * A_f * ((K_F + 1) / 2) ** (-(K_F + 1) / (2 * (K_F - 1)))))

ofs = mdot_oxs / mdot_fs

mdot_totals = mdot_oxs + mdot_fs
mdot_totals = np.array(mdot_totals)

c_ratio_ox = (2 / (K_OX+1)) ** (K_OX / (K_OX - 1)) 
c_ratio_f = (2 / (K_F+1)) ** (K_F / (K_F - 1)) 

###############################
# CEA C-STAR & PC CALCULATION
###############################

# Reactants and solver setup
reac_names = ["CH4", "O2"]
#T_reactant = np.array([FEED_TEMP, FEED_TEMP])  # Reactant temperatures (K)
#fuel_weights = np.array([1.0, 0.0])
#ox_weights = np.array([0.0, 1.0])

#reac = cea.Mixture(reac_names)
#prod = cea.Mixture(reac_names, products_from_reactants=True)
#solver = cea.RocketSolver(prod, reactants=reac)
#solution = cea.RocketSolution(solver)

A_throat = A_torch
Pc_set_Pa = 200 * psi2Pa
cstars = []
T_combs = []
T_exits = []
P_c_arr = []
for m_dot, of in zip(mdot_totals, ofs):
    #weights = reac.of_ratio_to_weights(ox_weights, fuel_weights, of_ratio=of)
    #hc = reac.calc_property(cea.ENTHALPY, weights, T_reactant) / cea.R

    # Evaluate c* at the nominal 200 psia pressure baseline
    #solver.solve(solution, weights, Pc_set_Pa / bar2pa, hc=hc, iac=True)
    #cstar_val = solution.c_star[-1] * CSTAR_EFF

    cstar = cea.get_Cstar(Pc=Pc_set_Pa/psi2Pa, MR=of) * CSTAR_EFF * ft2m # [m/s] (returns ft/s)

    # Back-calculate Pc from total mass flow, c*, and fixed sized throat
    pc_val_Pa = (m_dot * cstar) / A_throat
    cstar = cea.get_Cstar(Pc=pc_val_Pa/psi2Pa, MR=of) * CSTAR_EFF * ft2m # [m/s] (returns ft/s)
    #Tc = cea.get_Tcomb(Pc=Pc_set_Pa/psi2Pa, MR=OF_set) + r2f
    Ts = cea.get_Temperatures(Pc=pc_val_Pa/psi2Pa, MR=of, eps = 1)
    Tc = Ts[0]  + r2f
    Te = Ts[-1]  + r2f

    cstars.append(cstar)
    T_combs.append(Tc)
    T_exits.append(Te)
    P_c_arr.append(pc_val_Pa)

P_c_arr = np.array(P_c_arr)
cstars = np.array(cstars)

cstars = np.array(cstars)
T_combs = np.array(T_combs)
T_exit = np.array(T_exits)

P_c_arr = np.array(P_c_arr)
#P_c_arr = (mdot_oxs + mdot_fs) * cstars / A_throat

###############################
# EXPORT RESULTS TO CSV
###############################

df_results = pd.DataFrame({
    'Pinj_ox': pox_arr / psi2Pa,   # Or keep in Pa: pox_arr
    'Pinj_f': pf_arr / psi2Pa,     # Or keep in Pa: pf_arr
    'Pc': P_c_arr / psi2Pa,        # Or keep in Pa: P_c_arr
    'Cstar': cstars,               # m/s
    'O/F': ofs,
    'Tc': T_combs,                 # K
    'Te': T_exit,                   # K
    'mdot': mdot_totals
})

#output_csv_filename = r"Outputs/torch_test_matrix_results.csv"
#df_results.to_csv(output_csv_filename, index=False)
#print(f"Results successfully exported to {output_csv_filename}")
print(df_results)

###############################
#DEFINE SUBFUNCTIONS
###############################






###############################
#MAIN FUNCTION
###############################
'''
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
'''