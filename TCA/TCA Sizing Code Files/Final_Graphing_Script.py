#*****note**** ALL OUTPUTS USE NO EFFICIENCY FACTOR AND ARE DIRECTLY FROM NASA CEA

##################################
#Importing Necessary Libraries
##################################
import numpy as np
import matplotlib.pyplot as plt
from rocketcea.cea_obj import CEA_Obj

import yaml

np.set_printoptions(legacy='1.25')   #Fix for some numpy float printing isssues that happened

##################################
#Define Global Conversion Factors
##################################

#Rankine to Kelvin conversion factor
rankineToKelvin = 5.0 / 9.0

##################################
#Define NASA Chemical Equilibrium with Applications (CEA) Object
##################################

#Define CEA Object as C
C = CEA_Obj( oxName='LOX', fuelName='RP1')

##################################
#Define Function Inputs
##################################

#ARRAY BOUNDS
#Lower bound for O/F ratio graphs
mrLo = 0.75
#Upper bound for O/F ratio graphs
mrHi = 4

#Lower bound for chamber pressure graphs 
pLo = 45   #Lower bound of pc divided by 10: (pc / 10)
#Upper bound for chamber pressure graphs
pHi = 55   #Lower bound of pc divided by 10 plus 1:  (pc / 10) + 1

#VARIABLES

#Importing yaml file containing TCA parameters
with open(r'/Users/dl/Documents/GitHub/Turbopump/TCA/TCA_params.yaml') as file:
	tca_params = yaml.safe_load(file)

#Expansion Ratio of Nozzle
EPS = tca_params['tca_expansion_ratio']
#Exit Diameter of Nozzle
De = tca_params['tca_exit_diameter']
#Exit Area of Nozzle
Ae = np.pi * ((De / 2.0) ** 2)
#Target chamber pressure (psi)
pc = tca_params['tca_chamber_pressure']
#Ambient pressure (psi)
pamb= tca_params['ambient_pressure']
#Target exit pressure (psi), equal to pamb
pe = tca_params['tca_exit_pressure']
#O/F ratio
of = tca_params['oxidizer_fuel_ratio']
#Efficieny factor of engine (estimate)
ef_cstar = tca_params['c_star_efficiency']
#Efficieny factor of engine (estimate)
ef_cf = tca_params['thrust_coefficient_efficiency']
#Gravitational Acceleration (ft/s)
g = 32.174

psi_to_pa = 6894.76

N_to_lbf = 0.224809

#ARRAY DECLARATION
mr_range = np.arange(mrLo, mrHi, .05)  #range of O/F ratios, from mrLo to mrHi with spacing 0.05 in between to get proper curves
Comb_Ts = [0] * len(mr_range)          #Empty array for combustion temperatures
Isp_vals = [0] * len(mr_range)         #Empty array for specific impulse values
Cstars = [0] * len(mr_range)           #Empty array for characteristic velocity values
Cf_vals = [0] * len(mr_range)          #Empty array for coefficient of thrust values
T_vals = [0] * len(mr_range)          #Empty array for thrust values

##################################
#All Plotting Figures
##################################

#This figure plots specific impulse values (in seconds) across a set range of O/F ratios at the chosen chamber pressure
plt.figure(1)
count = 0
for mr in mr_range:
    Eps = C.get_eps_at_PcOvPe(Pc = pc, MR = mr, PcOvPe= (pc / pe))  #Optimal Nozzle Expansion Ratio for inputs
    Isp_vals[count] = C.estimate_Ambient_Isp(Pc = pc, MR = mr, eps = Eps, Pamb = pamb)
    count = count + 1
plt.plot(mr_range, Isp_vals, 'b')
Isp_chosen = C.estimate_Ambient_Isp(Pc = pc, MR = of, eps = EPS, Pamb = pamb)
plt.plot(of, Isp_chosen, 'ko')
plt.text(2.15, 244, f'({of : .1f}, {Isp_chosen : .1f}[s])', fontsize = 12)
plt.xlabel("O/F ratio by mass", fontsize = 14)
plt.ylabel("Specific Impulse [s]", fontsize = 14)
plt.title(f"Spec. Impulse vs O/F ratio @{pc} psia", fontsize = 17)
plt.tick_params(axis='x', labelsize=12)
plt.tick_params(axis='y', labelsize=12)
plt.grid()

plt.savefig(r"C:\Users\dl\Documents\GitHub\Turbopump\TCA\CEA Graphs\Specific_Impulse_Target_psi.png")

#This figure plots specific impulse values (in seconds) across a set range of O/F ratios and a range of chamber pressures
plt.figure(2)
for p in range(pLo, pHi):
    count = 0
    for mr in mr_range:
      po = p * 10.0     #Multiplies p values by 10 to reach actual chamber pressure to test
      Eps = C.get_eps_at_PcOvPe(Pc = po, MR = mr, PcOvPe= (po / pe))  #Optimal Nozzle Expansion Ratio for inputs
      Isp_vals[count] = C.estimate_Ambient_Isp(Pc = po, MR = mr, eps = Eps, Pamb = pamb)
      count = count + 1
    plt.plot(mr_range, Isp_vals, label = f'Chamber pressure = {po} psia')
plt.xlabel("O/F ratio by mass")
plt.ylabel("Specific Impulse [s]")
plt.title("Specific Impulse vs O/F ratio")
plt.grid()
plt.legend()

plt.savefig(r"C:\Users\dl\Documents\GitHub\Turbopump\TCA\CEA Graphs\Specific_Impulse_P_Range.png")

#This figure plots combustion temperatures (in Fahrenhiet) across a set range of O/F ratios at the chosen chamber pressure
plt.figure(3)
count = 0
for mr in mr_range:
    Comb_Ts[count] = C.get_Tcomb(Pc = pc, MR = mr) * rankineToKelvin   #Outputs combustion temperature in Rankine, converts to Fahrenheit
    count = count + 1
plt.plot(mr_range, Comb_Ts, 'r')
T_chosen = C.get_Tcomb(Pc = pc, MR = of) * rankineToKelvin
plt.plot(of, T_chosen, 'ko')
plt.text(2.15, 3250, f'({of : .1f}, {T_chosen : .1f}[$^\circ K$])', fontsize = 12)
plt.xlabel("O/F ratio by mass", fontsize = 14)
plt.ylabel("Combustion Temperature [K]", fontsize = 14)
plt.title(f"Combustion Temperature vs O/F ratio @{pc} psia", fontsize = 17)
plt.tick_params(axis='x', labelsize=12)
plt.tick_params(axis='y', labelsize=12)
plt.grid()

plt.savefig(r"C:\Users\dl\Documents\GitHub\Turbopump\TCA\CEA Graphs\Combustion_Temps_Target_psi.png")

#This figure plots combustion temperatures (in Fahrenhiet) across a set range of O/F ratios and a range of chamber pressures
plt.figure(4)
count = 0
for p in range(pLo, pHi):
    count = 0
    for mr in mr_range:
        po = p * 10.0      #Multiplies p values by 10 to reach actual chamber pressure to test
        Comb_Ts[count] = C.get_Tcomb(Pc = po, MR = mr) * rankineToKelvin   #Outputs combustion temperature in Rankine, converts to Fahrenheit
        count = count + 1
    plt.plot(mr_range, Comb_Ts, label = f'Chamber pressure = {po} psia')
plt.xlabel("O/F ratio by mass")
plt.ylabel("Combustion Temperature [K]")
plt.title("Combustion Temps vs O/F ratio")
plt.grid()
plt.legend()

plt.savefig(r"C:\Users\dl\Documents\GitHub\Turbopump\TCA\CEA Graphs\Combustion_Temps_P_Range.png")

#This figure plots characteristic velocity values (in ft/s) across a set range of O/F ratios at the chosen chamber pressure
plt.figure(5)
count = 0
for mr in mr_range:
    Cstars[count] = C.get_Cstar(Pc = pc, MR = mr) * ef_cstar
    count = count + 1
plt.plot(mr_range, Cstars)
plt.xlabel("O/F ratio by mass")
plt.ylabel("Characteristic Velocity [ft/s]")
plt.title(f"Characteristic Velocity vs O/F ratio @{pc} psia")
plt.grid()

plt.savefig(r"C:\Users\dl\Documents\GitHub\Turbopump\TCA\CEA Graphs\Characteristic_Velocity_Target_psi.png")

#This figure plots characteristic velocity values (in ft/s) across a set range of O/F ratios and a range of chamber pressures
plt.figure(6)
for p in range(pLo, pHi):
    count = 0
    for mr in mr_range:
        po = p * 10.0     #Multiplies p values by 10 to reach actual chamber pressure to test
        Cstars[count] = C.get_Cstar(Pc = po, MR = mr) * ef_cstar
        count = count + 1
    plt.plot(mr_range, Cstars, label = f'Chamber pressure = {po} psia')
plt.xlabel("O/F ratio by mass")
plt.ylabel("Characteristic Velocity [ft/s]")
plt.title("Characteristic Velocity vs O/F ratio")
plt.grid()
plt.legend()

plt.savefig(r"C:\Users\dl\Documents\GitHub\Turbopump\TCA\CEA Graphs\Characteristic_Velocity_P_Range.png")

#This figure plots coefficient of thrust values across a set range of O/F ratios at the chosen chamber pressure
plt.figure(7)
count = 0
for mr in mr_range:
    Eps = C.get_eps_at_PcOvPe(Pc = po, MR = mr, PcOvPe= (po / pe))   #Optimal Nozzle Expansion Ratio for inputs
    cfo = C.get_PambCf(Pamb = pamb, Pc = pc, MR = mr, eps = Eps)     #Calculates the coefficient of thrust at given set points
    Cf_vals[count] = cfo[0] * ef_cf   
    count = count + 1
plt.plot(mr_range, Cf_vals)
plt.xlabel("O/F ratio by mass")
plt.ylabel("Coefficient of Thrust")
plt.title(f"Coefficient of Thrust vs O/F ratio @{pc} psia")
plt.grid()

plt.savefig(r"C:\Users\dl\Documents\GitHub\Turbopump\TCA\CEA Graphs\Thrust_Coefficient_Target_psi.png")

#This figure plots coefficient of thrust values across a set range of O/F ratios and a range of chamber pressures
plt.figure(8)
for p in range(pLo, pHi):
    count = 0
    for mr in mr_range:
      po = p * 10.0    #Multiplies p values by 10 to reach actual chamber pressure to test
      Eps = C.get_eps_at_PcOvPe(Pc = po, MR = mr, PcOvPe= (po / pe))   #Optimal Nozzle Expansion Ratio for inputs
      cfo = C.get_PambCf(Pamb = pamb, Pc = po, MR = mr, eps = Eps)     #Calculates the coefficient of thrust at given set points
      Cf_vals[count] = cfo[0] * ef_cf
      count = count + 1
    plt.plot(mr_range, Cf_vals, label = f'Chamber pressure = {po} psia')
plt.xlabel("O/F ratio by mass")
plt.ylabel("Coefficient of Thrust")
plt.title("Coefficient of Thrust vs O/F ratio")
plt.grid()
plt.legend()

plt.savefig(r"C:\Users\dl\Documents\GitHub\Turbopump\TCA\CEA Graphs\Thrust_Coefficient_P_Range.png")

plt.figure(9)
mdot = 9.02
count = 0
for mr in mr_range:
    Ve = C.get_MachNumber(Pc = pc, MR = mr, eps = EPS) * 343   #Exit velocity in m/s
    pcpe = C.get_PcOvPe(Pc = pc, MR = mr, eps = EPS)
    pet = (1.0 / pcpe) * pc * psi_to_pa        #exit pressure for case
    T_vals[count] = (mdot * Ve) - ((pet - (pamb * psi_to_pa)) * Ae) * N_to_lbf
    count = count + 1
plt.plot(mr_range, T_vals)
plt.xlabel("O/F ratio by mass")
plt.ylabel("Thrust Force")
plt.title(f"Coefficient of Thrust vs O/F ratio @{pc} psia")
plt.grid()

plt.savefig(r"C:\Users\dl\Documents\GitHub\Turbopump\TCA\CEA Graphs\Thrust_Generation_Range.png")

plt.show()




