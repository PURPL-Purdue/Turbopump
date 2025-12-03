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
mrLo = 1
#Upper bound for O/F ratio graphs
mrHi = 3

#Lower bound for chamber pressure graphs 
pLo = 45   #Lower bound of pc divided by 10: (pc / 10)
#Upper bound for chamber pressure graphs
pHi = 55   #Lower bound of pc divided by 10 plus 1:  (pc / 10) + 1

#VARIABLES

#Importing yaml file containing TCA parameters
with open(r'C:\Users\igoto\Documents\GitHub\Turbopump\TCA\TCA_params.yaml') as file:
	tca_params = yaml.safe_load(file)

#Target chamber pressure (psi)
pc = tca_params['tca_chamber_pressure']
#Ambient pressure (psi)
pamb= tca_params['ambient_pressure']
#Target exit pressure (psi), equal to pamb
pe = tca_params['tca_exit_pressure']

#ARRAY DECLARATION
mr_range = np.arange(mrLo, mrHi, .05)  #range of O/F ratios, from mrLo to mrHi with spacing 0.05 in between to get proper curves
Comb_Ts = [0] * len(mr_range)          #Empty array for combustion temperatures
Isp_vals = [0] * len(mr_range)         #Empty array for specific impulse values
Cstars = [0] * len(mr_range)           #Empty array for characteristic velocity values
Cf_vals = [0] * len(mr_range)          #Empty array for coefficient of thrust values

##################################
#All Plotting Figures
##################################

#This figure plots specific impulse values (in seconds) across a set range of O/F ratios at the chosen chamber pressure
plt.figure(1)
count = 0
for mr in mr_range:
    Eps = C.get_eps_at_PcOvPe(Pc = pc, MR = mr, PcOvPe= (pc / pe))  #Optimal Nozzle Expansion Ratio for inputs
    Isp_vals[count] = C.get_Isp(Pc = pc, MR = mr, eps = Eps) #Calculates specific impulse (s) for given inputs
    count = count + 1
plt.plot(mr_range, Isp_vals)
plt.xlabel("O/F ratio by mass")
plt.ylabel("Specific Impulse (s)")
plt.title(f"Spec. Impulse vs O/F ratio @{pc} psi")
plt.grid()

plt.savefig(r"TCA\CEA Graphs\Specific_Impulse_Target_psi.png")

#This figure plots specific impulse values (in seconds) across a set range of O/F ratios and a range of chamber pressures
plt.figure(2)
for p in range(pLo, pHi):
    count = 0
    for mr in mr_range:
      po = p * 10.0     #Multiplies p values by 10 to reach actual chamber pressure to test
      Eps = C.get_eps_at_PcOvPe(Pc = po, MR = mr, PcOvPe= (po / pe))  #Optimal Nozzle Expansion Ratio for inputs
      Isp_vals[count] = C.get_Isp(Pc = po, MR = mr, eps = Eps)     #Calculates specific impulse (s) for given inputs
      count = count + 1
    plt.plot(mr_range, Isp_vals, label = f'Chamber pressure = {po} psia')
plt.xlabel("O/F ratio by mass")
plt.ylabel("Specific Impulse (s)")
plt.title("Specific Impulse vs O/F ratio")
plt.grid()
plt.legend()

plt.savefig(r"TCA\CEA Graphs\Specific_Impulse_P_Range.png")

#This figure plots combustion temperatures (in Fahrenhiet) across a set range of O/F ratios at the chosen chamber pressure
plt.figure(3)
count = 0
for mr in mr_range:
    Comb_Ts[count] = C.get_Tcomb(Pc = pc, MR = mr) * rankineToKelvin   #Outputs combustion temperature in Rankine, converts to Fahrenheit
    count = count + 1
plt.plot(mr_range, Comb_Ts)
plt.xlabel("O/F ratio by mass")
plt.ylabel("Combustion Temperature (F)")
plt.title(f"Combustion Temperature vs O/F ratio @{pc} psi")
plt.grid()

plt.savefig(r"TCA\CEA Graphs\Combustion_Temps_Target_psi.png")

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
plt.ylabel("Combustion Temperature (F)")
plt.title("Combustion Temps vs O/F ratio")
plt.grid()
plt.legend()

plt.savefig(r"TCA\CEA Graphs\Combustion_Temps_P_Range.png")

#This figure plots characteristic velocity values (in ft/s) across a set range of O/F ratios at the chosen chamber pressure
plt.figure(5)
count = 0
for mr in mr_range:
    Cstars[count] = C.get_Cstar(Pc = pc, MR = mr)
    count = count + 1
plt.plot(mr_range, Cstars)
plt.xlabel("O/F ratio by mass")
plt.ylabel("Characteristic Velocity (ft/s)")
plt.title(f"Characteristic Velocity vs O/F ratio @{pc} psi")
plt.grid()

plt.savefig(r"TCA\CEA Graphs\Characteristic_Velocity_Target_psi.png")

#This figure plots characteristic velocity values (in ft/s) across a set range of O/F ratios and a range of chamber pressures
plt.figure(6)
for p in range(pLo, pHi):
    count = 0
    for mr in mr_range:
        po = p * 10.0     #Multiplies p values by 10 to reach actual chamber pressure to test
        Cstars[count] = C.get_Cstar(Pc = po, MR = mr)
        count = count + 1
    plt.plot(mr_range, Cstars, label = f'Chamber pressure = {po} psia')
plt.xlabel("O/F ratio by mass")
plt.ylabel("Characteristic Velocity (ft/s)")
plt.title("Characteristic Velocity vs O/F ratio")
plt.grid()
plt.legend()

plt.savefig(r"TCA\CEA Graphs\Characteristic_Velocity_P_Range.png")

#This figure plots coefficient of thrust values across a set range of O/F ratios at the chosen chamber pressure
plt.figure(7)
count = 0
for mr in mr_range:
    Eps = C.get_eps_at_PcOvPe(Pc = po, MR = mr, PcOvPe= (po / pe))   #Optimal Nozzle Expansion Ratio for inputs
    cfo = C.get_PambCf(Pamb = pamb, Pc = pc, MR = mr, eps = Eps)     #Calculates the coefficient of thrust at given set points
    Cf_vals[count] = cfo[0]    
    count = count + 1
plt.plot(mr_range, Cf_vals)
plt.xlabel("O/F ratio by mass")
plt.ylabel("Coefficient of Thrust")
plt.title(f"Coefficient of Thrust vs O/F ratio @{pc}psia")
plt.grid()

plt.savefig(r"TCA\CEA Graphs\Thrust_Coefficient_Target_psi.png")

#This figure plots coefficient of thrust values across a set range of O/F ratios and a range of chamber pressures
plt.figure(8)
for p in range(pLo, pHi):
    count = 0
    for mr in mr_range:
      po = p * 10.0    #Multiplies p values by 10 to reach actual chamber pressure to test
      Eps = C.get_eps_at_PcOvPe(Pc = po, MR = mr, PcOvPe= (po / pe))   #Optimal Nozzle Expansion Ratio for inputs
      cfo = C.get_PambCf(Pamb = pamb, Pc = po, MR = mr, eps = Eps)     #Calculates the coefficient of thrust at given set points
      Cf_vals[count] = cfo[0]
      count = count + 1
    plt.plot(mr_range, Cf_vals, label = f'Chamber pressure = {po} psia')
plt.xlabel("O/F ratio by mass")
plt.ylabel("Coefficient of Thrust")
plt.title("Coefficient of Thrust vs O/F ratio")
plt.grid()
plt.legend()

plt.savefig(r"TCA\CEA Graphs\Thrust_Coefficient_P_Range.png")

plt.show()




