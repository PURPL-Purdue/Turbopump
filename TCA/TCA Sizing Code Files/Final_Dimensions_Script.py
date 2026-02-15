#*****note**** Change Output Path for CVS (line ) to your own computer and specific folder
#EXAMPLE: output_file_path = r'C:\Users\igoto\OneDrive - purdue.edu\Turbopump\Final_TCA_Scripts\cea_dimensions.csv'
#Make sure the path name ends with \cea_dimensions.csv

##################################
#Importing Necessary Libraries
##################################
import numpy as np
import csv
import yaml
from rocketcea.cea_obj import CEA_Obj

np.set_printoptions(legacy='1.25')    #Fix for some numpy float printing isssues that happened

##################################
#Define NASA Chemical Equilibrium with Applications (CEA) Object
##################################

#Define CEA Object as C
C = CEA_Obj( oxName='LOX', fuelName='RP1')

##################################
#Define Global Conversion Factors
##################################

#psia to pascals conversion
psi_to_pa = 6894.76
#feet to meters conversion
ft_to_m = 0.3048
#m to cm
mcm = 100
#Rankine to Kelvin conversion factor
rankineToKelvin = 5.0 / 9.0
#Rankine to Fahrenheit conversion
rankineToF = -459.67
#lbf to N conversion, N = kg-m/s^2
lbf_N = 4.4482216
#universal gas constant (J / kmol-K)
Ru_m = 8314.462618
#m to in
m_in = 39.3701
#gravity in m/s^2
g = 9.81
#kg to lbm conversion
kg_to_lbm = 2.20462

##################################
#Define Function Inputs
##################################

##################################
#REPLACE PATH WITH PATH YOU NEED FOR YOUR OWN COMPUTER
#Devin Path: 'C:\Users\igoto\Downloads\GH\Turbopump\TCA\TCA_params.yaml'
#Dani Path: '/Users/dl/Documents/GitHub/Turbopump/TCA/TCA_params.yaml'
#Other Path: 
#################################

#Importing yaml file containing TCA parameters
with open(r'C:\Users\igoto\Downloads\GH\Turbopump\TCA\TCA_params.yaml') as file:
	tca_params = yaml.safe_load(file)

#VARIABLES
#Target thrust (lbf)
F = tca_params['thrust']
#Target chamber pressure (psi)
pc = tca_params['tca_chamber_pressure']
#Efficieny factor of engine (estimate)
ef_cstar = tca_params['c_star_efficiency']
#Efficieny factor of engine (estimate)
ef_cf = tca_params['thrust_coefficient_efficiency']
#Ambient pressure (psi)
pamb= tca_params['ambient_pressure']
#Target exit pressure (psi), equal to pamb
pe = tca_params['tca_exit_pressure']
#O/F ratio
mr = tca_params['oxidizer_fuel_ratio']
#Characteristic Length (inches)
L_star_in = tca_params['characteristic_length']
#Characteristic Length calc converted to centimeters
L_star_cm = L_star_in / 12 * ft_to_m * mcm
#Convergent Half-Angle (degrees)
a = tca_params['tca_convergent_half_angle']
#System Mass Flow Rate (lbm/s)
mdot = tca_params['turbopump_mdot']
#Chosen chamber diameter to match pipe (inches)
Dc_in = tca_params['tca_chamber_diameter']

##################################
#Calculations
##################################

Ft = F * lbf_N                                                      #Converts force of thrust to Newtons
Eps = C.get_eps_at_PcOvPe(Pc = pc, MR = mr, PcOvPe= (pc / pe))      #Calculates optimal expansion ratio
Tc_F = C.get_Tcomb(Pc = pc, MR = mr) + rankineToF                   #Calculates combustion temperature in Fahrenheit
Cstar = C.get_Cstar(Pc = pc, MR = mr) * ft_to_m * ef_cstar          #Calculates characteristic velocity in m/s, with efficiency factor
cf_arr = C.get_PambCf(Pamb = pamb, Pc = pc, MR = mr, eps = Eps)     
cf = cf_arr[0] * ef_cf                                                 #Calculates coefficient of thrust, with efficiency factor
isp = C.estimate_Ambient_Isp(Pc = pc, MR = mr, eps = Eps, Pamb = pamb, frozen=0, frozenAtThroat=0)   #calculates isp in seconds

At = ((mdot / kg_to_lbm) * Cstar) / (pc * psi_to_pa)           #Calculates area of throat in m^2
At_in = At * ((m_in) ** 2)             #Converts area of throat to in^2
At_cm = At * (mcm ** 2)                #Converts area of throat to cm^2
Dt_in = 2 * np.sqrt(At_in / np.pi)     #Calculates throat diameter in inches

Ae_in = At_in * Eps                    #Calculates exit area in in^2
De_in = 2 * np.sqrt(Ae_in / np.pi)     #Calculates exit diameter in inches

MW_t, gam_t = C.get_Throat_MolWt_gamma(Pc = pc, MR = mr, eps = Eps, frozen=0)    #Outputs Molecular Weight at the throat
R_comb = Ru_m / MW_t                                                             #Calculates gas constant for our combination
Tc_1, Tt, Te = C.get_Temperatures(Pc = pc, MR = mr, eps = Eps, frozen=0, frozenAtThroat=0) #Outputs throat temps @ chamber, throat, exit
Tt_K = Tt * rankineToKelvin                #Converts the throat temperature to Kelvin
Tt_F = Tt + rankineToF                     #Converts the throat temperature to Fahrenheit

#Calculates the mass flow rate for choked flow at the throat
#choked_mdot = (At * pc * psi_to_pa / (Tt_K**0.5)) * ((gam_t / R_comb)**0.5) * (((gam_t + 1.0)/2.0) ** (-(gam_t + 1.0)/(2.0*(gam_t - 1.0)))) * kg_to_lbm

con_r = (Dc_in / Dt_in) ** 2

Ac_in = At_in * con_r                  #Calculates chamber area in in^2
Dc_in = 2 * np.sqrt(Ac_in / np.pi)     #Calculates chamber diameter in inches

Lc = (L_star_cm - (1.0/3.0) * np.sqrt(At_cm / np.pi) * (1 / np.tan(np.deg2rad(a))) * (con_r **(1.0/3.0) - 1)) / con_r  #Calculates chamber length in cm
Lc_in = Lc / 100 / ft_to_m * 12     #Converts chamber length to inches

#CEA Calcs for extra parameters
k = C.get_HeatCapacities
MW, gamma = C.get_Chamber_MolWt_gamma(Pc = pc, MR = mr, eps = Eps)
gamma, viscosity, k, Pnum = C.get_Chamber_Transport(Pc = pc, MR = 2.1, eps = Eps, frozen=0)

#Printing our outputs
print('\nGiven these inputs:')
print(f'Theoretical specific impulse is {isp[0]: .3f} seconds')
print(f'The temperature of combustion is {Tc_F: .3f} degrees Fahrenheit')
print(f'\nTheoretically maximum mass flow rate (choked at throat) is {mdot: .3f} lmb/s')
print(f'\nThe throat diameter should be {Dt_in: .3f} inches')
print(f'The throat radius should be {(Dt_in / 2): .3f} inches')
print(f'The throat area should be {(At_in): .3f} inches squared')
print(f'\nThe expansion ratio should be {Eps: .3f}')
print(f'\nThe exit diameter should be {De_in: .3f} inches')
print(f'The exit radius should be {(De_in / 2): .3f} inches')
print(f'\nThe chamber diameter should be {Dc_in: .3f} inches')
print(f'The contraction ratio is {con_r: .3f}')
print(f'\nThe combustion chamber length should be {Lc_in : .3f} inches')
print(f'\nThe molecular weight is {MW : .3f} lbm/lbmole')
print(f'The ratio of specific heats is {gamma : .3f}')
print(f'The thermal conductivity is {k : .3f} mcal/cm-K-s')

#Construction of data list to export to csv, includes input values and output values
data = [
['Force of thrust(lbf)', 'Chamber Pressure (psia)', 'Contraction Ratio', 'O/F Ratio', 'L* (in)', 'Convergent Half-Angle (deg)',
'Chamber Diameter (in)', 'Throat Diameter (in)', 'Exit Diameter (in)', 'Chamber Length (in)', 'Expansion Ratio', 'Mass Flow Rate (lbm/s)'],
[F, pc, round(con_r, 1), round(mr, 1), round(L_star_in, 1), round(a, 1), round(Dc_in, 3), round(Dt_in, 3), round(De_in, 3), round(Lc_in, 3),
 round(Eps, 3), round(mdot,3)]
]

##################################
#REPLACE PATH WITH PATH YOU NEED FOR YOUR OWN COMPUTER

#Devin Path: 'C:\Users\igoto\Downloads\GH\Turbopump\TCA\TCA Sizing Code Files\Dimensions CSV/dimensions.csv'
#Dani Path: '/Users/dl/Documents/GitHub/Turbopump/TCA/Dimensions CSV/dimensions.csv'
#Other Path: 
##################################
output_file_path = r'C:\Users\igoto\Downloads\GH\Turbopump\TCA\TCA Sizing Code Files\Dimensions CSV/dimensions.csv'

with open(output_file_path, 'w', newline='') as csvfile:
  csv_writer = csv.writer(csvfile)
  csv_writer.writerows(data)

