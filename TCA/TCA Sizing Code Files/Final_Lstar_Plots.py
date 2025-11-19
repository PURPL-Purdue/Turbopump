##################################
#Importing Necessary Libraries
##################################
import numpy as np
import matplotlib.pyplot as plt
from rocketcea.cea_obj import CEA_Obj

np.set_printoptions(legacy='1.25')    #Fix for some numpy float printing isssues that happened

##################################
#Define Global Conversion Factors
##################################

#psia to pascals conversion
psi_to_pa = 6894.76
#feet to meters conversion
ft_to_m = 0.3048
#m to cm
mcm = 100
#lbf to N conversion, N = kg-m/s^2
lbf_N = 4.4482216
#Efficieny factor of engine (guessed)
ef = 0.975

#Define CEA Object as C
C = CEA_Obj( oxName='LOX', fuelName='RP1')

##################################
#Define Function Inputs
##################################

#ARRAY BOUNDS
#Lower bound for characteristic length (inches)
L_lo = 40
#Upper bound for characteristic length (inches)
L_hi = 50 

#Lower bound for contraction ratio
cr_lo = 1.5
#Upper bound for contraction ratio
cr_hi = 5.0

#VARIABLES
#Force of thurst (lbf)
F = 5000
#Chamber pressure (psi)
pc = 500
#O/F ratio
mr = 2.1
#Exit pressure (psi), equal to ambient pressure
pe = 14.7
#Convergent half-angle (degrees)
a = 45
#Final Chosen L* (inches)
L_star_in = 40.0
#Final Chosen L* (centimeters conversion)
L_star_cm = 40 / 12.0 * ft_to_m * mcm

#ARRAY DECLARATION
cs = np.linspace(cr_lo, cr_hi, 100)    #Array of contraction ratios, with 100 values between the bounds
l_stars = np.linspace(L_lo,L_hi, 11)   #Array of L* values, every integer between 40 and 50, inclusive
Lc_ins = [0] * len(cs)                 #Empty array to store L* outputs for both graphs

##################################
#Calculations
##################################

Ft = F * lbf_N                                                    #Converts thrust to Newtons
Eps = C.get_eps_at_PcOvPe(Pc = pc, MR = mr, PcOvPe= (pc / pe))    #Calculates optimal expansion ratio
cf_arr = C.get_PambCf(Pamb = 14.7, Pc = pc, MR = mr, eps = Eps)   
cf = cf_arr[0] * ef                                               #Calculates the coefficient of thrust with efficiency factor

At = Ft / (cf * pc * psi_to_pa)  #Calculates area of throat (m^2)
At_cm = At * (mcm ** 2)          #Converts throat area to cm^2

##################################
#Graphing Chamber Lengths and Contraction Ratios
##################################

#Plots Chamber Length (inches) vs Contraction Ratio for given range derived from Huzel and Huang chamber volume 
# function and L* = Vc/At relation, plotting range of L* values
plt.figure(1)
for L in l_stars:
    idx = 0
    for cr in cs:
       L1 = L / 12 * ft_to_m * mcm #converts each L* to inches
       #calculates chambe length in inches
       Lc_ins[idx] = (L1 - (1.0/3.0) * np.sqrt(At_cm / np.pi) * (1 / np.tan(np.deg2rad(a))) * (cr **(1.0/3.0) - 1)) / cr / 100 / ft_to_m * 12
       idx = idx + 1
    plt.plot(cs, Lc_ins, label = f'L* = {L} in')
plt.xlabel('Contraction Ratio (Ac/At)')
plt.ylabel('Chamber Length (in)')
plt.title('Chamber Length vs Contraction Ratio')
plt.grid()
plt.legend()

#Plots Chamber Length (inches) vs Contraction Ratio for given range derived from Huzel and Huang chamber volume 
#function and L* = Vc/At relation @chosen L*
plt.figure(2)
idx = 0
for cr in cs:
   #calculates chambe length in inches
   Lc_ins[idx] = (L_star_cm - (1.0/3.0) * np.sqrt(At_cm / np.pi) * (1 / np.tan(np.deg2rad(a))) * (cr **(1.0/3.0) - 1)) / cr / 100 / ft_to_m * 12
   idx = idx + 1
plt.plot(cs, Lc_ins)
plt.xlabel('Contraction Ratio (Ac/At)')
plt.ylabel('Chamber Length (in)')
plt.title(f'Chamber Length vs Contraction Ratio @L* = {L_star_in}in')
plt.grid()
plt.show()

