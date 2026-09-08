from rocketcea.cea_obj import CEA_Obj
import numpy as np
from pyfluids import Fluid, FluidsList, Input
import matplotlib.pyplot as plt
import CEA_Wrap as CEA

def AtoD(Area):
    return 2 * np.sqrt(Area / np.pi)

def DtoA(Diameter):
    return np.pi * (Diameter/2)**2

# constants
n2lbf = 4.44822
psi2Pa = 6894.76
ft2m = 0.3048
lbm2kg = 0.453592
m2in = 39.3701
msq2insq = 1550
g0 = 9.81 # [m/s^2]
T_AMB_CELSIUS = 20 # [deg C]
k_f = 1.41  # gh2 specific heat ratio
k_ox = 1.40 # gox specific heat ratio
p_amb = 14.7 # ambient pressure
R2K = 0.555556 # rankine to Kelvin

# desgin parameters
mdot_main = 1.884 # main engine mDot [lbm/s]
p_chamber = 150 # chamber pressure [psi]
#thrust = 15 # [N]
OF = 2.9
Cf = 0.8 # efficiency factor

########## LINE PRESSURES ############

# calculate line pressure using stiffness
# chamber is in psi, convert to Pa for pyfluids

# calculate critical pressure ratios (minimum stiffness)
fuel_critRatio = ((k_f + 1) / 2) ** (k_f / (k_f - 1)) 
ox_critRatio = ((k_ox + 1) / 2) ** (k_ox / (k_ox - 1))

stiffness_f = fuel_critRatio + 0.2 # injector stiffness 
stiffness_ox = ox_critRatio + 0.2 # injector stiffness 
p_line_f = (1+stiffness_f) * p_chamber * psi2Pa # [Pa]
p_line_ox = (1+stiffness_ox) * p_chamber * psi2Pa # [Pa]
print(f'Pressures [psi]:\n Fuel Line: {p_line_f/psi2Pa:0.3f}\n Oxidizer Line: {p_line_ox/psi2Pa:0.3f}')
print(f' Chamber: {p_chamber}')

# make pyfluids instances of propellants and save parameters
fuel = Fluid(FluidsList.Hydrogen).with_state(Input.pressure(p_line_f), Input.temperature(T_AMB_CELSIUS))
ox = Fluid(FluidsList.Oxygen).with_state(Input.pressure(p_line_ox), Input.temperature(T_AMB_CELSIUS))

rho_f = fuel.density
rho_ox = ox.density 
########## END LINE PRESSURES ############
Cf = 1
# get CEA outputs
ispObj = CEA_Obj( oxName='GOX', fuelName='GH2')
[k, MolWt_t] = ispObj.get_Throat_MolWt_gamma(Pc= p_chamber, MR= OF, eps= 0)
cstar = Cf * ispObj.get_Cstar(Pc= p_chamber, MR= OF) * ft2m # [m/s]
rho_chamber = ispObj.get_Chamber_Density(Pc= p_chamber, MR= OF) * 16.0185 # [lb/ft3 -> kg/m3]
T_comb = ispObj.get_Tcomb(Pc= p_chamber, MR= OF) * R2K #[K]

# Matches CEA Website outputs, units are different but noted
print(ispObj.get_full_cea_output(Pc= p_chamber, MR= OF, eps=1))

# gives vaccuum Isp
isp = Cf * ispObj.get_Throat_Isp(Pc= p_chamber, MR= OF, frozen= 1) #[s]

print(ispObj.estimate_Ambient_Isp(Pc=p_chamber, MR=OF, eps=1, Pamb=14.7, frozen=1))

#calculate ambient isp using CF
x = ispObj.get_PambCf(Pamb=14.7, Pc=p_chamber, MR=OF, eps=1)
CF = x[0]
print(cstar*CF / g0)