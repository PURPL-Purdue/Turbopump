from rocketcea.cea_obj import CEA_Obj
import numpy as np
from pyfluids import Fluid, FluidsList, Input
import matplotlib.pyplot as plt

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
mdot_main = 9.5 / lbm2kg # main engine mDot [lbm/s]
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

# get CEA outputs
ispObj = CEA_Obj( oxName='GOX', fuelName='GH2')
[k, MolWt_t] = ispObj.get_Throat_MolWt_gamma(Pc= p_chamber, MR= OF, eps= 0)
cstar = Cf * ispObj.get_Cstar(Pc= p_chamber, MR= OF) * ft2m # [m/s]
isp = Cf * ispObj.get_Throat_Isp(Pc= p_chamber, MR= OF, frozen= 1) #[s]
rho_chamber = ispObj.get_Chamber_Density(Pc= p_chamber, MR= OF) * 16.0185 # [lb/ft3 -> kg/m3]
T_comb = ispObj.get_Tcomb(Pc= p_chamber, MR= OF) * R2K #[K]

#calculate ambient isp using CF
x = ispObj.get_PambCf(Pamb=14.7, Pc=p_chamber, MR=OF, eps=1)
CF = x[0]
isp = cstar * CF / g0


######### ORIFICE/NOZZLE SIZING #########

# calculate total mass flow rate:
# prior art indicates (main mass flow / 200) provides reliable main engine ignition
# increase mass flow to compensate for losses
mDot_torch = mdot_main * (1/200)  * lbm2kg # [kg/s]

# calculate ox and fuel mass flow rates
mDot_f = mDot_torch / (OF + 1) # [kg/s]
mDot_ox = mDot_torch - mDot_f  # [kg/s]

# calculate injection area that ensures choked flow [m^2]
A_f = mDot_f / (Cf * np.sqrt(k_f * rho_f * p_line_f * (2 / (k_f + 1))**( (k_f + 1) / (k_f - 1)))) #[m^2]
A_ox = mDot_ox / (Cf * np.sqrt(k_ox * rho_ox * p_line_ox * (2 / (k_ox + 1))**( (k_ox + 1) / (k_ox -1 )))) #[m^2]

# ensure choked flow at throat
throat_critRatio = ((k + 1) / 2) ** (k / (k - 1))
A_throat = mDot_torch * cstar / (p_chamber * psi2Pa) #[m^2]

# diameter calcs [in]
d_ox = AtoD(A_ox * msq2insq)
d_f = AtoD(A_f * msq2insq)
d_throat = AtoD(A_throat * msq2insq)

print('\nMass Flows')
print(f' Torch Mass Flow [lbm/s]: {mDot_torch/lbm2kg:.3f}')
print(f' Oxidizer Mass Flow [lbm/s]: {mDot_ox/lbm2kg:0.3f}')
print(f' Fuel Mass Flow [lbm/s]: {mDot_f/lbm2kg:0.3f}')

print('\nMisc')
print(f' Thrust [N] {isp * mDot_torch * g0:0.3f}')
print(f' Combustion Temperature [K]: {T_comb:0.3f}\n')
print(f'Orifice/Nozzle Diameters [in] \n GOX: {d_ox:.8f}\n GH2: {d_f:.8f}\n Throat: {d_throat:.4f}')
print(f'Choked Flow at Throat: {(p_chamber / p_amb) > throat_critRatio}')
######### END ORIFICE/NOZZLE SIZING #########

######### CHAMBER SIZING ########
"""
# chamber volume and stay time 
Lstar = 30 / m2in  # [in->m] 
V_chamber = Lstar * A_throat # [m^3]
t_stay = (V_chamber * rho_chamber) / mDot_torch # [s]
"""
#get chamber volume by defining stay time
t_stay = 0.001 # [s]
V_chamber = t_stay * mDot_torch / rho_chamber # [m^3]
Lstar = V_chamber / A_throat # [m]

#chamber volume -> dimensions
conv_angle = 45 # convergent angle [deg]
r_contraction = 30 # contraction ratio

A1 = r_contraction * A_throat #chamber area [m^2]
d_chamber = AtoD(A1)# chamber diameter [m]
L_conv = (d_chamber/2 - (d_throat / m2in)/2) / np.sin(np.deg2rad(conv_angle))# convergent length[m]
L1 = ( V_chamber - A1*L_conv * (1 + np.sqrt(A_throat/A1) + A_throat/A1) ) / A1 #chamber length [m]

# check chamber volume
#V_chamber_geo = (A1*L1 + A1*L_conv * (1 + np.sqrt(A_throat/A1) + A_throat/A1)) * (m2in**3) #[m^3 -> in^3]
#print(f' Chamber Volume From Geometry [in^3]: {V_chamber_geo:0.5f}\n')

print(f'\nChamber Sizing')
print(f' Stay Time [s]: {t_stay:0.4f}')
print(f' L Star [in]: {Lstar*m2in:0.3f}')
print(f" Chamber Volume [in^3]: {V_chamber * (m2in**3):0.3f}")
print(f' Contraction Ratio: {A1/A_throat:0.2f}')
print(f' Chamber Diameter [in^3]: {d_chamber * m2in:0.3f}')
print(f' Chamber Length [in]: {L1* m2in:0.3f}')
print(f' Convergent Length [in]: {L_conv* m2in:0.3f}')
print(f' Convergent Angle [deg] {conv_angle}\n')

######### END CHAMBER SIZING ########




######### STRESS ANALYSIS #########
# Hoop Stress for chamber wall thickness
# 7075 Aluminum
sigma_yield = 435 * (10**6) # [MPa->Pa] 

FoS = 10
sigma_FoS = sigma_yield / FoS

r_chamber = d_chamber / 2 #[m]
t_chamber =  ((p_chamber*psi2Pa) * r_chamber) / sigma_FoS #[m]

print(f'Structural')
print(f' Chamber Thickness: {t_chamber*m2in:0.5f}')
print(f' Factor of Safety: {FoS:0.3f}')

######### END STRESS ANALYSIS #########

######### ENSURE SONIC FLOW #########

# assume isentropic flow to define state at orifice throat
T_t_f = 2 * (T_AMB_CELSIUS + 273) / (k_f + 1) - 273
T_t_ox = 2 * (T_AMB_CELSIUS + 273) / (k_ox + 1) - 273

rho_t_f = rho_f * (2 / (k_f + 1)) ** (1 / (k_f - 1))
rho_t_ox = rho_ox * (2 / (k_ox + 1)) ** (1 / (k_ox - 1))

fuel_throat = Fluid(FluidsList.Hydrogen).with_state(Input.density(rho_t_f), Input.temperature(T_t_f))
ox_throat = Fluid(FluidsList.Oxygen).with_state(Input.density(rho_t_ox), Input.temperature(T_t_ox))

# find flow velocity and compute Mach
#v_f = mDot_f / (rho_t_f * A_f)
fuel_sonic = fuel_throat.sound_speed
A_f = mDot_f / (rho_t_f * fuel_sonic)
v_f = mDot_f / (rho_t_f * A_f)
fuel_Mach = v_f / fuel_sonic

#v_ox = mDot_ox / (rho_t_ox * A_ox)
ox_sonic = ox_throat.sound_speed
A_ox = mDot_ox / (rho_t_ox * ox_sonic)
v_ox = mDot_ox / (rho_t_ox * A_ox)
ox_Mach = v_ox / ox_sonic

print('\nMomentum Ratios')
print(f' Fuel Velocity: {v_f:0.3f}')
print(f' Fuel Sonic: {fuel_sonic:0.3f}')
print(f' Fuel Mach: {fuel_Mach:0.3f}')

print(f' Ox Velocity: {v_ox:0.3f}')
print(f' Ox Sonic: {ox_sonic:0.3f}')
print(f' Ox Mach: {ox_Mach:0.3f}')

# diameter calcs [in]
d_ox = AtoD(A_ox * msq2insq)
d_f = AtoD(A_f * msq2insq)
d_throat = AtoD(A_throat * msq2insq)
print(f'\nOrifice/Nozzle Diameters [in] \n GOX: {d_ox:.8f}\n GH2: {d_f:.8f}\n Throat: {d_throat:.4f}')

######### END SONIC FLOW #########

######### INJECTOR DIMENSIONS #########

# momentum ratios

LD = 5 # impingement dist / Orifice Diameter

#confirm orifice diamters:
# ideal diameter ratio for max mixing efficiency (MME) 
# (VALIDATED FROM LIQUIDS) | may use liquid densities?
M = 1 # mixing coefficient for doublet
r_MME = np.sqrt((M * (rho_ox / rho_f) * (mDot_f / mDot_ox)**2)**0.7)
print(f'Optimal d_f [in]: {r_MME * d_ox:0.3f}, Actual d_f [in]: {d_f:0.3f}\n')
# orifice diameters are not near equal: going to keep it because it aligns with 
# this formula, additionally, density is lower for GH2

#impingement distance: use mean of ox and fuel diameters
dist_imp = LD * np.mean([d_ox, d_f])
mv_ox = mDot_ox * ox_sonic
mv_fuel = mDot_f * fuel_sonic

theta_f = np.arange(0, 360, 0.1) # [deg]
Dstar = (1 / np.tan(np.arcsin(mv_ox/mv_fuel * np.sin(theta_f)))) + (1 / np.tan(theta_f))
D_0 = np.mean([d_ox, d_f])
k = 5

PE = (Dstar - d_chamber/ (k*D_0)) / (Dstar - d_chamber/(k*D_0)) 
#for i in range(len(Dstar)):
    #print(Dstar[i])
print(d_chamber/ (k*D_0))
print(np.argwhere(PE < 0.1))


#thickness of plate: fluids perspective
tD = 4 # thickness/ orifice diameter: typically 4-10: chose lower bound
t_plate = tD * np.mean([d_ox, d_f])

"""
print('Injector Sizing')
print(f' Impingement Angle [deg]: {theta_imp}')
print(f' Impingement Distance [in]: {dist_imp:0.3f}')
print(f' Orifice Distance [in]: {dist_orifices:0.3f}')
print(f' Plate Thickness [in]: {t_plate:0.3f}')
"""

######### END INJECTOR SIZING #########