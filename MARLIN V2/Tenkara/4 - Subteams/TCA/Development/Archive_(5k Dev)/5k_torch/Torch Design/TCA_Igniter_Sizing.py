# TCA Igniter Sizing Code
# 11/28/2025 : Created - Louis DeSano

# units standard: every variable is stored as SI, convert when inputting to external libraries if needed
# rocketCEA: uses English Engineering
# pyfluids: uses SI with Celsius

import numpy as np
import matplotlib.pyplot as plt
from rocketcea.cea_obj import CEA_Obj
from pyfluids import Fluid, FluidsList, Input
import pandas as pd

## Unit Conversions & Constants ##
n2lbf = 4.44822     # [N/lbf]

psi2Pa = 6894.76    # [psi/Pa]

lbm2kg = 0.453592   #[lbm/kg]

ft2m = 0.3048       # [ft/m]
m2in = 39.3701      # [m/in]

R2K = 0.555556      # [R/K]

# Constants
G0 = 9.81           # [m/s^2]
P_ATM = 101325      # [Pa]
T_AMB_CELSIUS = 20 # [deg C]

# Function: Calculate Area from Diameter
def diameter_to_area(diameter):
# Input: 
    # Diameter
# Output: 
    # Area
    return np.pi * (diameter/2)**2

# Function: Calculate Diameter from Area
def area_to_diameter(area):
# Input: 
    # Area
# Output: 
    # Diameter
    return 2 * np.sqrt(area / np.pi)

# Function: Calculate back pressure for choked flow
def choked_backpressure(k, stiffness, p_c):
# Inputs
    # k: specific heat ratio [--]
    # stiffness: amount stiffer than choked minimum as a decimal [--]
    # p_c: chamber pressure [Pa]

# Outputs
    # p_line: line pressure sized for choked flow [Pa]

    # calculate critical pressure ratio (minimum stiffness for choked flow)
    critRatio = (2 / (k+1)) ** (k / (k - 1)) 
    
    # add stiffness and calculate line pressure [Pa]
    p_line = (1 + stiffness) * critRatio**-1 * p_c

    return p_line

def unchoked_area(mDot, k, p_line, p_ratio, rho_line, Cd):
# Inputs
    # k: specific heat ratio [--]
    # stiffness: pressure drop / chamber pressure [--]
    # p_c: chamber pressure [Pa]

# Outputs
    # A_inj: orifice injection Area [m^2]

    A_inj = mDot/Cd * ( np.sqrt(2 * rho_line * p_line * (k / (k-1)) * (p_ratio**(2/k)  - p_ratio**((k+1) / k))) )**-1

    return A_inj

# Function: Construct line of excel sheet and add to sheet
def add_cad_param(sheet, name, dim, unit):
#Inputs
    #name: dimension name string
    #dim: dimension value
    #unit: unit string
#Output
    #sheet: Sheet with print to be printed to excel

    excelLine = [name, unit, f'{dim:0.3f}{unit}']
    sheet.append(excelLine)

    return sheet

# Function: Choked Flow at Flow Constriction
def throat_area(mDot, Cd, k, rho , p0):
# Inputs: 
    # k: Specific Heat Ratio    [--]
    # Cd: Discharge Coefficient [--]
    # rho: density              [kg/m^3]
    # p0: stagnation pressure   [Pa]
# Outputs:
    # At: Throat Area           [m^3]

    At = mDot / (Cd * np.sqrt(k * rho * p0 * ( 2 / (k + 1))**((k + 1)/(k - 1))))
    return At

# Function: Size Line for a velocity
def size_line(mDot, rho, v_line_chosen):
    # mDot: mass flow rate [kg/s]
    # rho: density [kg/m^3]
    # A_inj: injector area [m^2]
    # v_line_chosen: target line velocity [m/s]

    # Nominal SAE ORB inner diameters
    ID = np.arange(1/16, 2, 1/16) / m2in   # [m]

    # Line velocity for each ID:
    v_line = (mDot / rho) * (4 / (np.pi * ID**2))

    # Find the index of the diameter that gives closest velocity
    idx = np.abs(v_line - v_line_chosen).argmin()

    # Extract chosen diameter and velocity
    best_ID = ID[idx] * m2in
    best_v = v_line[idx]

    return best_ID, best_v

def main():

    ### Design Setpoints ###

    mDot_main = 9.5                 # main chamber mass flow [kg/s]
    mDot_torch = mDot_main / (2* 100)    # torch mass flow [kg/s] (Huzel and Huang)
    mDot_torch = 0.018
    p_c = 300 * psi2Pa              # torch chamber pressure [psi->Pa]

    OF = 2.2                       # OF ratio

    Cd = 0.9                       # Discharge Coefficient
    choked = 1

    # Propellant Specific
    ox_CEA = 'GOX'
    fuel_CEA = 'GH2'
    k_f = 1.41  # gh2 specific heat ratio
    k_ox = 1.40 # gox specific heat ratio  

    ### End Design Setpoints ###

    """ Mass Flow Rates """
    # calculate ox and fuel mass flow rates
    mDot_f = mDot_torch / (OF + 1) # [kg/s]
    mDot_ox = mDot_torch - mDot_f  # [kg/s]

    print('\nMass Flows')
    print(f' Torch Mass Flow [kg/s]: {mDot_torch:.5f}')
    print(f' Oxidizer Mass Flow [kg/s]: {mDot_ox:0.5f}')
    print(f' Fuel Mass Flow [kg/s]: {mDot_f:0.5f}')
    print('\n')


    """ Orifice/Throat Areas"""
    # Calculate Line Pressures based on chosen stiffness (uniform stiffness)
    stiffness = 0.2                    # deltaP / p_c [--]

    #define pressure ratio and line pressure
    if (choked):
        p_line = choked_backpressure(k_ox, stiffness, p_c)

    else:
        p_ratio = 1 / (1 + stiffness)   # pc/pline [--]
        p_line = p_c / p_ratio          # line pressure [Pa]
    ## Get fluid properties to size orifice/throat areas ##
    #pyfluids for ox and fuel
    fuel_PYF = Fluid(FluidsList.Hydrogen).with_state(Input.pressure(p_line), Input.temperature(T_AMB_CELSIUS))
    ox_PYF = Fluid(FluidsList.Oxygen).with_state(Input.pressure(p_line), Input.temperature(T_AMB_CELSIUS))

    rho_f = fuel_PYF.density            # [kg/m^3]
    rho_ox = ox_PYF.density             # [kg/m^3]
    sonic_f = fuel_PYF.sound_speed      # [m/s]
    sonic_ox = ox_PYF.sound_speed       # [m/s]

    print("\nPressures")
    print(f" Stiffness: {stiffness * 100}%")
    print(f" Ox   [psi]: {p_line/psi2Pa:0.3f}")
    print(f" Fuel [psi]: {p_line/psi2Pa:0.3f}")
    print(f" Chamber [psi]: {p_c/psi2Pa:0.3f}")
    print("\n")

    # CEA for chamber exit (frozen gives lower performance bound)
    ispObj = CEA_Obj(oxName=ox_CEA, fuelName=fuel_CEA)
    isp = 0.8 * ispObj.estimate_Ambient_Isp(Pc= (p_c/psi2Pa), MR= OF, Pamb= (P_ATM/psi2Pa), eps=1, frozen= 1)[0] # [s]
    rho_c = ispObj.get_Densities(Pc=(p_c/psi2Pa), MR= OF, eps=1, frozen=1)[0] * lbm2kg / (ft2m**3)             #[lbm/ft^3 -> kg/m^3]
    k_c = ispObj.get_Chamber_MolWt_gamma(Pc= (p_c/psi2Pa), MR= OF, eps=1)[1]

    thrust = mDot_torch * isp * G0                           # [N]
    cstar = ispObj.get_Cstar(Pc=p_c, MR=OF) * ft2m           # [m/s]
    T_comb = ispObj.get_Tcomb(Pc= p_c, MR= OF) * R2K         # [K]
    ## End of Fluid Properties ##

    # Assume Cd and cstar efficiency are unity to provide margin for choked flow
    A_t = throat_area(mDot_torch, Cd, k_c, rho_c, p_c)                  # Throat Area [m^2]
    #A_t = mDot_torch * cstar / (p_c)                                   # Throat Area [m^2]

    if choked:
        A_f = throat_area(mDot_f, Cd, k_f, rho_f, p_line)
        A_ox = throat_area(mDot_ox, Cd, k_ox, rho_ox, p_line)
    else:
        A_f = unchoked_area(mDot_f, k_f, p_line, p_ratio, rho_f, Cd)        # Fuel Injection Area [m^2]
        A_ox = unchoked_area(mDot_ox, k_ox, p_line, p_ratio, rho_ox, Cd)    # Ox Injection Area   [m^2]

    D_f = area_to_diameter(A_f) 
    D_ox = area_to_diameter(A_ox)
    D_t = area_to_diameter(A_t)
    
    print("\nOrifice/Nozzle Sizing [mm]")
    print(f" Fuel Injection Diameter: {D_f*10**3:0.3f}")
    print(f" Ox Injection Diameter: {D_ox*10**3:0.3f}")
    print(f" Throat Diameter: {D_t*10**3:0.3f}")
    print("\n")


    # Additional Outputs
    print("\nMisc")
    print(f" Thrust [N]: {thrust:0.3f} ")
    print(f" Tcomb  [C]: {T_comb + 273:0.3f} ")
    print("\n")


    """Chamber Dimensions""" 
    #get chamber volume by defining stay time
    t_stay = 0.0005 # [s]
    V_chamber = t_stay * mDot_torch / rho_c # [m^3]
    Lstar = V_chamber / A_t # [m]

    #chamber volume -> dimensions
    conv_angle = 45 # convergent angle [deg]
    r_contraction = 8 # contraction ratio

    A1 = r_contraction * A_t                                                # chamber area [m^2]
    D_c = area_to_diameter(A1)                                              # chamber diameter [m]
    D_c = 21/64 / m2in
    L_conv = (D_c/2 - D_t/2) / np.sin(np.deg2rad(conv_angle))               # convergent length[m]
    L1 = ( V_chamber - A1*L_conv * (1 + np.sqrt(A_t/A1) + A_t/A1) ) / A1    # chamber length [m]

    print(f'\nChamber Sizing')
    print(f' Stay Time [s]: {t_stay:0.4f}')
    print(f' L Star [m]: {Lstar:0.3f}')
    print(f" Chamber Volume [m^3]: {V_chamber:0.5f}")
    print(f' Contraction Ratio: {A1/A_t:0.2f}')
    print(f' Convergent Angle [deg] {conv_angle}')
    print(f' Chamber Diameter [mm]: {D_c*10**3:0.3f}')
    print(f' Chamber Length [mm]: {L1*10**3:0.3f}')
    print(f' Convergent Length [mm]: {L_conv*10**3:0.3f}')
    print("\n")


    """Balance Momentum Ratios"""
    D_o = np.mean([D_ox, D_f]) #mean orifice diameter 
    if choked:
        rho_star_ox = rho_ox * (2 / (k_ox+1)) ** (1 / (k_ox-1))
        rho_star_f = rho_f * (2 / (k_f+1)) ** (1 / (k_f-1))
        v_inj_ox = mDot_ox/(rho_star_ox*A_ox)
        v_inj_f = mDot_f/(rho_star_f*A_f)
    else:
        v_inj_ox = mDot_ox/(rho_ox*A_ox)
        v_inj_f = mDot_f/(rho_f*A_f)

    theta_ox = 90 - (45) # oxidizer injection angle [deg]
    mom_ox = mDot_ox * v_inj_ox * np.sin(np.deg2rad(theta_ox))
    theta_f = np.rad2deg(np.asin( mom_ox / (mDot_f * v_inj_f) ))

    D = 0.1 * L1
    #LDo = 2
    #D = LDo * D_o * (np.tan(np.deg2rad(theta_ox)) + np.tan(np.deg2rad(theta_f)))**-1
    l_f = D * np.tan(np.deg2rad(theta_f))
    l_ox = D * np.tan(np.deg2rad(theta_ox))

    print("\nInjector Sizing")
    print(f" Oxidizer Angle: {theta_ox:0.2f}")
    print(f" Fuel Angle: {theta_f:0.2f}")
    print(f" Impingement Distance [mm]: {D*10**3:0.3f}")
    print(f" Ox Orifice Radial Distance [mm]: {l_ox*10**3:0.3f}")
    print(f" Fuel Orifice Radial Distance [mm]: {l_f*10**3:0.3f}")
    print("\n")


    """Manifold Sizing"""
    tD = 3 # thickness/orifice diameter: typically 4-10
    t_orifice = tD * D_o # orifice thickness [m]

    # size orifice volume to be 5 times that of orifice
    V_orifice = t_orifice * np.max([A_ox, A_f])
    V_manifold = 5 * V_orifice

    #circle shaped manifolds
    D_manifold = 1/4 * D_c
    A_manifold = diameter_to_area(D_manifold)    #area of each manifold
    h_manifold = V_manifold/A_manifold            #height of each manifold

    print("\nManifold Sizing")
    print(f" Orifice Thickness [mm]: {t_orifice * 10**3:0.3f}")
    print(f" Manifold Volume [mm^3]: {V_manifold * (10**3)**3:0.3f}")
    print(f" Manifold Area [mm^2]: {A_manifold * (10**3)**2:0.3f}")
    print(f" Manifold Height [mm]: {h_manifold * 10**3:0.3f}")
    print(f" Manifold Diameter [mm]: {D_manifold * 10**3:0.3f}")
    print("\n")

    ## Check Injector Sizing Fits Chamber Diameter ##
    used = np.max([l_ox, l_f]) + D_manifold/2
    avail =  (D_c - 5*10**-3) / 2
    print("\nCheck Injector Fits Chamber Diameter")
    print(f" Radial Distance Occupied [mm]: {used*10**3:0.3f} ")
    print(f" Radial Distance Available [mm]: {avail * 10**3:0.3f} ")
    if (used > avail):
        print(" Sizing INVALID")
    else:
        print(" Sizing VALID")
    print("\n")

    """Bolt Calcs"""

    # Thrust Force on mounted plate

    # Pressure Force on Injector Face


    """Line Sizing"""
    # find injeciton velocity
    D_line_ox, v_line_ox = size_line(mDot_ox, rho_ox, 50)
    D_line_f, v_line_f = size_line(mDot_f, rho_f, 160)
    print("\nLine Sizing")
    print(f" Ox Line: ORB -{D_line_ox * 16:0.0f}")
    print(f" Ox Line Velocity [m/s]: {v_line_ox:0.3f}")
    print(f" Fuel Line: ORB -{D_line_f * 16:0.0f}")
    print(f" Fuel Line Velocity [m/s]: {v_line_f:0.3f}")

   

    print(f" Ox Injection Velocity [m/s]: {v_inj_ox:0.3f}")
    print(f" Fuel Injection Velocity [m/s]: {v_inj_f:0.3f}")
    print("\n")
    print(f"Ox Sonic [m/s]: {sonic_ox:0.3f}")
    print(f"Fuel Sonic [m/s]: {sonic_f:0.3f}")



    """Wall Sizing"""
    """
    sigma_yield = []
    temps = []
    FOS = 1.5
    for i in range(len(sigma_yield)):
        t = p_c*(D_c) / (2*sigma_yield[i])
    plt.plot(temps, t*(10**3))
    plt.ylabel("Wall Thickness [m]")
    plt.xlabel("Wall Temperature [C]")
    plt.title(f"Required Thickness for {FOS} FOS vs. Wall Temp")"""




    """Make Cad Parameter Sheet"""
    sheet = []

    sheet = add_cad_param(sheet, "Conv_Angle", conv_angle, "deg")
    sheet = add_cad_param(sheet, "Ox_Angle", theta_ox, "deg")
    sheet = add_cad_param(sheet, "Fuel_Angle", theta_f, "deg")

    sheet = add_cad_param(sheet, "d_t", D_t, "m")
    sheet = add_cad_param(sheet, "d_ox", D_ox, "m")
    sheet = add_cad_param(sheet, "d_f", D_f, "m")
    sheet = add_cad_param(sheet, "d_c", D_c, "m")

    sheet = add_cad_param(sheet, "L_c", L1, "m")
    sheet = add_cad_param(sheet, "L_conv", L_conv, "m")

    sheet = add_cad_param(sheet, "d_man", D_manifold, "m")
    sheet = add_cad_param(sheet, "h_man", h_manifold, "m")

    sheet = add_cad_param(sheet, "t_wall", 0.35, "in")

    sheet = add_cad_param(sheet, "l_ox", l_ox, "m")
    sheet = add_cad_param(sheet, "l_f", l_f, "m")
    sheet = add_cad_param(sheet, "t_inj", t_orifice, "m")
    sheet = add_cad_param(sheet, "D_imp", D, "m")


    df = pd.DataFrame(sheet)
    df.to_csv("TCA_Igniter_dims.csv", index=False, header=False)
    
    return

if __name__== "__main__":
    main()
