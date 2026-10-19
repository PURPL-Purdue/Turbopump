import numpy as np
from rocketcea.cea_obj import CEA_Obj, add_new_fuel
from scipy.optimize import brentq
from pyfluids import Fluid, FluidsList, Input

PSI2PA = 6894.76
LBM2KG = 0.453592
IN2M = 0.0254
THOU2M = IN2M/1000
FT2M = 0.3048
T_AMB_CELSIUS = 20
T_AMB_KELVIN = 273
R2F = -458.67
R2K = 0.555556

def dittus_boelter(Re, Pr, fluid_being_cooled):
    if fluid_being_cooled:
        n = 0.3
    else:
        n = 0.4
    Nu = 0.023 * Re**0.8 * Pr**n
    if Re < 10000:
        print(f"dittus-boelter equation does not hold:  Re = {Re:0.0f}")
        return None
    elif Pr < 0.6 or Pr > 160:
        print(f"dittus-boelter equation does not hold:  Pr = {Pr:0.3f}")
        return None
    return Nu

def gnielinski(f, Re, Pr):
    Nu = (f/8)*(Re - 1000)*Pr / (1 + 12.7*(f/8)**0.5*(Pr**(2/3) - 1))
    if Re < 3000 or Re > 5 * 10**6:
        print(f"gnielinski equation does not hold:  Re = {Re:0.0f}")
        return None
    elif Pr < 0.5 or Pr > 2000:
        print(f"gnielinski equation does not hold:  e/D = {Pr:0.0f}")
        return None
    return Nu

def friction_factor(e, D, Re):
    """Swamee and Jain approximation of Colebrook Equation"""

    f = 0.25 / (np.log10((e/D)/3.7 + 5.74/Re**0.9))**2
    
    if Re < 5000 or Re > 10**8:
        print(f"friction factor equation does not hold:  Re = {Re:0.0f}")
        return None
    elif e/D < 0.000001 or e/D > 0.05:
        print(f"friction factor equation does not hold:  e/D = {e/D:0.0f}")
        return None
    return f

def mach_from_area_ratio(AR, gam):
    """Subsonic Mach number from isentropic A/A*."""
    g = gam
    f = lambda M: (1/M)*((2/(g+1))*(1 + (g-1)/2*M**2))**((g+1)/(2*(g-1))) - AR
    return brentq(f, 1e-6, 1 - 1e-9, xtol=1e-6, rtol=1e-6)

def torch():
    """
    
    Pressures
        Ox: 505.1 psi [Stiffness: 153%]
        Fuel: 490.8 psi [Stiffness: 145%]
        Chamber: 200.0 psi

    Mass Flows
        Torch: 0.0783 lbm/s [MIE Requirement: 0.0550 lbm/s]
        Ox: 0.0470 lbm/s
        Fuel: 0.0313 lbm/s

    Characteristics
        Cstar: 1413.4 m/s [Efficiency: 0.90]
        AFT: 2694.6 deg F
    """
    #hardware/setpoint definition
    
    #chamber-throat
    D_t = 0.268 #in 
    D_c = 0.457 # in
    e = 250 * 10**-6 
    Pc = 200 # psi
    mdot = 0.0783 #lbm/s
    OF = 1.5 # --
    fuel = 'GCH4'
    ox = 'GOX'
    frozen_flag = 0

    # orifices
    D_ox = 63 # thou
    D_f = 63 # thou
    p_line_f = 505.1 #psi
    p_line_ox = 490.8 #psi
    M_f = 1
    M_ox = 1
    gam_ox = 1.3
    gam_f = 1.4
        
    
    # convert to SI
    Pc *= PSI2PA
    mdot *= LBM2KG
    D_t *= IN2M
    D_c *= IN2M
    e *= IN2M
    D_ox *= THOU2M
    D_f *= THOU2M
    p_line_ox *= PSI2PA
    p_line_f *= PSI2PA
    
    # run cea to get combustion product properties
    cea = CEA_Obj(fuelName=fuel, oxName=ox)
    [T_c, T_t, _] = cea.get_Temperatures(Pc=Pc/PSI2PA, MR=OF, eps=1, frozen=frozen_flag, frozenAtThroat=frozen_flag)
    T_c *= R2K
    T_t *= R2K

    # get transport properties # [(BTU/lbm degR), millipoise, mcal/cm-K-s, -- ] -> [J/(kg-K), Pa-s, W/(m-K), --]
    [Cp_c, mu_c, k_c, Pr_c] = cea.get_Chamber_Transport(Pc=Pc/PSI2PA, MR=OF, frozen=frozen_flag) 
    [Cp_t, mu_t, k_t, Pr_t] = cea.get_Throat_Transport(Pc=Pc/PSI2PA, MR=OF, eps=1, frozen=frozen_flag)
    Cp_c *= 4186.8
    Cp_t *= 4186.8
    mu_c *= 0.0001
    mu_t *= 0.0001
    k_c *= 0.4184
    k_t *= 0.4184

    # determine reynolds number
    
    #flow velocity [m/s] (returned in ft/s)
    
    [_, gam_c] = cea.get_Chamber_MolWt_gamma(Pc=Pc/PSI2PA, MR=OF, eps=1)
    [_, gam_t] = cea.get_Throat_MolWt_gamma(Pc=Pc/PSI2PA, MR=OF, eps=1, frozen=frozen_flag)
    M_c = mach_from_area_ratio((D_c/D_t)**2, gam_c)

    [a_c, a_t, _] = cea.get_SonicVelocities(Pc=Pc/PSI2PA, MR=OF, eps=1, frozen=frozen_flag, frozenAtThroat=frozen_flag)
    a_t *= FT2M
    a_c *= FT2M

    M_t = 1
    v_c = a_c * M_c
    v_t = a_t * M_t

    #densities [kg/m^3] (returned in lbm/ft^3)
    [rho_c, rho_t, _] = cea.get_Densities(Pc=Pc/PSI2PA, MR=OF, eps=1, frozen=frozen_flag, frozenAtThroat=frozen_flag)
    rho_c *=  (LBM2KG / FT2M**3)
    rho_t *=  (LBM2KG / FT2M**3)

    #use mdot for reynolds number since we asusme a Cd in our sizing
    Re_c = 4 * mdot / (np.pi * D_c * mu_c)
    Re_t = 4 * mdot / (np.pi * D_t * mu_t)
    
    # friciton factor:
    f_c = friction_factor(e, D_c, Re_c)
    f_t = friction_factor(e, D_t, Re_t)

    #nusselt number
    Nu_c = gnielinski(f_c, Re_c, Pr_c)
    Nu_t = gnielinski(f_t, Re_t, Pr_t)

    #HTC
    HTC_c = Nu_c * k_c / D_c
    HTC_t = Nu_t * k_t / D_t


    # calculate recovery temperatures
   
    r_c = Pr_c**(1/3)
    r_t = Pr_t**(1/3)
    
    Tr_c = T_c * (1 + r_c*((gam_c-1)/2)*M_c**2) / (1 + ((gam_c-1)/2)*M_c**2)
    Tr_t = T_t*(1 + r_t * ((gam_t-1)/2)*(M_t**2))

    print()
    print("Heating HTCs")
    print(f"Adiabatic Flame Temp: {T_c:0.0f} K")
    print(f"Chamber: {HTC_c:0.0f} W/m^2-K @ {Tr_c:0.0f} K")
    print(f"Throat: {HTC_t:0.0f} W/m^2-K @ {Tr_t:0.0f} K")
    print()
    print("Non Dimenional Properties")
    print(f"Chamber: M = {M_c:0.3f}, Re = {Re_c:0.0f}, Pr = {Pr_c:0.3f}, e/D = {e/D_c:.2e}")
    print(f"Throat: M = {M_t}, Re = {Re_t:0.0f}, Pr = {Pr_t:0.3f}, e/D = {e/D_t:.2e}")
    print()

    # HTCs for cooling at the orifices
    mdot_ox = mdot * OF / (1 + OF)
    mdot_f = mdot * 1 / (1 + OF)

    #pyfluids for ox and fuel
    fuel_PYF = Fluid(FluidsList.Methane).with_state(Input.pressure(p_line_f), Input.temperature(T_AMB_CELSIUS))
    ox_PYF = Fluid(FluidsList.Oxygen).with_state(Input.pressure(p_line_ox), Input.temperature(T_AMB_CELSIUS))

    mu_f = fuel_PYF.dynamic_viscosity            # [Pa*s]
    mu_ox = ox_PYF.dynamic_viscosity             # [Pa*s]    

    Re_ox = 4 * mdot_ox / (np.pi * D_ox * mu_ox)
    Re_f = 4 * mdot_f / (np.pi * D_f * mu_f)

    Pr_ox = fuel_PYF.prandtl
    Pr_f = ox_PYF.prandtl

    Nu_ox = dittus_boelter(Re_ox, Pr_ox, fluid_being_cooled=False)
    Nu_f = dittus_boelter(Re_f, Pr_f, fluid_being_cooled=False)

    HTC_ox = Nu_ox * k_c / D_c
    HTC_f = Nu_f * k_c / D_c

    # determine recovery temepratures
    r_ox = Pr_ox**(1/3)
    r_f = Pr_f**(1/3)
    
    Tr_ox = T_AMB_KELVIN * (1 + r_ox*((gam_ox-1)/2)*M_ox**2) / (1 + ((gam_ox-1)/2)*M_ox**2)
    Tr_f = T_AMB_KELVIN * (1 + r_f*((gam_f-1)/2)*M_f**2) / (1 + ((gam_f-1)/2)*M_f**2)
    
    print("Cooling HTCs")
    print(f"Ox Orifice: {HTC_ox:0.0f} W/m^2-K @ {Tr_ox:0.0f} K")
    print(f"Fuel Orifice: {HTC_f:0.0f} W/m^2-K @ {Tr_f:0.0f} K")
    print()
    print("Non Dimenional Properties")
    print(f"Ox: M Re = {Re_ox:0.0f}, Pr = {Pr_ox:0.3f}")
    print(f"Fuel: Re = {Re_f:0.0f}, Pr = {Pr_f:0.3f}")

    return

def chamber():
    
    #hardware/setpoint definition
    card_str = """
    fuel IPA  C 3 H 8 O 1
    h,cal=-65180.1  t(k)=298.0  wt%=100.0
    """
    add_new_fuel('IPA', card_str)

    #chamber-throat
    D_c = 3.826 # in
    D_t = 1.8012 # in
    e = 250 * 10**-6 
    Pc = 440 # psi
    mdot = 3.8050 #lbm/s
    OF = 1 # --
    fuel = 'IPA'
    ox = 'LOX'
    frozen_flag = 0

    
    # convert to SI
    Pc *= PSI2PA
    mdot *= LBM2KG
    D_c *= IN2M
    D_t *= IN2M
    e *= IN2M

    
    # run cea to get combustion product properties
    cea = CEA_Obj(fuelName=fuel, oxName=ox)
    [T_c, _, _] = cea.get_Temperatures(Pc=Pc/PSI2PA, MR=OF, eps=1, frozen=frozen_flag, frozenAtThroat=frozen_flag)
    T_c *= R2K

    # get transport properties # [(BTU/lbm degR), millipoise, mcal/cm-K-s, -- ] -> [J/(kg-K), Pa-s, W/(m-K), --]
    [Cp_c, mu_c, k_c, Pr_c] = cea.get_Chamber_Transport(Pc=Pc/PSI2PA, MR=OF, frozen=frozen_flag) 
    Cp_c *= 4186.8
    mu_c *= 0.0001
    k_c *= 0.4184


    # determine reynolds number
    
    #flow velocity [m/s] (returned in ft/s)
    
    [_, gam_c] = cea.get_Chamber_MolWt_gamma(Pc=Pc/PSI2PA, MR=OF, eps=1)
    M_c = mach_from_area_ratio((D_c/D_t)**2, gam_c)

    [a_c, _, _] = cea.get_SonicVelocities(Pc=Pc/PSI2PA, MR=OF, eps=1, frozen=frozen_flag, frozenAtThroat=frozen_flag)
    a_c *= FT2M

    #densities [kg/m^3] (returned in lbm/ft^3)
    [rho_c, _, _] = cea.get_Densities(Pc=Pc/PSI2PA, MR=OF, eps=1, frozen=frozen_flag, frozenAtThroat=frozen_flag)
    rho_c *=  (LBM2KG / FT2M**3)

    #use mdot for reynolds number since we asusme a Cd in our sizing
    Re_c = 4 * mdot / (np.pi * D_c * mu_c)
    
    # friciton factor:
    f_c = friction_factor(e, D_c, Re_c)

    #nusselt number
    Nu_c = gnielinski(f_c, Re_c, Pr_c)

    #HTC
    HTC_c = Nu_c * k_c / D_c


    # calculate recovery temperatures
    
    r_c = Pr_c**(1/3)
    
    Tr_c = T_c * (1 + r_c*((gam_c-1)/2)*M_c**2) / (1 + ((gam_c-1)/2)*M_c**2)

    print()
    print("Main Chamber HTC")
    print(f"Adiabatic Flame Temp: {T_c:0.0f} K")
    print(f"Chamber: {HTC_c:0.0f} W/m^2-K @ {Tr_c:0.0f} K")
    print()
    print("Non Dimenional Properties")
    print(f"Chamber: M = {M_c:0.3f}, Re = {Re_c:0.0f}, Pr = {Pr_c:0.3f}, e/D = {e/D_c:.2e}")
    print()
    
def main():
    torch()
    chamber()

if __name__ == "__main__":
    main()