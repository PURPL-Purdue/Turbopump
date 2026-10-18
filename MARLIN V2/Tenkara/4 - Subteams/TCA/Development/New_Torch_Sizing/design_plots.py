from rocketcea.cea_obj import CEA_Obj
import numpy as np
import matplotlib.pyplot as plt

PSI2PA = 6894.76
LBM2KG = 0.453592
IN2M = 0.0254
THOU2M = IN2M/1000
FT2M = 0.3048
T_AMB_CELSIUS = 20
T_AMB_KELVIN = 273
R2F = -458.67

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

# Function: Upstream pressure that chokes flow over an orifice
def choked_pressure(mdot, area, k, T, R, Cd=1.0):
    """Upstream stagnation pressure that chokes mdot through an orifice."""
    factor = (2 / (k + 1)) ** ((k + 1) / (2 * (k - 1)))
    return mdot / (Cd * area * np.sqrt(k / (R * T) * factor))

def main():
    #pick Pc, At, solve for cstar
    #-> determine mdot, if mdot is greater than required MIE, good. (provide minimum mdot)
    # now we have Pc, At, mdot, OF
    # solve for line pressures

    ## inputs ##
    #hardware definition
    dt = 0.268 #in 
    df = 63 #thou
    dox = 63 #thou
    Cdf = 0.8
    Cdox = 0.8
    eta_cstar = 0.9

    #setpoint definition
    Pc = 200 # psi
    mdot_min = 0.055 #lbm/s
    OF = 1.5 # --
    fuel = 'CH4'
    ox = 'GOX'

    #propellant properties
    kf = 1.30
    kox = 1.4
    Rox = 259.84 # J/kg-K
    Rf = 518.28 # J/kg-K

    ##preliminary calcs##
    
    # convert to SI
    Pc *= PSI2PA
    dt *= IN2M
    df *= THOU2M
    dox *= THOU2M
    mdot_min *= LBM2KG


    #calculate orifice areaa
    At = diameter_to_area(dt)
    Af = diameter_to_area(df)
    Aox = diameter_to_area(dox)
    
    ##sizing calcs##

    # run cea to get cstar, determine actual mdot of system
    cea = CEA_Obj(fuelName=fuel, oxName=ox)
    cstar = cea.get_Cstar(Pc=Pc/PSI2PA, MR=OF) * eta_cstar * FT2M # [m/s] (returns ft/s)
    Tc = cea.get_Tcomb(Pc=Pc/PSI2PA, MR=OF) + R2F
    mdot = Pc * At / cstar # kg/s
    
    #split mass flows
    mdotox = mdot * OF / (OF + 1)
    mdotf = mdot * 1/(OF+1)
    if mdot < mdot_min:
        print("Mass flow does not meet requirement.")
        return

    # determine line pressures
    ## Get fluid properties to size orifice/throat areas ##
    #pyfluids for ox and fuel
    Pox = choked_pressure(mdotox, Aox, kox, T_AMB_KELVIN, Rox, Cdox)
    Pf = choked_pressure(mdotf, Af, kf, T_AMB_KELVIN, Rf, Cdf)

    print("Pressures")
    print(f"  Ox: {Pox/PSI2PA:0.1f} psi [Stiffness: {(Pox-Pc)/Pc*100:0.0f}%]")
    print(f"  Fuel: {Pf/PSI2PA:0.1f} psi [Stiffness: {(Pf-Pc)/Pc*100:0.0f}%]")
    print(f"  Chamber: {Pc/PSI2PA:0.1f} psi")
    print()

    print("Mass Flows")
    print(f"  Torch: {mdot/LBM2KG:0.4f} lbm/s [MIE Requirement: {mdot_min/LBM2KG:0.4f} lbm/s]")
    print(f"  Ox: {mdotox/LBM2KG:0.4f} lbm/s")
    print(f"  Fuel: {mdotf/LBM2KG:0.4f} lbm/s")
    print()
    
    print("Characteristics")
    print(f"  Cstar: {cstar:0.1f} m/s [Efficiency: {eta_cstar:0.2f}]")
    print(f"  AFT: {Tc:0.1f} deg F")
    print(f"{cea.get_Throat_Isp() * mdot * 9.81}")
    print()
    return

if __name__ == "__main__":
    main()