"""
This file contains functions pertaining to compressible flow.

Index:
 - mach_to_area
 - area_to_mach
 - stagnation_properties
 - normal_shock

Future implementations:


Author: Krish Kalmadi
Created: 08/23/2026
Last Modified: 09/08/2026
"""

from pint import UnitRegistry as ureg, Quantity as Q_
import cea
import CoolProp as CoolProp

from cea_setup import cea_setup
import numpy as np

def mach_to_area(chamberDiameter, machChamber, localMach, gam=1.4,):
    """
    Description: Determine area from mach.
    
    Methods: Isentropic Flow Area-Mach relation

    Arguments:
        chamberDiameter (float) Plenum Diameter [length] (unit impartial)
        machChamber (float) Mach number in the chamber (unit impartial)
        localMach (float) Local mach number (unit impartial)
        gam (float) Specific heat ratio of working fluid 

    Returns: 
        localArea (float) Local area [length^2]
        
    """
    # determine area from diameter
    chamberArea = (np.pi/4) * chamberDiameter**2

    # determine throat area from chamber area and mach in chamber
    throatArea = chamberArea * (1/machChamber) * ((2/(gam+1)) * (1 + ((gam-1)/2)*(machChamber**2)))**((gam+1)/(2*(gam-1)))
    
    # isentropic flow area ratio
    areaRatio = (1/localMach) * ((2/(gam+1)) * (1 + ((gam-1)/2)*(localMach**2)))**((gam+1)/(2*(gam-1)))
    
    # find local area
    localArea = areaRatio * throatArea

    return localArea

def area_to_mach(chamberDiameter, localDiameter, machChamber, gam=1.4, supersonic=True):
    """
    Description: Determine mach from area.

    Methods: Isentropic Flow Area-Mach relation, finding mach where both sides are equal

    Arguments:
        chamberDiameter (float) Plenum Diameter [length] (unit impartial)
        localDiameter (float) Plenum Diameter [length] (unit impartial)
        machChamber (float) Mach number in the chamber (unit impartial)
        gam (float) Specific heat ratio of working fluid
        supersonic (bool) Flag to indicate if the solution should be supersonic

    Returns: 
        mach (float) Local mach number
    """

    # determine chamber area from diameter
    chamberArea = (np.pi/4) * chamberDiameter**2

    # determine throat area from chamber area and mach in chamber
    throatArea = chamberArea * (1/machChamber) * ((2/(gam+1)) * (1 + ((gam-1)/2)*(machChamber**2)))**((gam+1)/(2*(gam-1)))

    # determine local area from diameter
    localArea = (np.pi/4) * localDiameter**2
    
    # isentropic flow area ratio
    areaRatio = localArea / throatArea

    # area-mach equation minus known area ratio
    def find_zero(mach):
        return ((1/mach) * ((2/(gam+1)) * (1 + ((gam-1)/2)*(mach**2)))**((gam+1)/(2*(gam-1)))- areaRatio)

    # choose bounds depending on desired solution
    if supersonic:
        low = 1.000001
        high = 25.0
    else:
        low = 0.000001
        high = 0.999999

    # bisection tolerance
    tolerance = 1e-6

    while (high - low) > tolerance:

        mach = (low + high) / 2

        if abs(find_zero(mach)) < tolerance:
            return mach

        elif find_zero(low) * find_zero(mach) < 0:
            high = mach

        else:
            low = mach

    # final mach estimate
    mach = (low + high) / 2

    return mach

def stagnation_properties(temperature_initial, pressure_initial, density_initial, mach, gam=1.4):
    """
    Description: Determine stagnation properties from mach and gamma

    Methods: Isentropic ideal-gas stagnation equations

    Arguments:
        temperature_initial (float) temperature of moving fluid [K or R] (unit impartial)
        pressure_initial (float) pressure of moving fluid [Pa, kPa, MPa, psi, bar] (unit impartial)
        density_initial (float) density of moving fluid [kg/m^3 or lbm/ft^3] (unit impartial)
        mach (float) Local mach number
        gam (float) Specific heat ratio of working fluid

    Returns: 
        temperature_stagnation (float) temperature of stagnation fluid [K or R] (unit impartial)
        pressure_stagnation (float) pressure of stagnation fluid [Pa, kPa, MPa, psi, bar] (unit impartial)
        density_stagnation (float) density of stagnation fluid [kg/m^3 or lbm/ft^3] (unit impartial)
        
    """

    # determine isentropic flow factor
    isentrop_ff = 1+((gam-1)/2) * mach**2

    # determine stag temperature
    temperature_stagnation = temperature_initial * isentrop_ff

    # determine stag pressure
    pressure_stagnation = pressure_initial * (isentrop_ff)**((gam)/(gam-1))

    # determine stag density
    density_stagnation = density_initial * (isentrop_ff)**((1)/(gam-1))

    return temperature_stagnation, pressure_stagnation, density_stagnation

def normal_shock(mach_initial, gam=1.4):
    """
    Description:
        Determine downstream Mach number and property ratios across a normal shock.

    Methods:
        Ideal-gas normal shock relations.

    Arguments:
        mach_initial (float): Upstream Mach number, M1 
        gam (float): Specific heat ratio of working fluid

    Returns:
        mach_final (float): Downstream Mach number, M2
        pressure_ratio (float): P2 / P1 (unit impartial)
        density_ratio (float): rho2 / rho1 (unit impartial)
        temperature_ratio (float): T2 / T1 (unit impartial)
        stagnation_pressure_ratio (float): P02 / P01 (unit impartial)
    """
    if mach_initial <= 1:
        raise ValueError("Upstream Mach number must be greater than 1.")

    # downstream Mach number
    mach_final = np.sqrt((2 + ((gam - 1)) * mach_initial**2) / (2 * gam * mach_initial**2 - ((gam - 1))))

    # static pressure ratio
    pressure_ratio = (2 * gam * mach_initial**2 - (gam - 1)) / (gam + 1)

    # density ratio
    density_ratio = (((gam + 1) * mach_initial**2) / ((gam - 1) * mach_initial**2 + 2))

    # temperature ratio
    temperature_ratio = pressure_ratio / density_ratio

    # stagnation pressure ratio
    stagnation_pressure_ratio = ((((gam + 1) * mach_initial**2) / ((gam - 1) * mach_initial**2 + 2)) ** (gam / (gam - 1)) * ((gam + 1)
     / (2 * gam * mach_initial**2 - (gam - 1))) ** (1 / (gam - 1)))

    return (mach_final, pressure_ratio, density_ratio, temperature_ratio, stagnation_pressure_ratio)

##def Bartz(reac_names, reac_temps, of_ratio, Pc, Dt, Twall):
    """
    Description: Determine gas side convection coefficient using Bartz correlation.

    Methods: Bartz correlation for convective heat transfer in rocket nozzles.

    Arguments:
        ureg (pint.UnitRegistry) Unit registry for handling units
        reac_names (list) List of reactant names
        reac_temps (list) List of reactant temperatures
        of_ratio (float) Oxidizer-to-fuel ratio
        Dt (float) Throat diameter
        Pc (float) Chamber pressure
        Twall (list) Wall temperature

    Returns: 
        h_gas (float) Gas side convection coefficient [W/m^2K]
        
    """

     # Attach units if plain numbers were supplied
    if not hasattr(Twall, "units"):
        Twall = Q_(Twall, 'K')

    if not hasattr(Pc, "units"):
        Pc = Q_(Pc, 'psi')

    if not hasattr(Dt, "units"):
        Dt = Q_(Dt, 'm')

    if not hasattr(reac_temps, "units"):
        reac_temps = np.array(reac_temps) * ureg.K

    # Establish CEA solution for given reactants, temperatures, and conditions
    fuel_weights = np.array([1.0, 0.0])
    oxidant_weights = np.array([0.0, 1.0])
    cont_ratio = 1.8
    Prat = (12 * ureg.psi / (1 * ureg.bar)).to(ureg.dimensionless).magnitude

    solution = cea_setup(reac_names, oxidant_weights, fuel_weights, of_ratio, reac_temps.magnitude, Pc.to(ureg.bar).magnitude, Prat, cont_ratio)
    num_pts = solution.num_pts
    T = solution.T * ureg.K
    rho = solution.density * ureg.kg / ureg.m**3

    visc = solution.viscosity * ureg.millipoise
    Cp = solution.cp_eq * ureg.kJ / (ureg.kg * ureg.K)
    Pr = solution.Pr_eq
    Ma = solution.Mach
    sonVel = solution.sonic_velocity * ureg.m / ureg.s
    V = Ma * sonVel

    # Convert units to SI
    visc = visc.to(ureg.Pa * ureg.s)
    Cp = Cp.to(ureg.J / (ureg.kg * ureg.K))

    # Calculate at stagnation temperature
    T0 = T[0]
    Mu0 = visc[0]
    Cp0 = Cp[0]
    Pr0 = Pr[0]

    hg_list = []

    # Calculated at static temperatures
    for i in range(num_pts):
    
        Tinf = T[i]
        Muinf = visc[i]
        Rhoinf = rho[i]
        Vinf = V[i]

        # Rhoam / Rhoinf = Tinf / Tam | ideal gas relation
        Tam = (Tinf + Twall[i]) / 2
        RhoamqRhoinf = Tinf / Tam

        # currently aproximating using frestream and stagnation for parity
        w = np.log(Muinf/Mu0) / np.log(Tinf / T0)
        # use this power law for Muam/Mu0
        MuamqMu0 = (Tam / T0)**w

        hg = (0.026 / Dt**0.2) * ((Mu0**0.2 * Cp0) / Pr0**0.6) * (Rhoinf * Vinf) ** 0.8 *(RhoamqRhoinf)**0.8 * (MuamqMu0)**0.2
        hg = hg.to(ureg.W / (ureg.m**2 * ureg.K))

        hg_list.append(hg)

    return hg_list

##def coolantconvection(temps, pressures, rho, mu, k, cp, mdot, D):
    """
    Description: Determine heat transfer coefficient for coolant convection.

    Methods: Dittus-Boelter correlation for turbulent flow in a pipe.

    Arguments:
        temps (array): Array of temperatures [K] (unit impartial)
        pressures (array): Array of pressures [Pa] (unit impartial)
        rho (array): Array of densities [kg/m^3] (unit impartial)
        mu (array): Array of viscosities [Pa·s] (unit impartial)
        k (array): Array of thermal conductivities [W/(m·K)] (unit impartial)
        cp (array): Array of specific heats [J/(kg·K)] (unit impartial)
        mdot (array): Array of mass flow rates [kg/s] (unit impartial)
        D (array): Array of hydraulic diameters [m] (unit impartial)
     Returns:
        h_coolant (array): Array of heat transfer coefficients [W/m^2-K] (unit impartial)
    """
    if not hasattr(temps, "units"):
        temps = Q_(temps, 'K')

    if not hasattr(pressures, "units"):
        pressures = Q_(pressures, 'Pa')

    if not hasattr(rho, "units"):
        rho = Q_(rho, 'kg/m^3')

    if not hasattr(mu, "units"):
        mu = Q_(mu, 'Pa*s')

    if not hasattr(k, "units"):
        k = Q_(k, 'W/(m*K)')

    if not hasattr(cp, "units"):
        cp = Q_(cp, 'J/(kg*K)')

    if not hasattr(mdot, "units"):
        mdot = Q_(mdot, 'kg/s')

    if not hasattr(D, "units"):
        D = Q_(D, 'm')

    h_coolant = []

    for T, P, r, m, k_val, cp_val, m_dot, d in zip(temps, pressures, rho, mu, k, cp, mdot, D):
        # Calculate Reynolds number
        Re = (4 * m_dot) / (np.pi * d * mu)

        # Calculate Prandtl number
        Pr = (cp_val * mu) / k_val

        # Dittus-Boelter correlation for turbulent flow in a pipe
        Nu = 0.023 * Re**0.8 * Pr**0.4

        # Calculate heat transfer coefficient
        h = (Nu * k_val) / d

        h_coolant.append(h)

    return h_coolant
