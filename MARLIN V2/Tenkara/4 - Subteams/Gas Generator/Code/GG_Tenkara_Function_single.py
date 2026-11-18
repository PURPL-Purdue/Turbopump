import numpy as np
import matplotlib.pyplot as plt
import sys
from pint import Quantity as Q_
from pathlib import Path
from rocketcea.cea_obj import CEA_Obj

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from Inputs.Gas_Generator_DataClass import hardware as GG

cea = CEA_Obj(oxName=GG.GG_ox_type, fuelName=GG.GG_fuel_type)

# Solver settings
WEIGHT = 0.1          # under-relaxation factor
LOOP_LIMIT = 100
TOLERANCE = 0.1        # convergence tolerance on chamber pressure

def Pc_initial(p1: Q_, p2: Q_) -> Q_:
    """Guesses initial Chamber Pressure based on inlet pressures."""
    out = 0.8 * min(p1.to('psi').m, p2.to('psi').m)
    return Q_(out, 'psi')

def series_cda(*cdas: Q_) -> Q_:
    """Effective Cda of orifices in series: (sum(1/c^2))^-1/2.
    Returns a zero with area units if any Cda is zero."""
    if any(c == 0 for c in cdas):
        return 0.0 * cdas[0]
    return sum(1.0 / c**2 for c in cdas) ** -0.5

def get_cstar(P_c: Q_, OF: Q_) -> Q_:
    """Return c*(OF, Pc) in ft/s (from NASA CEA via RocketCEA, scaled by efficiency)."""
    cstar = cea.get_Cstar(Pc=P_c.to('psi').m, MR=OF.m)
    return GG.Cstar_eff * Q_(cstar, 'ft/s')

def orifice_mdot(cda: Q_, P_up: Q_, P_c: Q_, rho: Q_) -> Q_:
    """Incompressible orifice flow: mdot = Cda * sqrt(2 * rho * dP), with dP floored at 0."""
    dp = max(P_up - P_c, Q_(0.0, 'psi'))
    return (cda * np.sqrt(2 * rho * dp)).to('lb/s')

def GasGen_state_single(P_ox: Q_, P_fu: Q_, insert_id: str) -> dict[str, Q_]:
    """Calculates Gas Generator Total State
 
    Inputs: OX and FU feed pressures (pint), insert id (str).
    Outputs: OX Mdot and FU Mdot dictionary (str) (pint).
 
    Method:
      Guess a chamber pressure, compute the mass flows that result from that dP,
      check whether c* of that combustion reproduces the guessed Pc, and iterate
      (under-relaxed) until the guessed and calculated Pc agree.
    """
        
    def find_mass_flows(P_c: Q_) -> tuple[Q_, Q_, Q_]:
        mdot_ox = orifice_mdot(GG.LOx_Cda, P_ox, P_c, GG.LOx_rho)
        mdot_fuel = orifice_mdot(GG.IPA_Cda, P_fu, P_c, GG.IPA_rho)
        film_Cda = series_cda(GG.Insert_Cda[insert_id], GG.IPA_Cda)
        mdot_film = orifice_mdot(film_Cda, P_fu, P_c, GG.IPA_rho)
        return mdot_ox, mdot_fuel, mdot_film
 
    def solve(Pc_guess: Q_) -> tuple[Q_, Q_, Q_]:
        mdot_ox, mdot_fuel, mdot_film = find_mass_flows(Pc_guess)
        mdot_fuel_total = mdot_fuel + mdot_film
        total_mdot = (mdot_ox + mdot_fuel_total)
        if mdot_fuel_total.m != 0:
            total_OF = (mdot_ox / mdot_fuel_total).to('dimensionless')
        else:
            total_OF = Q_(0.0, 'dimensionless')
        cstar = get_cstar(Pc_guess, total_OF)
        Pc_calc = (cstar * total_mdot / GG.Stator_Area).to('psi')
        return Pc_calc, mdot_ox, mdot_fuel_total

    Pc_guess = Pc_initial(P_ox, P_fu)
    for _ in range(LOOP_LIMIT):
        Pc_calc, mdot_ox, mdot_fu = solve(Pc_guess)
        if abs(Pc_calc.to('psi').m - Pc_guess.to('psi').m) < TOLERANCE:
            if Pc_calc.to('psi').m < min(P_ox.to('psi').m, P_fu.to('psi').m):
                return {"OX Mdot": mdot_ox.to('lb/s'), "FU Mdot": mdot_fu.to('lb/s')}
            else:
                return {"OX Mdot": Q_(np.nan, 'lb/s'), "FU Mdot": Q_(np.nan, 'lb/s')}
        Pc_guess = Pc_guess + WEIGHT * (Pc_calc - Pc_guess)
    return {"OX Mdot": Q_(np.nan, 'lb/s'), "FU Mdot": Q_(np.nan, 'lb/s')}