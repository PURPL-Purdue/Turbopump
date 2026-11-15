
import numpy as np
import sys
from pint import Quantity as Q_
from pathlib import Path
from rocketcea.cea_obj import CEA_Obj

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from Inputs.Gas_Generator_DataClass import hardware as GG

_cea = CEA_Obj(oxName=GG.GG_ox_type, fuelName=GG.GG_fuel_type)

# Solver settings
WEIGHT = 0.2          # under-relaxation factor
LOOP_LIMIT = 1000
TOLERANCE = 0.05        # convergence tolerance on chamber pressure

def Pc_guess(p1: Q_, p2: Q_) -> Q_:
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
    cstar = _cea.get_Cstar(Pc=P_c.to('psi').m, MR=float(OF))  # MR must be a plain float
    return GG.Cstar_eff * Q_(cstar, 'ft/s')

def orifice_mdot(cda: Q_, pup: Q_, Pc: Q_, rho: Q_) -> Q_:
    """Incompressible orifice flow: mdot = Cda * sqrt(2 * rho * dP), with dP floored at 0."""
    Pc = Q_(0.0, 'psi') if Pc is None else Pc
    dp = max(pup - Pc, Q_(0.0, 'psi'))
    return (cda * np.sqrt(2 * rho * dp)).to('lb/s')

def GasGen_state_single(P_ox: Q_, P_fu: Q_, insert_id: str) -> tuple[Q_, Q_, Q_]:
    """Calculates Gas Generator Total State
 
    Inputs: OX and FU feed pressures (pint), insert id (str).
    Outputs: Mdot, Pc, OF, or NaN if the solve does not converge.
 
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
 
    def solve(P_c: Q_) -> tuple[Q_, Q_, Q_]:
        mdot_ox, mdot_fuel, mdot_film = find_mass_flows(P_c)
        mdot_fuel_total = mdot_fuel + mdot_film
        total_mdot = (mdot_ox + mdot_fuel_total).to('lb/s')
        total_OF = (mdot_ox / mdot_fuel_total).to('dimensionless')
        cstar = get_cstar(P_c, total_OF)
        Pc_calc = (cstar * total_mdot / GG.Stator_Area).to('psi')
        return total_mdot, Pc_calc, total_OF

    Pc_init = Pc_guess(P_ox, P_fu)
    for _ in range(LOOP_LIMIT):
        total_mdot, Pc_calc, total_OF = solve(Pc_init)
        if abs(Pc_init.to('psi').m - Pc_calc.to('psi').m) < TOLERANCE:
            return total_mdot.to('lb/s'), Pc_calc.to('psi'), total_OF.to('dimensionless')
        Pc_init = Pc_init + WEIGHT * (Pc_calc - Pc_init)
 
    return Q_(np.nan, 'lb/s'), Q_(np.nan, 'psi'), Q_(np.nan, 'dimensionless')