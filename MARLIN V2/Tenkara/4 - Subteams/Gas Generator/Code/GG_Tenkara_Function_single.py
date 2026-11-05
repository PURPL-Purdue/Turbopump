
import numpy as np
from pint import Quantity as Q_
from pathlib import Path
import sys
sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from Inputs.Gas_Generator_DataClass import hardware as GG


def GasGen_mdot_single(P_ox: Q_, P_Fuel: Q_, insert_id: str) -> Q_:
    ## Map of iterative function
    # Inputs: Pressures, Insert id, and Hardware class
    # Outputs: Mass Flow
    # Convergence criteria: Mdot_new_fuel = Mdot_old_fuel && same for ox
    # Process of doing calculation:
    #   Guess a chamber pressure, calculate mass flow which results from this dp,
    #   Check if the cstar of this combustion yields the guessed chamber pressure,
    #   Iterate the guessed chamber pressure until loop is satisfied.

    weight = 0.01
    loop_limit = 1000
    
    def Pc_guess(p1: Q_, p2: Q_) -> Q_:
        out = 0.81 * np.min([p1.to('psi').m, p2.to('psi').m])
        return Q_(out, 'psi')
    
    def series_cda(*cdas):                # cdas: 1-D float array
        total = 0.0
        for c in cdas:
            total += 1.0 / c**2
        return total ** -0.5

    def get_cstar(of: float, pc_pa: Q_) -> Q_:
        """Return c*(OF, P_c[Pa]) -> ft/s (ideal, from NASA CEA via RocketCEA)."""
        from rocketcea.cea_obj import CEA_Obj
        cea = CEA_Obj(oxName=GG.GG_ox_type, fuelName=GG.GG_fuel_type)
        cstar = GG.CSTAR_eff * Q_(cea.get_Cstar(Pc=pc_pa.to('psi').m, MR=of), 'ft/s')
        return cstar

    
    def mdot(cda: Q_, pup: Q_, Pc: Q_, rho: Q_) -> Q_:
        Pc = Q_(0.0, "psi") if Pc is None else Pc
        if (pup - Pc) == 0:
            mdot_out = 0
        else:
            mdot_out = cda * np.sqrt(2 * rho * np.max((pup - Pc),0))
        return mdot_out.to('lb/s')
    
    def find_mass_flow(P_ox: Q_, P_fuel: Q_, Pc: Q_):
        mdot_ox = mdot(GG.LOx_Cda, P_ox, Pc, GG.LOx_rho)
        mdot_fuel = mdot(GG.Fuel_Cda, P_fuel, Pc, GG.IPA_rho)
        film_Cda = series_cda(GG.Insert_Cda[insert_id], GG.Fuel_Cda)
        mdot_film = mdot(film_Cda, P_fuel, Pc, GG.IPA_rho)
        return mdot_ox, mdot_fuel, mdot_film

    Pc_init = Pc_guess(P_ox, P_Fuel)
    mdot_ox, mdot_fuel, mdot_film = find_mass_flow(P_ox, P_Fuel, Pc_init)
    total_mdot = (mdot_ox + mdot_fuel+ mdot_film).to('lb/s')
    core_OF = (mdot_fuel / mdot_ox).to('dimensionless')
    total_OF = ((mdot_fuel + mdot_film) / mdot_ox).to('dimensionless')
    cstar_init = get_cstar(total_OF, Pc_init)

    Pc = Q_(cstar_init * total_mdot / GG.Stator_Area, 'psi')

    loop_num = 0
    while True:
        loop_num += 1

        if np.abs(Pc_init.to('psi').m - Pc.to('psi').m) < 0.05:
            break
        elif loop_num >= loop_limit:
            return np.nan
        
        error = Pc - Pc_init
        Pc_init = Pc_init + weight * error
        mdot_ox, mdot_fuel, mdot_film = find_mass_flow(P_ox, P_Fuel, Pc_init)
        total_mdot = (mdot_ox + mdot_fuel + mdot_film).to('lb/s')
        core_OF = (mdot_fuel / mdot_ox).to('dimensionless')
        total_OF = ((mdot_fuel + mdot_film) / mdot_ox).to('dimensionless')
        cstar_init = get_cstar(total_OF, Pc_init)
        Pc = Q_(cstar_init * total_mdot / GG.Stator_Area, 'psi')

    print(f"Final Pc: {round(Pc.to('psi').m, 1)}")
    return total_mdot.to('lb/s')