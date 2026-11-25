"""
TRBP Gas Generator V2.0
Injector Test Matrix Calculator

Inputs : Total mdot, core OF, insert label, Torch chamber pressure, Torch OF
Outputs: LOX and IPA feed pressures, chamber pressure, overall OF, film cooling %, 
            torch GOx and CH4 feed pressures, torch mass flow, and diagnostics

Usage : python TRBP_GasGenV2_TestMatrixCalculator.py
--input TestMatrix_Inputs.csv --output TestMatrix_Outputs.csv
"""

import sys
import csv
import argparse
from pint import Quantity as Q_
from math import sqrt
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from Inputs.Gas_Generator_DataClass import hardware as GG

from GG_filmcooling import film_core_mixed_temp

from rocketcea.cea_obj import CEA_Obj
ceaGG = CEA_Obj(oxName=GG.GG_ox_type, fuelName=GG.GG_fuel_type)
ceaTorch = CEA_Obj(oxName=GG.torch_ox_type, fuelName=GG.torch_fuel_type)

# --------------------------------------------------------------------------
# Settings
# --------------------------------------------------------------------------

MIN_STIFF = 0.20     # feed-system decoupling rule of thumb
VENTURI_DP = Q_(30, 'psi')  # Venturi pressure drop

# --------------------------------------------------------------------------
# Helper functions
# --------------------------------------------------------------------------

def series_cda(*cdas: Q_) -> Q_:
    """Effective Cda of orifices in series: (sum(1/c^2))^-1/2.
    Returns a zero with area units if any Cda is zero."""
    if any(c == 0 for c in cdas):
        return 0.0 * cdas[0]
    return Q_(sum(1.0 / c**2 for c in cdas) ** -0.5, cdas[0].units)

def orifice_dp(mdot: Q_, cda: Q_, rho: Q_) -> Q_:
    """Pressure drop needed to pass mdot through an orifice of effective area cda."""
    if cda == 0.0:
        return Q_(0.0, 'psi')
    else:
        return ((mdot / cda) ** 2 / (2.0 * rho)).to('psi')

def mass_flow_calculation(mdot_total: Q_, of_core: Q_, insert_id: str) -> dict:
    """Closed-form mass flows. The film flow is set by geometry, not chosen."""
    film_cda_total = series_cda(GG.Insert_Cda[insert_id], GG.IPA_Cda)
    r = film_cda_total / GG.IPA_Cda                  # mdot_film / mdot_fuel
    mdot_fuel = mdot_total / (of_core + 1.0 + r)
    return {
        "mdot_ox": of_core * mdot_fuel,
        "mdot_fuel": mdot_fuel,
        "mdot_film": r * mdot_fuel,
    }

def overall_of(mdot_ox: Q_, mdot_fuel: Q_, mdot_film: Q_) -> Q_:
    """Overall OF: all LOX over all fuel (core + film)."""
    return Q_(mdot_ox / (mdot_fuel + mdot_film), 'dimensionless')

def cooling_percent(mdot_film: Q_, mdot_total: Q_) -> float:
    """Film flow as a percent of total mass flow."""
    return float(100.0 * mdot_film / mdot_total)

def get_GG_cstar(of: Q_, pc: Q_) -> Q_:
    """Return c* from NASA CEA via RocketCEA)."""
    cstar = Q_(GG.Cstar_eff.m * ceaGG.get_Cstar(Pc = pc.to('psi').m, MR=of.m), 'ft/s')
    return cstar.to('m/s')

def get_torch_cstar(of: Q_, pc: Q_) -> Q_:
    """Return c*(OF, P_c[Pa]) -> m/s (ideal, from NASA CEA via RocketCEA)."""
    cstar = Q_(GG.torch_cstar_eff.m * ceaTorch.get_Cstar(Pc = pc.to('psi').m, MR=of.m), 'ft/s')
    return cstar.to('m/s')

def chamber_pressure(mdot_total: Q_, of_overall: Q_):
    """Fixed-point iteration on P_c = c*(OF, P_c) * eff * mdot / (Cd*A_t)."""
    max_iter = 300
    weight = 0.1
    tol = 0.1
    pc = Q_(100, 'psi')
    for _ in range(max_iter):
        cstar = get_GG_cstar(of_overall, pc)
        pc_new = cstar * mdot_total / GG.Stator_Area
        if abs(pc_new - pc) / pc_new < tol:
            return pc_new.to('psi'), cstar
        pc = pc + weight * (pc_new - pc)
    raise RuntimeError(f"P_c iteration did not converge.")

def stiffness(dp: Q_, p_c: Q_) -> float:
    """Stiffness: dp / p_c."""
    return float(dp / p_c)

def add_venturi_dp(p: Q_) -> Q_:
    """Adds venturi pressure drop to all feed pressures"""
    return (p + VENTURI_DP)

# --------------------------------------------------------------------------
# Solver functions
# --------------------------------------------------------------------------

def solve_torch(P_c_torch: Q_, of_torch: Q_) -> dict:
    cstar = get_torch_cstar(of_torch, P_c_torch)
    mdot_torch = (P_c_torch.to('Pa') * GG.torch_eff_throat_area / cstar).to('kg/s')
    mdot_gox = mdot_torch * (of_torch/(1 + of_torch))
    mdot_ch4 = mdot_torch * 1/(1 + of_torch)

    T_up = Q_(278.0, 'K')   # Kelvin

    CH4_term = (2 / (1 + GG.CH4_gamma)) ** ((GG.CH4_gamma + 1) / (GG.CH4_gamma - 1))
    GOx_term = (2 / (1 + GG.GOx_gamma)) ** ((GG.GOx_gamma + 1) / (GG.GOx_gamma - 1))
    
    P_CH4 = ( (mdot_ch4 / GG.CH4_Cda) ** 2 / GG.CH4_gamma * (GG.CH4_R * T_up) / CH4_term ) ** 0.5
    P_GOx = ( (mdot_gox / GG.GOx_Cda) ** 2 / GG.GOx_gamma * (GG.GOx_R * T_up) / GOx_term ) ** 0.5

    outputs = {
        "mdot_torch": mdot_torch,
        "mdot_gox": mdot_gox,
        "mdot_ch4": mdot_ch4,
        "P_CH4": P_CH4,
        "P_GOx": P_GOx,
    }
    return outputs

def solve_GG(mdot_total: Q_, of_core: Q_, insert_id: str) -> dict:
    flows = mass_flow_calculation(mdot_total, of_core, insert_id)
    of_tot = overall_of(flows["mdot_ox"], flows["mdot_fuel"], flows["mdot_film"])
    cool_pct = cooling_percent(flows["mdot_film"], mdot_total)
    pc, cstar = chamber_pressure(mdot_total, of_tot)

    film_cda_total = series_cda(GG.Insert_Cda[insert_id], GG.IPA_Cda)

    dp_ox = orifice_dp(flows["mdot_ox"], GG.LOx_Cda, GG.LOx_rho)
    dp_fuel = orifice_dp(flows["mdot_fuel"], GG.IPA_Cda, GG.IPA_rho)
    dp_film = orifice_dp(flows["mdot_film"], film_cda_total, GG.IPA_rho)

    p_ox = pc + dp_ox
    p_fuel = pc + dp_fuel

    stiff_ox = stiffness(dp_ox, pc)
    stiff_fuel = stiffness(dp_fuel, pc)
    stiff_film = stiffness(dp_film, pc)

    p_ox = add_venturi_dp(p_ox)
    p_fuel = add_venturi_dp(p_fuel)

    temp = Q_(film_core_mixed_temp(p_ox.to('Pa').m, p_fuel.to('Pa').m, insert_id), 'K')

    res = {
        **flows,
        "P_c": pc,
        "of_overall": of_tot,
        "cooling_pct": cool_pct,
        "P_ox": p_ox,
        "P_fuel": p_fuel,
        "stiff_ox": stiff_ox,
        "stiff_fuel": stiff_fuel,
        "stiff_film": stiff_film,
        "cstar": cstar,
        "temp": temp,
    }

    # Warnings
    warn = []
    if res["stiff_ox"] < MIN_STIFF:
        warn.append(f"LOX stiff={res['stiff_ox']:.2f}<{MIN_STIFF}")
    if res["stiff_fuel"] < MIN_STIFF:
        warn.append(f"IPA stiff={res['stiff_fuel']:.2f}<{MIN_STIFF}")
    res["warnings"] = "; ".join(warn)
    return res

# --------------------------------------------------------------------------
# CSV I/O
# --------------------------------------------------------------------------

IN_COLS = ["mdot_total_kg_s", "OF_core", "insert_id", "P_N2_purge_psi", "P_c_torch_psi", "OF_torch"]

OUT_COLS = [
    "P_IPA_psi", "P_LOx_psi", "P_N2_purge_psi", "P_CH4_psi", "P_GOx_psi",
    "P_c_psi", "temp_K", "OF_core", "OF_overall", "mdot_total_kg_s", 
    "mdot_torch_kg_s", "P_c_torch_psi", "OF_torch",
    "insert_id", "cooling_pct",
    "mdot_IPA_kg_s", "mdot_LOx_kg_s", "mdot_film_kg_s",
    "mdot_CH4_kg_s", "mdot_GOx_kg_s",
    "stiff_fuel", "stiff_ox", "stiff_film", "cstar_m_s", "warnings"
]

def resolve_script_path(raw_path: str, script_dir: Path) -> str:
    path = Path(raw_path)
    if path.is_absolute():
        return str(path)
    return str(script_dir / path)

def read_cases(path: str):
    with open(path, newline="") as f:
        reader = csv.DictReader(f)
        missing = [c for c in IN_COLS if c not in (reader.fieldnames or [])]
        if missing:
            sys.exit(f"Input CSV is missing column(s): {missing}. Expected: {IN_COLS}")
        for i, row in enumerate(reader, start=2):
            try:
                yield (Q_(row["mdot_total_kg_s"], 'kg/s'), Q_(row["OF_core"], 'dimensionless'), row["insert_id"].strip(), 
                       Q_(row["P_N2_purge_psi"], 'psi'), Q_(row["P_c_torch_psi"], 'psi'), Q_(row["OF_torch"], 'dimensionless'))
            except ValueError:
                sys.exit(f"Bad number in input CSV on line {i}: {row}")

def to_row(mdot: Q_, of_core: Q_, insert_id: str, P_purge: Q_, P_c_torch: Q_, of_torch: Q_, r_GG: dict, r_torch: dict) -> dict:
    return {
        "P_c_psi": r_GG["P_c"].to('psi').m,
        "P_IPA_psi": r_GG["P_fuel"].to('psi').m,
        "P_LOx_psi": r_GG["P_ox"].to('psi').m,
        "P_N2_purge_psi": P_purge.to('psi').m,
        "mdot_total_kg_s": mdot.to('kg/s').m,
        "OF_core": of_core.to('dimensionless').m,
        "OF_overall": r_GG["of_overall"].to('dimensionless').m,
        "insert_id": insert_id,
        "cooling_pct": r_GG["cooling_pct"],
        "temp_K": r_GG["temp"].to('K').m,
        "cstar_m_s": r_GG["cstar"].to('m/s').m,
        "mdot_IPA_kg_s": r_GG["mdot_fuel"].to('kg/s').m,
        "mdot_LOx_kg_s": r_GG["mdot_ox"].to('kg/s').m,
        "mdot_film_kg_s": r_GG["mdot_film"].to('kg/s').m,
        "stiff_fuel": r_GG["stiff_fuel"],
        "stiff_ox": r_GG["stiff_ox"],
        "stiff_film": r_GG["stiff_film"],
        "P_c_torch_psi": P_c_torch.to('psi').m,
        "P_CH4_psi": r_torch["P_CH4"].to('psi').m,
        "P_GOx_psi": r_torch["P_GOx"].to('psi').m,
        "mdot_torch_kg_s": r_torch["mdot_torch"].to('kg/s').m,
        "OF_torch": of_torch.to('dimensionless').m,
        "mdot_CH4_kg_s": r_torch["mdot_ch4"].to('kg/s').m,
        "mdot_GOx_kg_s": r_torch["mdot_gox"].to('kg/s').m,
        "warnings": r_GG["warnings"],
    }

# --------------------------------------------------------------------------
# Main
# --------------------------------------------------------------------------

def main():
    script_dir = Path(__file__).resolve().parent
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--input", default="TestMatrix_Inputs.csv")
    ap.add_argument("--output", default="TestMatrix_Outputs.csv")
    args = ap.parse_args()

    args.input = resolve_script_path(args.input, script_dir)
    args.output = resolve_script_path(args.output, script_dir)

    rows = []
    for mdot, of_core, insert_id, P_purge, P_c_torch, of_torch in read_cases(args.input):
        try:
            r_GG = solve_GG(mdot, of_core, insert_id)
            r_torch = solve_torch(P_c_torch, of_torch)
            rows.append(to_row(mdot, of_core, insert_id, P_purge, P_c_torch, of_torch, r_GG, r_torch))
        except Exception as e:   # keep going so one bad case doesn't kill the batch
            row = {c: None for c in OUT_COLS}
            row.update(warnings=f"FAILED: {e}")
            rows.append(row)

    with open(args.output, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=OUT_COLS)
        w.writeheader()
        w.writerows(rows)

    print(f"Wrote {len(rows)} case(s) to {args.output}")

if __name__ == "__main__":
    main()
