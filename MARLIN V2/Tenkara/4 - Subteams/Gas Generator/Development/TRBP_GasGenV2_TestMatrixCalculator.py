"""
TRBP Gas Generator V2.0
Injector Feed-Pressure Solver

Inputs : mdot_total_kg_s, OF_core
Outputs: upstream LOX and IPA pressure, chamber pressure, overall OF, film cooling %, and diagnostics.

Model summary
-------------
Unknowns: mdot_ox, mdot_fuel, mdot_film, P_ox, P_fuel, P_c.
  - Orifice:   mdot = CdA * sqrt(2 * rho * dP)
  - Core OF:   mdot_ox / mdot_fuel = OF_core
  - Mass bal.: mdot_total = mdot_ox + mdot_fuel + mdot_film
  - c*:        P_c = c* * eta * mdot_total / A_t

Usage
-----
    python gg_injector_solver.py --config GG_hardware_definition.yaml \
                                 --input TestMatrix_Inputs.csv --output TestMatrix_Outputs.csv
"""

import argparse
import csv
import sys
from dataclasses import dataclass
from typing import Callable
import yaml

from tc_vs_cstar_eff import get_temp_at_eta

# --------------------------------------------------------------------------
# Settings
# --------------------------------------------------------------------------
# Warning thresholds
MIN_STIFF = 0.20     # feed-system decoupling rule of thumb
CEA_PC_GUESS_PA = 1.0e6   # starting guess for the P_c iteration
# Unit conversions
PA_PER_PSI = 6894.757
FT_TO_M = 0.3048

# --------------------------------------------------------------------------
# Hardware definition
# --------------------------------------------------------------------------
@dataclass
class Hardware:
    cda_ox: float                 # m^2
    cda_fuel: float               # m^2 (core fuel orifices)
    cda_film: float               # m^2 (total Cda of film orifices + insert)
    rho_ox: float                 # kg/m^3
    rho_fuel: float               # kg/m^3
    throat_cda: float             # m^2 (total CdA of the stator nozzles)
    cstar_eff: float              # dimensionless (0 < eff <= 1)
    cea_ox: str                   # RocketCEA oxidizer name
    cea_fuel: str                 # RocketCEA fuel name


def series_cda(*cdas: float) -> float:
    """Effective CdA of orifices in series (same mdot and rho): [sum 1/CdA^2]^-1/2.
    If any orifice is fully blocked (CdA = 0), nothing flows, so the result is 0."""
    if any(c <= 0.0 for c in cdas):
        return 0.0
    return sum(1.0 / c**2 for c in cdas) ** -0.5


def load_hardware(path: str, insert_id: str) -> Hardware:
    """Load the YAML and select the film insert named by insert_id."""
    with open(path) as f:
        cfg = yaml.safe_load(f)

    cda_film_insert = cfg["film_insert_Cda"][insert_id]
    cda_film_orifice = cfg["film_Cda"]
    cda_film = series_cda(cda_film_insert, cda_film_orifice)

    return Hardware(
        cda_ox=cfg["LOx_Cda"],
        cda_fuel=cfg["fuel_Cda"],
        cda_film=cda_film,
        rho_ox=cfg["LOx_rho"],
        rho_fuel=cfg["IPA_rho"],
        throat_cda=cfg["stator_throat_area"],
        cstar_eff=cfg["CSTAR_eff"],
        cea_ox=cfg["GG_ox_type"],
        cea_fuel=cfg["GG_fuel_type"],
    )


# --------------------------------------------------------------------------
# Physics
# --------------------------------------------------------------------------
def orifice_dp(mdot: float, cda: float, rho: float) -> float:
    """Pressure drop (Pa) needed to pass mdot through an orifice of effective area cda."""
    return (mdot / cda) ** 2 / (2.0 * rho)


def mass_flow_calculation(mdot_total: float, of_core: float, hw: Hardware) -> dict:
    """Closed-form mass flows. The film flow is set by geometry, not chosen."""
    r = hw.cda_film / hw.cda_fuel                  # mdot_film / mdot_fuel
    mdot_fuel = mdot_total / (of_core + 1.0 + r)
    return {
        "mdot_ox": of_core * mdot_fuel,
        "mdot_fuel": mdot_fuel,
        "mdot_film": r * mdot_fuel,
    }


def overall_of(mdot_ox: float, mdot_fuel: float, mdot_film: float) -> float:
    """Overall OF: all LOX over all fuel (core + film)."""
    return mdot_ox / (mdot_fuel + mdot_film)


def cooling_percent(mdot_film: float, mdot_total: float) -> float:
    """Film flow as a percent of total mass flow."""
    return 100.0 * mdot_film / mdot_total


def get_cstar(of: float, pc_pa: float, hw: Hardware) -> float:
    """Return c*(OF, P_c[Pa]) -> m/s (ideal, from NASA CEA via RocketCEA)."""
    from rocketcea.cea_obj import CEA_Obj
    cea = CEA_Obj(oxName=hw.cea_ox, fuelName=hw.cea_fuel)
    cstar = cea.get_Cstar(Pc=pc_pa / PA_PER_PSI, MR=of) * FT_TO_M
    return cstar


def chamber_pressure(mdot_total: float, of_overall: float, hw: Hardware,
                     tol: float = 1e-4, max_iter: int = 1000):
    """Fixed-point iteration on P_c = c*(OF, P_c) * eff * mdot / (Cd*A_t)."""
    pc = CEA_PC_GUESS_PA
    for _ in range(max_iter):
        cstar = get_cstar(of_overall, pc, hw) * hw.cstar_eff
        pc_new = cstar * mdot_total / hw.throat_cda
        if abs(pc_new - pc) / pc_new < tol:
            return pc_new, cstar
        pc = pc_new
    raise RuntimeError(f"P_c iteration did not converge (last P_c = {pc:.0f} Pa)")


def stiffness(dp: float, p_in: float) -> float:
    """Stiffness: dp / p_in."""
    return dp / p_in


def solve(mdot_total: float, of_core: float, hw: Hardware) -> dict:
    flows = mass_flow_calculation(mdot_total, of_core, hw)
    of_tot = overall_of(flows["mdot_ox"], flows["mdot_fuel"], flows["mdot_film"])
    cool_pct = cooling_percent(flows["mdot_film"], mdot_total)
    pc, cstar = chamber_pressure(mdot_total, of_tot, hw)

    dp_ox = orifice_dp(flows["mdot_ox"], hw.cda_ox, hw.rho_ox)
    dp_fuel = orifice_dp(flows["mdot_fuel"], hw.cda_fuel, hw.rho_fuel)
    dp_film = orifice_dp(flows["mdot_film"], hw.cda_film, hw.rho_fuel)

    p_ox = pc + dp_ox
    p_fuel = pc + dp_fuel    # inlet pressure feeding both the core and film paths

    stiff_ox = stiffness(dp_ox, p_ox)
    stiff_fuel = stiffness(dp_fuel, p_fuel)
    stiff_film = stiffness(dp_film, p_fuel)

    temp = get_temp_at_eta(hw.cstar_eff, of_core, pc / PA_PER_PSI)

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
        warn.append(f"LOX dP/Pc={res['stiff_ox']:.2f}<{MIN_STIFF}")
    if res["stiff_fuel"] < MIN_STIFF:
        warn.append(f"IPA dP/Pc={res['stiff_fuel']:.2f}<{MIN_STIFF}")
    res["warnings"] = "; ".join(warn)
    return res


# --------------------------------------------------------------------------
# CSV I/O
# --------------------------------------------------------------------------
IN_COLS = ["mdot_total_kg_s", "OF_core", "insert_id"]

OUT_COLS = [
    "mdot_total_kg_s", "OF_core", "insert_id",
    "P_c_psia", "P_ox_psia", "P_fuel_psia",
    "cooling_pct", "OF_overall",
    "mdot_ox_kg_s", "mdot_fuel_kg_s", "mdot_film_kg_s",
    "stiff_ox", "stiff_fuel", "stiff_film",
    "cstar_ms", "temp",
    "warnings"
]


def read_cases(path: str):
    with open(path, newline="") as f:
        reader = csv.DictReader(f)
        missing = [c for c in IN_COLS if c not in (reader.fieldnames or [])]
        if missing:
            sys.exit(f"Input CSV is missing column(s): {missing}. Expected: {IN_COLS}")
        for i, row in enumerate(reader, start=2):
            try:
                yield (float(row["mdot_total_kg_s"]), float(row["OF_core"]),
                       row["insert_id"].strip())
            except ValueError:
                sys.exit(f"Bad number in input CSV on line {i}: {row}")


def to_row(mdot: float, of_core: float, insert_id: str, r: dict) -> dict:
    psi = lambda pa: None if pa is None else pa / PA_PER_PSI
    return {
        "mdot_total_kg_s": mdot,
        "OF_core": of_core,
        "insert_id": insert_id,
        "P_ox_psia": psi(r["P_ox"]),
        "P_fuel_psia": psi(r["P_fuel"]),
        "P_c_psia": psi(r["P_c"]),
        "cooling_pct": r["cooling_pct"],
        "OF_overall": r["of_overall"],
        "mdot_ox_kg_s": r["mdot_ox"],
        "mdot_fuel_kg_s": r["mdot_fuel"],
        "mdot_film_kg_s": r["mdot_film"],
        "stiff_ox": r["stiff_ox"],
        "stiff_fuel": r["stiff_fuel"],
        "stiff_film": r["stiff_film"],
        "cstar_ms": r["cstar"],
        "temp": r["temp"],
        "warnings": r["warnings"],
    }


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--config", default="../Inputs/GG_hardware_definition.yaml")
    ap.add_argument("--input", default="TestMatrix_Inputs.csv")
    ap.add_argument("--output", default="TestMatrix_Outputs.csv")
    args = ap.parse_args()

    hw_cache = {}        # insert_id -> Hardware (YAML is read once per insert)
    cstar_fn = None      # built once; propellants don't depend on the insert

    rows = []
    for mdot, of_core, insert_id in read_cases(args.input):
        try:
            if insert_id not in hw_cache:
                hw_cache[insert_id] = load_hardware(args.config, insert_id)
            hw = hw_cache[insert_id]
            
            rows.append(to_row(mdot, of_core, insert_id,
                               solve(mdot, of_core, hw)))
        except Exception as e:   # keep going so one bad case doesn't kill the batch
            row = {c: None for c in OUT_COLS}
            row.update(mdot_total_kg_s=mdot, OF_core=of_core, insert_id=insert_id,
                       warnings=f"FAILED: {e}")
            rows.append(row)

    with open(args.output, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=OUT_COLS)
        w.writeheader()
        w.writerows(rows)

    print(f"Wrote {len(rows)} case(s) to {args.output}")


if __name__ == "__main__":
    main()
