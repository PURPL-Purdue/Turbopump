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
  - c*:        P_c = c* * eff * mdot_total / A_t

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

# --------------------------------------------------------------------------
# Settings
# --------------------------------------------------------------------------
# Warning thresholds
MIN_DP_OVER_PC = 0.20     # feed-system decoupling rule of thumb
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
    cda_insert: float             # m^2 (selected film orifice insert; 0 = solid)
    cda_film: float               # m^2 (film orifices downstream of the manifold)
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

    # Insert IDs are compared as strings so 1, "1" and " 1 " all match
    inserts = {str(k).strip(): v for k, v in cfg["film_insert_Cda"].items()}
    key = str(insert_id).strip()
    if key not in inserts:
        raise ValueError(f"insert_id '{key}' not found in film_insert_Cda "
                         f"(available: {', '.join(inserts)})")
    if inserts[key] is None:
        raise ValueError(f"film_insert_Cda for insert '{key}' is empty in the YAML")

    return Hardware(
        cda_ox=cfg["LOx_Cda"],
        cda_fuel=cfg["fuel_Cda"],
        cda_insert=inserts[key],
        cda_film=cfg["film_Cda"],
        rho_ox=cfg["LOx_rho"],
        rho_fuel=cfg["IPA_rho"],
        throat_cda=cfg["stator_throat_area"],   # taken as the total nozzle CdA
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
    cda_film_path = series_cda(hw.cda_insert, hw.cda_film)
    r = cda_film_path / hw.cda_fuel                  # mdot_film / mdot_fuel
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


def make_cstar_fn(hw: Hardware) -> Callable[[float, float], float]:
    """Return c*(OF, P_c[Pa]) -> m/s (ideal, from NASA CEA via RocketCEA)."""
    try:
        from rocketcea.cea_obj import CEA_Obj
    except ImportError:
        sys.exit("rocketcea is not installed. Run: pip install rocketcea")

    cea = CEA_Obj(oxName=hw.cea_ox, fuelName=hw.cea_fuel)

    def cstar(of: float, pc_pa: float) -> float:
        # RocketCEA takes Pc in psia and returns c* in ft/s
        return cea.get_Cstar(Pc=pc_pa / PA_PER_PSI, MR=of) * FT_TO_M

    return cstar


def chamber_pressure(mdot_total: float, of_overall: float, hw: Hardware,
                     cstar_fn: Callable[[float, float], float],
                     tol: float = 1e-4, max_iter: int = 100):
    """Fixed-point iteration on P_c = c*(OF, P_c) * eff * mdot / (Cd*A_t)."""
    pc = CEA_PC_GUESS_PA
    for _ in range(max_iter):
        cstar = cstar_fn(of_overall, pc) * hw.cstar_eff
        pc_new = cstar * mdot_total / hw.throat_cda
        if abs(pc_new - pc) / pc_new < tol:
            return pc_new, cstar
        pc = pc_new
    raise RuntimeError(f"P_c iteration did not converge (last P_c = {pc:.0f} Pa)")


def solve(mdot_total: float, of_core: float, hw: Hardware,
          cstar_fn: Callable[[float, float], float]) -> dict:
    flows = mass_flow_calculation(mdot_total, of_core, hw)
    of_tot = overall_of(flows["mdot_ox"], flows["mdot_fuel"], flows["mdot_film"])
    cool_pct = cooling_percent(flows["mdot_film"], mdot_total)
    pc, cstar = chamber_pressure(mdot_total, of_tot, hw, cstar_fn)

    dp_ox = orifice_dp(flows["mdot_ox"], hw.cda_ox, hw.rho_ox)
    dp_fuel = orifice_dp(flows["mdot_fuel"], hw.cda_fuel, hw.rho_fuel)
    p_ox = pc + dp_ox
    p_fuel = pc + dp_fuel    # inlet pressure feeding both the core and film paths

    # Film path internals: inlet -> insert -> manifold -> film orifices -> chamber
    # (with a solid insert there is no flow: manifold sits at P_c and the whole
    #  P_fuel - P_c drop is across the insert)
    dp_film_orifice = orifice_dp(flows["mdot_film"], hw.cda_film, hw.rho_fuel)
    p_man = pc + dp_film_orifice
    dp_insert = p_fuel - p_man

    res = {
        **flows,
        "cooling_pct": cool_pct,
        "of_overall": of_tot,
        "P_c": pc,
        "P_ox": p_ox,
        "P_fuel": p_fuel,
        "cstar": cstar,
        "dp_ox_over_pc": dp_ox / pc,
        "dp_fuel_over_pc": dp_fuel / pc,
        "P_man": p_man,
        "dp_insert": dp_insert,
        "dp_film_orifice": dp_film_orifice,
    }

    # Warnings
    warn = []
    if res["dp_ox_over_pc"] < MIN_DP_OVER_PC:
        warn.append(f"LOX dP/Pc={res['dp_ox_over_pc']:.2f}<{MIN_DP_OVER_PC}")
    if res["dp_fuel_over_pc"] < MIN_DP_OVER_PC:
        warn.append(f"IPA dP/Pc={res['dp_fuel_over_pc']:.2f}<{MIN_DP_OVER_PC}")
    res["warnings"] = "; ".join(warn)
    return res


# --------------------------------------------------------------------------
# CSV I/O
# --------------------------------------------------------------------------
IN_COLS = ["mdot_total_kg_s", "OF_core", "insert_id"]

OUT_COLS = [
    "mdot_total_kg_s", "OF_core", "insert_id",
    "P_ox_psia", "P_fuel_psia", "P_c_psia",
    "cooling_pct", "OF_overall",
    "mdot_ox_kg_s", "mdot_fuel_kg_s", "mdot_film_kg_s",
    "P_ox_Pa", "P_fuel_Pa", "P_c_Pa",
    "cstar_eff_ms", "dP_ox_over_Pc", "dP_fuel_over_Pc",
    "P_man_psia", "dP_insert_psi",
    "warnings",
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
        "P_ox_Pa": r["P_ox"],
        "P_fuel_Pa": r["P_fuel"],
        "P_c_Pa": r["P_c"],
        "cstar_eff_ms": r["cstar"],
        "dP_ox_over_Pc": r["dp_ox_over_pc"],
        "dP_fuel_over_Pc": r["dp_fuel_over_pc"],
        "P_man_psia": psi(r["P_man"]),
        "dP_insert_psi": psi(r["dp_insert"]),
        "warnings": r["warnings"],
    }


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--config", default="GG_hardware_definition.yaml")
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
            if cstar_fn is None:
                cstar_fn = make_cstar_fn(hw)
            rows.append(to_row(mdot, of_core, insert_id,
                               solve(mdot, of_core, hw, cstar_fn)))
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
