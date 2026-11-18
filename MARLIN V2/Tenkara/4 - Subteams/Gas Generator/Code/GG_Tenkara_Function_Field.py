import numpy as np
from pint import Quantity as Q_
from pathlib import Path
from scipy.interpolate import RegularGridInterpolator
from GG_Tenkara_Function_single import solve_GasGen_state

# Global field parameters
P_MIN = Q_(100.0, 'psi')
P_MAX = Q_(1000.0, 'psi')
NUM_POINTS = 181
SPACING = (P_MAX - P_MIN).m / (NUM_POINTS - 1)
P_AXIS = np.linspace(P_MIN, P_MAX, NUM_POINTS)
VALID_INSERTS = ["0%", "5%", "10%", "15%", "20%", "25%", "30%", "35%", "40%"]

field_cache = {}

def check_insert(insert_id: str):
    """Check if insert is valid, raise error if not."""
    if insert_id not in VALID_INSERTS:
        raise ValueError(f"Invalid insert_id. Must be one of: {VALID_INSERTS}")

def field_path(insert_id: str) -> Path:
    """Return the path for saved npy field."""
    name = f"GasGenerator_Field_{insert_id.replace('%', 'pct')}.npy"
    return Path(__file__).with_name(name)

def get_field(insert_id: str) -> np.ndarray:
    """Return the field for an insert: cache -> npy -> generate and save."""

    if insert_id in field_cache:
        return field_cache[insert_id]
    
    path = field_path(insert_id)
    field = None
    if path.exists():
        field = np.load(path)

    if field is None:
        field = generate_field(insert_id)
        np.save(path, field)

    field_cache[insert_id] = field
    return field

def generate_field(insert_id: str) -> np.ndarray:
    """Builds NxN field for Mdot, Pc and OF for one insert.
    field[name, row, col]: name = Mdot, Pc, OF ; row = P_fu index ; col = P_ox index."""

    check_insert(insert_id)
    field = np.full((3, NUM_POINTS, NUM_POINTS), np.nan)

    for row, p_fu in enumerate(P_AXIS):
        print(f"Row {row+1}/{NUM_POINTS}")
        for col, p_ox in enumerate(P_AXIS):
            if col % 10 == 0:
                print(f"Column {col+1}/{NUM_POINTS}")
            state = solve_GasGen_state(Q_(p_ox, 'psi'), Q_(p_fu, 'psi'), insert_id)
            field[:, row, col] = (state["Mdot"].to('lb/s').m, state["Pc"].to('psi').m, state["OF"].to('dimensionless').m)

    return field

def interpolation_functions(field: np.ndarray) -> tuple:
    """Return interpolation functions for Mdot, Pc and OF for one insert."""
    interp_Mdot = RegularGridInterpolator((P_AXIS.m, P_AXIS.m), field[0], method='linear', bounds_error=False, fill_value=np.nan)
    interp_Pc = RegularGridInterpolator((P_AXIS.m, P_AXIS.m), field[1], method='linear', bounds_error=False, fill_value=np.nan)
    interp_OF = RegularGridInterpolator((P_AXIS.m, P_AXIS.m), field[2], method='linear', bounds_error=False, fill_value=np.nan)
    return interp_Mdot, interp_Pc, interp_OF

def get_GasGen_state(P_ox: Q_, P_fu: Q_, insert_id: str) -> dict[str, Q_]:
    """Look up total state for the given feed pressures and insert.
    Returns pint quantities for Mdot (lb/s), Pc (psi), and OF (dimensionless)."""

    check_insert(insert_id)
    if not (P_MIN <= P_ox <= P_MAX and P_MIN <= P_fu <= P_MAX):
        raise ValueError("Gas Generator feed pressure values are out of bounds.")
    
    field = get_field(insert_id)
    interp_Mdot, interp_Pc, interp_OF = interpolation_functions(field)
    Mdot = Q_(interp_Mdot([P_ox.m, P_fu.m]).item(), 'lb/s')
    Pc = Q_(interp_Pc([P_ox.m, P_fu.m]).item(), 'psi')
    OF = Q_(interp_OF([P_ox.m, P_fu.m]).item(), 'dimensionless')

    return {"Mdot": Mdot, "Pc": Pc, "OF": OF}

def get_GasGen_Mdot_total(P_ox: Q_, P_fu: Q_, insert_id: str) -> Q_:
    """Look up total mass flow rates for the given feed pressures and insert.
    Returns a pint Quantity in lb/s."""

    check_insert(insert_id)
    if not (P_MIN <= P_ox <= P_MAX and P_MIN <= P_fu <= P_MAX):
        raise ValueError("Gas Generator feed pressure values are out of bounds.")
    
    field = get_field(insert_id)
    interp_mdot = interpolation_functions(field)[0]

    return Q_(interp_mdot([P_ox.m, P_fu.m]).item(), 'lb/s')

def get_GasGen_Pc(P_ox: Q_, P_fu: Q_, insert_id: str) -> Q_:
    """Look up chamber pressure for the given feed pressures and insert.
    Returns a pint Quantity in psi."""

    check_insert(insert_id)
    if not (P_MIN <= P_ox <= P_MAX and P_MIN <= P_fu <= P_MAX):
        raise ValueError("Gas Generator feed pressure values are out of bounds.")
    
    field = get_field(insert_id)
    interp_Pc = interpolation_functions(field)[1]

    return Q_(interp_Pc([P_ox.m, P_fu.m]).item(), 'psi')

def get_GasGen_OF(P_ox: Q_, P_fu: Q_, insert_id: str) -> Q_:
    """Look up OF ratio for the given feed pressures and insert.
    Returns a pint Quantity in dimensionless."""

    check_insert(insert_id)
    if not (P_MIN <= P_ox <= P_MAX and P_MIN <= P_fu <= P_MAX):
        raise ValueError("Gas Generator feed pressure values are out of bounds.")
    
    field = get_field(insert_id)
    interp_OF = interpolation_functions(field)[2]

    return Q_(interp_OF([P_ox.m, P_fu.m]).item(), 'dimensionless')

def get_GasGen_Mdot_separate(P_ox: Q_, P_fu: Q_, insert_id: str) -> dict[str, Q_]:
    """Look up OX and FU mass flow rates for the given feed pressures and insert.
    Returns dict of pint Quantity in lb/s."""
    check_insert(insert_id)
    if not (P_MIN <= P_ox <= P_MAX and P_MIN <= P_fu <= P_MAX):
        raise ValueError("Gas Generator feed pressure values are out of bounds.")

    field = get_field(insert_id)
    interp_Mdot, _, interp_OF = interpolation_functions(field)

    total_mdot = Q_(interp_Mdot([P_ox.m, P_fu.m]).item(), 'lb/s')
    total_OF = Q_(interp_OF([P_ox.m, P_fu.m]).item(), 'dimensionless')
    mdot_ox = total_mdot * total_OF / (1 + total_OF)
    mdot_fu = total_mdot / (1 + total_OF)

    return {"OX Mdot": mdot_ox, "FU Mdot": mdot_fu}

def main():
    for insert in VALID_INSERTS:
        print(f"Generating field for insert: {insert}")
        out = get_GasGen_state(Q_(550, 'psi'), Q_(550, 'psi'), insert)
        print(f"Total Mdot: {out["Mdot"].to('lb/s'):.4f~P}")
        print(f"Chamber Pressure: {out["Pc"].to('psi'):.4f~P}")
        print(f"OF Ratio: {out["OF"]:.4f~P}")

if __name__ == "__main__":
    main()