import numpy as np
from pint import Quantity as Q_
from pathlib import Path
from scipy.interpolate import RegularGridInterpolator
from GG_Tenkara_Function_single import GasGen_state_single

# Global field parameters
P_MIN = Q_(100.0, 'psi')
P_MAX = Q_(1000.0, 'psi')
NUM_POINTS = 91
SPACING = (P_MAX - P_MIN).m / (NUM_POINTS - 1)
P_AXIS = np.linspace(P_MIN, P_MAX, NUM_POINTS)

VALID_INSERTS = ["0%", "5%", "10%", "15%", "20%", "25%", "30%", "35%", "40%"]

field_cache = {}

def check_insert(insert_id: str):
    if insert_id not in VALID_INSERTS:
        raise ValueError(f"Invalid insert_id. Must be one of: {VALID_INSERTS}")

def field_path(insert_id: str) -> Path:
    name = f"GasGenerator_Field_{insert_id.replace('%', 'pct')}.npy"
    return Path(__file__).with_name(name)

def get_field(insert_id: str) -> np.ndarray:
    """Return the field for an insert: memory cache -> npy -> generate and save."""
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
    field[name, row, col]: name = OX or FU Mdot, row = fuel pressure index, col = ox pressure index."""
    check_insert(insert_id)
    field = np.full((2, NUM_POINTS, NUM_POINTS), np.nan)
    for row, p_fu in enumerate(P_AXIS):
        print(f"Row {row+1}/{NUM_POINTS}")
        for col, p_ox in enumerate(P_AXIS):
            if col % 10 == 0:
                print(f"Column {col+1}/{NUM_POINTS}")
            state = GasGen_state_single(Q_(p_ox, 'psi'), Q_(p_fu, 'psi'), insert_id)
            field[:, row, col] = (state["OX Mdot"].to('lb/s').m, state["FU Mdot"].to('lb/s').m)
    return field

def get_GasGen_mdot(P_ox: Q_, P_fu: Q_, insert_id: str) -> dict[str, Q_]:
    """Look up mass flow rates for the given pressures and insert.
    Returns a pint Quantity in lb/s."""
    check_insert(insert_id)
    if not (P_MIN <= P_ox <= P_MAX and P_MIN <= P_fu <= P_MAX):
        raise ValueError("Gas Generator feed pressure values are out of bounds.")

    field = get_field(insert_id)
    interp_ox_mdot = RegularGridInterpolator((P_AXIS.m, P_AXIS.m), field[0], method='linear', bounds_error=False, fill_value=np.nan)
    interp_fu_mdot = RegularGridInterpolator((P_AXIS.m, P_AXIS.m), field[1], method='linear', bounds_error=False, fill_value=np.nan)

    mdot_ox = Q_(interp_ox_mdot([P_ox.m, P_fu.m]).item(), 'lb/s')
    mdot_fu = Q_(interp_fu_mdot([P_ox.m, P_fu.m]).item(), 'lb/s')

    return {"OX Mdot": mdot_ox, "FU Mdot": mdot_fu}

def main():
    out = get_GasGen_mdot(Q_(557, 'psi'), Q_(575, 'psi'), '10%')
    # out2 = GasGen_state_single(Q_(557, 'psi'), Q_(575, 'psi'), '10%')
    print(f"OX Mdot: {out["OX Mdot"].to('lb/s'):.4f~P}")
    print(f"FU Mdot: {out["FU Mdot"].to('lb/s'):.4f~P}")

if __name__ == "__main__":
    main()