import numpy as np
from pint import Quantity as Q_
from pathlib import Path
from GG_Tenkara_Function_single import GasGen_mdot_single

# ---- Global field parameters ----
FIELD_MIN_PRESSURE = Q_(0.0, 'psi')
FIELD_MAX_PRESSURE = Q_(800.0, 'psi')
NUM_POINTS = 100

VALID_INSERTS = ["0%", "5%", "10%", "15%", "20%", "25%", "30%", "35%", "40%"]

_P_MIN = FIELD_MIN_PRESSURE.to('psi').m
_P_MAX = FIELD_MAX_PRESSURE.to('psi').m
_SPACING = (_P_MAX - _P_MIN) / (NUM_POINTS - 1)

_field_cache = {}  # insert_id -> field, so each CSV is read at most once per run

def _check_insert(insert_id: str):
    if insert_id not in VALID_INSERTS:
        raise ValueError(f"Invalid insert_id. Must be one of: {VALID_INSERTS}")

def _csv_path(insert_id: str) -> Path:
    name = f"Gas_Generator_Field_{insert_id.replace('%', 'pct')}.csv"
    return Path(__file__).with_name(name)

def _get_field(insert_id: str) -> np.ndarray:
    """Return the field for an insert: memory cache -> CSV -> generate and save."""
    if insert_id in _field_cache:
        return _field_cache[insert_id]

    path = _csv_path(insert_id)
    field = None
    if path.exists():
        field = np.loadtxt(path, delimiter=",", ndmin=2)
        if field.shape != (NUM_POINTS, NUM_POINTS):  # stale file, rebuild
            field = None
    if field is None:
        field = generate_GasGen_field(insert_id)
        np.savetxt(path, field, delimiter=",")

    _field_cache[insert_id] = field
    return field

def _snap_index(p: Q_) -> int:
    idx = round((p.to('psi').m - _P_MIN) / _SPACING)
    return max(0, min(NUM_POINTS - 1, idx))

def generate_GasGen_field(insert_id: str) -> np.ndarray:
    """Build the NUM_POINTS x NUM_POINTS mdot field (lb/s) for one insert.
    field[row, col]: row = fuel pressure index, col = ox pressure index."""
    _check_insert(insert_id)
    p_axis = np.linspace(_P_MIN, _P_MAX, NUM_POINTS)  # same axis for ox and fuel
    field = np.full((NUM_POINTS, NUM_POINTS), np.nan)
    for row, p_fu in enumerate(p_axis):
        print(f"{row}")
        for col, p_ox in enumerate(p_axis):
            result = GasGen_mdot_single(Q_(float(p_ox), 'psi'), Q_(float(p_fu), 'psi'), insert_id)
            field[row, col] = result.to('lb/s').m if isinstance(result, Q_) else float(result)
    return field

def get_GasGen_mdot(P_ox: Q_, P_fu: Q_, insert_id: str) -> Q_:
    """Look up mdot for the given pressures and insert by snapping to the nearest grid point.
    Returns a pint Quantity in lb/s."""
    _check_insert(insert_id)
    if not (FIELD_MIN_PRESSURE <= P_ox <= FIELD_MAX_PRESSURE
            and FIELD_MIN_PRESSURE <= P_fu <= FIELD_MAX_PRESSURE):
        raise ValueError("Gas Generator feed pressure values are out of bounds")

    field = _get_field(insert_id)
    return Q_(field[_snap_index(P_fu), _snap_index(P_ox)], 'lb/s')

def main():
    out = get_GasGen_mdot(Q_(557, 'psi'), Q_(575, 'psi'), '10%')
    print(f"Result: {out.to('kg/s'):.4f~P}")

if __name__ == "__main__":
    main()