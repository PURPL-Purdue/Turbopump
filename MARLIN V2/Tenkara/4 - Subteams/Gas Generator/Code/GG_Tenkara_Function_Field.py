import numpy as np
from pint import Quantity as Q_
from pathlib import Path
from GG_Tenkara_Function_single import GasGen_mdot_single
from functools import lru_cache


field_min_pressure = Q_(170.0, 'psi')
field_max_pressure = Q_(600.0, 'psi')



def calculate_field(func, x_values, y_values, insert_id):
    """Return func(x, y) on a grid, with y indexing rows and x indexing columns."""
    field = np.full((y_values.size, x_values.size), np.nan)
    for row, y in enumerate(y_values):
        for col, x in enumerate(x_values):
            result = func(Q_(float(x), 'psi'), Q_(float(y), 'psi'), insert_id)
            val = result.to('lb/s').m if isinstance(result, Q_) else float(result)
            field[row, col] = val
    return field


@lru_cache(maxsize=2)
def _calculate_cached(min_psi: Q_, max_psi: Q_, insert_id: str):
    points = (max_psi.to("psi").m - min_psi.to("psi").m) / 43  # Adjust the step size as needed
    xs = np.linspace(min_psi.to("psi").m, max_psi.to("psi").m, int(points))
    ys = np.linspace(min_psi.to("psi").m, max_psi.to("psi").m, int(points))
    return calculate_field(GasGen_mdot_single, xs, ys, insert_id)


def get_gas_generator_field(min_psi: Q_ = field_min_pressure, max_psi: Q_ = field_max_pressure, insert_id: str = "0%"):
    return _calculate_cached(min_psi, max_psi, insert_id).copy()

def get_GasGen_mdot(P_ox: Q_, P_fuel: Q_, insert_id: str = "0%"):
    """This is the Gas Generator mass flow rate Tenkara function. Input 2 pint objects for oxygen and fuel pressures plus the desired film cooling insert id"""
    if P_ox < field_min_pressure or P_ox > field_max_pressure or P_fuel < field_min_pressure or P_fuel > field_max_pressure:
        raise ValueError("Gas Generator Feed pressure values are out of bounds")
    if insert_id not in ["0%", "5%", "10%", "15%", "20%", "25%", "30%", "35%", "40%"]:
        raise ValueError("Invalid Gas Generator insert_id. Must be one of: '0%', '5%', '10%', '15%', '20%', '25%', '30%', '35%', '40%'")
    
    field = get_gas_generator_field(insert_id=insert_id)
    csv_path = Path(__file__).with_name(
        f"Gas_Generator_Field_{insert_id.replace('%', 'pct')}.csv"
    )
    np.savetxt(csv_path, field, delimiter=",")

    col = int(round(P_ox.to("psi").m - field_min_pressure.to("psi").m)/43)
    row = int(round(P_fuel.to("psi").m - field_min_pressure.to("psi").m)/43)

    return field[row, col]


def __main__():
    pin = Q_(510, "psi")
    ins = '10%'
    #Test Case
    print(_calculate_cached.cache_info())
    out = get_GasGen_mdot(pin, pin, ins)
    out = get_GasGen_mdot(pin, pin, ins)
    print(_calculate_cached.cache_info())
    print(f"Result: {round(Q_(out, 'lb/s'), 2)}")

if __name__ == "__main__":
	__main__()