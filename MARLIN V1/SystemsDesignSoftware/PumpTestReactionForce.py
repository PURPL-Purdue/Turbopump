import numpy as np
import matplotlib.pyplot as mpl
from CoolProp.CoolProp import PropsSI
import CoolProp
import yaml


## Conversion factors ##
psi2Pa = 6894.76
Ntolbf = 0.224809
lbtokg = 0.453592
kgpm3tolbpft3 = 0.062428
ft3togal = 7.48052
R2K = 0.555556
fttom = 0.3048
ft2tom2 = 0.092903
ft3tom3 = 0.0283168


## YAML import ##
with open('pumps_params.yaml') as file:
    pumps_params = yaml.safe_load(file)


## Var initialize ##
lox_temp = pumps_params['lox_inlet_design_temperature'] #(R)
kero_temp = pumps_params['kero_inlet_design_temperature'] #(R)

lox_inlet_press = pumps_params['lox_inlet_design_pressure'] #(psi)
kero_inlet_press = pumps_params['kero_inlet_design_pressure'] #(psi)

lox_outlet_press = pumps_params['lox_outlet_design_pressure'] #(psi)
kero_outlet_press = pumps_params['kero_outlet_design_pressure'] #(psi)

lox_mdot = pumps_params['lox_design_mdot'] #(lbm/s)
kero_mdot = pumps_params['kero_design_mdot'] #(lbm/s)

lox_vapor_P = PropsSI('P','T', lox_temp*R2K,'Q',0,'oxygen')/psi2Pa #(psi)
kero_vapor_P = 0.031 #(psi)

lox_rho = kgpm3tolbpft3 * PropsSI('D','T', lox_temp*R2K,'P', 14.7*psi2Pa, 'oxygen') #(lbm/ft^3)
kero_rho = 50.3 #(lbm/ft^3)

kero_vdot = kero_mdot / (kero_rho) #(ft^3/s)
lox_vdot = lox_mdot / (lox_rho) #(ft^3/s)

## Generate equivalent fluid values ##

# Water for kero
water_rho = kgpm3tolbpft3 * PropsSI('D','T', kero_temp*R2K,'P', 14.7*psi2Pa, 'water') #(lbm/ft^3)
water_vdot = kero_vdot
water_mdot = water_vdot * water_rho #(lbm/s)

water_headrise = 144 * (kero_outlet_press - kero_inlet_press) / kero_rho #(ft)
water_outlet_head = 144 * (kero_inlet_press - kero_vapor_P) / kero_rho + water_headrise #(ft)
water_outlet_press = water_outlet_head * water_rho / 144 #(psi)

# N2 for LOX
N2_temp = (PropsSI('T','P', 14.7*psi2Pa,'Q',0,'Nitrogen') - 10) / R2K #(R) Enforcing 10K subcool to make it liquid
N2_rho = kgpm3tolbpft3 * PropsSI('D','T', N2_temp*R2K,'P', 14.7*psi2Pa, 'Nitrogen') #(lbm/ft^3)
N2_vdot = kero_vdot
N2_mdot = N2_vdot * N2_rho #(lbm/s)

N2_headrise = 144 * (lox_outlet_press - lox_inlet_press) / lox_rho
N2_outlet_head = 144 * (lox_inlet_press - lox_vapor_P) / lox_rho + N2_headrise #(ft)
N2_outlet_press = N2_outlet_head * N2_rho / 144 #(psi)


## Orifice sizing ##
# Sorry, the unit conversions get pretty nasty past this point...

C_d = 0.8 #0.8 is pretty valid for most thick plate orifices

# Kero and kero simulator (water)
kero_orifice_A = kero_vdot * ft3tom3 / (C_d * (2 * (kero_outlet_press - 14.7) * psi2Pa / (kero_rho / kgpm3tolbpft3))**0.5) #(m^3)
water_orifice_A = water_vdot * ft3tom3 / (C_d * (2 * (water_outlet_press - 14.7) * psi2Pa / (water_rho / kgpm3tolbpft3))**0.5) #(m^3)

kero_velocity = kero_vdot * ft3tom3 / kero_orifice_A #(m/s)
water_velocity = water_vdot * ft3tom3 / water_orifice_A #(m/s)

kero_orifice_diam = 12 / fttom * 2 * (kero_orifice_A / np.pi) ** 0.5 #(in)
water_orifice_diam = 12 / fttom * 2 * (water_orifice_A / np.pi) ** 0.5 #(in)

# LOX and LOX simulator (LN2)
lox_orifice_A = lox_vdot * ft3tom3 / (C_d * (2 * (lox_outlet_press - 14.7) * psi2Pa / (lox_rho / kgpm3tolbpft3))**0.5) #(m^3)
N2_orifice_A = N2_vdot * ft3tom3 / (C_d * (2 * (N2_outlet_press - 14.7) * psi2Pa / (N2_rho / kgpm3tolbpft3))**0.5) #(m^3)

lox_velocity = lox_vdot * ft3tom3 / lox_orifice_A #(m/s)
N2_velocity = N2_vdot * ft3tom3 / N2_orifice_A #(m/s)

lox_orifice_diam = 12 / fttom * 2 * (lox_orifice_A / np.pi) ** 0.5 #(in)
N2_orifice_diam = 12 / fttom * 2 * (N2_orifice_A / np.pi) ** 0.5 #(in)


## Reaction forces ##

# Kero and kero simulator (water)
kero_force = kero_mdot * lbtokg * (kero_velocity) * Ntolbf #(lbf)
water_force = water_mdot * lbtokg * (water_velocity) * Ntolbf #(lbf)

# LOX and LOX simulator (LN2)
lox_force = lox_mdot * lbtokg * (lox_velocity) * Ntolbf #(lbf)
N2_force = N2_mdot * lbtokg * (N2_velocity) * Ntolbf #(lbf)


## Outputs ##
table = [
    [kero_outlet_press, kero_velocity / fttom, kero_force, 0, 0, kero_orifice_diam],
    [water_outlet_press, water_velocity / fttom, water_force, 0, 0, water_orifice_diam],
    [lox_outlet_press, lox_velocity / fttom, lox_force, 0, 0, lox_orifice_diam],
    [N2_outlet_press, N2_velocity / fttom, N2_force, 0, 0, N2_orifice_diam]
]

# Column/row names
headers = ["Commodity", "Press (psi)", "Velocity (ft/s)", "Force (lbf)", "Extra column", "Extra column", "Orifice diam (in)"]
row_labels = ["Kero", "Water", "LOX", "LN2"]

# Build string versions with correct decimal formatting first
rows_str = []
for label, row in zip(row_labels, table):
    row_str = [label] \
              + [f"{val:.2f}" for val in row[:-1]] \
              + [f"{row[-1]:.4f}"]
    rows_str.append(row_str)

# Determine column widths (max of header or any row value)
col_widths = []
for col in range(len(headers)):
    max_len = max(len(headers[col]), max(len(r[col]) for r in rows_str))
    col_widths.append(max_len + 2)  # add padding

# Print header (left-align headers, right-align numbers)
header_fmt = "".join(
    f"{{:<{col_widths[0]}}}" if i == 0 else f"{{:>{w}}}" 
    for i, w in enumerate(col_widths)
)
row_fmt = "".join([f"{{:<{col_widths[0]}}}"] + [f"{{:>{w}}}" for w in col_widths[1:]])

print(header_fmt.format(*headers))
print("-" * sum(col_widths))

for r in rows_str:
    print(row_fmt.format(*r))