import numpy as np
import yaml
from CoolProp.CoolProp import PropsSI
import pint

with open(r'Inputs/TCA_params.yaml') as file:
    tca_params = yaml.safe_load(file)

ureg = pint.UnitRegistry()

mdot = tca_params['tp_mdot'] * ureg.kg/ureg.s  # Total propellant mass flow rate [kg/s]
OF_Ratio = tca_params['of_ratio']  # Mixture ratio O/F

rho_fuel = tca_params['densities']['fuel'] * ureg.kg/ureg.m**3  # Ethanol Density [kg/m^3] at injector conditions
rho_lox = tca_params['densities']['oxidizer'] * ureg.kg/ureg.m**3  # LOX Density [kg/m^3] at injector conditions

Pc = tca_params['chamber_pressure'] * ureg.psi  # Chamber pressure 

Cd_ox = tca_params['discharge_coeff']['oxidizer']  # Discharge coefficient for LOX
Cd_fuel = tca_params['discharge_coeff']['fuel']  # Discharge coefficient for fuel

Length_chamber = tca_params['chamber_length'] * ureg.m  # Combustion chamber length [m]
CombDiam = tca_params['chamber_diameter'] * ureg.inch  # chamber inner diameter [in]
throat_diameter = tca_params['throat_diameter'] * ureg.inch  # throat diameter [in]

fuel_diameter = tca_params['hole_diameters']['fuel'] * ureg.mm  # fuel orifice diameter [mm]
ox_diameter = tca_params['hole_diameters']['oxidizer'] * ureg.mm  # LOX orifice diameter [mm]

num_holes_fuel_inj = tca_params['hole_number']['fuel']  # number of fuel holes in injector faceplate
num_holes_ox_inj = tca_params['hole_number']['oxidizer']  # number of LOX holes in injector faceplate

# == Injector Geometry Parameters ==
wall_thickness_fuel = tca_params['wall_thickness']['fuel'] * ureg.mm  # [mm]
wall_thickness_ox = tca_params['wall_thickness']['oxidizer'] * ureg.mm  # [mm]

impinge_distance = tca_params['impinge_percentage'] * Length_chamber  # streams will impinge at the specified percentage of chamber length
ox_hole_distance = tca_params['ox_distance'] * ureg.mm  # distance between LOX holes to fuel hole

# Surface tension properties
sigma_fuel= PropsSI('SURFACE_TENSION', 'T', 298, 'Q', 0, 'Ethanol') * ureg.N/ureg.m  # Surface tension of Ethanol [N/m] at ~298K (literature value)
sigma_lox = PropsSI('SURFACE_TENSION', 'T', 90, 'Q', 0, 'Oxygen') * ureg.N/ureg.m  # Surface tension of LOX [N/m]

### Massflow Distribution ###

mdot_fuel = mdot / (1 + OF_Ratio)  # Fuel mass flow rate
mdot_ox = mdot * OF_Ratio / (1 + OF_Ratio)  # Oxidizer mass flow rate

def mdot_prop(mdot_propellant, number_holes):
    mdot_per_hole = mdot_propellant / number_holes
    return mdot_propellant, mdot_per_hole

mdot_fuel_inj, mdot_fuel_inj_per_hole = mdot_prop(mdot_fuel, num_holes_fuel_inj)
mdot_ox_inj, mdot_lox_inj_per_hole = mdot_prop(mdot_ox, num_holes_ox_inj)

print(f"#### Mass Flow Distribution ####")
print(f"Fuel total mass flow: {mdot_fuel_inj:.3f}")
print(f"Oxidizer total mass flow: {mdot_ox_inj:.3f}")
print(f"\nPer-hole mass flows:")
print(f"Fuel per-hole mass flow: {mdot_fuel_inj_per_hole:.5f}")
print(f"Oxidizer per-hole mass flow: {mdot_lox_inj_per_hole:.5f}")

### Orifice Sizing ###

# Calculate orifice areas
area_LOX = np.pi * (ox_diameter.to(ureg('m')) / 2)**2   # m^2
area_fuel = np.pi * (fuel_diameter.to(ureg('m'))/ 2)**2   # m^2

# Calculate jet velocities from continuity equation: v = mdot / (rho * A)
velocity_oxidizer = mdot_lox_inj_per_hole / (rho_lox * area_LOX)
velocity_fuel = mdot_fuel_inj_per_hole / (rho_fuel * area_fuel)

print(f"#### Orifice Sizing ####")
print(f"\n=== Jet Velocities ===")
print(f"LOX velocity: {velocity_oxidizer:.2f}")
print(f"Fuel velocity: {velocity_fuel:.2f}")
print(f"\n=== Orifice Diameters ===")
print(f"LOX orifice diameter: {ox_diameter:.4f} - ({ox_diameter.to('mm'):.3f})")
print(f"Fuel orifice diameter: {fuel_diameter:.4f} - ({fuel_diameter.to('mm'):.3f})")

### Pressure Drop Calculations ###

# Calculate pressure drop
pressure_drop_oxidizer = (velocity_oxidizer**2 * rho_lox) / (Cd_ox**2 * 2)
pressure_drop_fuel = (velocity_fuel**2 * rho_fuel) / (Cd_fuel**2 * 2)

# Calculate required inlet pressures
pressure_inlet_LOX = Pc + pressure_drop_oxidizer.to(ureg.psi)
pressure_inlet_fuel = Pc + pressure_drop_fuel.to(ureg.psi)

print(f"\n=== Pressure Drop Analysis ===")
print(f"LOX pressure drop: {pressure_drop_oxidizer.to(ureg.psi):.2f} - ({pressure_drop_oxidizer.to(ureg.kPa):.2f})")
print(f"Fuel pressure drop: {pressure_drop_fuel.to(ureg.psi):.2f} psi ({pressure_drop_fuel.to(ureg.kPa):.2f})")
print(f"\n=== Required Inlet Pressures ===")
print(f"LOX inlet pressure: {pressure_inlet_LOX:.2f}")
print(f"Fuel inlet pressure: {pressure_inlet_fuel:.2f}")
print(f"\nChamber pressure: {Pc:.2f}")
print(f"Pressure drop as % of chamber pressure:")
print(f"  LOX: {(pressure_drop_oxidizer.to(ureg.psi)/Pc).to('%'):.1f}")
print(f"  Fuel: {(pressure_drop_fuel.to(ureg.psi)/Pc).to('%'):.1f}")

### Angle of Impingement ###

def ox_impinge_angle(impinge_distance, ox_hole_distance):
    # Calculate angle of impingement for LOX jets
    theta = np.arctan(ox_hole_distance / impinge_distance).to(ureg.radian)
    return theta

def ox_offset(theta_ox, ox_hole_thickness):
    # Calculate offset of LOX jet from centerline at impingement point
    offset = np.tan(theta_ox) * ox_hole_thickness
    return offset

theta_ox = ox_impinge_angle(impinge_distance, ox_hole_distance)
ox_offset = ox_offset(theta_ox, wall_thickness_ox)

print(f"\n=== Angle of Impingement ===")
print(f"Ox angle of impingement: {theta_ox.to(ureg.degrees):.2f}")
print(f"Ox hole distance: {ox_hole_distance.to(ureg.mm):.2f}")
print(f"Impinge distance: {impinge_distance.to(ureg.mm):.2f}")


print(f"\n=== Ox Offset ===")
print(f"Ox offset: {ox_offset.to(ureg.mm):.2f}")

### Webber Number Calculation ###

v_ox_transverse = velocity_oxidizer * np.sin(theta_ox)
v_fuel_transverse = velocity_fuel 
    
We_lox = ((rho_lox.to(ureg.kg/ureg.m**3) * v_ox_transverse.to(ureg.m/ureg.s)**2 * ox_diameter.to(ureg.m)) / sigma_lox.to(ureg.N/ureg.m)).to_base_units()
We_fuel = ((rho_fuel.to(ureg.kg/ureg.m**3) * v_fuel_transverse.to(ureg.m/ureg.s)**2 * fuel_diameter.to(ureg.m)) / sigma_fuel.to(ureg.N/ureg.m)).to_base_units()
    
print(f"\n=== Weber Number Check ===")
print(f"Weber Number for LOX jets: {We_lox:.2f}")
print(f"Weber Number for fuel jets: {We_fuel:.2f}")
print(f"\nAtomization quality assessment:")
print(f"  We < 100: Poor atomization")
print(f"  100 < We < 1000: Good atomization")
print(f"  We > 1000: Excellent atomization")
    
if We_lox > 1000 and We_fuel > 1000:
    print(f"\n✓ Both propellants show excellent atomization characteristics")
elif We_lox > 100 and We_fuel > 100:
    print(f"\n✓ Both propellants show good atomization characteristics")
else:
    print(f"\n⚠ Warning: Atomization may be poor (We < 100)")


#Diameter ratio Max Mixing Condition

M = 1.6 # Mixing ratio for 2:1 impingement

ratio1 = (fuel_diameter.to(ureg.m) / ox_diameter.to(ureg.m))**2
ratio2 = M*((rho_lox/rho_fuel)*(mdot_fuel/mdot_ox)**2)**0.7

print(f"\n=== Diameter Ratio for Max Mixing Condition ===")
print(f"(Df/Dox)^2_Real = {ratio1:.2f}")
print(f"(Df/Dox)^2_MMC = {ratio2:.2f}")
print(f"(Df/Dox)^2_Real / (Df/Dox)^2_MMC = {(ratio1/ratio2):.2f}")

#Ratio of Diameters
print(f"\n=== Ratio of Diameters - Needs to be in the order of 1.22 for max mixing ===")
print(f"(Df/Dox) {(fuel_diameter.to(ureg.m) / ox_diameter.to(ureg.m)):.2f}")

