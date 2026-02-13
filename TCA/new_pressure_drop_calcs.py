import numpy as np

inches_into_meters = 0.0254 
psi_to_pa = 6894.76
lbm_to_kg = 0.453592
degrees_into_rad = np.pi / 180

Cd = 0.8

rho_lox = 1141 #[kg/m^3]
rho_rp1 = 810 #[kg/m^3]

pressure_out = 500 #[psi]
mass_flow_RP1 = 2.909/80 #[kg/s]
mass_flow_LOX = 6.110/80 #[kg/s]
diameter_RP1 = 0.0430 #0.0430 #[in]
diameter_LOX = 0.0595 #0.0595 #[in]

diameter_LOX_m = diameter_LOX * inches_into_meters
diameter_RP1_m = diameter_RP1 * inches_into_meters

area_LOX = np.pi*(diameter_LOX_m/2)**2
area_RP1 = np.pi*(diameter_RP1_m/2)**2

velocity_LOX = mass_flow_LOX / (rho_lox * area_LOX)
velocity_RP1 = mass_flow_RP1 / (rho_rp1 * area_RP1)

pressure_drop_LOX = ((velocity_LOX**2 * rho_lox) / (Cd**2 * 2)) / psi_to_pa
pressure_drop_RP1 = ((velocity_RP1**2 * rho_rp1) / (Cd**2 * 2)) / psi_to_pa
pressure_inlet_LOX = pressure_out + pressure_drop_LOX
pressure_inlet_RP1 = pressure_out + pressure_drop_RP1

print(velocity_LOX)
print(velocity_RP1)

print("\n=== Injector Inlet Conditions ===")
print(f"LOX Inlet Pressure: {pressure_inlet_LOX:.2f} psi")
print(f"RP-1 Inlet Pressure: {pressure_inlet_RP1:.2f} psi")

def find_optimal_theta_ox(theta_fuel_deg, fuel_mass_flow, ox_mass_flow, fuel_velocity, ox_velocity):
    theta_f_rad = theta_fuel_deg*degrees_into_rad
    #transverse component factor for fuel (toward center)
    sin_theta_f = np.sin(theta_f_rad)
    required_sin_theta_ox = (fuel_mass_flow * fuel_velocity * sin_theta_f) / (ox_mass_flow * ox_velocity)
    theta_ox_rad = np.arcsin(required_sin_theta_ox)
    theta_ox_deg = theta_ox_rad/degrees_into_rad
    return theta_ox_deg

#computes the LOX injection angle rquired so that the resultant 
#post-impingement sheet angle matches the desired beta.
print("\n=== Injection Angle Results ===")
#user provides fuel angles and target sheet angle beta.
theta_rp1_deg = float(input("Fuel angle (deg): "))
theta_lox_deg = find_optimal_theta_ox(theta_rp1_deg, mass_flow_RP1, mass_flow_LOX, velocity_RP1, velocity_LOX)
print(f"Optimal LOX injection angle (deg): {theta_lox_deg:.3f}")