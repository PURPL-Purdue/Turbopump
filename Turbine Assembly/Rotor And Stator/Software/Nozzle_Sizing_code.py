# Hello everybody welcome to my super awesome code

import math

def main():

# Conversions
    lb_to_kg = 0.45359237
    bar_to_pa = 100000
    psi_to_pa = 6894.75729


# Nozzle Constants
    chamber_pressure = 350 * psi_to_pa # psi => N/m^2 (pa)
    exit_pressure = 30 * psi_to_pa # psi => N/m^2 (pa)
    mass_flow = 0.81 * lb_to_kg # lbm/s => kg/s
    inlet_temp = 876.05 # K
    gamma = 1.1201 # Specific heat ratio (unitless)
    gas_constant = 8.3145 # J/(mol*K)
    molar_mass = 11.328 # g/mol
    n = 8 # Number of nozzles (unitless)
    c_star = 1056.9 # Characteristic velocity (m/s)


# Isentropic Equations
    specific_gas_constant = (gas_constant / molar_mass) * 1000
        # J/(kg*K)

    area_throat = (c_star * mass_flow) / chamber_pressure
        # m^2

    inlet_density = chamber_pressure / (specific_gas_constant * inlet_temp)
        # kg/m^3

    exhaust_velocity = (inlet_temp * specific_gas_constant * ( (2 * gamma) / (gamma - 1) ) 
                        * (1 - (exit_pressure / chamber_pressure) ** ( (gamma - 1) / gamma ) ) ) ** 0.5
        # m/s

    temp_throat = inlet_temp / (1 + ( (gamma - 1) / 2 ) )
        # K

    pressure_throat = chamber_pressure / (1 + ( (gamma - 1) / 2 ) ) ** (gamma / (gamma - 1) )
        # pa

    density_throat = pressure_throat / (specific_gas_constant * temp_throat)
        # kg/m^3

    mach_exit = ( (2 * (math.e ** ( -( (math.log(exit_pressure / chamber_pressure) ) * (gamma - 1) ) / gamma) - 1) ) / (gamma - 1) ) ** 0.5
        # unitless

    density_exit = inlet_density / (1 + ( (gamma - 1) / 2 ) * (mach_exit ** 2) ) ** (1 / (gamma - 1))
        # kg/m^3

    area_exit = mass_flow / (density_exit * exhaust_velocity)
        # m^2

    temp_exit = inlet_temp / (1 + ( (gamma - 1) / 2 ) * (mach_exit ** 2) )
        # K

    thrust_force = mass_flow * exhaust_velocity
        # N


# Rao Calculations
    expansion_ratio = area_exit / area_throat
        # unitless

    radius_throat = (area_throat / math.pi) ** 0.5
        # m

    diameter_throat = 2 * radius_throat
        # m

    radius_exit = (expansion_ratio ** 0.5) * radius_throat
        # m

    diameter_exit = 2 * radius_exit
        # m

    angle_in_radians = math.radians(15)
        # unitless

    rao_length = ( ( (expansion_ratio ** 0.5) - 1) * radius_throat) / math.tan(angle_in_radians)
        # m

    nozzle_length = 0.8 * rao_length # length from throat to exit
        # m


# Sizing for Multiple nozzles
    area_throat_n = area_throat / n # This is the circular area
        # m^2

    radius_throat_n = (area_throat_n / math.pi) ** 0.5
        # m

    diameter_throat_n = radius_throat_n * 2 # Technically the length of one side of the throat square area
        # m

    square_area_throat = diameter_throat_n ** 2
        # m^2

    area_exit_n = area_exit / n # This is the circular area
        # m^2

    radius_exit_n = (area_exit_n / math.pi) ** 0.5
        # m

    diameter_exit_n = radius_exit_n * 2 # Technically the length of one side of the exit square area
        # m

    square_area_exit = diameter_exit_n ** 2
        # m^2

    nozzle_length_n = nozzle_length / n
        # m


# Prints (feel free to add more prints if you want to see other values)
    print(f"Mass Flow: 0.81 kg/s or {mass_flow:.4f} lbm/s")
    print(f"Exhaust Velocity: {exhaust_velocity:.4f} m/s")
    print(f"Exit Mach: {mach_exit:.4f}")
    print(f"Exit Pressure: 30 psi")
    print(f"Expansion Ratio: {expansion_ratio:.4f}")
    print(f"Square Throat Area: {square_area_throat:.4e} m^2")
    print(f"Square Exit Area: {square_area_exit:.4e} m^2")
    print(f"Number of Nozzles: {n}")
    print(f"Nozzle length from throat to exit: {nozzle_length_n:.4f} m")

# Inlet and exit angles (given from NASA SP-8120)
    print("Throat Angle: ~21 degrees")
    print("Exit Angle: ~15 degrees")

if __name__ == "__main__":
    main()