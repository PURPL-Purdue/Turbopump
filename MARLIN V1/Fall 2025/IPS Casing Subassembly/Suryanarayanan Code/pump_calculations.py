import math

def pump_thrust_empricial(lbs_per_second, psi_delta, fluid_grams_per_ml, discharge_coefficient):
    kg_per_second = lbs_per_second * 0.453592
    newtons_per_meter_squared = psi_delta * 6894.76
    kg_per_meter_cubed = fluid_grams_per_ml * 1000  # g/mL * 1000 mL/L * 1 kg/1000 g * 1000 L/m^3

    newton_thrust = pump_thrust_metric(kg_per_second, newtons_per_meter_squared, kg_per_meter_cubed, discharge_coefficient)
    return newton_thrust * 0.224809

def pump_thrust_metric(mass_rate, delta_pressure, fluid_density, discharge_coefficient):
    return mass_rate * math.sqrt(delta_pressure/fluid_density) * (math.sqrt(2) + 1/discharge_coefficient/math.sqrt(2))

def main():
    flow_input = float(input("Enter the mass flow (lb/s): "))
    pressure_delta = float(input("Enter the pressure delta (PSI): "))
    fluid_density = float(input("Enter fluid density (g/mL): "))
    discharge_coefficient = float(input("Enter the discharge coefficient: "))
    thrust = pump_thrust_empricial(flow_input, pressure_delta, fluid_density, discharge_coefficient)
    print("Estimated pump thrust (lbf): "+str(thrust))

if __name__ == '__main__':
    main()