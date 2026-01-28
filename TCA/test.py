import CoolProp.CoolProp as CP

fluid = "CO2"
T = 1000.0   # K
P = 101325   # Pa

# Absolute enthalpy [J/kg]
h_T = CP.PropsSI("H", "T", T, "P", P, fluid)
h_ref = CP.PropsSI("H", "T", 298.15, "P", P, fluid)  # reference T = 298.15 K

# Sensible enthalpy relative to 298 K
h_sensible = h_T - h_ref

print(f"Sensible enthalpy of CO2 relative to 298 K: {h_sensible:.2f} J/kg")

