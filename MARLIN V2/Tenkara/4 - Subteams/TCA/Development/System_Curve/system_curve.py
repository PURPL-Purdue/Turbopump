import numpy as np
import yaml
import matplotlib.pyplot as plt

IN2M = 0.0254
PSI2PA = 6894.76

# Import TCA Parameters
with open(r'Inputs/TCA_params.yaml') as file:
    tca_params = yaml.safe_load(file)

# inputs
mdot_design = tca_params["tp_mdot"]           # nominal mass flow [kg/s]

rho_f = tca_params["densities"]["fuel"]        # [kg/m^3]
rho_ox = tca_params["densities"]["oxidizer"]   # [kg/m^3]

dt = tca_params["throat_diameter"]  * IN2M     # in -> [m]
At = np.pi/4 * dt**2                           # [m^2]

cstar = tca_params["cstar"]                    # [m/s]

# System curve driving function

# i want this flowrate at this OF
mdot = np.linspace(0,6,500)
OF = 1

# that requires this chamber pressure
Pc = mdot * cstar / At

# given the size of my injector holes and the mass flow through my orifices, this is the required upstream pressure
mdot_ox = mdot * (OF / (OF + 1))
mdot_f = mdot * (1 / (OF + 1))

ox_count = tca_params["hole_number"]["oxidizer"]
fuel_count = tca_params["hole_number"]["fuel"]

ox_element_mdot = mdot_ox / ox_count
fuel_element_mdot = mdot_f / fuel_count


d_ox = tca_params["hole_diameters"]["oxidizer"] / 1000 # mm -> [m]
d_f = tca_params["hole_diameters"]["fuel"] / 1000 # mm -> [m]

A_ox = np.pi / 4 * d_ox**2 # [m^2]
A_f = np.pi / 4 * d_f**2   # [m^2]

Cd_ox = 0.6
Cd_f = 0.6

dP_inj_ox =  (ox_element_mdot / (Cd_ox * A_ox))**2 / (2 * rho_ox)
dP_inj_f =  (fuel_element_mdot / (Cd_f * A_f))**2 / (2 * rho_f)

# assume a pressure drop through the lines that does not change with velocity
dP_line_f = 134 * PSI2PA # [Pa]
dP_line_ox = 10 * PSI2PA # [Pa]

P_pump_f = Pc + dP_line_ox + dP_inj_f 
P_pump_ox = Pc + dP_line_ox + dP_inj_ox

#design mass flow rate
design_idx = np.abs(mdot - mdot_design).argmin()
ox_design = P_pump_ox[design_idx]
fuel_design = P_pump_f[design_idx]


plt.plot(mdot, P_pump_f * 10**-5)
plt.plot(mdot, P_pump_ox * 10**-5)
plt.plot(mdot_design, ox_design * 10**-5, "ro")
plt.plot(mdot_design, fuel_design * 10**-5, "ro")

plt.legend(["Fuel Pump", "Ox Pump", "TCA Design Point"])
plt.xlabel("Mass Flow Rate [kg/s]")
plt.ylabel("Required Pressure [Bar]")
plt.title("TCA Head-Flow Curve")

plt.grid(True)
plt.show()

