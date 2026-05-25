"""
Calculates LOX/RP-1 chamber pressure and required fuel/oxidizer feed pressures for a given total mass flow rate
and O/F ratio, using RocketCEA for c* and an incompressible orifice (CdA) model for the injector pressure drops.

Author: Joaquin Alarcon
"""

import numpy as np
import matplotlib.pyplot as plt
from rocketcea.cea_obj import CEA_Obj
from CoolProp.CoolProp import PropsSI

# Unit conversion constants
psi_to_pa = 6894.76 
ft_to_m = 0.3048
meters_to_inches = 39.3701
kg_to_lbm = 2.20462

# Parameters for the calculator (inputs)
mdot_total_lbm = 23.630 # Total propellant mass flow rate [lbm/s]
mdot_total_kg = mdot_total_lbm / kg_to_lbm # Convert mass flow rate to kg/s for calculations

OF_low_end = 1.3 # Lower end of O/F ratio range to analyze
OF_high_end = 2.0 # Upper end of O/F ratio range to analyze
number_of_OF_points = 8 # Number of O/F points to evaluate between low and high end
OF_ratio = np.linspace(OF_low_end, OF_high_end, number_of_OF_points) # Oxidizer to Fuel mass ratio

tca_throat_diameter_inches = 3.0746 # Throat diameter in inches, used to calculate throat area. Adjust this value based on your actual nozzle throat size.
tca_throat_diameter_meters = tca_throat_diameter_inches / meters_to_inches
throat_area = np.pi * (tca_throat_diameter_meters / 2) ** 2 # Throat area in m^2

cstar_eff = 0.875 # C* efficiency to apply to the ideal C* calculation, accounting for non-idealities in the combustion process. 

number_of_fuel_holes = 80 # number of injector holes per side (fuel and oxidizer)
Cd = 0.8 # Discharge coefficient for the injector holes, accounting for non-ideal flow. 
diameter_fuel_hole_meters = 0.0015113 # Diameter of each fuel injector hole in meters. [1.2 mm]
diameter_ox_hole_meters = 0.0017018 # Diameter of each oxidizer injector hole in meters. [1.7018 mm]
fuel_CdA = number_of_fuel_holes * Cd * np.pi * (diameter_fuel_hole_meters / 2) ** 2 # Effective flow area for fuel side (m^2)
ox_CdA = number_of_fuel_holes * Cd * np.pi * (diameter_ox_hole_meters / 2) ** 2 # Effective flow area for oxidizer side (m^2)

# Fluid properties (density). Using approximation inlet conditions. 
P_ref = 700*psi_to_pa # Reference pressure for density calculation (Pa)
rho_fuel = 789 # PropsSI('D', 'T', 298, 'P', P_ref, 'Isopropanol') fuel density at 298 K [kg/m^3]
rho_ox = PropsSI('D', 'T', 90, 'P', P_ref, 'Oxygen') # oxidizer density at 90 K [kg/m^3]

def feed_pressures(mdot_total_kg, OF_ratio, throat_area, cstar_eff, fuel_CdA, ox_CdA, rho_fuel, rho_ox, ):
    cea = CEA_Obj(oxName='LOX', fuelName='Isopropanol') # RocketCEA object to compute c* for LOX/Isopropanol given Pc and MR
    Pc_psi = 500 # initial guess for chamber pressure (psi)
    print("-----------------------------------------------------------------------")
    print("| OF Ratio | Chamber Pressure (psi) | Fuel Feed (psi) | Ox Feed (psi) |")
    for OF in OF_ratio:
        for _ in range(100): 
            mdot_fuel = mdot_total_kg / (1 + OF) # split total mass flow into fuel [kg/s]
            mdot_ox = mdot_total_kg - mdot_fuel # split total mass flow into oxidizer [kg/s]
            cstar_ft_s = cea.get_Cstar(Pc=Pc_psi, MR=OF) # ideal characteristic velocity from CEA [ft/s] 
            cstar_m_s = cstar_ft_s * ft_to_m # Convert c* to [m/s]
            Pc_Pa = mdot_total_kg * (cstar_eff * cstar_m_s) / throat_area # Chocked throat relation
            Pc_psi = Pc_Pa / psi_to_pa # Convert chamber pressure back to psi for the next iteration
        Pc_Pa = Pc_psi * psi_to_pa # final chamber pressure in Pa after convergence
        dp_fuel = (mdot_fuel / fuel_CdA) ** 2 / (2.0 * rho_fuel) # Pressure drop across fuel side (Pa)
        dp_ox = (mdot_ox / ox_CdA) ** 2 / (2.0 * rho_ox) # Pressure drop across oxidizer side (Pa)
        P_fuel_feed_psi = (Pc_Pa + dp_fuel) / psi_to_pa # Pressure feed for fuel [psi]
        P_ox_feed_psi = (Pc_Pa + dp_ox) / psi_to_pa # Pressure feed for oxidizer [psi]
        print("-----------------------------------------------------------------------")
        print(f"| {OF:.2f}     | {Pc_psi:.4f}               | {P_fuel_feed_psi:.4f}        | {P_ox_feed_psi:.4f}      |")
    print("-----------------------------------------------------------------------")

feed_pressures(mdot_total_kg, OF_ratio, throat_area, cstar_eff, fuel_CdA, ox_CdA, rho_fuel, rho_ox)