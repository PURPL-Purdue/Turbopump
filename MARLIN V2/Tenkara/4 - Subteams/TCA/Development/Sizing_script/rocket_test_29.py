import os
import numpy as np
import cea
from pint import UnitRegistry
import yaml
import Contour_Script as cs
import pandas as pd


ureg = UnitRegistry()
with open('TCA_params_test.yaml') as f:
    p = yaml.safe_load(f)

### Engine Parameters/Inputs

output_filename = 'CSV_DXF_OUTPUTS/nozzle_contour2600-600psi' # Output file for nozzle contour data

F     = p['Thrust_target']  * ureg.lbf      # Target thrust             [lbf]
pc    = p['Chamber_pressure'] * ureg.psi    # Chamber pressure          [psia]
pe    = p['Exit_pressure'] * ureg.psi       # Exit pressure             [psia]

pi_p = [pc / pe]                            # Pressure ratio chamber to exit

eta_cstar = p['cstar_efficiency']           # C* efficiency
eta_cf = p['cf_efficiency']                 # Cf efficiency

L_star = p['L_star'] * ureg.m             # Characteristic length (m)
alpha = p['alpha_divergence'] * ureg.deg     # Divergence half-angle (degrees)
Dc = p['chamber_diameter'] * ureg.inch       # Chamber diameter (in)


of_ratio = p['of_ratio']

### CEA Setup ###

reac_names = ["C2H5OH(L)", "O2(L)"]
T_reactant = np.array([298.15, 90.17]) * ureg.K
fuel_weights = np.array([1.0, 0.0])
oxidant_weights = np.array([0.0, 1.0])

reac = cea.Mixture(reac_names)
prod = cea.Mixture(reac_names, products_from_reactants=True)

solver = cea.RocketSolver(prod, reactants=reac, transport=True)
solution = cea.RocketSolution(solver)

weights = reac.of_ratio_to_weights(oxidant_weights, fuel_weights, of_ratio) #Convert OF to weights.

#Compute chamber enthalpy. Normalized.
hc = reac.calc_property(cea.ENTHALPY, weights, T_reactant.to(ureg.kelvin).magnitude)/cea.R
#Solve the rocket problem for given inputs. Normalized.

# IAC run to get C* and Cf for this pc
solver.solve(solution, weights, pc.to(ureg.bar).magnitude, pi_p, iac=True, hc=hc)
Cf    = solution.coefficient_of_thrust[-1] * eta_cf
cstar = solution.c_star[-1] * eta_cstar * ureg.m/ureg.s # m/s

#Obtain throat area from mass flow rate, C*, and chamber pressure.
mdot  = F.to(ureg.N) / (Cf * cstar)                  # kg/s
At    = (mdot * cstar) / pc.to(ureg.Pa)              # m²

#Obtain the required contraction ratio.
A_c   = np.pi / 4 * Dc.to(ureg.m)**2              # m²
ac_at = A_c.magnitude / At.magnitude

#Running the solver with the obtained contraction ratio.
solver.solve(solution, weights, pc.to(ureg.bar).magnitude, pi_p, ac_at=ac_at ,iac=False, hc=hc)

#### CEA OUTPUTS ####

num_pts = solution.num_pts

ae_at = solution.ae_at                  # Expansion ratio (Ae/At)
Cf_i = solution.coefficient_of_thrust   # Coefficient of thrust
cstar_i  = solution.c_star              # Characteristic velocity (m/s)
gamma = solution.gamma_s                # Specific heat ratio Frozen?
Cp = solution.cp * ureg.kJ / (ureg.kg * ureg.K)  # Specific heat at constant pressure (KJ/kg-K)
MW = solution.MW * ureg.kg / ureg.kmol           # Molecular weight (kg/kmol)
k = solution.conductivity_eq * (ureg.uW / (ureg.cm * ureg.K))  # Thermal conductivity (uW/cm-K).
R = 8.31446261815324 * ureg.kJ / (ureg.kmol * ureg.K) / MW.to(ureg.kg/ureg.kmol)  # Specific gas constant (kJ/kg-K)
gamma2 = Cp / (Cp - R)                                # Specific heat ratio from Cp and R
mass_fractions = solution.mass_fractions           # Mass fractions of species in the flow


T = solution.T #Temperature (K)
P = solution.P #Pressure (bar)

### Efficiency Calculations ###

cstar = cstar_i[-1] * eta_cstar *ureg.m/ureg.s  # Adjusted characteristic velocity (m/s)
Cf = Cf_i[-1] * eta_cf                          # Adjusted coefficient of thrust
Isp = Cf * cstar.to(ureg.m/ureg.s) / ureg.g0         # Specific impulse at sea level (s)

### Output Results ###

mdot    = F.to(ureg.N) / (Cf * cstar.to(ureg.m/ureg.s))     # Mass flow rate (kg/s)
At   = (mdot * cstar.to(ureg.m/ureg.s)) / pc.to(ureg.Pa)    # throat area  [m²]
dt = 2 * np.sqrt(At / np.pi)                                # throat diameter (m)
Ae = ae_at[-1] * At                                         # Exit area (m²)
de = 2 * np.sqrt(Ae / np.pi)                                # Exit diameter (m)

### Chamber Geometry ###

Ac = ac_at * At                                             # Chamber area (m²)
dc = 2 * np.sqrt(Ac / np.pi)                                # Chamber diameter (m)

# ── Chamber cylindrical length from L* ───────────────────────────────────────
# L* = Vc / At  where Vc includes the cylindrical section + converging cone
# Lc = (L*_cm - V_cone/At_cm) / con_r
# V_cone/At_cm = (1/3) * Rt_cm * (1/tan(α)) * (con_r^(1/3) - 1)

L_cone   = (dc.to(ureg.m)/2 - dt.to(ureg.m)/2) / np.tan(alpha.to(ureg.rad))
V_cone  = (np.pi / 3) * ((dc.to(ureg.m)/2)**2 + dc.to(ureg.m)/2 *dt.to(ureg.m)/2 + (dt.to(ureg.m)/2)**2) * L_cone

Vc_total = L_star.to(ureg.m) * At.to(ureg.m**2)            # total required volume
V_cyl    = Vc_total - V_cone                               # subtract cone

Lc = V_cyl / Ac.to(ureg.m**2)                              # cylindrical length (m)
Ltotal = Lc + L_cone                                       # Total chamber length (m)

#### Nozzle Contour Generation and Plotting ###
angles, contour, R2 = cs.bell_nozzle(ae_at[-1], dt.to(ureg.mm).magnitude/2, 80, ac_at, alpha.to(ureg.deg).magnitude, Lc.to(ureg.mm).magnitude)

# Plot the engine overview contour #
title = (f'Bell Nozzle\n'
        f'[ε = {round(ae_at[-1], 1)}, '
        f'Rt = {round(dt.to(ureg.mm).magnitude/2, 2)} mm, '
        f'L% = {80}%]')    
#cs.plot_overview(title, dt.to(ureg.mm).magnitude/2, angles, contour,output_filename) 

# Export the contour in csv and dxf #
if output_filename:
    base     = os.path.splitext(output_filename)[0]
    csv_file = base + '.csv'
    dxf_file = base + '.dxf'
else:
    csv_file = None
    dxf_file = None

#cs.export_nozzle_csv(contour, filename=csv_file)
#cs.export_nozzle_dxf(contour, filename=dxf_file)

### Output performance parameters ###

print()
print("PERFORMANCE PARAMETERS")
print()


def format_values(values, skip_index=0, width=10, precision=5):
    return " ".join(
        f"{float(values[i]):{width}.{precision}f}"
        for i in range(len(values))
        if i != skip_index
    )

print(f"{'Ae/At':<15}{format_values(ae_at)}")
print(f"{'Cf':<15}{format_values(Cf_i)}")
print(f"{'C*[m/s]':<15}{format_values(cstar_i)}")
print(f"{'T[K]':<15}{format_values(T)}")
print(f"{'P[bar]':<15}{format_values(P)}")
print(f"{'gamma':<15}{format_values(gamma)}")
print(f"{'gamma2':<15}{format_values(gamma2)}")
print(f"{'Cp[J/kg-K]':<15}{format_values(Cp.to(ureg.J / (ureg.kg * ureg.K)).magnitude)}")
print(f"{'k[uW/cm-K]':<15}{format_values(k.to(ureg.uW / (ureg.cm * ureg.K)).magnitude)}") 
print(f"{'MW[kg/kmol]':<15}{format_values(MW.magnitude)}")


print()
print(f"Ac/At: {ac_at:.4f}")
print(f"Ae/At: {ae_at[-1]:.4f}")
print(f"Cstar (m/s): {cstar.to(ureg.m/ureg.s):.2f}")
print(f"Cf: {Cf:.4f}")
print(f"Isp (s): {Isp.to(ureg.s):.2f}")
print(f"Mass flow rate (kg/s): {mdot.to(ureg.kg/ureg.s):.4f}")
print(f"Throat area (m²): {At.to(ureg.m**2):.6f}")
print(f"Throat diameter (m): {dt.to(ureg.inch):.4f}")
print(f"Exit area (m²): {Ae.to(ureg.m**2):.6f}")
print(f"Exit diameter (in): {de.to(ureg.inch):.4f}")
print()
print("CHAMBER GEOMETRY")
print(f"Chamber diameter (in): {dc.to(ureg.inch):.4f}")
print(f"Cylindrical length (m): {Lc.to(ureg.m):.4f}")
print(f"Conical length (in): {L_cone.to(ureg.inch):.4f}")
print(f"Total chamber length (in): {Ltotal.to(ureg.inch):.4f}")

ox, fuel = weights
total = ox + fuel
x_ox = ox / total
x_fuel = fuel / total

print(f"Oxidizer fraction: {x_ox:.6f}")
print(f"Fuel fraction: {x_fuel:.6f}")

mf_df = pd.DataFrame(mass_fractions)
print("\nMass Fractions of Species in the Flow:")
print(mf_df)
# Save to Excel
mf_df.to_excel("mass_fractions.xlsx", index=True)

print( T)





