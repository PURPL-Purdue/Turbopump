import numpy as np
import cea
from pint import UnitRegistry
import yaml
import contour_script_revised as cs

ureg = UnitRegistry()
with open('TCA_params_test.yaml') as f:
    p = yaml.safe_load(f)

### Engine Parameters/Inputs

csv_dxf_output = 'CSV_DXF_OUTPUTS/nozzle_contour3.csv ' # Output file for nozzle contour data (CSV)

F     = p['Thrust_target']  * ureg.lbf      # Target thrust             [lbf]
pc    = p['Chamber_pressure'] * ureg.psi    # Chamber pressure          [psia]
pe    = p['Exit_pressure'] * ureg.psi       # Exit pressure             [psia]

pi_p = [pc / pe]                            # Pressure ratio chamber to exit

eta_cstar = p['cstar_efficiency']           # C* efficiency
eta_cf = p['cf_efficiency']                 # Cf efficiency

L_star = p['L_star'] * ureg.inch             # Characteristic length (in)
alpha = p['alpha_divergence'] * ureg.deg     # Divergence half-angle (degrees)
ac_at = p['contraction_ratio']               #Contraction ratio (Ac/At)

of_ratio = p['of_ratio']

### CEA Setup ###

reac_names = ["C3H8O,2propanol", "O2(L)"]
T_reactant = np.array([298.15, 90.17]) * ureg.K
fuel_weights = np.array([1.0, 0.0])
oxidant_weights = np.array([0.0, 1.0])

reac = cea.Mixture(reac_names)
prod = cea.Mixture(reac_names, products_from_reactants=True)

solver = cea.RocketSolver(prod, reactants=reac)
solution = cea.RocketSolution(solver)

weights = reac.of_ratio_to_weights(oxidant_weights, fuel_weights, of_ratio) #Convert OF to weights.

#Compute chamber enthalpy. Normalized.
hc = reac.calc_property(cea.ENTHALPY, weights, T_reactant.to(ureg.kelvin).magnitude)/cea.R
#Solve the rocket problem for given inputs. Normalized.
solver.solve(solution, weights, pc.to(ureg.bar).magnitude, pi_p, ac_at=ac_at ,iac=False, hc=hc)

#### CEA OUTPUTS ####

num_pts = solution.num_pts

ae_at = solution.ae_at                  # Expansion ratio (Ae/At)
Isp_i = solution.Isp                    # Specific impulse at sea level (s)
Isp_vac = solution.Isp_vacuum           # Specific impulse in vacuum (s)
Cf_i = solution.coefficient_of_thrust   # Coefficient of thrust
cstar_i  = solution.c_star              # Characteristic velocity (m/s)
gamma = solution.gamma_s                # Specific heat ratio
Cp = solution.cp                        # Specific heat at constant pressure (J/kg-K)
MW = solution.MW                        # Molecular weight (kg/kmol)
k = solution.conductivity_eq            # Thermal conductivity (W/m-K). ####NEEEDS FIXING

T = solution.T #Temperature (K)
P = solution.P #Pressure (bar)

### Efficiency Calculations ###

cstar = cstar_i[-1] * eta_cstar *ureg.m/ureg.s  # Adjusted characteristic velocity (m/s)
Cf = Cf_i[-1] * eta_cf                          # Adjusted coefficient of thrust

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


##Run contour script - Export csv and dxf files.
angles, contour = cs.contour_script(80, ae_at[-1], ac_at, alpha, Lc.to(ureg.mm).magnitude, dt.to(ureg.mm).magnitude/2, csv_dxf_output)



print()
print("PERFORMANCE PARAMETERS")
print()


def format_values(values, skip_index=1, width=10, precision=3):
    return " ".join(
        f"{float(values[i]):{width}.{precision}f}"
        for i in range(len(values))
        if i != skip_index
    )

print(f"{'Ae/At':<15}{format_values(ae_at)}")
print(f"{'Isp_SL[m/s]':<15}{format_values(Isp_i)}")
print(f"{'Isp_vac[m/s]':<15}{format_values(Isp_vac)}")
print(f"{'Cf':<15}{format_values(Cf_i)}")
print(f"{'C*[m/s]':<15}{format_values(cstar_i)}")
print(f"{'T[K]':<15}{format_values(T)}")
print(f"{'P[bar]':<15}{format_values(P)}")
print(f"{'gamma':<15}{format_values(gamma)}")
print(f"{'Cp[J/kg-K]':<15}{format_values(Cp)}")
print(f"{'k[mW/m-K]':<15}{format_values(k)}") ###NEEDS FIXING.
print(f"{'MW[kg/kmol]':<15}{format_values(MW)}")


print()
print(f"Ac/At: {ac_at:.4f}")
print(f"Ae/At: {ae_at[-1]:.4f}")
print(f"Mass flow rate (lbm/s): {mdot.to(ureg.lb/ureg.s):.4f}")
print(f"Throat area (m²): {At.to(ureg.m**2):.6f}")
print(f"Throat diameter (m): {dt.to(ureg.inch):.4f}")
print(f"Exit area (m²): {Ae.to(ureg.m**2):.6f}")
print(f"Exit diameter (m): {de.to(ureg.inch):.4f}")
print()
print("CHAMBER GEOMETRY")
print(f"Chamber diameter (in): {dc.to(ureg.inch):.4f}")
print(f"Cylindrical length (in): {Lc.to(ureg.inch):.4f}")
print(f"Conical length (in): {L_cone.to(ureg.inch):.4f}")
print(f"Total chamber length (in): {Ltotal.to(ureg.inch):.4f}")







