import os
from matplotlib import contour
import numpy as np
import cea
from pint import UnitRegistry
import yaml
import pandas as pd
from Contour_Script import contour_length_mm
from Contour_Script import bell_nozzle

ureg = UnitRegistry()
with open('TCA_params_test.yaml') as f:
    p = yaml.safe_load(f)

### Engine Parameters/Inputs
F                  = np.arange(1000, 6000, 200) * ureg.lbf
pc_array           = np.arange(300, 900, 20) * ureg.psi

### https://usaindustries.com/piping-isolation-testing-products/pipe-schedule-chart/ ###
chamber_diameters  = np.array([4.813,4.563,4.313,4.063,5.761,5.501,5.189,4.897,3.826]) * ureg.inch
chamber_thickness = np.array([.375,.500,.625,.750,.432,.562,.719,.864,.337]) * ureg.inch 

pe = p['Exit_pressure'] * ureg.psi #Target exit pressure

eta_cstar = p['cstar_efficiency']
eta_cf    = p['cf_efficiency']
L_star    = p['L_star'] * ureg.m
alpha     = p['alpha_divergence'] * ureg.deg
of_ratio  = p['of_ratio']

### CEA Setup ###
# "C2H5OH(L)" -> Ethanol 
reac_names      = ["C2H5OH(L)", "O2(L)"]
T_reactant      = np.array([298.15, 90.17])   # plain floats, K
fuel_weights    = np.array([1.0, 0.0])
oxidant_weights = np.array([0.0, 1.0])

reac = cea.Mixture(reac_names)
prod = cea.Mixture(reac_names, products_from_reactants=True)
solver   = cea.RocketSolver(prod, reactants=reac, transport=True)
solution = cea.RocketSolution(solver)

weights = reac.of_ratio_to_weights(oxidant_weights, fuel_weights, of_ratio)
hc      = reac.calc_property(cea.ENTHALPY, weights, T_reactant) / cea.R  # plain float

### Build design matrix ###
designs = []

for pc in pc_array:
    pi_p    = [pc.to(ureg.psi).magnitude / pe.to(ureg.psi).magnitude]
    pc  = pc.to(ureg.bar)

    # IAC run to get C* and Cf for this pc
    solver.solve(solution, weights, pc.magnitude, pi_p, iac=True, hc=hc)

    Cf    = solution.coefficient_of_thrust[-1] * eta_cf
    cstar = solution.c_star[-1] * eta_cstar * ureg.m/ureg.s # m/s
    Isp  = cstar * Cf / ureg.g0

    for F_target in F:
        F_N   = F_target.to(ureg.N)
        mdot  = F_N / (Cf * cstar)                  # kg/s
        At    = (mdot * cstar) / pc.to(ureg.Pa)              # m²

        for D_c, t_c in zip(chamber_diameters, chamber_thickness):
            A_c   = np.pi / 4 * D_c.to(ureg.m)**2              # m²
            ac_at = A_c / At

            #Running the solver with the contraction ratio.
            solver.solve(solution, weights, pc.magnitude, pi_p,
                         ac_at=ac_at, iac=False, hc=hc)

            Dt = 2 * np.sqrt(At / np.pi)          # throat diameter, m
            Ae = solution.ae_at[-1] * At.to(ureg.m**2) # exit area, m²
            De = 2 * np.sqrt(Ae / np.pi)          # exit diameter,

            L_cone   = (D_c.to(ureg.m)/2 - Dt.to(ureg.m)/2) / np.tan(alpha.to(ureg.rad))
            V_cone  = (np.pi / 3) * ((D_c.to(ureg.m)/2)**2 + D_c.to(ureg.m)/2 *Dt.to(ureg.m)/2 + (Dt.to(ureg.m)/2)**2) * L_cone

            Vc_total = L_star.to(ureg.m) * At.to(ureg.m**2)            # total required volume
            V_cyl    = Vc_total - V_cone                               # subtract cone

            Lc = V_cyl / A_c.to(ureg.m**2)                              # cylindrical length (m)
            Ltotal = Lc + L_cone                                       # Total chamber length (m)
            
            #Von Mises stress calculations#
            sigma_th = (pc.to(ureg.Pa) * D_c.to(ureg.m)) / (2 * t_c.to(ureg.m))  # hoop stress - Seamless Pipe
            sigma_ax = (pc.to(ureg.Pa) * D_c.to(ureg.m)) / (4 * t_c.to(ureg.m)*0.6)  # axial stress - Weld coefficient of 0.6 for welded joints
            sigma_vM = np.sqrt(sigma_th**2 + sigma_ax**2 - sigma_th * sigma_ax)  # von Mises stress

            angles, contour, R2 = bell_nozzle(solution.ae_at[-1],Dt.to(ureg.mm).magnitude/2, 80,ac_at,alpha.to(ureg.deg).magnitude,Lc.to(ureg.mm).magnitude)
            total_length_mm = contour_length_mm(contour)


            designs.append({
                'F_kN'        : F_target.to(ureg.kN).magnitude,
                'pc_bar'       : pc.to(ureg.bar).magnitude,
                'D_chamber_in' : D_c.to(ureg.inch).magnitude,
                't_chamber_in' : t_c.to(ureg.inch).magnitude,
                'Lc_m'         : Lc.to(ureg.m).magnitude,
                'ac_at'        : ac_at.magnitude,
                'ae_at'        : solution.ae_at[-1],
                'mdot_kg/s'    : mdot.to(ureg.kg/ureg.s).magnitude,     # kg/s
                'Dt_in'        : Dt.to(ureg.inch).magnitude,     # in
                'De_in'        : De.to(ureg.inch).magnitude,     # in
                'C_star_ms'    : cstar.to(ureg.m/ureg.s).magnitude, # m/s
                'Cf'           : Cf,
                'Isp_s'        : Isp.to(ureg.s).magnitude,       # s
                'T_chamber_K'  : solution.T[0],
                'sigma_vM_MPa'  : sigma_vM.to(ureg.MPa).magnitude, # MPa
                'total_length_mm' : total_length_mm,
            })

df = pd.DataFrame(designs)

# Size thresholds
max_L_m   = 0.5 * ureg.m   # meters
min_L_m   = 0.3 * ureg.m   # meters
max_De_in  = 0.4 * ureg.m   # inches

# Filter out designs exceeding thresholds
df = df[(df['total_length_mm'] <= max_L_m.to(ureg.mm).magnitude) & (df['total_length_mm'] >= min_L_m.to(ureg.mm).magnitude) & (df['De_in'] <= max_De_in.to(ureg.inch).magnitude)]

print(f"Designs after filtering: {len(df)}")
print(df.to_string())
df.to_excel('designs.xlsx', index=False)

