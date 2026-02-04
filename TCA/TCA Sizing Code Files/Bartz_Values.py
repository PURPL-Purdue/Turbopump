"""BartZ_Values.py
Purpose: compute thermophysical properties at key engine stations across an O/F range
and export two CSVs:
  - `turbopump_properties.csv` : rows of [MR, Station, Temperature_K, Pressure_Pa, Density_kg_m3,
                                      Mach, Speed_of_Sound_m_s, Gamma, Cp_J_kgK, Thermal_Conductivity, Prandtl]
  - `turbopump_poly_coeffs.csv` : polynomial coefficients (degree up to 3) for each property vs Temperature (K)

Notes:
  - Reads optional overrides from `TCA_params.yaml` when available (cp, tc, visc conversions, contraction ratio)
  - Uses rocketcea `CEA_Obj` for gas properties (LOX / RP1)
  - Code is written for clarity and reproducibility; small dataset checks are handled gracefully.
"""

import os
import csv
import yaml
import numpy as np
from rocketcea.cea_obj import CEA_Obj

np.set_printoptions(legacy='1.25')

# Initialize CEA for LOX / RP-1
C = CEA_Obj(oxName='LOX', fuelName='RP1')

# --------------------------
# Global conversion factors (defaults; can be overridden by TCA_params.yaml)
# --------------------------
psi_to_pa = 6894.76
rankine_to_kelvin = 5.0 / 9.0
Ru_m = 8314.462618   # J / kmol-K

# Conversion defaults (can be overridden in YAML)
c_siconv = 4186.8    # cp conversion to J/(kg-K)
tc_siconv = 418.4    # thermal conductivity conversion (kept consistent with repo conventions)
vis_siconv = 0.0001  # viscosity conversion to Pa-s

# --------------------------
# Load YAML config (look in repo root or cwd)
# --------------------------
repo_root = os.path.abspath(os.path.join(os.path.dirname(__file__), '..', '..'))
cfg_path = os.path.join(repo_root, 'TCA_params.yaml')
if not os.path.exists(cfg_path):
    cfg_path = os.path.join(os.getcwd(), 'TCA_params.yaml')

try:
    with open(cfg_path) as f:
        tca_params = yaml.safe_load(f) or {}
except FileNotFoundError:
    tca_params = {}

c_siconv = tca_params.get('cp_conversion', c_siconv)
tc_siconv = tca_params.get('tc_conversion', tc_siconv)
vis_siconv = tca_params.get('visc_conversion', vis_siconv)
contraction_ratio = tca_params.get('tca_contraction_ratio', 10.0)

# --------------------------
# Helper functions
# --------------------------

def compute_gamma(cp_JkgK, MW_kg_per_kmol):
    """Estimate specific heat ratio gamma from cp and molecular weight.
    gamma = cp / (cp - Rspec), where Rspec = Ru / MW (J/kg-K)
    Returns nan if invalid values encountered.
    """
    if cp_JkgK is None or MW_kg_per_kmol is None:
        return np.nan
    R_spec = Ru_m / MW_kg_per_kmol
    denom = cp_JkgK - R_spec
    if denom <= 0:
        return np.nan
    return cp_JkgK / denom

# --------------------------
# User inputs (single O/F evaluation)
# --------------------------
mr = float(input('O/F ratio to evaluate: '))
pc = float(input('Chamber pressure (psia): '))
pe = 14.7  # exit pressure (psia) used for expansion calculations (kept as default)

# Create a one-element array so the existing loop logic is reused without an O/F sweep
MRs = np.array([mr])

# --------------------------
# Output containers and headers
# --------------------------
# Note: Added `Molecular_Weight_kg_kmol` to properties and to fitted properties list
properties_header = ['MR', 'Station', 'Temperature_K', 'Pressure_Pa', 'Density_kg_m3', 'Mach',
                     'Speed_of_Sound_m_s', 'Gamma', 'Cp_J_kgK', 'Molecular_Weight_kg_kmol', 'Thermal_Conductivity', 'Viscosity_Pa_s', 'Prandtl']
properties_rows = []

stations = ['Injector', 'Combustor', 'Throat', 'Exit']
# Properties kept for general output; polynomial fittings will be restricted to a subset below
prop_keys = ['Temperature_K', 'Pressure_Pa', 'Density_kg_m3', 'Mach',
             'Speed_of_Sound_m_s', 'Gamma', 'Cp_J_kgK', 'Molecular_Weight_kg_kmol', 'Thermal_Conductivity', 'Viscosity_Pa_s', 'Prandtl']

# Structure for polynomial data: poly_data[station][prop] = list of values (aligned by Temperature_K)
poly_data = {st: {pk: [] for pk in prop_keys} for st in stations}

# --------------------------
# Main loop: sample properties for each O/F
# --------------------------
for mr in MRs:
    # Compute expansion (eps) and temperatures from CEA
    Eps = C.get_eps_at_PcOvPe(Pc=pc, MR=mr, PcOvPe=(pc / pe), frozen=0)
    Tc_R, Tt_R, Te_R = C.get_Temperatures(Pc=pc, MR=mr, eps=Eps, frozen=0, frozenAtThroat=0)

    # Convert temperatures to Kelvin for consistent fitting and outputs
    Tc_K = Tc_R * rankine_to_kelvin
    Tt_K = Tt_R * rankine_to_kelvin
    Te_K = Te_R * rankine_to_kelvin

    # Chamber (used for Injector and Combustor values)
    heat_cap_c, viscosity_c, therm_con_c, pr_c = C.get_Chamber_Transport(Pc=pc, MR=mr, eps=Eps, frozen=0)
    cp_c = heat_cap_c * c_siconv
    k_c = therm_con_c * tc_siconv
    mu_c = viscosity_c * vis_siconv
    prandtl_c = pr_c
    MW_C = C.get_Chamber_MolWt_gamma(Pc=pc, MR=mr, eps=Eps)
    MW_C_val = MW_C[0] if isinstance(MW_C, (list, tuple, np.ndarray)) else MW_C

    # Chamber Mach number at combustor end
    M_comb = C.get_Chamber_MachNumber(Pc=pc, MR=mr, fac_CR=contraction_ratio)

    # Throat
    heat_cap_t, viscosity_t, therm_con_t, pr_t = C.get_Throat_Transport(Pc=pc, MR=mr, eps=Eps, frozen=0)
    cp_t = heat_cap_t * c_siconv
    k_t = therm_con_t * tc_siconv
    mu_t = viscosity_t * vis_siconv
    prandtl_t = pr_t
    MW_T = C.get_Throat_MolWt_gamma(Pc=pc, MR=mr, eps=Eps)
    MW_T_val = MW_T[0] if isinstance(MW_T, (list, tuple, np.ndarray)) else MW_T
    M_th = 1.0

    # Exit
    heat_cap_e, viscosity_e, therm_con_e, pr_e = C.get_Exit_Transport(Pc=pc, MR=mr, eps=Eps, frozen=0)
    cp_e = heat_cap_e * c_siconv
    k_e = therm_con_e * tc_siconv
    mu_e = viscosity_e * vis_siconv
    prandtl_e = pr_e
    MW_E = C.get_exit_MolWt_gamma(Pc=pc, MR=mr, eps=Eps, frozen=0)
    MW_E_val = MW_E[0] if isinstance(MW_E, (list, tuple, np.ndarray)) else MW_E
    M_ex = C.get_MachNumber(Pc=pc, MR=mr, eps=Eps, frozen=0)

    # Total chamber pressure in Pa
    Pc_Pa = pc * psi_to_pa

    # ----- Injector (approx stagnation conditions; Mach ~ 0) -----
    M_inj = 0.0
    T_inj_K = Tc_K
    gamma_inj = compute_gamma(cp_c, MW_C_val)
    p_inj = Pc_Pa  # static approx = total when M~0
    R_spec_inj = Ru_m / MW_C_val
    rho_inj = p_inj / (R_spec_inj * T_inj_K) if R_spec_inj > 0 else np.nan
    a_inj = np.sqrt(gamma_inj * R_spec_inj * T_inj_K) if not np.isnan(gamma_inj) else np.nan

    # Injector: include molecular weight from chamber composition (MW_C_val)
    properties_rows.append([mr, 'Injector', T_inj_K, p_inj, rho_inj, M_inj, a_inj, gamma_inj, cp_c, MW_C_val, k_c, mu_c, prandtl_c])
    poly_data['Injector']['Temperature_K'].append(T_inj_K)
    poly_data['Injector']['Pressure_Pa'].append(p_inj)
    poly_data['Injector']['Density_kg_m3'].append(rho_inj)
    poly_data['Injector']['Mach'].append(M_inj)
    poly_data['Injector']['Speed_of_Sound_m_s'].append(a_inj)
    poly_data['Injector']['Gamma'].append(gamma_inj)
    poly_data['Injector']['Cp_J_kgK'].append(cp_c)
    poly_data['Injector']['Molecular_Weight_kg_kmol'].append(MW_C_val)
    poly_data['Injector']['Thermal_Conductivity'].append(k_c)
    poly_data['Injector']['Viscosity_Pa_s'].append(mu_c)
    poly_data['Injector']['Prandtl'].append(prandtl_c) 

    # ----- Combustor (use chamber static at Mach = M_comb) -----
    gamma_comb = compute_gamma(cp_c, MW_C_val)
    # Isentropic relation: p = p0 / (1 + 0.5*(gamma-1)*M^2)^(gamma/(gamma-1))
    p_comb = Pc_Pa / (1 + 0.5 * (gamma_comb - 1) * M_comb ** 2) ** (gamma_comb / (gamma_comb - 1)) if not np.isnan(gamma_comb) else np.nan
    T_comb_K = Tc_K / (1 + 0.5 * (gamma_comb - 1) * M_comb ** 2) if not np.isnan(gamma_comb) else Tc_K
    R_spec_comb = Ru_m / MW_C_val
    rho_comb = p_comb / (R_spec_comb * T_comb_K) if R_spec_comb > 0 else np.nan
    a_comb = np.sqrt(gamma_comb * R_spec_comb * T_comb_K) if not np.isnan(gamma_comb) else np.nan

    # Combustor end: molecular weight same as chamber composition (MW_C_val)
    properties_rows.append([mr, 'Combustor', T_comb_K, p_comb, rho_comb, M_comb, a_comb, gamma_comb, cp_c, MW_C_val, k_c, mu_c, prandtl_c])
    poly_data['Combustor']['Temperature_K'].append(T_comb_K)
    poly_data['Combustor']['Pressure_Pa'].append(p_comb)
    poly_data['Combustor']['Density_kg_m3'].append(rho_comb)
    poly_data['Combustor']['Mach'].append(M_comb)
    poly_data['Combustor']['Speed_of_Sound_m_s'].append(a_comb)
    poly_data['Combustor']['Gamma'].append(gamma_comb)
    poly_data['Combustor']['Cp_J_kgK'].append(cp_c)
    poly_data['Combustor']['Molecular_Weight_kg_kmol'].append(MW_C_val)
    poly_data['Combustor']['Thermal_Conductivity'].append(k_c)
    poly_data['Combustor']['Viscosity_Pa_s'].append(mu_c)
    poly_data['Combustor']['Prandtl'].append(prandtl_c) 

    # ----- Throat -----
    gamma_th = compute_gamma(cp_t, MW_T_val)
    T_th_K = Tt_K
    # Pressure at throat using isentropic relation from chamber total (approx)
    p_th = Pc_Pa / (1 + 0.5 * (gamma_th - 1) * M_th ** 2) ** (gamma_th / (gamma_th - 1)) if not np.isnan(gamma_th) else np.nan
    R_spec_th = Ru_m / MW_T_val
    rho_th = p_th / (R_spec_th * T_th_K) if R_spec_th > 0 else np.nan
    a_th = np.sqrt(gamma_th * R_spec_th * T_th_K) if not np.isnan(gamma_th) else np.nan

    # Throat: molecular weight from throat composition (MW_T_val)
    properties_rows.append([mr, 'Throat', T_th_K, p_th, rho_th, M_th, a_th, gamma_th, cp_t, MW_T_val, k_t, mu_t, prandtl_t])
    poly_data['Throat']['Temperature_K'].append(T_th_K)
    poly_data['Throat']['Pressure_Pa'].append(p_th)
    poly_data['Throat']['Density_kg_m3'].append(rho_th)
    poly_data['Throat']['Mach'].append(M_th)
    poly_data['Throat']['Speed_of_Sound_m_s'].append(a_th)
    poly_data['Throat']['Gamma'].append(gamma_th)
    poly_data['Throat']['Cp_J_kgK'].append(cp_t)
    poly_data['Throat']['Molecular_Weight_kg_kmol'].append(MW_T_val)
    poly_data['Throat']['Thermal_Conductivity'].append(k_t)
    poly_data['Throat']['Viscosity_Pa_s'].append(mu_t)
    poly_data['Throat']['Prandtl'].append(prandtl_t) 

    # ----- Exit -----
    gamma_ex = compute_gamma(cp_e, MW_E_val)
    T_ex_K = Te_K / (1 + 0.5 * (gamma_ex - 1) * M_ex ** 2) if not np.isnan(gamma_ex) else Te_K
    p_ex = Pc_Pa / (1 + 0.5 * (gamma_ex - 1) * M_ex ** 2) ** (gamma_ex / (gamma_ex - 1)) if not np.isnan(gamma_ex) else np.nan
    R_spec_ex = Ru_m / MW_E_val
    rho_ex = p_ex / (R_spec_ex * T_ex_K) if R_spec_ex > 0 else np.nan
    a_ex = np.sqrt(gamma_ex * R_spec_ex * T_ex_K) if not np.isnan(gamma_ex) else np.nan

    # Exit: molecular weight from exit composition (MW_E_val)
    properties_rows.append([mr, 'Exit', T_ex_K, p_ex, rho_ex, M_ex, a_ex, gamma_ex, cp_e, MW_E_val, k_e, mu_e, prandtl_e])
    poly_data['Exit']['Temperature_K'].append(T_ex_K)
    poly_data['Exit']['Pressure_Pa'].append(p_ex)
    poly_data['Exit']['Density_kg_m3'].append(rho_ex)
    poly_data['Exit']['Mach'].append(M_ex)
    poly_data['Exit']['Speed_of_Sound_m_s'].append(a_ex)
    poly_data['Exit']['Gamma'].append(gamma_ex)
    poly_data['Exit']['Cp_J_kgK'].append(cp_e)
    poly_data['Exit']['Molecular_Weight_kg_kmol'].append(MW_E_val)
    poly_data['Exit']['Thermal_Conductivity'].append(k_e)
    poly_data['Exit']['Viscosity_Pa_s'].append(mu_e)
    poly_data['Exit']['Prandtl'].append(prandtl_e) 

# --------------------------
# Write properties CSV
# --------------------------
properties_csv = os.path.join(os.getcwd(), 'TCA/turbopump_properties.csv')
with open(properties_csv, 'w', newline='') as f:
    writer = csv.writer(f)
    writer.writerow(properties_header)
    for r in properties_rows:
        writer.writerow(r)

print(f"Saved properties CSV: {properties_csv}")

# --------------------------
# Polynomial fitting: fit Gamma, Cp, Thermal Conductivity, and Molecular Weight
# as functions of Temperature (Kelvin) using three station points: Combustor, Throat, Exit
# Use quadratic (degree=2) when three valid points are available; otherwise degree = n-1
# --------------------------
poly_csv = os.path.join(os.getcwd(), 'TCA/turbopump_poly_coeffs.csv')
with open(poly_csv, 'w', newline='') as f:
    writer = csv.writer(f)
    writer.writerow(['O/F', 'Pc_psia', 'Property', 'Degree', 'Coefficients_high->low', 'T_points_K', 'Y_points'])

    fit_props = ['Gamma', 'Cp_J_kgK', 'Thermal_Conductivity', 'Viscosity_Pa_s', 'Molecular_Weight_kg_kmol']

    # Temperatures at the three stations (should each have one entry for the chosen O/F)
    T_comb = poly_data['Combustor']['Temperature_K'][0] if len(poly_data['Combustor']['Temperature_K']) > 0 else np.nan
    T_th = poly_data['Throat']['Temperature_K'][0] if len(poly_data['Throat']['Temperature_K']) > 0 else np.nan
    T_ex = poly_data['Exit']['Temperature_K'][0] if len(poly_data['Exit']['Temperature_K']) > 0 else np.nan
    T_points = np.array([T_comb, T_th, T_ex])

    for pk in fit_props:
        y_comb = poly_data['Combustor'].get(pk, [np.nan])[0]
        y_th = poly_data['Throat'].get(pk, [np.nan])[0]
        y_ex = poly_data['Exit'].get(pk, [np.nan])[0]
        Y_points = np.array([y_comb, y_th, y_ex])

        # filter out any invalid points
        mask = (~np.isnan(T_points)) & (~np.isnan(Y_points))
        Tfit = T_points[mask]
        Yfit = Y_points[mask]

        if Tfit.size < 2:
            deg = 0
            coeffs = []
        else:
            deg = min(2, Tfit.size - 1)
            coeffs = np.polyfit(Tfit, Yfit, deg).tolist()

        writer.writerow([mr, pc, pk, deg, coeffs, Tfit.tolist(), Yfit.tolist()])

print(f"Saved polynomial coefficients CSV: {poly_csv}")

# End


