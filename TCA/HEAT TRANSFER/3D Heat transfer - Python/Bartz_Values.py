"""Bartz_Values.py
Purpose: compute thermophysical properties at key engine stations and export:
  - turbopump_properties.csv : rows of [MR, Station, Temperature_K, ...]
  - turbopump_poly_coeffs.csv : polynomial coefficients for each property vs T

Inputs are read from TCA_params.yaml (O/F ratio, chamber pressure, etc.).
Can be run standalone or imported and called via run().

Usage (standalone):
    python Bartz_Values.py
    # reads TCA_params.yaml from TCA/ relative to repo root, writes CSVs next to this script

Usage (from run.py):
    import Bartz_Values
    Bartz_Values.run(tca_params, output_dir)
"""

import os
import csv
import yaml
import numpy as np
from rocketcea.cea_obj import CEA_Obj

np.set_printoptions(legacy='1.25')

# --------------------------
# Module-level constants
# --------------------------
PSI_TO_PA        = 6894.76
RANKINE_TO_KELVIN = 5.0 / 9.0
RU_M             = 8314.462618   # J / kmol-K

# Default unit-conversion factors (can be overridden in TCA_params.yaml)
_DEFAULT_C_SICONV   = 4186.8    # cp:   raw CEA units -> J/(kg-K)
_DEFAULT_TC_SICONV  = 418.4     # k:    raw CEA units -> W/(m-K)
_DEFAULT_VIS_SICONV = 0.0001    # mu:   raw CEA units -> Pa-s


def compute_gamma(cp_JkgK, MW_kg_per_kmol):
    """Estimate gamma from cp and molecular weight using ideal-gas relation."""
    if cp_JkgK is None or MW_kg_per_kmol is None:
        return np.nan
    R_spec = RU_M / MW_kg_per_kmol
    denom  = cp_JkgK - R_spec
    if denom <= 0:
        return np.nan
    return cp_JkgK / denom


def run(tca_params, output_dir):
    """
    Compute gas properties via CEA and write CSVs to output_dir.

    Parameters
    ----------
    tca_params : dict
        Contents of TCA_params.yaml.  Required keys:
            oxidizer_fuel_ratio   (dimensionless)
            tca_chamber_pressure  (psia)
        Optional keys (with defaults):
            tca_exit_pressure     (psia, default 14.7)
            tca_contraction_ratio (dimensionless, default 10.0)
            cp_conversion         (default 4186.8)
            tc_conversion         (default 418.4)
            visc_conversion       (default 0.0001)
    output_dir : str
        Directory where turbopump_properties.csv and
        turbopump_poly_coeffs.csv will be written.

    Returns
    -------
    str
        Absolute path to the written turbopump_properties.csv.
    """
    # --- Unpack inputs ---
    mr = float(tca_params['oxidizer_fuel_ratio'])
    pc = float(tca_params['tca_chamber_pressure'])   # psia
    pe = float(tca_params.get('tca_exit_pressure', 14.7))   # psia

    contraction_ratio = float(tca_params.get('tca_contraction_ratio', 10.0))
    c_siconv  = float(tca_params.get('cp_conversion',   _DEFAULT_C_SICONV))
    tc_siconv = float(tca_params.get('tc_conversion',   _DEFAULT_TC_SICONV))
    vis_siconv = float(tca_params.get('visc_conversion', _DEFAULT_VIS_SICONV))

    # --- CEA object (LOX / RP-1) ---
    C = CEA_Obj(oxName='LOX', fuelName='RP1')

    # --- Output containers ---
    properties_header = [
        'MR', 'Station', 'Temperature_K', 'Pressure_Pa', 'Density_kg_m3',
        'Mach', 'Speed_of_Sound_m_s', 'Gamma', 'Cp_J_kgK',
        'Molecular_Weight_kg_kmol', 'Thermal_Conductivity',
        'Viscosity_Pa_s', 'Prandtl'
    ]
    properties_rows = []
    stations  = ['Injector', 'Combustor', 'Throat', 'Exit']
    prop_keys = [
        'Temperature_K', 'Pressure_Pa', 'Density_kg_m3', 'Mach',
        'Speed_of_Sound_m_s', 'Gamma', 'Cp_J_kgK',
        'Molecular_Weight_kg_kmol', 'Thermal_Conductivity',
        'Viscosity_Pa_s', 'Prandtl'
    ]
    poly_data = {st: {pk: [] for pk in prop_keys} for st in stations}

    # --- CEA evaluation (single O/F point) ---
    Eps = C.get_eps_at_PcOvPe(Pc=pc, MR=mr, PcOvPe=(pc / pe), frozen=0)
    Tc_R, Tt_R, Te_R = C.get_Temperatures(Pc=pc, MR=mr, eps=Eps,
                                           frozen=0, frozenAtThroat=0)
    Tc_K = Tc_R * RANKINE_TO_KELVIN
    Tt_K = Tt_R * RANKINE_TO_KELVIN
    Te_K = Te_R * RANKINE_TO_KELVIN

    # Chamber transport
    heat_cap_c, viscosity_c, therm_con_c, pr_c = C.get_Chamber_Transport(
        Pc=pc, MR=mr, eps=Eps, frozen=0)
    cp_c  = heat_cap_c  * c_siconv
    k_c   = therm_con_c * tc_siconv
    mu_c  = viscosity_c * vis_siconv
    MW_C  = C.get_Chamber_MolWt_gamma(Pc=pc, MR=mr, eps=Eps)
    MW_C_val = MW_C[0] if isinstance(MW_C, (list, tuple, np.ndarray)) else MW_C
    M_comb   = C.get_Chamber_MachNumber(Pc=pc, MR=mr, fac_CR=contraction_ratio)

    # Throat transport
    heat_cap_t, viscosity_t, therm_con_t, pr_t = C.get_Throat_Transport(
        Pc=pc, MR=mr, eps=Eps, frozen=0)
    cp_t  = heat_cap_t  * c_siconv
    k_t   = therm_con_t * tc_siconv
    mu_t  = viscosity_t * vis_siconv
    MW_T  = C.get_Throat_MolWt_gamma(Pc=pc, MR=mr, eps=Eps)
    MW_T_val = MW_T[0] if isinstance(MW_T, (list, tuple, np.ndarray)) else MW_T
    M_th  = 1.0

    # Exit transport
    heat_cap_e, viscosity_e, therm_con_e, pr_e = C.get_Exit_Transport(
        Pc=pc, MR=mr, eps=Eps, frozen=0)
    cp_e  = heat_cap_e  * c_siconv
    k_e   = therm_con_e * tc_siconv
    mu_e  = viscosity_e * vis_siconv
    MW_E  = C.get_exit_MolWt_gamma(Pc=pc, MR=mr, eps=Eps, frozen=0)
    MW_E_val = MW_E[0] if isinstance(MW_E, (list, tuple, np.ndarray)) else MW_E
    M_ex  = C.get_MachNumber(Pc=pc, MR=mr, eps=Eps, frozen=0)

    Pc_Pa = pc * PSI_TO_PA

    # ----- Injector (stagnation; M ~ 0) -----
    gamma_inj  = compute_gamma(cp_c, MW_C_val)
    R_spec_inj = RU_M / MW_C_val
    rho_inj    = Pc_Pa / (R_spec_inj * Tc_K) if R_spec_inj > 0 else np.nan
    a_inj      = np.sqrt(gamma_inj * R_spec_inj * Tc_K) if not np.isnan(gamma_inj) else np.nan
    properties_rows.append([mr, 'Injector', Tc_K, Pc_Pa, rho_inj,
                             0.0, a_inj, gamma_inj, cp_c, MW_C_val, k_c, mu_c, pr_c])
    for pk, val in zip(prop_keys, [Tc_K, Pc_Pa, rho_inj, 0.0, a_inj,
                                    gamma_inj, cp_c, MW_C_val, k_c, mu_c, pr_c]):
        poly_data['Injector'][pk].append(val)

    # ----- Combustor (chamber static at M = M_comb) -----
    gamma_comb = compute_gamma(cp_c, MW_C_val)
    R_spec_comb = RU_M / MW_C_val
    if not np.isnan(gamma_comb):
        factor_comb = (1.0 + 0.5 * (gamma_comb - 1.0) * M_comb**2)
        T_comb_K = Tc_K / factor_comb
        p_comb   = Pc_Pa / factor_comb ** (gamma_comb / (gamma_comb - 1.0))
    else:
        T_comb_K, p_comb = Tc_K, np.nan
    rho_comb = p_comb / (R_spec_comb * T_comb_K) if R_spec_comb > 0 else np.nan
    a_comb   = np.sqrt(gamma_comb * R_spec_comb * T_comb_K) if not np.isnan(gamma_comb) else np.nan
    properties_rows.append([mr, 'Combustor', T_comb_K, p_comb, rho_comb,
                             M_comb, a_comb, gamma_comb, cp_c, MW_C_val, k_c, mu_c, pr_c])
    for pk, val in zip(prop_keys, [T_comb_K, p_comb, rho_comb, M_comb, a_comb,
                                    gamma_comb, cp_c, MW_C_val, k_c, mu_c, pr_c]):
        poly_data['Combustor'][pk].append(val)

    # ----- Throat (M = 1) -----
    gamma_th   = compute_gamma(cp_t, MW_T_val)
    R_spec_th  = RU_M / MW_T_val
    if not np.isnan(gamma_th):
        p_th = Pc_Pa / (1.0 + 0.5 * (gamma_th - 1.0) * M_th**2) ** (gamma_th / (gamma_th - 1.0))
    else:
        p_th = np.nan
    rho_th = p_th / (R_spec_th * Tt_K) if R_spec_th > 0 else np.nan
    a_th   = np.sqrt(gamma_th * R_spec_th * Tt_K) if not np.isnan(gamma_th) else np.nan
    properties_rows.append([mr, 'Throat', Tt_K, p_th, rho_th,
                             M_th, a_th, gamma_th, cp_t, MW_T_val, k_t, mu_t, pr_t])
    for pk, val in zip(prop_keys, [Tt_K, p_th, rho_th, M_th, a_th,
                                    gamma_th, cp_t, MW_T_val, k_t, mu_t, pr_t]):
        poly_data['Throat'][pk].append(val)

    # ----- Exit -----
    gamma_ex  = compute_gamma(cp_e, MW_E_val)
    R_spec_ex = RU_M / MW_E_val
    if not np.isnan(gamma_ex):
        factor_ex = (1.0 + 0.5 * (gamma_ex - 1.0) * M_ex**2)
        T_ex_K = Te_K / factor_ex
        p_ex   = Pc_Pa / factor_ex ** (gamma_ex / (gamma_ex - 1.0))
    else:
        T_ex_K, p_ex = Te_K, np.nan
    rho_ex = p_ex / (R_spec_ex * T_ex_K) if R_spec_ex > 0 else np.nan
    a_ex   = np.sqrt(gamma_ex * R_spec_ex * T_ex_K) if not np.isnan(gamma_ex) else np.nan
    properties_rows.append([mr, 'Exit', T_ex_K, p_ex, rho_ex,
                             M_ex, a_ex, gamma_ex, cp_e, MW_E_val, k_e, mu_e, pr_e])
    for pk, val in zip(prop_keys, [T_ex_K, p_ex, rho_ex, M_ex, a_ex,
                                    gamma_ex, cp_e, MW_E_val, k_e, mu_e, pr_e]):
        poly_data['Exit'][pk].append(val)

    # --- Write turbopump_properties.csv ---
    os.makedirs(output_dir, exist_ok=True)
    properties_csv = os.path.join(output_dir, 'turbopump_properties.csv')
    with open(properties_csv, 'w', newline='') as f:
        writer = csv.writer(f)
        writer.writerow(properties_header)
        for r in properties_rows:
            writer.writerow(r)
    print(f'Bartz: saved {properties_csv}')

    # --- Write turbopump_poly_coeffs.csv ---
    poly_csv = os.path.join(output_dir, 'turbopump_poly_coeffs.csv')
    T_comb_pt = poly_data['Combustor']['Temperature_K'][0] if poly_data['Combustor']['Temperature_K'] else np.nan
    T_th_pt   = poly_data['Throat']['Temperature_K'][0]   if poly_data['Throat']['Temperature_K']   else np.nan
    T_ex_pt   = poly_data['Exit']['Temperature_K'][0]     if poly_data['Exit']['Temperature_K']     else np.nan
    T_points  = np.array([T_comb_pt, T_th_pt, T_ex_pt])
    fit_props = ['Gamma', 'Cp_J_kgK', 'Thermal_Conductivity',
                 'Viscosity_Pa_s', 'Molecular_Weight_kg_kmol']
    with open(poly_csv, 'w', newline='') as f:
        writer = csv.writer(f)
        writer.writerow(['O/F', 'Pc_psia', 'Property', 'Degree',
                         'Coefficients_high->low', 'T_points_K', 'Y_points'])
        for pk in fit_props:
            y_pts = np.array([
                poly_data['Combustor'].get(pk, [np.nan])[0],
                poly_data['Throat'].get(pk,    [np.nan])[0],
                poly_data['Exit'].get(pk,      [np.nan])[0],
            ])
            mask = (~np.isnan(T_points)) & (~np.isnan(y_pts))
            Tfit, Yfit = T_points[mask], y_pts[mask]
            if Tfit.size < 2:
                deg, coeffs = 0, []
            else:
                deg    = min(2, Tfit.size - 1)
                coeffs = np.polyfit(Tfit, Yfit, deg).tolist()
            writer.writerow([mr, pc, pk, deg, coeffs,
                             Tfit.tolist(), Yfit.tolist()])
    print(f'Bartz: saved {poly_csv}')

    return properties_csv


# =============================================================================
# Standalone entry point
# =============================================================================

if __name__ == '__main__':
    # Locate TCA_params.yaml: two levels up from this file (TCA/)
    # Path: .../TCA/HEAT TRANSFER/3D Heat transfer - Python/ -> ../../ -> .../TCA/
    _script_dir = os.path.dirname(os.path.abspath(__file__))
    _tca_yaml   = os.path.normpath(os.path.join(_script_dir, '..', '..', 'TCA_params.yaml'))
    if not os.path.exists(_tca_yaml):
        _tca_yaml = os.path.join(os.getcwd(), 'TCA_params.yaml')

    with open(_tca_yaml) as _f:
        _tca_params = yaml.safe_load(_f) or {}

    # When run standalone, write CSVs next to this script
    run(_tca_params, _script_dir)
