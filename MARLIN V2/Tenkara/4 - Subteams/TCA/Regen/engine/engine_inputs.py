"""
Engine definition, read from regen.yaml (engine, chamber, contour sections).

Every script imports this so they all describe the same engine. Values are
set in regen.yaml, not here; this file only names them and converts to SI.
"""
import math
import numpy as np
from util import config as c

# ============================== UNIT CONVERSIONS ===========================
PSI_TO_PA = 6894.76
PSI_TO_BAR = 0.0689476
N_TO_LBF = 0.224809
M_TO_IN = 39.37
G0 = 9.80665

_e, _ch, _co = c.engine, c.chamber, c.contour

# ============================== PROPELLANTS ================================
reac_names = list(_e["reactants"])
T_reactant = np.array(_e["T_reactant_K"], dtype=float)
fuel_weights = np.array(_e["fuel_weights"], dtype=float)
oxidant_weights = np.array(_e["oxidant_weights"], dtype=float)
fuel_dH_vap = _e["fuel_dH_vap_kJ_per_mol"] * 1e3 / (_e["fuel_molar_mass_g_per_mol"] * 1e-3)  # J/kg

# ============================== OPERATING POINT ============================
F_target_lbf = _e["F_target_lbf"]
p0_psi = _e["p0_psi"]
pa_psi = _e["pa_psi"]
OF = _e["OF"]
m_dot_pump_max = _e["m_dot_pump_max_kg_s"]
c_star_eff = _e["c_star_eff"]
c_thrust_eff = _e["c_thrust_eff"]

# ============================== CHAMBER ====================================
contraction_ratio = _ch["contraction_ratio"]
lstar = _ch["lstar_m"]

# ============================== CONTOUR ====================================
theta_conv_deg = _co["theta_conv_deg"]
R_up_over_Rt = _co["R_up_over_Rt"]
R_down_over_Rt = _co["R_down_over_Rt"]
bell_percent = int(_co["bell_percent"])

# ============================== SI VERSIONS ================================
p0_Pa = p0_psi * PSI_TO_PA
p0_bar = p0_psi * PSI_TO_BAR
pa_Pa = pa_psi * PSI_TO_PA
pa_bar = pa_psi * PSI_TO_BAR
theta_conv = math.radians(theta_conv_deg)
