"""
Jacket design inputs shared by both jacket solvers, read from regen.yaml
(jacket section, plus contour.stations). Step 1 of Heister 6.5.2.

JACKET is what the physics (jacket/march.py, jacket/sizing.py) takes. Channel
count / width / height are per solver (wall_temp_solver and channel_solver
sections of regen.yaml), read by the scripts in scripts/.
"""
import coolant
from heat_transfer import coolant_side
from jacket.march import Jacket
from util import config as c
from engine import engine_inputs as ei

_j = c.jacket

# --- system design ---
counterflow  = bool(_j["counterflow"])
m_stations   = int(c.contour["stations"])

# --- coolant inlet (fuel, all of it through the jacket - Eq. 6.31) ---
p_l_in       = _j["p_l_in_psi"] * ei.PSI_TO_PA   # Pa
T_l_in       = _j["T_l_in_K"]                    # K

# --- hot-gas wall ---
t_w          = _j["t_w_mm"] * 1e-3               # m
k_wall       = _j["k_wall_W_mK"]                 # W/m-K
T_wg_allow   = _j["T_wg_allow_K"]                # K

# --- channel count rule (step 2): most channels that fit at the throat ---
land_min     = _j["land_min_mm"] * 1e-3          # m

# --- roughness for the friction factor (Eq. 6.36) ---
eps_rough    = _j["eps_rough_mm"] * 1e-3         # m

# --- models: coolant properties and coolant-side correlation ---
coolant_name = _j["coolant"]
htc_name     = _j["coolant_htc"]
fins         = bool(_j["fin_conduction"])         # lands as fins, heat_transfer/fins.py
entrance     = bool(_j["entrance_effects"])       # developing flow near the coolant inlet

JACKET = Jacket(t_w=t_w, k_w=k_wall, p_l_in=p_l_in, T_l_in=T_l_in, eps_rough=eps_rough,
                counterflow=counterflow, fluid=coolant.get(coolant_name),
                hl=coolant_side.get(htc_name), hl_name=htc_name, fins=fins, entrance=entrance)
