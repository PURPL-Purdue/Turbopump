"""
Regen jacket march - Heister et al., Rocket Propulsion, 6.5.2. Forward problem:
given the channel geometry, find the wall and coolant temperatures and the
coolant pressure. At every station the heat leaving the gas equals the heat
through the wall equals the heat into the coolant (Eq. 6.29):

    q = h_g (T_r - T_wg) = k_w (T_wg - T_wl) / t_w = h_l (T_wl - T_l)

Procedure 6.5.2, step by step:
    1   system design, inflow, geometry      regen.yaml -> Jacket, build_case()
    2   N, m_dot_chan (Eq. 6.31)              the run script
    3   step to the next segment              march() loop
    4   guess T_wg                            brentq on residual()
    5   h_g (Bartz, Eq. 6.13), q = h_g (T_r - T_wg)          station_fluxes()
    6   T_wl from wall conduction                            station_fluxes()
    7   h_l (jk.hl, the regen.yaml correlation)              station_fluxes()
        times the fin factor phi if fins are on             coolant_area_factor()
    8   q = h_l phi (T_wl - T_l), compare with step 6        residual()
    9   iterate until they agree                             brentq
    10  T_l update (Eq. 6.39), p_l update (Eq. 6.36/6.37)    advance_coolant()
    11  next segment                                         march() loop

Only the book's equations are used, plus one optional extension: fin conduction
through the lands (jacket.fin_conduction, heat_transfer/fins.py, Blong 2024
Eqs. 21-27), which the book mentions but leaves out of the procedure. When fin_conduction is False, phi = 1
and the procedure is the book's. Not included: the momentum pressure loss from
coolant density change, film cooling, and combustion efficiency on the gas
temperature. Step 12 (structural check) is left to FEA.

Indexing: every array is indexed by geometry station, i = 0 at the injector
face, i = m-1 at the nozzle exit. The march walks the stations in coolant order
(m-1 -> 0 for counterflow) and stores results at the geometry index.

Pure physics: no config reads. jacket/jacket_inputs.py builds the Jacket from regen.yaml.
"""
import math
from dataclasses import dataclass
from types import ModuleType
from typing import Callable

import numpy as np
from scipy.optimize import brentq

from engine.cea_interface import size_engine
from engine.contour import build_contour
from heat_transfer.gas_side import gas_environment, bartz_at_station
from heat_transfer.friction import friction_factor
from heat_transfer import fins


@dataclass(frozen=True)
class Jacket:
    """Everything the march needs besides the geometry and the hot gas."""
    t_w: float            # m      hot wall thickness
    k_w: float            # W/m-K  wall conductivity
    p_l_in: float         # Pa     coolant pressure entering the jacket
    T_l_in: float         # K      coolant temperature entering the jacket
    eps_rough: float      # m      channel wall roughness (0 = smooth-wall fits)
    counterflow: bool     # coolant enters at the nozzle exit
    fluid: ModuleType     # coolant properties (coolant/ipa_props.py)
    hl: Callable          # coolant-side correlation (heat_transfer/coolant_side.py)
    hl_name: str          # its name in regen.yaml, for reports
    fins: bool = False    # lands as fins (heat_transfer/fins.py); False = book procedure
    entrance: bool = False  # h_l depends on distance from the coolant inlet; False = fully developed


def build_case(cfg, jk, m_stations):
    """Step 1: coolant property table, engine sizing, contour, gas conditions."""
    jk.fluid.init(jk.p_l_in)    # check this
    eng = size_engine(cfg)
    x, r, dims = build_contour(eng, m_stations, cfg)  # check this
    gas = gas_environment(eng, x, r, dims, cfg)
    return eng, x, r, dims, gas


def coolant_order(m, jk):
    """Station indices in the order the coolant meets them."""
    return np.arange(m)[::-1] if jk.counterflow else np.arange(m)


# ============================== STATION BALANCE (steps 4-9) ================

def inlet_distance(i, chan, jk):
    """
    Distance the coolant has travelled along the channel to station i [m], for
    the entrance terms of the coolant correlation. inf with entrance effects
    off (fully developed). The inlet station gets half its own segment length,
    not 0, where the developing-flow terms would be infinite (our choice).
    """
    if not jk.entrance:
        return math.inf
    s = chan["s"]
    dist = (s[-1] - s[i]) if jk.counterflow else (s[i] - s[0])
    return max(dist, 0.5 * chan["ds"][i])


def coolant_area_factor(h_l, w, h, i, chan, jk):
    """
    Fin factor phi at station i for a channel w x h, and the fin efficiency eta:
    q_l = h_l phi (T_wl - T_l). Fins off: (1, nan), the book's procedure.
    w and h are passed so the channel solver can try heights not yet in chan.
    """
    if not jk.fins:
        return 1.0, float("nan")
    return fins.fin_factor(h_l, w, h, chan["r"][i], jk.t_w, chan["N"], jk.k_w, chan["ds"][i])


def station_fluxes(T_wg, i, T_l, gas, chan, mdot_chan, jk):
    """
    For a guessed T_wg at station i, the heat flux two ways:
        step 5  q_g = h_g (T_r - T_wg)
        step 6  T_wl = T_wg - q_g t_w / k_w
        step 7  h_l at the coolant temperature of the previous segment,
                fin factor phi (1 with fins off)
        step 8  q_l = h_l phi (T_wl - T_l)
    Returns (q_g, q_l, h_g, h_l, T_wl, phi, eta).
    """
    h_g = bartz_at_station(gas, i, T_wg)
    q_g = h_g * (gas["T_r"][i] - T_wg)
    T_wl = T_wg - q_g * jk.t_w / jk.k_w
    h_l = jk.hl(T_l, T_wl, mdot_chan, chan["A_c"][i], chan["D_h"][i], jk.fluid,
                inlet_distance(i, chan, jk))
    phi, eta = coolant_area_factor(h_l, chan["w"][i], chan["h"][i], i, chan, jk)
    q_l = h_l * phi * (T_wl - T_l)
    return q_g, q_l, h_g, h_l, T_wl, phi, eta


def residual(T_wg, i, T_l, gas, chan, mdot_chan, jk):
    """Step 8-9: zero when the gas-side and coolant-side fluxes agree."""
    q_g, q_l, *_ = station_fluxes(T_wg, i, T_l, gas, chan, mdot_chan, jk)
    return q_g - q_l


def solve_station(i, T_l, gas, chan, mdot_chan, jk):
    """
    Steps 4-9. The root is bracketed: at T_wg = T_l the gas side gives q > 0 and
    the coolant side q < 0; at T_wg = T_r the reverse (phi > 0 keeps the signs).
    Returns (T_wg, q, h_g, h_l, T_wl, phi, eta).
    """
    args = (i, T_l, gas, chan, mdot_chan, jk)
    T_wg = brentq(residual, T_l, gas["T_r"][i], args=args)
    q, _, h_g, h_l, T_wl, phi, eta = station_fluxes(T_wg, *args)
    return T_wg, q, h_g, h_l, T_wl, phi, eta


def advance_coolant(i, q, T_l, p_l, chan, mdot_chan, jk):
    """
    Step 10 over segment i:
        Eq. 6.39  T_l,i+1 = T_l,i + q dA / (m_dot cp),  dA = b ds
        Eq. 6.36  p_l,i+1 = p_l,i - 4 cf (ds / d) rho v^2 / 2
    Liquid properties updated at the new temperature before the pressure step.
    Returns (T_l_new, p_l_new, v, Re).
    """
    fluid = jk.fluid
    T_new = T_l + q * chan["b"][i] * chan["ds"][i] / (mdot_chan * fluid.cp(T_l))
    rho = fluid.rho(T_new)
    v = mdot_chan / (rho * chan["A_c"][i])
    Re = rho * v * chan["D_h"][i] / fluid.mu(T_new)
    cf = friction_factor(Re, jk.eps_rough / chan["D_h"][i])
    p_new = p_l - 4 * cf * (chan["ds"][i] / chan["D_h"][i]) * 0.5 * rho * v**2
    return T_new, p_new, v, Re


# ============================== MARCH (steps 3-11) =========================

def march(x, r, gas, chan, mdot_chan, jk):
    """
    Walk the passage in coolant order. Returns a dict of (m,) arrays indexed by
    geometry station (T_l, p_l are the values entering each segment), plus the
    jacket outlet state T_l_out, p_l_out.
    """
    m = len(x)
    out: dict = {k: np.full(m, np.nan) for k in
           ("q", "h_g", "h_l", "phi", "eta_fin", "T_wg", "T_wl", "T_l", "p_l", "v_l", "Re_l")}
    T_l, p_l = jk.T_l_in, jk.p_l_in

    for i in coolant_order(m, jk):
        T_wg, q, h_g, h_l, T_wl, phi, eta = solve_station(i, T_l, gas, chan, mdot_chan, jk)
        out["q"][i], out["h_g"][i], out["h_l"][i] = q, h_g, h_l
        out["phi"][i], out["eta_fin"][i] = phi, eta
        out["T_wg"][i], out["T_wl"][i] = T_wg, T_wl
        out["T_l"][i], out["p_l"][i] = T_l, p_l
        T_l, p_l, out["v_l"][i], out["Re_l"][i] = advance_coolant(
            i, q, T_l, p_l, chan, mdot_chan, jk)

    # Step 11: at the outlet manifold the dynamic pressure is lost, so the
    # manifold pressure is the static pressure of the last segment - p_l as is.
    out["T_l_out"], out["p_l_out"] = T_l, p_l
    return out
