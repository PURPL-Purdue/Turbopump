"""
Channel height sizing - procedure 6.5.2 run with T_wg as the input (inverse of
jacket/march.py).

The march takes the channel (N, w, h) and finds T_wg. This takes N, w and a
target T_wg and finds the channel height h at every station. The book allows
d = d(x) ("there is a tendency to want to decrease d in the throat region"),
and this picks it station by station.

At each station, with T_wg fixed at the target, steps 5-6 are explicit:
    step 5  h_g at T_wg (Bartz, Eq. 6.13),  q = h_g (T_r - T_wg)
    step 6  T_wl = T_wg - q t_w / k_w
and steps 7-8 become the condition the channel has to meet:
    h_l(h) phi(h) = h_l_req = q / (T_wl - T_l)
h_l from the regen.yaml correlation, phi the fin factor (1 with fins off,
heat_transfer/fins.py with jacket.fin_conduction on). The root-find varies h
only; q, T_wl and h_l_req don't depend on it.

Without fins a shorter channel is faster and has a smaller D_h, so h_l rises
monotonically as h falls. With fins a shorter channel also loses fin area, so
h_l phi can peak inside the range. The search scans N_SCAN heights, takes the
tallest one that meets the target and refines it with brentq. That is the
TALLEST channel that holds T_wg at the target, i.e. the one with the least
pressure drop.

If the target can't be met inside [h_min, h_max] the height is clamped:
    h_max  station runs cooler than the target even with the tallest channel
    h_min  even the shortest channel can't hold the target (T_wg runs hot)
    best   no height holds the target; the coolest height (the peak of
           h_l phi) is used and T_wg runs hot. Only possible with fins on.
Step 10 (coolant temperature update) then uses the actual heat flux at the
chosen height, so the next station sees the right coolant temperature.

Run the result through jacket.march.march() so every temperature and pressure
reported comes from the standard procedure.
"""
import numpy as np
from scipy.optimize import brentq

from heat_transfer.gas_side import bartz_at_station
from jacket.geometry import channel_geometry, set_station_height
from jacket.march import solve_station, advance_coolant, coolant_order, coolant_area_factor, inlet_distance

N_SCAN = 60   # heights scanned per station, log-spaced over [h_min, h_max]


def height_for_target(i, T_l, gas, chan, mdot_chan, w, T_target, h_min, h_max, jk):
    """
    Channel height at station i that makes T_wg = T_target. chan supplies the
    station's r, N and ds for the fin factor.
    Returns (h, status), status 'ok', 'h_max', 'h_min' or 'best'.
    """
    T_r = gas["T_r"][i]
    if T_target >= T_r:
        return h_max, "h_max"                        # gas can't heat the wall that far

    h_g = bartz_at_station(gas, i, T_target)
    q = h_g * (T_r - T_target)                       # step 5
    T_wl = T_target - q * jk.t_w / jk.k_w            # step 6
    if T_wl <= T_l:
        return h_min, "h_min"                        # wall conduction alone uses up the budget
    h_l_req = q / (T_wl - T_l)                       # steps 7-8 as a requirement
    x_in = inlet_distance(i, chan, jk)               # entrance terms, inf if off

    def excess(h):
        """Coolant side's effective h_l (h_l phi) minus the required one, at height h."""
        A_c = w * h
        D_h = 4 * A_c / (2 * (w + h))                # Eq. 6.18
        h_l = jk.hl(T_l, T_wl, mdot_chan, A_c, D_h, jk.fluid, x_in)
        phi, _ = coolant_area_factor(h_l, w, h, i, chan, jk)
        return h_l * phi - h_l_req

    # Scan the range first. Without fins excess only falls as h rises, but with
    # fins a taller channel also adds fin area, so it can peak inside the range.
    # The scan finds the tallest height that still meets the target, and
    # brentq refines it between the two scan points around it.
    hs = np.geomspace(h_min, h_max, N_SCAN)
    ex = np.array([excess(h) for h in hs])
    if ex[-1] >= 0:
        return h_max, "h_max"
    ok = np.nonzero(ex >= 0)[0]
    if ok.size == 0:                                 # no height meets the target
        k = int(np.argmax(ex))                       # the coolest one instead
        return hs[k], ("h_min" if k == 0 else "best")
    j = ok[-1]                                       # tallest scan point that meets it
    h: float = brentq(excess, hs[j], hs[j + 1], xtol=1e-7)   # type: ignore[assignment]  (scipy stub: may be a tuple with full_output)
    return h, "ok"


def size_heights(x, r, gas, N, mdot_chan, w, T_target, h_min, h_max, jk):
    """
    March in coolant order picking h at each station, updating the coolant
    temperature with the actual heat flux at that height (step 10, Eq. 6.39).
    Returns (h, status), both (m,) indexed by geometry station.
    """
    m = len(x)
    chan = channel_geometry(x, r, N, w, h_max, jk.t_w)
    status = np.empty(m, dtype=object)
    T_l, p_l = jk.T_l_in, jk.p_l_in

    for i in coolant_order(m, jk):
        h, status[i] = height_for_target(i, T_l, gas, chan, mdot_chan, w, T_target, h_min, h_max, jk)
        set_station_height(chan, i, h)
        _, q, *_ = solve_station(i, T_l, gas, chan, mdot_chan, jk)
        T_l, p_l, _, _ = advance_coolant(i, q, T_l, p_l, chan, mdot_chan, jk)

    return chan["h"].copy(), status
