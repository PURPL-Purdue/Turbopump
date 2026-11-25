"""
Channel solver (design): give the channel width and a target gas-side wall
temperature, get the channel height h(x). Heister 6.5.2 run inverse; the method
is in jacket/sizing.py. The result is then run through the standard march, so
every temperature and pressure reported comes from the forward procedure.

Inputs: regen.yaml, channel_solver section.
Writes to output/:
    channel_geometry.csv   always: N, w(x), h(x) on this wall. Input for
                           scripts/wall_temp_solver.py (geometry_file).
    contour.dxf            if write_dxf: wall, channel floor, channel top for NX
    channel_heights.csv    if write_csv: plus NX point files (radial offsets)
    inputs_used.yaml       what made all of the above

Run:  python scripts/channel_solver.py
"""
import sys, pathlib
sys.path.insert(0, str(pathlib.Path(__file__).resolve().parents[1]))

import numpy as np

from engine import engine_inputs as ei
from jacket import jacket_inputs as ji
from jacket.geometry import channel_geometry, channels_from_land
from jacket.march import build_case, march
from jacket.sizing import size_heights
from util import config as c
from util import export
from util import provenance as prov
from util import report as rep

_cs = c.channel_solver
w_chan       = _cs["w_chan_mm"] * 1e-3      # m  channel width (constant along the length)
N_channels   = _cs["N_channels"]            # None = as many as fit at the throat with land >= land_min
T_wg_target  = _cs["T_wg_target_K"] or ji.T_wg_allow   # K  gas-side wall temperature to design to
h_min        = _cs["h_min_mm"] * 1e-3       # m  height search range
h_max        = _cs["h_max_mm"] * 1e-3       # m


def run(cfg=ei, jk=ji.JACKET):
    eng, x, r, dims, gas = build_case(cfg, jk, ji.m_stations)        # step 1
    rep.print_contour_dims(dims)
    N = N_channels or channels_from_land(eng["Rt"], w_chan, jk.t_w, ji.land_min)   # step 2
    mdot_chan = eng["m_dot_fuel"] / N                                 # Eq. 6.31

    h, status = size_heights(x, r, gas, N, mdot_chan, w_chan, T_wg_target, h_min, h_max, jk)
    chan = channel_geometry(x, r, N, w_chan, h, jk.t_w)
    out = march(x, r, gas, chan, mdot_chan, jk)                       # standard procedure

    out_dir = c.ROOT / "output"
    print()
    export.write_channel_geometry(out_dir / "channel_geometry.csv", x, r, N, chan)
    if _cs["write_csv"]:
        export.write_height_csv(out_dir / "channel_heights.csv", x, r, chan, jk.t_w)
    if _cs.get("write_dxf"):
        export.write_contour_dxf(out_dir / "contour.dxf", dims, t_w=jk.t_w, x=x, r=r, h=chan["h"])
    prov.write_inputs_used(out_dir)

    mm = 1e3
    i_t = rep.throat_index(x, dims)
    print(f"\n=============== CHANNEL SIZING ({jk.fluid.NAME}, O/F {cfg.OF}) ===============")
    print(f"  Engine Isp {eng['Isp']:.1f} s")
    print(f"  Thrust {eng['F_N']:.1f} N ({eng['F_lbf']:.1f} lbf)")
    print(f"  {N} channels x {w_chan*mm:.2f} mm wide, t_w {jk.t_w*mm:.2f} mm, "
          f"target T_wg {T_wg_target:.0f} K, h search {h_min*mm:.2f}-{h_max*mm:.2f} mm, h_l: {jk.hl_name}, fins {'on' if jk.fins else 'off'}")
    print(f"  h {np.min(chan['h'])*mm:.3f} - {np.max(chan['h'])*mm:.3f} mm   "
          f"(throat {chan['h'][i_t]*mm:.3f} mm), land {chan['land'].min()*mm:.2f} - "
          f"{chan['land'].max()*mm:.2f} mm")
    rep.print_summary(x, dims, out, jk, T_wg_target, "target")
    n_lo, n_hi = np.sum(status == "h_min"), np.sum(status == "h_max")
    if n_lo:
        print(f"  ** {n_lo} stations at h_min: target not reachable there, T_wg runs hot **")
    n_best = np.sum(status == "best")
    if n_best:
        print(f"  ** {n_best} stations at the coolest height (peak of h_l*phi): target not "
              f"reachable there, T_wg runs hot **")
    if n_hi:
        print(f"  {n_hi} stations at h_max: cooler than the target even at max height")
    rep.print_checks(out, jk, cfg)
    rep.print_station_table(x, r, dims, chan, gas, out, status)
    rep.plot_run(x, r, dims, chan, out,
                 f"Channel sizing - {N} ch, w {w_chan*mm:.2f} mm, T_wg target {T_wg_target:.0f} K",
                 T_wg_target, "target")
    return h, status, out


if __name__ == "__main__":
    run()
