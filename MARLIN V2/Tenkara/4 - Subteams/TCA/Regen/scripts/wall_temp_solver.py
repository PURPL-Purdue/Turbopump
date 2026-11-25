"""
Wall temperature solver (analysis): give the channels, get the wall and coolant
temperatures and the coolant pressure along the jacket. Heister 6.5.2, forward.

Channels come from regen.yaml, wall_temp_solver section:
    geometry_file: null   constant N, w, h from the same section
    geometry_file: path   N, w(x), h(x) from a channel geometry file, e.g.
                          output/channel_geometry.csv from scripts/channel_solver.py,
                          or a copy of it with edited heights. It has to match
                          the current contour (checked).

Use it for any geometry the channel solver didn't produce directly: rounded or
smoothed heights, printer limits, roughness or inlet-temperature what-ifs.

Run:  python scripts/wall_temp_solver.py
"""
import sys, pathlib
sys.path.insert(0, str(pathlib.Path(__file__).resolve().parents[1]))

from engine import engine_inputs as ei
from jacket import jacket_inputs as ji
from jacket.geometry import channel_geometry, read_channel_geometry
from jacket.march import build_case, march
from util import config as c
from util import report as rep


def load_channels(x, r, sec):
    """(N, w, h, description) from the wall_temp_solver section of regen.yaml."""
    if sec.get("geometry_file"):
        path = pathlib.Path(sec["geometry_file"])
        if not path.is_absolute():
            path = c.ROOT / path
        N, w, h = read_channel_geometry(path, x, r)
        return N, w, h, f"from {sec['geometry_file']}"
    return (int(sec["N_channels"]), sec["w_chan_mm"] * 1e-3, sec["h_chan_mm"] * 1e-3,
            "constant, from regen.yaml")


def run(cfg=ei, jk=ji.JACKET):
    eng, x, r, dims, gas = build_case(cfg, jk, ji.m_stations)        # step 1
    rep.print_contour_dims(dims)

    N, w, h, source = load_channels(x, r, c.wall_temp_solver)
    chan = channel_geometry(x, r, N, w, h, jk.t_w)
    mdot_chan = eng["m_dot_fuel"] / N                                 # step 2, Eq. 6.31
    mm = 1e3
    print(f"\n--- Channels ({source}) ---\n  N = {N}, w {chan['w'].min()*mm:.3f}-{chan['w'].max()*mm:.3f} mm, "
          f"h {chan['h'].min()*mm:.3f}-{chan['h'].max()*mm:.3f} mm, t_w = {jk.t_w*mm:.2f} mm, "
          f"land {chan['land'].min()*mm:.2f}-{chan['land'].max()*mm:.2f} mm, "
          f"mdot/channel = {mdot_chan*1e3:.1f} g/s\n  {jk.fluid.NAME} coolant, h_l: {jk.hl_name}, fins {'on' if jk.fins else 'off'}")

    out = march(x, r, gas, chan, mdot_chan, jk)                       # steps 3-11

    print(f"\n=============== WALL TEMPERATURES ({jk.fluid.NAME}, O/F {cfg.OF}) ===============")
    rep.print_summary(x, dims, out, jk, ji.T_wg_allow, "allowable")
    rep.print_checks(out, jk, cfg)
    rep.print_station_table(x, r, dims, chan, gas, out)
    rep.plot_run(x, r, dims, chan, out,
                 f"Wall temperatures ({jk.fluid.NAME}, O/F {cfg.OF}) - {N} ch, channels {source}, "
                 f"t_w {jk.t_w*mm:.2f} mm, " + ("counterflow" if jk.counterflow else "coflow"),
                 ji.T_wg_allow, "allowable")
    return out


if __name__ == "__main__":
    run()
