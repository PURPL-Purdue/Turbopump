"""
Engine sizing sanity check: CEA sizing results and the hot-gas state at a few
area ratios.

Run:  python scripts/engine_sizing.py
"""
import sys, pathlib
sys.path.insert(0, str(pathlib.Path(__file__).resolve().parents[1]))

import numpy as np

from engine import engine_inputs as ei
from engine.cea_interface import size_engine, gas_side_states


def run(cfg=ei):
    eng = size_engine(cfg)
    print(f"\n{cfg.reac_names[0]} / {cfg.reac_names[1]}   Pc = {cfg.p0_psi} psi   O/F = {cfg.OF}")
    for k in ["m_dot_fuel", "m_dot", "c_star", "Cf", "eps", "Rt", "Rc", "Re", "F_lbf", "Isp", "T0"]:
        print(f"  {k:<12}{eng[k]:>12.4g}")

    # Coarse station sweep: chamber 5 -> throat -> exit eps
    ar = np.concatenate([np.linspace(5, 1.2, 5), [1.0], np.linspace(1.2, eng["eps"], 5)])
    sup = np.concatenate([np.zeros(5, bool), [False], np.ones(5, bool)])
    g = gas_side_states(ar, sup, cfg=cfg)
    print(f"\n  {'A/At':>7}{'M':>8}{'T [K]':>9}{'p [bar]':>9}{'mu [Pa.s]':>11}{'Pr':>7}")
    for i in range(len(ar)):
        print(f"  {ar[i]:7.3f}{g['M'][i]:8.3f}{g['T'][i]:9.1f}{g['p'][i]/1e5:9.2f}"
              f"{g['mu'][i]:11.3e}{g['Pr'][i]:7.3f}")
    return eng


if __name__ == "__main__":
    run()
