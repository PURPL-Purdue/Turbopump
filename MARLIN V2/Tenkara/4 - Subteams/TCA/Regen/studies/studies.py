"""
Parameter studies on the engine definition in engine/engine_inputs.py.

Currently: CEA chamber temperature vs O/F at the chamber pressure in engine_inputs.
Every other input (propellants, injection temps, Pc) is taken from engine_inputs.
"""
import sys, pathlib, types
import numpy as np
import matplotlib.pyplot as plt

sys.path.insert(0, str(pathlib.Path(__file__).resolve().parents[1]))
from engine import engine_inputs as ei
from engine import cea_interface as ci
from util import provenance as prov

# ============================== SWEEP ======================================
OF_min, OF_max, n_pts = 0.8, 2.5, 60      # O/F range to plot


def with_OF(OF):
    """Copy of engine_inputs with only OF changed."""
    cfg = types.SimpleNamespace(**{k: v for k, v in vars(ei).items() if not k.startswith("__")})
    cfg.OF = OF
    return cfg


def chamber_temp_vs_OF():
    OFs = np.linspace(OF_min, OF_max, n_pts)
    T0 = np.array([ci._solve(with_OF(of)).T[0] for of in OFs])   # CEA station 0 = chamber
    return OFs, T0


if __name__ == "__main__":
    OFs, T0 = chamber_temp_vs_OF()
    i_max = np.argmax(T0)

    fig, ax = plt.subplots(figsize=(8, 5))
    ax.plot(OFs, T0, lw=2)
    ax.plot(OFs[i_max], T0[i_max], "o", color="C3",
            label=f"max {T0[i_max]:.0f} K at O/F = {OFs[i_max]:.2f}")
    ax.axvline(ei.OF, color="gray", ls="--", label=f"current O/F = {ei.OF}")
    ax.set_xlabel("O/F (mass)")
    ax.set_ylabel("Chamber temperature [K]")
    ax.set_title(f"{ei.reac_names[0]} / {ei.reac_names[1]}   Pc = {ei.p0_psi} psi")
    ax.grid(alpha=0.3)
    ax.legend()
    fig.tight_layout(rect=(0, 0.02, 1, 1))
    prov.stamp_figure(fig)
    plt.show()
