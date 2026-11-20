"""Plot a Gas Generator field .npy file as two contour plots:

  1. Total mass flow rate (filled) with O/F ratio lines (black).
  2. Chamber pressure (filled) with minimum-stiffness lines (black), and red diagonal
     hatching over the region where minimum stiffness < STIFF_THRESHOLD.

File format: array of shape (3, n_ox, n_fuel) = field[layer, ox, fuel] (rows are ox feed
pressure, columns are fuel feed pressure), layers in the order mdot, Pc, OF.
"""
import argparse
from pathlib import Path
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.patches import Patch

# Select .npy to plot
NPY_PATH = "GasGenerator_Field_10pct.npy"

# Must match the field generator
P_MIN = 100.0    # psi
P_MAX = 1000.0   # psi

# O/F line spacing
OF_STEP_LOW = 0.2
OF_STEP_HIGH = 0.4
OF_BREAK = 1.0
OF_TRIM_PERCENTILES = (0, 99)

# Stiffness lines and mask
STIFF_STEP = 0.1
STIFF_THRESHOLD = 0.2

DEFAULT_LAYERS = ("mdot", "Pc", "OF")

def load_field(path: Path) -> np.ndarray:
    field = np.load(path, allow_pickle=False)
    if field.ndim != 3:
        raise ValueError(f"{path.name}: expected shape (layers, N, N), got {field.shape}")
    return field

def style(ax, title):
    ax.set_xlabel("Oxidizer pressure (psi)")
    ax.set_ylabel("Fuel pressure (psi)")
    ax.set_title(title)
    ax.set_xlim(P_MIN, P_MAX)
    ax.set_ylim(P_MIN, P_MAX)

def no_data(ax):
    ax.text(0.5, 0.5, "Not enough valid points to contour", ha="center",
            va="center", transform=ax.transAxes)

def of_levels(lo, hi):
    """Adaptive O/F contour levels: OF_STEP_LOW up to OF_BREAK, then OF_STEP_HIGH above it.
    Only levels inside [lo, hi] are returned."""
    low = OF_STEP_LOW * np.arange(0, int(round(OF_BREAK / OF_STEP_LOW)) + 1)       # 0, 0.2, ... 1.0
    n_high = int(np.ceil((hi - OF_BREAK) / OF_STEP_HIGH)) if hi > OF_BREAK else 0
    high = OF_BREAK + OF_STEP_HIGH * np.arange(1, n_high + 1)                      # 1.4, 1.8, ...
    levels = np.concatenate([low, high])
    return levels[(levels >= lo) & (levels <= hi)]

def plot_mdot_with_of(fig, ax, X, Y, mdot, of):
    """Total mdot (filled) with black O/F lines."""
    good = ~np.isnan(mdot)
    if good.sum() < 4 or np.ptp(mdot[good]) == 0:
        no_data(ax)
        style(ax, "Total mdot (lb/s) with O/F lines")
        return

    cf = ax.contourf(X, Y, np.ma.masked_invalid(mdot), levels=12, cmap="viridis")
    fig.colorbar(cf, ax=ax, label="Total mdot (lb/s)")

    good_of = ~np.isnan(of)
    if good_of.sum() >= 4 and np.ptp(of[good_of]) > 0:
        print(f"O/F range: {np.nanmin(of):.2f} to {np.nanmax(of):.2f}")
        # Always trim to the chosen percentile band so a few extreme O/F points
        # (e.g. near-zero fuel flow) can't flood the plot with lines.
        p_lo, p_hi = np.nanpercentile(of, OF_TRIM_PERCENTILES)
        levels = of_levels(p_lo, p_hi)
        if len(levels):
            ocs = ax.contour(X, Y, np.ma.masked_invalid(of), levels=levels,
                             colors="k", linestyles="--", linewidths=0.9)
            ax.clabel(ocs, fontsize=7, fmt="%.1f")
            ax.legend(handles=[Line2D([], [], color="k", ls="--", lw=1.0, label="O/F")],
                      loc="upper left", fontsize=8)
    style(ax, "Total mdot (lb/s) with O/F lines")

def plot_pc_with_stiffness(fig, ax, X, Y, pc):
    """Chamber pressure (filled) with dashed black min-stiffness lines, a solid red line at
    the threshold, and red hatching below it."""
    good = ~np.isnan(pc)
    if good.sum() < 4 or np.ptp(pc[good]) == 0:
        no_data(ax)
        style(ax, "Chamber pressure (psi) with min stiffness")
        return

    cf = ax.contourf(X, Y, np.ma.masked_invalid(pc), levels=12, cmap="viridis")
    fig.colorbar(cf, ax=ax, label="Chamber pressure (psi)")

    # Minimum stiffness: the lower feed pressure against the shared chamber pressure
    p_up_min = np.minimum(X, Y)
    with np.errstate(divide="ignore", invalid="ignore"):
        stiff = np.where(pc > 0, (p_up_min - pc) / pc, np.nan)
    stiff = np.ma.masked_invalid(stiff)
    s_lo, s_hi = np.nanmin(stiff), np.nanmax(stiff)
    print(f"Min stiffness range: {s_lo:.3f} to {s_hi:.3f}")

    handles = []

    # All stiffness lines at STIFF_STEP spacing (dashed black)
    levels = np.arange(STIFF_STEP, s_hi + STIFF_STEP, STIFF_STEP)
    levels = levels[(levels >= s_lo) & (levels <= s_hi)]
    if len(levels):
        scs = ax.contour(X, Y, stiff, levels=levels,
                         colors="k", linestyles="dashed", linewidths=0.9)
        # the threshold level is labeled on its own red line below
        ax.clabel(scs, levels=[lv for lv in levels if not np.isclose(lv, STIFF_THRESHOLD)],
                  fontsize=7, fmt="%.1f")
        handles.append(Line2D([], [], color="k", ls="--", lw=0.9, label="Min stiffness"))

    # Threshold line (solid red) on top of the dashed lines
    if s_lo < STIFF_THRESHOLD < s_hi:
        ths = ax.contour(X, Y, stiff, levels=[STIFF_THRESHOLD],
                         colors="red", linestyles="solid", linewidths=0.9)
        ax.clabel(ths, fontsize=7, fmt="%.1f")
        handles.append(Line2D([], [], color="red", lw=0.9, label=f"Stiffness = {STIFF_THRESHOLD:g}"))

    # Red diagonal hatching where stiffness < threshold; the fill is left clear so Pc stays readable
    if s_lo < STIFF_THRESHOLD:
        with plt.rc_context({"hatch.color": "red", "hatch.linewidth": 0.9}):
            ax.contourf(X, Y, stiff, levels=[s_lo - 1.0, STIFF_THRESHOLD],
                        colors="none", hatches=["///"])
        handles.append(Patch(facecolor="none", edgecolor="red", hatch="///",
                             label=f"Stiffness < {STIFF_THRESHOLD:g}"))

    if handles:
        ax.legend(handles=handles, loc="upper left", fontsize=8)
    style(ax, "Chamber pressure (psi) with min stiffness")

def plot_field(path: Path, field: np.ndarray, names):
    n_layers, n_ox, n_fuel = field.shape       # field[layer, ox, fuel]: rows = ox, cols = fuel
    p_ox = np.linspace(P_MIN, P_MAX, n_ox)     # rows -> ox pressure
    p_fu = np.linspace(P_MIN, P_MAX, n_fuel)   # cols -> fuel pressure
    # indexing="ij" keeps X and Y the same shape/order as each layer: X[i, j] = ox, Y[i, j] = fuel
    X, Y = np.meshgrid(p_ox, p_fu, indexing="ij")

    mdot, pc, of = (field[names.index(n)] for n in DEFAULT_LAYERS)

    fig, axes = plt.subplots(1, 2, figsize=(12, 5), squeeze=False)
    plot_mdot_with_of(fig, axes[0][0], X, Y, mdot, of)
    plot_pc_with_stiffness(fig, axes[0][1], X, Y, pc)
    fig.suptitle(path.stem)
    fig.tight_layout()

def main():
    parser = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    parser.add_argument("file", nargs="?", type=Path, help="the .npy field to plot")
    parser.add_argument("--layers", help="comma-separated layer names in file order "
                                         "(must be mdot, Pc, OF in some order)")
    args = parser.parse_args()

    path = args.file
    if path is None and NPY_PATH:
        path = Path(NPY_PATH)
        if not path.is_absolute():
            path = Path(__file__).parent / path   # relative paths are relative to this script
    if path is None:
        found = sorted(Path(__file__).parent.glob("*Field*.npy"))
        if len(found) != 1:
            listing = "\n  ".join(p.name for p in found) or "(none found)"
            raise SystemExit(f"Set NPY_PATH in the script or pass a .npy path. "
                             f"Candidates next to this script:\n  {listing}")
        path = found[0]
    if not path.exists():
        raise SystemExit(f"File not found: {path}")

    field = load_field(path)
    names = [s.strip() for s in args.layers.split(",")] if args.layers else list(DEFAULT_LAYERS)
    if len(names) != field.shape[0] or set(names) != set(DEFAULT_LAYERS):
        raise SystemExit(f"{path.name} has {field.shape[0]} layers; expected 3 layers named "
                         f"{DEFAULT_LAYERS} (got {names}).")
    plot_field(path, field, names)
    plt.show()

if __name__ == "__main__":
    main()
