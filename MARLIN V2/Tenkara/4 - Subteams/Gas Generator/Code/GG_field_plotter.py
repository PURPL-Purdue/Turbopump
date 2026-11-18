"""Plot a Gas Generator field .npy file as two contour plots:

  1. Total mass flow rate (filled) with O/F ratio lines (black).
  2. Chamber pressure (filled) with minimum-stiffness lines (black), and red diagonal
     hatching over the region where minimum stiffness < STIFF_THRESHOLD.

Usage:
    python plot_GasGen_field.py file.npy                # plots that one file
    python plot_GasGen_field.py                         # uses NPY_PATH below; if that is None and exactly
                                                        # one *Field*.npy sits next to this script, plots it

File format: array of shape (3, n_rows, n_cols) = field[layer, fuel, ox], layers in the order
    0: total mdot (lb/s), 1: chamber pressure Pc (psi), 2: O/F ratio (ox / fuel)
Rows are fuel feed pressure, columns are ox feed pressure, both linspace(P_MIN, P_MAX, n).

Stiffness of each feed = (P_up - Pc) / Pc. The minimum stiffness at a grid point is the
smaller of the ox and fuel values, i.e. (min(P_ox, P_fuel) - Pc) / Pc, since both share Pc.
"""
import argparse
from pathlib import Path
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.patches import Patch

# Which .npy to plot. Set the path here (absolute, or relative to this script)
NPY_PATH = "GasGenerator_Field_10pct.npy"

# Must match the field generator (the .npy file does not store its pressure axis)
P_MIN = 100.0    # psi
P_MAX = 1000.0   # psi

# O/F line spacing: OF_STEP_LOW up to OF_BREAK, then OF_STEP_HIGH above it
OF_STEP_LOW = 0.2
OF_STEP_HIGH = 0.4
OF_BREAK = 1.0
OF_MAX_LINES = 40            # if more lines than this are needed, trim to the bulk of the data
OF_TRIM_PERCENTILES = (0, 99)

# Stiffness lines and mask
STIFF_STEP = 0.1             # spacing between stiffness lines
STIFF_THRESHOLD = 0.2        # region below this is masked in red; its line is drawn solid red
STIFF_LINE_WIDTH = 0.9       # width of the dashed black stiffness lines
THRESH_LINE_WIDTH = 0.9      # width of the red threshold line
MASK_HATCH = "///"           # diagonal red hatching: "//" = sparser, "////" = denser, "xx" = cross-hatch
MASK_LINEWIDTH = 0.9         # thickness of the hatch lines

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
        lo, hi = np.nanmin(of), np.nanmax(of)
        print(f"O/F range: {lo:.2f} to {hi:.2f}")
        levels = of_levels(lo, hi)
        if len(levels) > OF_MAX_LINES:
            # A few extreme O/F points (e.g. near-zero fuel flow) stretch the range;
            # keep the lines to the bulk of the data so they stay readable.
            p_lo, p_hi = np.nanpercentile(of, OF_TRIM_PERCENTILES)
            levels = of_levels(p_lo, p_hi)
            # print(f"O/F lines limited to {p_lo:.2f}-{p_hi:.2f} "
            #       f"(more than {OF_MAX_LINES} lines needed)")
        if len(levels):
            ocs = ax.contour(X, Y, np.ma.masked_invalid(of), levels=levels,
                             colors="k", linestyles="--", linewidths=0.9)
            ax.clabel(ocs, fontsize=7, fmt="%.1f")
            ax.legend(handles=[Line2D([], [], color="k", ls="--", lw=1.0, label="O/F")],
                      loc="upper left", fontsize=8)
    style(ax, "Total mdot (lb/s) with O/F lines")


def plot_pc_with_stiffness(fig, ax, X, Y, pc):
    """Chamber pressure (filled) with black min-stiffness lines and a red mask below threshold."""
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
    s_lo, s_hi = np.nanmin(stiff), np.nanmax(stiff)
    print(f"Min stiffness range: {s_lo:.3f} to {s_hi:.3f}")

    handles = []
    # Red grid (hatching) where stiffness < threshold; the fill is left clear so Pc stays readable
    if s_lo < STIFF_THRESHOLD:
        with plt.rc_context({"hatch.color": "red", "hatch.linewidth": MASK_LINEWIDTH}):
            ax.contourf(X, Y, np.ma.masked_invalid(stiff), levels=[s_lo - 1.0, STIFF_THRESHOLD],
                        colors="none", hatches=[MASK_HATCH])
            handles.append(Patch(facecolor="none", edgecolor="red", hatch=MASK_HATCH,
                                 label=f"Stiffness < {STIFF_THRESHOLD:g}"))

    # Dashed black stiffness lines; the threshold (separating) line is solid red
    levels = np.arange(STIFF_STEP, s_hi + STIFF_STEP, STIFF_STEP)
    levels = levels[(levels >= s_lo) & (levels <= s_hi)]
    if len(levels):
        is_thr = np.isclose(levels, STIFF_THRESHOLD)
        scs = ax.contour(X, Y, np.ma.masked_invalid(stiff), levels=levels,
                         colors=["red" if t else "k" for t in is_thr],
                         linestyles=["solid" if t else "dashed" for t in is_thr],
                         linewidths=[THRESH_LINE_WIDTH if t else STIFF_LINE_WIDTH for t in is_thr])
        ax.clabel(scs, fontsize=7, fmt="%.1f")
        handles.append(Line2D([], [], color="k", ls="--", lw=1.0, label="Min stiffness"))
        if is_thr.any():
            handles.append(Line2D([], [], color="red", lw=THRESH_LINE_WIDTH,
                                  label=f"Stiffness = {STIFF_THRESHOLD:g}"))

    if handles:
        ax.legend(handles=handles, loc="upper left", fontsize=8)
    style(ax, "Chamber pressure (psi) with min stiffness")


def plot_field(path: Path, field: np.ndarray, names):
    n_layers, n_rows, n_cols = field.shape
    p_fu = np.linspace(P_MIN, P_MAX, n_rows)   # rows -> fuel pressure
    p_ox = np.linspace(P_MIN, P_MAX, n_cols)   # cols -> ox pressure
    X, Y = np.meshgrid(p_ox, p_fu)             # X = ox (horizontal), Y = fuel (vertical)

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
