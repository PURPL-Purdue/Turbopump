"""Plot Gas Generator field .npy file(s) as filled contour plots, one panel per layer.

Usage:
    python plot_GasGen_field.py                         # plots every *Field*.npy next to this script
    python plot_GasGen_field.py file.npy ...            # plots the files you name

File format: array of shape (n_layers, n_rows, n_cols) = field[layer, fuel, ox],
i.e. rows are fuel pressure, columns are ox pressure, both linspace(P_MIN, P_MAX, n).

When the layers are mdot_ox and mdot_fuel, a third panel shows total mdot (ox + fuel)
as filled contours with O/F (ox / fuel) contour lines on top.
"""
import argparse
from pathlib import Path
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D

# Must match the field generator (the .npy file does not store its pressure axis)
P_MIN = 100.0   # psi
P_MAX = 1000.0   # psi

OF_STEP_LOW = 0.2            # O/F line spacing from 0 up to OF_BREAK
OF_STEP_HIGH = 0.4           # O/F line spacing above OF_BREAK
OF_MAX_LINES = 20            # if more than this many lines are needed, trim to percentiles
OF_BREAK = 1.0
OF_TRIM_PERCENTILES = (0, 99)

DEFAULT_2_LAYERS = ("mdot_ox", "mdot_fuel")
LABELS = {
    "mdot_ox": "Ox mdot (lb/s)",
    "mdot_fuel": "Fuel mdot (lb/s)",
}

def load_field(path: Path) -> np.ndarray:
    field = np.load(path, allow_pickle=False)
    if field.ndim == 2:                      # a single NxN layer
        field = field[None]
    return field

def style(ax, title):
    ax.set_xlabel("Oxidizer pressure (psi)")
    ax.set_ylabel("Fuel pressure (psi)")
    ax.set_title(title)
    ax.set_xlim(P_MIN, P_MAX)
    ax.set_ylim(P_MIN, P_MAX)

def of_levels(lo, hi):
    """Adaptive O/F contour levels: OF_STEP_LOW up to OF_BREAK, then OF_STEP_HIGH above it.
    Only levels inside [lo, hi] are returned."""
    low = OF_STEP_LOW * np.arange(0, int(round(OF_BREAK / OF_STEP_LOW)) + 1)       # 0, 0.2, ... 1.0
    n_high = int(np.ceil((hi - OF_BREAK) / OF_STEP_HIGH)) if hi > OF_BREAK else 0
    high = OF_BREAK + OF_STEP_HIGH * np.arange(1, n_high + 1)                      # 1.4, 1.8, ...
    levels = np.concatenate([low, high])
    return levels[(levels >= lo) & (levels <= hi)]

def plot_total_with_of(fig, ax, X, Y, mdot_ox, mdot_fuel):
    """Total mdot (filled) with O/F = ox / fuel contour lines."""
    total = mdot_ox + mdot_fuel
    with np.errstate(divide="ignore", invalid="ignore"):
        of = np.where(mdot_fuel > 0, mdot_ox / mdot_fuel, np.nan)

    good = ~np.isnan(total)
    if good.sum() >= 4 and np.ptp(total[good]) > 0:
        cf = ax.contourf(X, Y, np.ma.masked_invalid(total), levels=12, cmap="viridis")
        fig.colorbar(cf, ax=ax, label="Total mdot (lb/s)")
        handles = []
        good_of = ~np.isnan(of)
        if good_of.sum() >= 4 and np.ptp(of[good_of]) > 0:
            lo, hi = np.nanmin(of), np.nanmax(of)
            # print(f"O/F range: {lo:.2f} to {hi:.2f}")
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
                                 colors="black", linestyles="--", linewidths=1.0)
                ax.clabel(ocs, fontsize=7, fmt="%.1f")
            handles.append(Line2D([], [], color="black", ls="--", lw=1.2, label="O/F"))
        if handles:
            ax.legend(handles=handles, loc="upper right", fontsize=8)
    else:
        ax.text(0.5, 0.5, "Not enough valid points to contour", ha="center",
                va="center", transform=ax.transAxes)
    style(ax, "Total mdot (lb/s) with O/F lines")

def plot_layers(path: Path, field: np.ndarray, names):
    n_layers, n_rows, n_cols = field.shape
    p_fu = np.linspace(P_MIN, P_MAX, n_rows)   # rows -> fuel pressure
    p_ox = np.linspace(P_MIN, P_MAX, n_cols)   # cols -> ox pressure
    X, Y = np.meshgrid(p_ox, p_fu)             # X = ox (horizontal), Y = fuel (vertical)

    has_total = "mdot_ox" in names and "mdot_fuel" in names
    n_panels = n_layers + (1 if has_total else 0)

    fig, axes = plt.subplots(1, n_panels, figsize=(5.4 * n_panels, 4.8), squeeze=False)
    for ax, layer, name in zip(axes[0], field, names):
        bad = np.isnan(layer)
        label = LABELS.get(name, name)
        if (~bad).sum() >= 4 and np.ptp(layer[~bad]) > 0:
            cf = ax.contourf(X, Y, np.ma.masked_invalid(layer), levels=12, cmap="viridis")
            cs = ax.contour(X, Y, np.ma.masked_invalid(layer), levels=cf.levels,
                            colors="k", linewidths=0.4)
            ax.clabel(cs, fontsize=7, fmt="%.2f")
            fig.colorbar(cf, ax=ax, label=label)
        else:
            ax.text(0.5, 0.5, "Not enough valid points to contour", ha="center",
                    va="center", transform=ax.transAxes)
        style(ax, label)

    if has_total:
        plot_total_with_of(fig, axes[0][n_layers], X, Y,
                           field[names.index("mdot_ox")], field[names.index("mdot_fuel")])
    fig.suptitle(f"{path.stem}")
    fig.tight_layout()

def main():
    parser = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    parser.add_argument("files", nargs="*", type=Path)
    parser.add_argument("--layers", help="comma-separated layer names in file order")
    args = parser.parse_args()

    paths = args.files or sorted(Path(__file__).parent.glob("*Field*.npy"))
    if not paths:
        raise SystemExit("No .npy files found.")

    for p in paths:
        field = load_field(p)
        n = field.shape[0]
        if args.layers:
            names = [s.strip() for s in args.layers.split(",")]
        elif n == 2:
            names = list(DEFAULT_2_LAYERS)
        else:
            names = [f"layer {i}" for i in range(n)]
        if len(names) != n:
            raise SystemExit(f"{p.name} has {n} layers but {len(names)} names were given.")
        plot_layers(p, field, names)
    plt.show()


if __name__ == "__main__":
    main()
