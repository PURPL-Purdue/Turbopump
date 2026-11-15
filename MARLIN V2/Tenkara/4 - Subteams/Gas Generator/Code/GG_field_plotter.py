"""Plot Gas Generator field CSV(s) as filled contour plots.

Usage:
    python plot_GasGen_field.py                  # plots every Gas_Generator_Field_*.csv next to this script
    python plot_GasGen_field.py file1.csv ...    # plots the files you name
"""
import sys
from pathlib import Path
import numpy as np
import matplotlib.pyplot as plt

# Must match the field generator
P_MIN = 0.0   # psi
P_MAX = 800.0   # psi


def plot_field(csv_path: Path):
    field = np.loadtxt(csv_path, delimiter=",", ndmin=2)  # field[row, col] = [fuel, ox]
    n_rows, n_cols = field.shape
    p_fu = np.linspace(P_MIN, P_MAX, n_rows)  # rows -> fuel pressure
    p_ox = np.linspace(P_MIN, P_MAX, n_cols)  # cols -> ox pressure
    X, Y = np.meshgrid(p_ox, p_fu)            # X = ox (horizontal), Y = fuel (vertical)

    bad = np.isnan(field)
    data = np.ma.masked_invalid(field)

    fig, ax = plt.subplots(figsize=(7, 5.5))
    if (~bad).sum() >= 4:
        cf = ax.contourf(X, Y, data, levels=15, cmap="viridis")
        cs = ax.contour(X, Y, data, levels=cf.levels, colors="k", linewidths=0.4)
        ax.clabel(cs, fontsize=7, fmt="%.2f")
        fig.colorbar(cf, ax=ax, label="mdot (lb/s)")
    else:
        ax.text(0.5, 0.5, "Not enough valid points to contour", ha="center", va="center",
                transform=ax.transAxes)

    # # Mark grid points with no solution
    # if bad.any():
    #     ax.scatter(X[bad], Y[bad], marker="x", color="red", s=18, label=f"NaN ({bad.sum()})")
    #     ax.legend(loc="upper right")

    ax.set_xlabel("Oxidizer pressure (psi)")
    ax.set_ylabel("Fuel pressure (psi)")
    ax.set_title(csv_path.stem.replace("_", " "))
    ax.set_xlim(P_MIN, P_MAX)
    ax.set_ylim(P_MIN, P_MAX)
    fig.tight_layout()


def main():
    if len(sys.argv) > 1:
        paths = [Path(a) for a in sys.argv[1:]]
    else:
        paths = sorted(Path(__file__).parent.glob("Gas_Generator_Field_*.csv"))
    if not paths:
        sys.exit("No CSV files found.")
    for p in paths:
        plot_field(p)
    plt.show()


if __name__ == "__main__":
    main()