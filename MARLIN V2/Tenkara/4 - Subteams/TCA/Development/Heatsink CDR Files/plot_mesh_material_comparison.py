"""
Chamber / Nozzle (convergent) / Gasket maximum wall temperature vs. time,
comparing Carbon Steel and SS316, each at coarse (0.1in) and fine (0.03in)
mesh density. One subplot per portion, shared legend, right-hand Kelvin axis,
engine-cutoff marker, and per-material melting point reference lines.
"""

import re
from pathlib import Path

import pandas as pd
import matplotlib.pyplot as plt

BASE = Path(__file__).resolve().parent
CARBON_DATA = BASE / "Carbon Steel" / "Data"
SS_DATA = BASE / "SS" / "Data"

CUTOFF_TIME_S = 4.0

# ---------------------------------------------------------------------------
# Materials: color, melting point [F] (skipped on the gasket subplot), and
# fine/coarse file per portion.
# ---------------------------------------------------------------------------
MATERIALS = {
    "Carbon Steel": {
        "color": "tab:blue",
        "melt_F": 2600,  # A36-ish generic carbon steel, ~2570-2800F depending on grade
        "files": {
            "chamber": {"Fine": CARBON_DATA / "carbon_chamber_0.03in.txt", "Coarse": CARBON_DATA / "carbon_chamber_0.1in.txt"},
            "convergent": {"Fine": CARBON_DATA / "carbon_convergent_0.03in.txt", "Coarse": CARBON_DATA / "carbon_convergent_0.1in.txt"},
            "gasket": {"Fine": CARBON_DATA / "carbon_gasket_0.03in.txt", "Coarse": CARBON_DATA / "carbon_gasket_0.1in.txt"},
        },
    },
    "SS316": {
        "color": "tab:orange",
        "melt_F": 2530,  # 316 SS melting range ~2507-2552F
        "files": {
            "chamber": {"Fine": SS_DATA / "SS_chamber_0.03.txt", "Coarse": SS_DATA / "SS_chamber_0.1.txt"},
            "convergent": {"Fine": SS_DATA / "SS_Convergent_0.03.txt", "Coarse": SS_DATA / "SS_Convergent_0.1.txt"},
            "gasket": {"Fine": SS_DATA / "SS_gasket_0.03.txt", "Coarse": SS_DATA / "SS_gasket_0.1.txt"},
        },
    },
}

MESH_STYLE = {
    "Fine": {"marker": "o", "markersize": 4, "markevery": 3},
    "Coarse": {"marker": "^", "markersize": 4, "markevery": 3},
}

PORTIONS = ["chamber", "convergent", "gasket"]
PORTION_TITLES = {
    "chamber": "Chamber Maximum Temperature",
    "convergent": "Nozzle Maximum Temperature",
    "gasket": "Gasket Maximum Temperature",
}
PORTIONS_WITHOUT_MELT_LINE = {"gasket"}


def load_series(path):
    """Return (time_s, temp, temp_label) from one exported temperature file."""
    df = pd.read_csv(path, sep="\t", encoding="utf-8-sig")
    df = df.dropna(axis=1, how="all")
    time_col = next(c for c in df.columns if "Time" in c)
    temp_col = next(c for c in df.columns if "Temperature" in c)
    df = df.dropna(subset=[time_col, temp_col])
    return (
        df[time_col].to_numpy(dtype=float),
        df[temp_col].to_numpy(dtype=float),
        temp_col,
    )


def unit_of(label):
    match = re.search(r"\[(.*?)\]", label)
    return match.group(1) if match else ""


def f_to_k(f):
    return (f - 32) * 5.0 / 9.0 + 273.15


def k_to_f(k):
    return (k - 273.15) * 9.0 / 5.0 + 32.0


fig, axes = plt.subplots(len(PORTIONS), 1, figsize=(11, 13), sharex=True)

for ax, portion in zip(axes, PORTIONS):
    temp_unit = None

    for material, spec in MATERIALS.items():
        color = spec["color"]
        for mesh, path in spec["files"][portion].items():
            if not path.exists():
                print(f"Skipping missing file for {material} / {portion} / {mesh}: {path}")
                continue

            time_s, temp, temp_label = load_series(path)
            temp_unit = temp_unit or unit_of(temp_label)
            style = MESH_STYLE[mesh]
            ax.plot(
                time_s, temp,
                color=color,
                marker=style["marker"],
                markersize=style["markersize"],
                markevery=style["markevery"],
                linewidth=1.3,
                label=f"{material} ({mesh})",
            )

        melt_F = spec.get("melt_F")
        if melt_F is not None and portion not in PORTIONS_WITHOUT_MELT_LINE:
            ax.axhline(
                melt_F,
                color=color,
                linestyle="--",
                linewidth=1.2,
                label=f"{material} melting pt.",
            )

    ax.axvline(
        CUTOFF_TIME_S,
        color="red",
        linestyle="--",
        linewidth=1.4,
        label="Engine Cutoff",
    )

    ax.set_title(PORTION_TITLES[portion])
    ax.set_ylabel(f"Temperature [{temp_unit or '?'}]")
    ax.grid(True, alpha=0.4)

    secax = ax.secondary_yaxis("right", functions=(f_to_k, k_to_f))
    secax.set_ylabel("Temperature [K]")

axes[-1].set_xlabel("Time [s]")
fig.suptitle("Wall Temperature vs. Time — Material & Mesh Comparison")

# One shared legend for the whole figure.
handles_by_label = {}
for ax in axes:
    for handle, label in zip(*ax.get_legend_handles_labels()):
        handles_by_label.setdefault(label, handle)
fig.legend(
    handles_by_label.values(),
    handles_by_label.keys(),
    loc="center left",
    bbox_to_anchor=(0.90, 0.5),
    borderaxespad=0.0,
)

fig.tight_layout(rect=(0, 0, 0.86, 1))

out_path = BASE / "mesh_material_comparison.png"
fig.savefig(out_path, dpi=150, bbox_inches="tight")
print(f"Saved plot to {out_path}")

plt.show()
