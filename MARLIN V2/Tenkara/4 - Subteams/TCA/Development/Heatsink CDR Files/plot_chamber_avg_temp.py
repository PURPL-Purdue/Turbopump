"""
Chamber wall space-weighted average temperature vs. time, Carbon Steel vs.
SS316, both at 0.03in mesh. Source files carry Minimum/Maximum/Average
columns (through-wall-thickness values); this plots the Average column,
which is a space-weighted average across the wall thickness at each time.
"""

from pathlib import Path

import pandas as pd
import matplotlib.pyplot as plt

BASE = Path(__file__).resolve().parent
CARBON_DATA = BASE / "Carbon Steel" / "Data"
SS_DATA = BASE / "SS" / "Data"

CUTOFF_TIME_S = 4.0

FILES = {
    "Carbon Steel": {"path": CARBON_DATA / "carbon_chamber_with_avg_0.03.txt", "color": "tab:blue"},
    "SS316": {"path": SS_DATA / "SS_chamber_with_avg_0.03.txt", "color": "tab:orange"},
}


def load_average(path):
    df = pd.read_csv(path, sep="\t", encoding="utf-8-sig")
    df = df.dropna(axis=1, how="all")
    time_col = next(c for c in df.columns if "Time" in c)
    avg_col = next(c for c in df.columns if "Average" in c)
    df = df.dropna(subset=[time_col, avg_col])
    return df[time_col].to_numpy(dtype=float), df[avg_col].to_numpy(dtype=float)


fig, ax = plt.subplots(figsize=(9, 6))

for material, spec in FILES.items():
    path = spec["path"]
    if not path.exists():
        print(f"Skipping missing file for {material}: {path}")
        continue
    time_s, avg_temp = load_average(path)
    ax.plot(time_s, avg_temp, color=spec["color"], marker="o", markersize=3, linewidth=1.4, label=material)

ax.axvline(CUTOFF_TIME_S, color="red", linestyle="--", linewidth=1.4, label="Engine Cutoff")

ax.set_title("Chamber Wall Space-Weighted Average Temperature vs. Time (0.03in Mesh)")
ax.set_xlabel("Time [s]")
ax.set_ylabel("Temperature [°F]")
ax.grid(True, alpha=0.4)
ax.legend()

fig.tight_layout()

out_path = BASE / "chamber_avg_temp.png"
fig.savefig(out_path, dpi=150, bbox_inches="tight")
print(f"Saved plot to {out_path}")

plt.show()
