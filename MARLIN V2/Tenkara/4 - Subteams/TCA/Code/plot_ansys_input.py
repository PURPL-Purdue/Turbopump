"""
Plots column 2 (hg, heat transfer coefficient) vs. column 1 (x, axial
position) from Outputs/ansys_input2.csv. That file is headerless -- it's the
spline-resampled (x_ansys, hg_ansys, Tr_ansys) output written by bartz.py.
"""

import numpy as np
import matplotlib.pyplot as plt

data = np.loadtxt("Outputs/ansys_input2.csv", delimiter=",")
x = data[:, 0]
hg = data[:, 1]

fig, ax = plt.subplots(figsize=(9, 6))
ax.plot(x, hg, linewidth=1.5)

ax.set_title("Gas-Side Heat Transfer Coefficient vs. Axial Position")
ax.set_xlabel("Axial position, x [m]")
ax.set_ylabel("Heat transfer coefficient, $h_g$ [W/m²K]")
ax.grid(True, alpha=0.4)

fig.tight_layout()

out_path = "Outputs/ansys_input2_hg_vs_x.png"
fig.savefig(out_path, dpi=150, bbox_inches="tight")
print(f"Saved plot to {out_path}")

plt.show()
