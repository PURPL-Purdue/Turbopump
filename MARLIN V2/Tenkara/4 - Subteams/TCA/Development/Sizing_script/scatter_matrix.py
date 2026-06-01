import pandas as pd
import matplotlib.pyplot as plt
from pandas.plotting import scatter_matrix

# ── Load data ──────────────────────────────────────────────────────────────────
df = pd.read_excel('designs.xlsx')

# ── Select columns to plot ─────────────────────────────────────────────────────
# Inputs
input_cols  = ['F_lbf', 'pc_psi', 'D_chamber_in', 't_chamber_in']

# Outputs
output_cols = ['Isp_s', 'Lc_m', 'De_in']

# Combined
cols = input_cols + output_cols

# Rename for nicer labels on the plot
rename = {
    'F_lbf':        'Thrust (lbf)',
    'pc_psi':       'Pc (psi)',
    'D_chamber_in': 'D chamber (in)',
    't_chamber_in': 'thickness (in)',
    'Isp_s':        'Isp (s)',
    'Lc_m':         'L chamber (m)',
    'De_in':        'D exit (in)',
}

df_plot = df[cols].rename(columns=rename)

# ── Plot ───────────────────────────────────────────────────────────────────────
scatter_matrix(
    df_plot,
    alpha=0.15,       # point transparency
    figsize=(16, 10),
    diagonal='kde',   # show smooth distribution on diagonal instead of histogram
    marker='.',
    color='steelblue',
)

plt.suptitle('Engine Design — Scatter Matrix', fontsize=12, y=0.95)
plt.tight_layout(rect=[0, 0, 1, 0.96])
plt.savefig('scatter_matrix.png', dpi=150, bbox_inches='tight')
plt.show()
print("Saved → scatter_matrix.png")