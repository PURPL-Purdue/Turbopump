import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

# ── Load data ──────────────────────────────────────────────────────────────────
df = pd.read_excel('designs.xlsx')

# ── Pareto front function ──────────────────────────────────────────────────────
# objectives: dict of {column: 'max' or 'min'}
# Returns boolean mask: True = Pareto optimal
def is_pareto(df, objectives):
    # Convert everything to "higher is better" by flipping minimization columns
    values = np.zeros((len(df), len(objectives)))
    for i, (col, direction) in enumerate(objectives.items()):
        values[:, i] = df[col].values if direction == 'max' else -df[col].values

    mask = np.ones(len(df), dtype=bool)
    for i in range(len(df)):
        # Dominated if another design is >= in all AND > in at least one
        dominated = (
            np.all(values >= values[i], axis=1) &
            np.any(values >  values[i], axis=1)
        )
        dominated[i] = False  # don't compare against itself
        if dominated.any():
            mask[i] = False
    return mask

# ── Define objectives ──────────────────────────────────────────────────────────
objectives = {
    'sigma_vM_MPa': 'min',   # minimize von Mises stress
    'Isp_s':        'max',   # maximize Isp
    'F_lbf':        'max',   # maximize thrust
    'pc_psi':       'min',   # minimize chamber pressure
}

# ── Compute Pareto front ───────────────────────────────────────────────────────
mask   = is_pareto(df, objectives)
pareto = df[mask]
rest   = df[~mask]

print(f"Total designs:  {len(df)}")
print(f"Pareto optimal: {mask.sum()}")
print(pareto[['F_lbf', 'pc_psi', 'sigma_vM_MPa', 'Isp_s', 'D_chamber_in', 'Lc_m']].sort_values('F_lbf').to_string(index=False))

# ── Plot: 2D projections of the Pareto front ──────────────────────────────────
fig, axes = plt.subplots(1, 4, figsize=(15, 5))

pairs = [
    ('F_lbf',  'Isp_s',        'Thrust (lbf)', 'Isp (s)'),
    ('pc_psi', 'Isp_s',        'Pc (psi)',      'Isp (s)'),
    ('pc_psi', 'sigma_vM_MPa', 'Pc (psi)',      'Von Mises Stress (MPa)'),
    ('pc_psi', 'D_chamber_in', 'Pc (psi)',      'D_chamber_in (in)'),
]

for ax, (x_col, y_col, x_label, y_label) in zip(axes, pairs):
    ax.scatter(rest[x_col],   rest[y_col],   color='lightgrey', s=8,  alpha=0.5, label='All designs')
    ax.scatter(pareto[x_col], pareto[y_col], color='red',       s=40, zorder=3,  label='Pareto optimal')
    ax.set_xlabel(x_label)
    ax.set_ylabel(y_label)
    ax.set_title(f'{y_label} vs {x_label}')
    ax.legend(fontsize=7)
    ax.grid(True, alpha=0.3)

plt.suptitle('Pareto Front — Minimize Von Mises Stress, Maximize Isp, Thrust / Minimize Pc', fontsize=11, y=1.02)
plt.tight_layout()
plt.savefig('pareto_front.png', dpi=150, bbox_inches='tight')
plt.show()