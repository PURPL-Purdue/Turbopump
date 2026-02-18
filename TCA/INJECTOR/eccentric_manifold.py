import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import minimize_scalar

# ─────────────────────────────────────────────
# UNIT CONVERSIONS
# ─────────────────────────────────────────────
psi_into_pa         = 6894.76
pa_to_psi           = 1 / psi_into_pa
meters_into_inches  = 39.37
inches_into_meters  = 0.0254
kg_into_lbm         = 2.20462

# ─────────────────────────────────────────────
# DESIGN PARAMETERS  (edit these)
# ─────────────────────────────────────────────
mdot_total   = 9.0          # Total propellant mass flow [kg/s]
OF_Ratio     = 2.1          # O/F ratio
mdot_kero    = mdot_total / (1 + OF_Ratio)   # Fuel mass flow [kg/s]

rho_rp1      = 810.0        # RP-1 density [kg/m³]
P_in_psi     = 710.12       # Inlet pressure [psi]

# Injector element layout
n_type1      = 8            # mini-manifolds with 6 holes each
n_type2      = 8            # mini-manifolds with 4 holes each
holes_type1  = 6
holes_type2  = 4
n_holes_total = n_type1 * holes_type1 + n_type2 * holes_type2   # 80 holes
mdot_per_hole = mdot_kero / n_holes_total

# Manifold geometry  ← FIXED OUTER, ECCENTRIC INNER
R_outer      = 76.06*10**-3        # Outer circle radius [m]  ← fixed
h_manifold   = 19*10**-3    # Manifold height (axial depth) [m]  ← constant

# Dynamic pressure assumption (sets velocity)
dyn_pressure_fraction = 0.01

# ─────────────────────────────────────────────
# DERIVED FLOW QUANTITIES
# ─────────────────────────────────────────────
dyn_pressure = dyn_pressure_fraction * P_in_psi * psi_into_pa
v_manifold   = np.sqrt(2 * dyn_pressure / rho_rp1)

print(v_manifold)

mdot_branch  = mdot_kero / 2.0
A_inlet      = mdot_branch / (rho_rp1 * v_manifold)
width_inlet  = A_inlet / h_manifold

print("=" * 55)
print("ECCENTRIC MANIFOLD SIZING  (fixed outer, eccentric inner)")
print("=" * 55)
print(f"Fuel mass flow (total):   {mdot_kero*kg_into_lbm:.4f} lbm/s  |  {mdot_kero:.4f} kg/s")
print(f"Mass flow per branch:     {mdot_branch*kg_into_lbm:.4f} lbm/s  |  {mdot_branch:.4f} kg/s")
print(f"Manifold velocity:        {v_manifold:.3f} m/s  ({v_manifold*3.28084:.3f} ft/s)")
print(f"Inlet area (per branch):  {A_inlet*1e6:.3f} mm²  ({A_inlet*meters_into_inches**2:.4f} in²)")
print(f"Inlet width:              {width_inlet*1e3:.3f} mm  ({width_inlet*meters_into_inches:.4f} in)")

# ─────────────────────────────────────────────
# ANGULAR STATIONS
# ─────────────────────────────────────────────
holes_per_step = []
for _ in range(4):
    holes_per_step.append(holes_type2)   # 4 holes
    holes_per_step.append(holes_type1)   # 6 holes
n_steps = len(holes_per_step)            # 8 stations per branch

theta_elements = np.linspace(0, np.pi, n_steps + 2)[1:-1]

mdot_remaining  = mdot_branch
areas_required  = []
mdot_at_station = []

for i, n_holes in enumerate(holes_per_step):
    areas_required.append(mdot_remaining / (rho_rp1 * v_manifold))
    mdot_at_station.append(mdot_remaining)
    mdot_remaining -= n_holes * mdot_per_hole

areas_required  = np.array(areas_required)
mdot_at_station = np.array(mdot_at_station)
widths_required = areas_required / h_manifold

print(f"\n{'Station':>8} {'θ (°)':>8} {'ṁ (kg/s)':>12} {'Area (mm²)':>12} {'Width (mm)':>12}")
print("-" * 55)
for i in range(n_steps):
    print(f"{i+1:>8} {np.degrees(theta_elements[i]):>8.1f} "
          f"{mdot_at_station[i]:>12.5f} "
          f"{areas_required[i]*1e6:>12.3f} "
          f"{widths_required[i]*1e3:>12.3f}")

# ─────────────────────────────────────────────
# ECCENTRIC FIT  — fixed outer, moving inner
#
# The inner circle has radius R_inner and its center is shifted
# by eccentricity e toward θ=180° (away from inlet).
#
# Gap at angle θ:
#   width(θ) = R_outer - (R_inner - e·cos(θ))
#            = (R_outer - R_inner) + e·cos(θ)
#
# This is identical in form to before, but now e describes
# the shift of the INNER circle (toward the narrow end).
#
# Two free params:  mean_gap = R_outer - R_inner,  and  e
# Inlet constraint: mean_gap + e = width_inlet
#                   → mean_gap = width_inlet - e
# Outlet gap:       mean_gap - e  (must stay > 0)
#
# Solve for e that minimises RMS width error at all stations.
# ─────────────────────────────────────────────

def gap_model(theta, mean_gap, e):
    """Channel width at angle theta for eccentric inner circle."""
    return mean_gap + e * np.cos(theta)

def fit_error(e):
    mean_gap = width_inlet - e
    if mean_gap - e <= 0:       # gap would close at 180°
        return 1e10
    if mean_gap - e > R_outer:  # inner radius would be negative
        return 1e10
    predicted = gap_model(theta_elements, mean_gap, e)
    return np.sqrt(np.mean((predicted - widths_required)**2))

result   = minimize_scalar(fit_error, bounds=(0, width_inlet * 0.99), method='bounded')
e_best   = result.x
mean_gap = width_inlet - e_best

# Derived radii
R_inner  = R_outer - mean_gap   # inner circle radius (from its own center)
e_offset = e_best               # inner circle center shifted by this toward θ=180°

# Sanity check
assert R_inner > 0, "Inner radius went negative — reduce R_outer or h_manifold"
assert mean_gap - e_best > 0,  "Gap closes before 180° — increase R_outer"

print(f"\n{'='*55}")
print("ECCENTRIC GEOMETRY SOLUTION")
print(f"{'='*55}")
print(f"Outer circle radius R_outer:  {R_outer*1e3:.3f} mm  (fixed)")
print(f"Inner circle radius R_inner:  {R_inner*1e3:.3f} mm")
print(f"Eccentricity e:               {e_offset*1e3:.3f} mm  (inner center shift)")
print(f"Mean gap:                     {mean_gap*1e3:.3f} mm")
print(f"Gap at inlet  (θ=0°):         {gap_model(0,      mean_gap, e_best)*1e3:.3f} mm  [target: {width_inlet*1e3:.3f} mm]")
print(f"Gap at outlet (θ=180°):       {gap_model(np.pi,  mean_gap, e_best)*1e3:.3f} mm")
print(f"RMS fit error:                {result.fun*1e3:.4f} mm")

print(f"\n{'Station':>8} {'θ (°)':>8} {'Required (mm)':>15} {'Model (mm)':>12} {'Error (mm)':>12}")
print("-" * 55)
for i in range(n_steps):
    w_model = gap_model(theta_elements[i], mean_gap, e_best)
    err     = (w_model - widths_required[i]) * 1e3
    print(f"{i+1:>8} {np.degrees(theta_elements[i]):>8.1f} "
          f"{widths_required[i]*1e3:>15.3f} "
          f"{w_model*1e3:>12.3f} "
          f"{err:>12.4f}")

# ─────────────────────────────────────────────
# PLOTS
# ─────────────────────────────────────────────
theta_full   = np.linspace(0, np.pi, 300)
widths_model = gap_model(theta_full, mean_gap, e_best)

fig, axes = plt.subplots(1, 2, figsize=(14, 6))
fig.suptitle("Eccentric Manifold Design — Fixed Outer, Eccentric Inner", fontsize=13, fontweight='bold')

# ── Plot 1: Width profile ──────────────────────────────────────────────────
ax1 = axes[0]
ax1.plot(np.degrees(theta_full), widths_model * 1e3, 'b-', linewidth=2, label='Eccentric model')
ax1.scatter(np.degrees(theta_elements), widths_required * 1e3,
            color='red', zorder=5, s=80, label='Required (discrete stations)')
ax1.set_xlabel("Angular position θ [°]", fontsize=12)
ax1.set_ylabel("Manifold channel width [mm]", fontsize=12)
ax1.set_title("Channel Width vs. Angular Position", fontsize=12)
ax1.legend()
ax1.grid(True, alpha=0.3)
ax1.set_xlim(0, 180)

# ── Plot 2: Top-view cross section ────────────────────────────────────────
ax2 = axes[1]
ax2.set_aspect('equal')

theta_circle = np.linspace(0, 2 * np.pi, 500)

# Outer circle — fixed, centered at origin
ax2.plot(R_outer * np.cos(theta_circle) * 1e3,
         R_outer * np.sin(theta_circle) * 1e3,
         'b-', linewidth=2.5, label='Outer wall (fixed)')

# Inner circle — center shifted by e toward θ=180° (negative x direction)
inner_cx = -e_offset * 1e3   # mm
inner_cy = 0.0
ax2.plot(inner_cx + R_inner * np.cos(theta_circle) * 1e3,
         inner_cy + R_inner * np.sin(theta_circle) * 1e3,
         'k-', linewidth=2.5, label='Inner wall (eccentric)')

# Mark centers
ax2.plot(0, 0, 'b+', markersize=12, markeredgewidth=2, label='Outer center')
ax2.plot(inner_cx, inner_cy, 'k+', markersize=12, markeredgewidth=2, label='Inner center')

# Eccentricity arrow
ax2.annotate('', xy=(inner_cx, -R_outer*1e3*0.15), xytext=(0, -R_outer*1e3*0.15),
             arrowprops=dict(arrowstyle='<->', color='green', lw=1.5))
ax2.text(inner_cx / 2, -R_outer*1e3*0.22, f'e = {e_offset*1e3:.2f} mm',
         ha='center', color='green', fontsize=9)

# Element stations (on the inner circle surface, upper half)
for th in theta_elements:
    x_pt = inner_cx + R_inner * np.cos(th) * 1e3
    y_pt = inner_cy + R_inner * np.sin(th) * 1e3
    ax2.plot(x_pt, y_pt, 'ro', markersize=5)

# Shade manifold channel (upper half)
theta_shade  = np.linspace(0, np.pi, 300)
x_outer_half = R_outer * np.cos(theta_shade) * 1e3
y_outer_half = R_outer * np.sin(theta_shade) * 1e3
x_inner_half = inner_cx + R_inner * np.cos(theta_shade) * 1e3
y_inner_half = inner_cy + R_inner * np.sin(theta_shade) * 1e3

ax2.fill_between(np.degrees(theta_shade),
                 y_inner_half, y_outer_half,
                 alpha=0.0)  # skip fill_between in polar sense; use patches

# Proper shading: fill between the two arcs in x-y space
for j in range(len(theta_shade) - 1):
    xs = [x_outer_half[j], x_outer_half[j+1], x_inner_half[j+1], x_inner_half[j]]
    ys = [y_outer_half[j], y_outer_half[j+1], y_inner_half[j+1], y_inner_half[j]]
    ax2.fill(xs, ys, color='steelblue', alpha=0.2)

# Inlet marker
ax2.annotate('Inlet (max width)', xy=(R_outer*1e3, 0), xytext=(R_outer*1e3 + 5, 8),
             fontsize=8, color='darkred',
             arrowprops=dict(arrowstyle='->', color='darkred'))

ax2.set_xlabel("x [mm]", fontsize=12)
ax2.set_ylabel("y [mm]", fontsize=12)
ax2.set_title("Top View — Eccentric Circle Geometry", fontsize=12)
ax2.legend(fontsize=8, loc='lower right')
ax2.grid(True, alpha=0.3)

plt.tight_layout()
plt.show()
print("\nPlot saved.")