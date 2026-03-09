"""
LOX Toroidal Manifold Sizing Tool
==================================
Toroid cross-section: D-shape (rectangle + optimized curved top)
The curved top is an ellipse arc fitted to minimize centrifugal ΔP
while smoothly connecting the rectangle top corners.

Coordinate system:
  - R_mean : mean radius of the toroid (centerline of cross-section)
  - rect_width  (w): width of the rectangular part  [m]
  - rect_height (h): height of the rectangular part [m]
  - Curved top: ellipse arc from (-w/2, h) to (+w/2, h) relative to cross-section center
"""

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib.patches import Circle
from scipy.optimize import minimize_scalar

# =============================================================================
# ── USER INPUTS ──────────────────────────────────────────────────────────────
# =============================================================================

# --- Fluid ---
mdot        = 6.145          # [kg/s]  total LOX mass flow into manifold
rho         = 1141.0       # [kg/m³] LOX density at operating conditions (~1 bar, ~90 K)

# --- Inlet pipe ---
D_inlet     = 0.021        # [m]     inlet pipe inner diameter
theta_swirl = 55.0         # [deg]   tangential inlet angle (0 = radial, 90 = fully tangential)

# --- Toroid geometry (rectangular part) ---
rect_width  = 0.044        # [m]     w — width of rectangular cross-section
rect_height = 0.015        # [m]     h — height of rectangular cross-section
R_mean      = 0.047        # [m]     mean toroid radius (centerline)

# --- Injector holes ---
n_holes     = 80           # [-]     number of outlet holes
d_hole      = 0.003048        # [m]     hole diameter
Cd_hole     = 0.8          # [-]     discharge coefficient

# --- Pressure budget ---
dP_injector = 780969.159          # [Pa]    injector ΔP (pump outlet - chamber pressure)
                           #         e.g. 5 bar = 5e5 Pa

b_optimal = 0.02                   # [m]  fixed ellipse cap height = 2 mm

# =============================================================================
# ── CALCULATIONS ─────────────────────────────────────────────────────────────
# =============================================================================

# --- Inlet geometry ---
A_inlet     = np.pi * (D_inlet / 2)**2          # [m²]
V_inlet     = mdot / (rho * A_inlet)            # [m/s]

theta_rad   = np.radians(theta_swirl)
V_axial     = V_inlet * np.cos(theta_rad)       # component along toroid axis (into dome)
V_swirl     = V_inlet * np.sin(theta_rad)       # tangential component

# --- Rectangular cross-section area ---
A_rect      = rect_width * rect_height          # [m²]

# --- Curved top: ellipse arc optimization ---
# Ellipse: x²/a² + y²/b² = 1, arc from x=-w/2 to x=+w/2 at height h above rect top
# a = w/2 (fixed, must pass through corners), b = free parameter (ellipse height)
# Objective: minimize centrifugal ΔP ~ ½ρ V_swirl_dome²
# Constraint: b > 0 (must form a closed curve)
# Area of ellipse half = π*a*b/2

a_ellipse   = rect_width / 2                    # semi-axis in x = w/2

def cross_section_area(b):
    """Total D-shape area for a given ellipse semi-axis b."""
    A_ellipse_half = np.pi * a_ellipse * b / 2
    return A_rect + A_ellipse_half

def dome_velocity(b):
    """Mean axial velocity inside dome cross-section."""
    return mdot / (rho * cross_section_area(b))

def centrifugal_dP(b):
    """Centrifugal pressure gradient across dome radius due to swirl.
    ΔP_cent ~ ½ρ(V_swirl_dome)² — conservative upper bound.
    V_swirl_dome estimated by angular momentum conservation:
    R_mean * V_swirl_inlet = R_mean * V_swirl_dome  (same radius, just decays)
    We scale by area ratio as swirl decays proportionally."""
    V_s_dome = V_swirl * (A_inlet / cross_section_area(b))
    return 0.5 * rho * V_s_dome**2

# Minimize centrifugal ΔP → maximize cross-section area → maximize b
# But b can't be arbitrarily large (manufacturability: b ≤ w)
# So we minimize centrifugal ΔP subject to b ≤ rect_width
result = minimize_scalar(
    centrifugal_dP,
    bounds=(0.001, rect_width),
    method='bounded'
)
# centrifugal_dP is monotonically decreasing in b → optimal b = rect_width (upper bound)
# But we also check the constraint that ΔP_cent / ΔP_injector < 5%
b_max = rect_width  # start with max allowed

# Find minimum b that satisfies ΔP_cent < 5% of dP_injector
target_fraction = 0.05
b_values = np.linspace(0.001, rect_width * 2, 5000)
fractions = np.array([centrifugal_dP(b) / dP_injector for b in b_values])
valid = b_values[fractions < target_fraction]

if len(valid) > 0:
    b_min_required = valid[0]   # smallest b that satisfies constraint
else:
    b_min_required = None



# --- Final cross-section with optimal b ---
A_curve     = np.pi * a_ellipse * b_optimal / 2
A_total     = A_rect + A_curve

# --- Dome velocities ---
V_dome_mean = mdot / (rho * A_total)
V_swirl_dome = V_swirl * (A_inlet / A_total)    # scaled swirl in dome

# --- Pressure checks ---
dP_cent     = 0.5 * rho * V_swirl_dome**2
dP_dynamic  = 0.5 * rho * V_dome_mean**2
frac_cent   = dP_cent   / dP_injector * 100     # %
frac_dyn    = dP_dynamic / dP_injector * 100    # %

# --- Toroid volume ---
# V = 2π * R_mean * A_cross  (Pappus theorem)
V_torus     = 2 * np.pi * R_mean * A_total

# --- Per-hole flow ---
mdot_per_hole = mdot / n_holes
A_hole      = np.pi * (d_hole / 2)**2
V_hole      = mdot_per_hole / (rho * A_hole)

# --- Orifice equation check: what ΔP do holes actually need? ---
# mdot_hole = Cd * A_hole * sqrt(2 * rho * dP_hole)
dP_hole_required = (mdot_per_hole / (Cd_hole * A_hole))**2 / (2 * rho)

# --- Velocity ratio ---
V_ratio = V_dome_mean / V_inlet

# =============================================================================
# ── PRINT RESULTS ────────────────────────────────────────────────────────────
# =============================================================================

sep = "=" * 60
print(sep)
print("  LOX TOROIDAL MANIFOLD — SIZING RESULTS")
print(sep)

print("\n── INLET ──────────────────────────────────────────────────")
print(f"  Mass flow rate         : {mdot:.3f} kg/s")
print(f"  Inlet pipe diameter    : {D_inlet*1000:.1f} mm")
print(f"  Inlet area             : {A_inlet*1e4:.3f} cm²")
print(f"  Inlet velocity         : {V_inlet:.2f} m/s")
print(f"  Swirl angle            : {theta_swirl:.1f} °")
print(f"  Axial component        : {V_axial:.2f} m/s")
print(f"  Tangential component   : {V_swirl:.2f} m/s")

print("\n── TOROID CROSS-SECTION (D-SHAPE) ─────────────────────────")
print(f"  Rectangle width        : {rect_width*1000:.1f} mm")
print(f"  Rectangle height       : {rect_height*1000:.1f} mm")
print(f"  Rectangle area         : {A_rect*1e4:.3f} cm²")
print(f"  Curved top (ellipse)   :")
print(f"    Semi-axis a (x)      : {a_ellipse*1000:.1f} mm  (= w/2, fixed)")
print(f"    Semi-axis b (y)      : {b_optimal*1000:.1f} mm  (optimized)")
print(f"    Curve area           : {A_curve*1e4:.3f} cm²")
print(f"  TOTAL cross-section    : {A_total*1e4:.3f} cm²")
print(f"  Area ratio A/A_inlet   : {A_total/A_inlet:.1f}×")

print("\n── TOROID GEOMETRY ────────────────────────────────────────")
print(f"  Mean toroid radius     : {R_mean*1000:.1f} mm")
print(f"  Inner radius           : {(R_mean - rect_width/2)*1000:.1f} mm")
print(f"  Outer radius           : {(R_mean + rect_width/2)*1000:.1f} mm")
print(f"  Toroid volume          : {V_torus*1e6:.2f} cm³")

print("\n── VELOCITIES INSIDE DOME ─────────────────────────────────")
print(f"  Mean axial velocity    : {V_dome_mean:.3f} m/s")
print(f"  Swirl velocity (dome)  : {V_swirl_dome:.3f} m/s")
print(f"  Velocity ratio V/V_in  : {V_ratio:.3f}  (target < 0.2)")
if V_ratio < 0.2:
    print(f"  ✓ Velocity ratio OK")
else:
    print(f"  ✗ Velocity ratio too HIGH — increase cross-section area")

print("\n── PRESSURE CHECK ─────────────────────────────────────────")
print(f"  Injector ΔP            : {dP_injector/1e5:.2f} bar")
print(f"  Centrifugal ΔP (dome)  : {dP_cent/1e5:.4f} bar  ({frac_cent:.2f}% of inj. ΔP)")
print(f"  Dynamic pressure(dome) : {dP_dynamic/1e5:.4f} bar  ({frac_dyn:.2f}% of inj. ΔP)")
if frac_cent < 5.0:
    print(f"  ✓ Centrifugal ΔP within 5% limit")
else:
    print(f"  ✗ Centrifugal ΔP exceeds 5% — reduce swirl angle or increase dome area")
if b_min_required is not None:
    print(f"  Min ellipse b for <5%  : {b_min_required*1000:.1f} mm")

print("\n── INJECTOR HOLES ─────────────────────────────────────────")
print(f"  Number of holes        : {n_holes}")
print(f"  Hole diameter          : {d_hole*1000:.2f} mm")
print(f"  Cd                     : {Cd_hole}")
print(f"  Flow per hole          : {mdot_per_hole*1000:.3f} g/s")
print(f"  Velocity in hole       : {V_hole:.2f} m/s")
print(f"  ΔP required by holes   : {dP_hole_required/1e5:.3f} bar")
print(f"  Available ΔP (inj.)    : {dP_injector/1e5:.2f} bar")
if dP_hole_required <= dP_injector:
    print(f"  ✓ Holes can pass required flow")
else:
    print(f"  ✗ Holes too small — increase d_hole or Cd")

print(f"\n── DESIGN SCORECARD ───────────────────────────────────────")
checks = {
    "Velocity ratio < 0.2"        : V_ratio < 0.2,
    "Centrifugal ΔP < 5% inj. ΔP" : frac_cent < 5.0,
    "Dynamic ΔP < 5% inj. ΔP"     : frac_dyn  < 5.0,
    "Holes can pass flow"          : dP_hole_required <= dP_injector,
    "b ≥ b_min_required"           : (b_min_required is None) or (b_optimal >= b_min_required),
}
for label, ok in checks.items():
    print(f"  {'✓' if ok else '✗'}  {label}")
n_pass = sum(checks.values())
print(f"\n  Score: {n_pass}/{len(checks)} checks passed", "— GOOD TO GO" if n_pass == len(checks) else "— NEEDS ATTENTION")
print(f"{sep}")

# =============================================================================
# ── VISUALIZATION ─────────────────────────────────────────────────────────────
# =============================================================================

fig, axes = plt.subplots(1, 2, figsize=(14, 6))
fig.suptitle("LOX Toroidal Manifold — Sizing", fontsize=14, fontweight='bold')

# ── Plot 1: D-shape cross-section ──
ax1 = axes[0]
ax1.set_title("Toroid Cross-Section (D-shape)", fontsize=12)
ax1.set_aspect('equal')

w = rect_width
h = rect_height
b = b_optimal
a = a_ellipse

# Rectangle corners
rect_x = [-w/2, w/2, w/2, -w/2, -w/2]
rect_y = [0,    0,   h,    h,    0   ]
ax1.fill(rect_x, rect_y, color='#AED6F1', alpha=0.7, label='Rectangle')
ax1.plot(rect_x, rect_y, 'b-', linewidth=1.5)

# Ellipse top arc
theta_arc = np.linspace(np.pi, 0, 300)    # from left to right across top
ex = a * np.cos(theta_arc)
ey = h + b * np.sin(theta_arc)            # offset up by h (rect height)
ax1.fill_between(ex, h, ey, color='#85C1E9', alpha=0.7, label='Ellipse cap')
ax1.plot(ex, ey, 'b-', linewidth=1.5)

# Dimensions
ax1.annotate('', xy=(w/2, -0.005), xytext=(-w/2, -0.005),
             arrowprops=dict(arrowstyle='<->', color='black'))
ax1.text(0, -0.008, f'w = {w*1000:.0f} mm', ha='center', va='top', fontsize=9)

ax1.annotate('', xy=(w/2+0.005, h), xytext=(w/2+0.005, 0),
             arrowprops=dict(arrowstyle='<->', color='black'))
ax1.text(w/2+0.008, h/2, f'h = {h*1000:.0f} mm', ha='left', va='center', fontsize=9)

ax1.annotate('', xy=(w/2+0.005, h+b), xytext=(w/2+0.005, h),
             arrowprops=dict(arrowstyle='<->', color='#1a6ea8'))
ax1.text(w/2+0.008, h+b/2, f'b = {b*1000:.0f} mm', ha='left', va='center',
         fontsize=9, color='#1a6ea8')

ax1.set_xlim(-w/2 - 0.025, w/2 + 0.04)
ax1.set_ylim(-0.015, h + b + 0.015)
ax1.set_xlabel('Width [m]')
ax1.set_ylabel('Height [m]')
ax1.legend(loc='upper right', fontsize=8)
ax1.grid(True, alpha=0.3)

# ── Plot 2: Top view of toroid with inlet swirl angle ──
ax2 = axes[1]
ax2.set_title("Toroid Top View — Inlet Swirl Angle", fontsize=12)
ax2.set_aspect('equal')

# Draw toroid as annulus
R_inner = R_mean - rect_width / 2
R_outer = R_mean + rect_width / 2
theta_full = np.linspace(0, 2 * np.pi, 500)

outer_patch = Circle((0, 0), R_outer, color='#AED6F1', alpha=0.5)
inner_patch = Circle((0, 0), R_inner, color='white')
ax2.add_patch(outer_patch)
ax2.add_patch(inner_patch)
# Draw rings
for R in [R_inner, R_mean, R_outer]:
    lw = 1.0 if R != R_mean else 1.5
    ls = '-' if R != R_mean else '--'
    col = 'blue' if R != R_mean else 'gray'
    ax2.plot(R * np.cos(theta_full), R * np.sin(theta_full),
             color=col, linewidth=lw, linestyle=ls)

# Inlet position: place at angle 0 (right side), tangent to toroid
inlet_angle = 0.0  # angle on toroid where inlet is placed
ix = R_mean * np.cos(inlet_angle)
iy = R_mean * np.sin(inlet_angle)

# Radial direction at inlet point
radial = np.array([np.cos(inlet_angle), np.sin(inlet_angle)])
tangential = np.array([-np.sin(inlet_angle), np.cos(inlet_angle)])  # CCW tangent

# Inlet pipe direction vector = cos(theta)*radial + sin(theta)*tangential  (into toroid)
# theta_swirl = angle from radial toward tangential
inlet_dir = -np.cos(theta_rad) * radial + np.sin(theta_rad) * tangential
inlet_dir = inlet_dir / np.linalg.norm(inlet_dir)

pipe_len = 0.04
ax2.annotate('', xy=(ix, iy),
             xytext=(ix - inlet_dir[0]*pipe_len, iy - inlet_dir[1]*pipe_len),
             arrowprops=dict(arrowstyle='->', color='red', lw=2))
ax2.plot([ix - inlet_dir[0]*pipe_len, ix], [iy - inlet_dir[1]*pipe_len, iy],
         'r-', linewidth=2, label='Inlet flow direction')

# Mark inlet point
ax2.plot(ix, iy, 'ro', markersize=8)

# Draw angle arc
arc_r = 0.025
arc_theta = np.linspace(np.pi + inlet_angle, np.pi + inlet_angle + theta_rad, 50)
ax2.plot(ix + arc_r * np.cos(arc_theta), iy + arc_r * np.sin(arc_theta), 'k-', lw=1)
mid_arc = np.pi + inlet_angle + theta_rad / 2
ax2.text(ix + (arc_r+0.008) * np.cos(mid_arc), iy + (arc_r+0.008) * np.sin(mid_arc),
         f'{theta_swirl:.0f}°', fontsize=9, ha='center', color='black')

# Annotations
ax2.text(ix + 0.01, iy - 0.012, 'Inlet', color='red', fontsize=9)
ax2.text(0, 0, f'R = {R_mean*1000:.0f} mm', ha='center', va='center',
         fontsize=9, color='gray')

ax2.set_xlim(-(R_outer + 0.06), R_outer + 0.06)
ax2.set_ylim(-(R_outer + 0.06), R_outer + 0.06)
ax2.set_xlabel('X [m]')
ax2.set_ylabel('Y [m]')
ax2.legend(loc='upper right', fontsize=8)
ax2.grid(True, alpha=0.3)

plt.tight_layout()
plt.savefig('lox_manifold_sizing.png', dpi=150, bbox_inches='tight')
plt.show()
print("\nPlot saved to lox_manifold_sizing.png")
