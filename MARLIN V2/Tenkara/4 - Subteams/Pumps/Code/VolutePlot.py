"""
Volute plan view + nested cross-sections, from A(theta).

Constant-velocity volute: A(theta) = (theta/360) * Q / c3, so the spiral is
algebraic -- no integration. For a cross-section that scales self-similarly,
one constant turns area into a radial extent:

    diamond, aspect k = b_axial / a_radial:   A = 2*a*b = 2*k*a^2
        ->  a(theta) = sqrt( A(theta) / (2k) )
        ->  r_outer(theta) = r_3 + 2*a(theta)

Swap `SectionHalfWidth` for another shape and nothing else changes.
"""
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.path import Path
from matplotlib.patches import PathPatch
from matplotlib.colors import LinearSegmentedColormap

# ---------------------------------------------------------------- inputs (mm, m^2)
r_2       = 26.67       # impeller outlet radius
gap       = 0.12        # (D3 - D2)/D2, Stepanoff Fig 7.5
A_throat  = 54.92       # mm^2, = Q / c3
k_aspect  = 1.0         # b_axial / a_radial; 1.0 = square rotated 45 deg
e_sp_deg  = 360.0       # wrap angle of one channel
n_sections = 8          # how many sections to draw in the right-hand panel

r_3 = r_2 * (1 + gap)

def AreaAtAngle(th_deg):
    """(theta/360) * A at the throat station -- the constant-velocity rule."""
    return A_throat * (th_deg / e_sp_deg)

def SectionHalfWidth(A):
    """Radial half-extent of a diamond of area A. Shape lives here and nowhere else."""
    return np.sqrt(A / (2 * k_aspect))

# ---------------------------------------------------------------- style
INK, MUTED, FAINT = '#1c1b1a', '#6b6966', '#c9c5be'
RAMP = LinearSegmentedColormap.from_list('vol', ['#cfe0f5', '#5b8dd9', '#1b3f78'])

fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(13.0, 6.0),
                               gridspec_kw={'width_ratios': [1.15, 1]})
fig.patch.set_facecolor('#fbfaf8')
for ax in (ax1, ax2):
    ax.set_facecolor('#fbfaf8')
    for s in ax.spines.values():
        s.set_color(FAINT)
    ax.tick_params(colors=MUTED, labelsize=9)

# ================================================================ PLAN VIEW
th = np.linspace(0, e_sp_deg, 721)
a  = SectionHalfWidth(AreaAtAngle(th))
r_out = r_3 + 2 * a
rad = np.radians(th)

# filled band between base circle and outer wall, shaded by theta
for i in range(len(th) - 1):
    s = slice(i, i + 2)
    verts = np.concatenate([
        np.column_stack([r_3 * np.cos(rad[s]),    r_3 * np.sin(rad[s])]),
        np.column_stack([r_out[s][::-1] * np.cos(rad[s][::-1]),
                         r_out[s][::-1] * np.sin(rad[s][::-1])])])
    ax1.add_patch(PathPatch(Path(verts), facecolor=RAMP(th[i] / e_sp_deg),
                            edgecolor='none', alpha=0.9, zorder=1))

ax1.plot(r_out * np.cos(rad), r_out * np.sin(rad), color=INK, lw=2, zorder=4,
         label='volute outer wall')
ax1.plot(r_3 * np.cos(rad), r_3 * np.sin(rad), color=MUTED, lw=1.6, zorder=4,
         label=f'base circle  ⌀{2*r_3:.1f} mm')
ax1.plot(r_2 * np.cos(rad), r_2 * np.sin(rad), color=FAINT, lw=1.6, ls='--',
         zorder=4, label=f'impeller OD  ⌀{2*r_2:.1f} mm')

# tongue: the step between the throat (theta=360) and the base circle (theta=0)
ax1.plot([r_3, r_out[-1]], [0, 0], color=INK, lw=3.5, solid_capstyle='round', zorder=5)
ax1.annotate('tongue', xy=(r_3 + a[-1], 0), xytext=(r_3 + 2*a[-1] + 13, -13),
             color=INK, fontsize=10, ha='left',
             arrowprops=dict(arrowstyle='-', color=INK, lw=1))
ax1.plot(0, 0, marker='+', color=MUTED, ms=11, mew=1.4)

ax1.set_aspect('equal'); ax1.axis('off')
ax1.set_title('Plan view — spiral grows with collected flow',
              color=INK, fontsize=12, pad=12, loc='left')
ax1.legend(loc='lower left', frameon=False, fontsize=9, labelcolor=MUTED,
           bbox_to_anchor=(-0.06, -0.04))

# ================================================================ SECTIONS
for i in range(n_sections, 0, -1):
    t = e_sp_deg * i / n_sections
    ai = SectionHalfWidth(AreaAtAngle(t)); bi = k_aspect * ai
    cx = r_3 + ai
    ax2.add_patch(plt.Polygon([[r_3, 0], [cx, bi], [r_3 + 2*ai, 0], [cx, -bi]],
                              closed=True, facecolor=RAMP(t / e_sp_deg),
                              edgecolor=INK, lw=1.1, alpha=0.55, zorder=2))
    ax2.annotate(f'{t:.0f}°', xy=(cx, bi), xytext=(5, 3),      # top vertex: they fan out
                 textcoords='offset points', fontsize=8.5, color=MUTED)

ax2.axvline(r_3, color=MUTED, lw=1.6, zorder=3)
ax2.axvline(r_2, color=FAINT, lw=1.6, ls='--', zorder=3)
ax2.axhline(0, color=FAINT, lw=1, zorder=1)
ax2.text(r_3 - 0.7, -k_aspect * SectionHalfWidth(A_throat) * 1.35, 'base circle',
         rotation=90, va='bottom', ha='right', fontsize=9, color=MUTED)
ax2.text(r_2 - 0.7, -k_aspect * SectionHalfWidth(A_throat) * 1.35, 'impeller OD',
         rotation=90, va='bottom', ha='right', fontsize=9, color=FAINT)

a_max = SectionHalfWidth(A_throat)
ax2.set_xlim(r_2 - 2, r_3 + 2 * a_max + 3)
ax2.set_ylim(-k_aspect * a_max * 1.30, k_aspect * a_max * 1.55)
ax2.set_xlabel('radial position  (mm)', color=MUTED, fontsize=10)
ax2.set_ylabel('axial position  (mm)', color=MUTED, fontsize=10)
ax2.set_aspect('equal')
ax2.grid(True, color=FAINT, lw=0.6, alpha=0.5); ax2.set_axisbelow(True)
ax2.set_title(f'Cross-sections — throat {A_throat/100:.3f} cm² at {e_sp_deg:.0f}°',
              color=INK, fontsize=12, pad=12, loc='left')

fig.suptitle('Constant-velocity volute geometry', color=INK, fontsize=15,
             fontweight='bold', x=0.065, ha='left', y=0.975)
fig.tight_layout(rect=[0, 0.02, 1, 0.93])
fig.savefig('volute_geometry.png', dpi=170, facecolor=fig.get_facecolor())
plt.show()

r"""
Volute cross-sections only.  White background, no axes, no grid.

        apex            diamond: square on corner, side a(theta),
         /\             lower vertex truncated at the neck width b_3
        /  \
       /    \           A(theta) = a^2 + b_3*h_gap - b_3^2/4
       \    /                    = Throat_Area * theta/360
       |    |           neck: b_3 wide (axial) x h_gap tall (radial)
       |____|
"""
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.patches import Polygon
from matplotlib.colors import LinearSegmentedColormap

# ----------------------------------------------------------------- parameters
impeller_radius = 26.65      # mm   (d2 / 2)
clearance       = 0.2032     # mm   (0.008 in)
b_3             = 2.4        # mm   axial width of the gap/neck
h_gap           = 1.0        # mm   radial height of the gap/neck
Throat_Area     = 50.71      # mm^2 at theta = 360 deg     <- LOX: use 54.9

SHOW_LABELS = False          # the theta tags on each section
SHOW_REFS   = False          # impeller OD / neck floor reference lines

# ------------------------------------------------------------- derived radii
R0   = impeller_radius + clearance     # neck floor
Rtop = R0 + h_gap                      # neck ceiling = where the diamond is b_3 wide
Rbv  = Rtop - b_3 / 2                  # diamond's virtual bottom vertex (CONSTANT)

theta_min = ((b_3*h_gap - b_3**2/4 + b_3**2/2) / Throat_Area) * 360.0

def side_length(theta_deg):
    """Rhombus side a as a function of wrap angle."""
    return np.sqrt(Throat_Area*theta_deg/360.0 + b_3**2/4 - b_3*h_gap)

def profile(theta_deg):
    """Closed cross-section outline as (radial, axial) vertices."""
    hw = side_length(theta_deg) / np.sqrt(2)      # diamond half-width (axial)
    return np.array([
        (R0,        -b_3/2),          # neck floor, -axial
        (R0,         b_3/2),          # neck floor, +axial
        (Rtop,       b_3/2),          # neck ceiling, +axial
        (Rbv + hw,   hw),             # diamond side vertex, +axial
        (Rbv + 2*hw, 0.0),            # apex
        (Rbv + hw,  -hw),             # diamond side vertex, -axial
        (Rtop,      -b_3/2),          # neck ceiling, -axial
    ])

# --------------------------------------------------------------------- plot
shown = [theta_min, 90, 180, 270, 360]        # wrap angles to draw
cmap  = LinearSegmentedColormap.from_list("v", ["#c3d7ef", "#2c5f8a"])

fig, ax = plt.subplots(figsize=(7.5, 7), facecolor="white")
ax.set_facecolor("white")

for i, t in enumerate(shown):
    pr = profile(t)
    ax.add_patch(Polygon(pr, closed=True, facecolor=cmap(i/(len(shown)-1)),
                         edgecolor="#5a6671", lw=.9, zorder=len(shown)-i))
    if SHOW_LABELS:
        mid = (pr[3] + pr[4]) / 2                 # outer diamond edge midpoint
        ax.annotate(f"{t:.1f}°" if t == theta_min else f"{t:.0f}°",
                    xy=mid, xytext=(4, 4), textcoords="offset points",
                    fontsize=10, fontweight="bold", color=cmap(i/(len(shown)-1)),
                    ha="left", va="bottom", zorder=50)

if SHOW_REFS:
    ax.axvline(impeller_radius, ls="--", color="#b4bcc4", lw=1.3)
    ax.axvline(R0,              ls="-",  color="#6b7680", lw=1.3)

allp = np.vstack([profile(t) for t in shown])
ax.set_xlim(allp[:, 0].min()-1.5, allp[:, 0].max()+1.5)
ax.set_ylim(allp[:, 1].min()-1.5, allp[:, 1].max()+1.5)
ax.set_aspect("equal")
ax.axis("off")                                    # no grid, no ticks, no frame

plt.tight_layout(pad=0.2)
#plt.savefig("volute_sections.png", dpi=200, facecolor="white",
#            bbox_inches="tight", pad_inches=0.1)
plt.show()
print("saved volute_sections.png")