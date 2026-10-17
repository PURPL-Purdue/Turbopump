#   TODO: 
#   Integrate this with the impeller sizing (automatically wooooo)
#   integrate auto cross section dimension calcs

#   Assumed single volute
#   Diamond cross section?
#   Consider joined discharge double volute to reduce radial forces on rotor
#   Gulich p. 403:
#   If nevertheless a 6-bladed impeller (for example) is to be employed 
#   it is advisable to reduce the wrap angle of the inner volute εi to
#   165° or 170° so that two blades never pass the cutwaters simultaneously

from pint import Quantity as Q_
from dataclasses import dataclass
from matplotlib import pyplot as plt
import numpy as np
import math
import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parents[2]))	# Pumps folder, so the Code package imports when run directly
from Code.Impeller import LOX, IPA
# Constants
g = Q_(9.81, 'm/s**2')

@dataclass
class VoluteDesignPoint:
    N_shaft: Q_[float]
    Q: Q_[float]
    H: Q_[float]

@dataclass
class VoluteGeometry:
    e_sp: Q_[float] # Wrap angle
    n_volutes: int # Number of volutes

def getVoluteVelocity(K_v: float, H: Q_[float]) -> Q_[float]:   # Uses constant velocity Stepanoff method
    avg_V = K_v * np.sqrt(2 * g.to('ft/s**2') * H.to('feet'))
    return avg_V

def getAreaAtAngle(Q: Q_[float], avg_V: Q_[float], e: Q_[float], e_sp: Q_[float], n_volutes: int) -> Q_[float]:  # Get areas as a function of angle
    A = (e / e_sp) * (Q / (n_volutes * avg_V))
    return A

class Volute:
    def __init__(self, geometry: VoluteGeometry, design_point: VoluteDesignPoint):
        pass

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
# System inputs
#Q_design = Q_(1.78, 'liters/second')
#H_design = Q_(294.8, 'meters')

# Impeller inputs
#C_u2 =  Q_(34.516, 'm/s')   # tangential velocity component
#r2 = Q_(1.1, 'inches')      # impeller outlet radius
#C_2 = Q_(68.2, 'm/s')       # Absolute velocity

# Volute inputs
#e_sp = Q_(2 * math.pi, 'radians')          # volute wrap angle (2pi for single volutes)
#K_v = 0.55                                 # experimental average velocity coefficient from h&h equation 6-70 pg 220
                                           # and from https://ijsea.com/archive/volume8/issue8/IJSEA08081016.pdf 

# Throat area calculations
#avg_V = K_v * np.sqrt(2 * g.to('ft/s**2') * H_design.to('feet'))    # use empirical relation from h&h pg 220
#A_throat = Q_design / avg_V
#throat_diameter = np.sqrt((4 / math.pi) * A_throat)                 # if it was a circle
#HH_throatV = Q_design.to('gallons per minute') / (A_throat.to('in**2') * 3.12) # flow velocity at throat outlet according to H&H pg 222

#print(f"Throat area: {A_throat.to('cm**2'):.4f}")
#print(f"Throat diameter: {throat_diameter.to('cm'):.4f}")
#print(f"Average velocity: {avg_V.to('m/s'):.4f}")
#print(f"Throat velocity: {HH_throatV.to('m/s'):.4f}")

# From gulich pg 420:
#X_sp = (Q_design / (math.pi * C_u2 * r2)) * ((e_sp / Q_(2 * math.pi, 'radians')).to('dimensionless'))
