"""
Bell Nozzle Contour Generator
==============================
Thrust-optimised parabolic (TOP) nozzle contour.

Reference:
    "The Thrust Optimised Parabolic Nozzle"
    http://www.aspirespace.org.uk/downloads/Thrust%20optimised%20parabolic%20nozzle.pdf

Unit convention (internal):
    All geometry is computed in MILLIMETRES.
    Conversion to inches happens only at plot / export time.

Fixes applied vs. original script
-----------------------------------
  [F1] export_nozzle_dxf divided by 1000 (treated mm as m) → changed to /25.4 (mm→in)
  [F2] find_wall_angles used len(x_index) on a scalar     → changed to len(aratio)
  [F3] bell_nozzle called find_wall_angles(... throat_radius ...) using the
       global variable instead of the local parameter Rt  → now passes Rt
  [F4] export_nozzle_csv had a hardcoded Windows path as default → uses _here
  [F5] DXF segments reversed but CSV segments not reversed → both now consistently
       oriented throat→exit (chamber wall reversed, others forward)
  [F6] con_len labelled as "chamber length" but is the cylindrical-only length
       → label updated in the plot
  [F7] l_percent fallback was silent → now prints a warning
"""

import math
import os
import csv

import numpy as np
import matplotlib.pyplot as plt
import yaml
import ezdxf

from bisect import bisect_left
from matplotlib.patches import Arc


# ──────────────────────────────────────────────────────────────────────────────
# Wall-angle empirical tables (Rao, 1958)
# ──────────────────────────────────────────────────────────────────────────────

_ARATIO    = [  4,    5,   10,   20,   30,   40,   50,  100]
_TN_60     = [26.5, 28.0, 32.0, 35.0, 36.2, 37.1, 35.0, 40.0]
_TN_80     = [21.5, 23.0, 26.3, 28.8, 30.0, 31.0, 31.5, 33.5]
_TN_90     = [20.0, 21.0, 24.0, 27.0, 28.5, 29.5, 30.2, 32.0]
_TE_60     = [20.5, 20.5, 16.0, 14.5, 14.0, 13.5, 13.0, 11.2]
_TE_80     = [14.0, 13.0, 11.0,  9.0,  8.5,  8.0,  7.5,  7.0]
_TE_90     = [11.5, 10.5,  8.0,  7.0,  6.5,  6.0,  6.0,  6.0]

_SUPPORTED_L_PERCENT = {60, 80, 90}

N_PTS = 100  # number of points per contour segment


# ──────────────────────────────────────────────────────────────────────────────
# Helper: linear interpolation
# ──────────────────────────────────────────────────────────────────────────────

def _interpolate(x_list, y_list, x):
    """1-D linear interpolation / extrapolation clamp."""
    if any(y - xv <= 0 for xv, y in zip(x_list, x_list[1:])):
        raise ValueError("x_list must be strictly ascending.")
    if x <= x_list[0]:
        return y_list[0]
    if x >= x_list[-1]:
        return y_list[-1]
    i = bisect_left(x_list, x) - 1
    slope = (y_list[i + 1] - y_list[i]) / (x_list[i + 1] - x_list[i])
    return y_list[i] + slope * (x - x_list[i])


def _find_nearest(array, value):
    """Return (index, value) of the nearest element in array."""
    arr = np.asarray(array)
    idx = int((np.abs(arr - value)).argmin())
    return idx, arr[idx]


# ──────────────────────────────────────────────────────────────────────────────
# Wall-angle lookup
# ──────────────────────────────────────────────────────────────────────────────

def find_wall_angles(ar, Rt, l_percent=80):
    """
    Return (nozzle_length_mm, theta_n_rad, theta_e_rad) for a given
    area ratio ar and throat radius Rt [mm].

    [F2] fixed: was len(x_index) on a scalar — now correctly len(aratio)
    [F7] fixed: unsupported l_percent now prints a warning before defaulting
    """
    if l_percent not in _SUPPORTED_L_PERCENT:
        print(f"Warning: l_percent={l_percent} is not supported "
              f"(valid: {_SUPPORTED_L_PERCENT}). Defaulting to 80 %.")
        l_percent = 80

    if   l_percent == 60:  theta_n_deg, theta_e_deg, Lnp = _TN_60, _TE_60, 0.6
    elif l_percent == 80:  theta_n_deg, theta_e_deg, Lnp = _TN_80, _TE_80, 0.8
    else:                  theta_n_deg, theta_e_deg, Lnp = _TN_90, _TE_90, 0.9

    f1 = ((math.sqrt(ar) - 1) * Rt) / math.tan(math.radians(15))
    Ln = Lnp * f1

    x_index, _ = _find_nearest(_ARATIO, ar)

    # exact match
    if round(_ARATIO[x_index], 1) == round(ar, 1):
        return Ln, math.radians(theta_n_deg[x_index]), math.radians(theta_e_deg[x_index])

    # interpolation — [F2] len(aratio) not len(x_index)
    n = len(_ARATIO)
    if x_index >= 2:
        sl = slice(x_index - 2, min(x_index + 2, n))
    elif n - x_index <= 1:
        sl = slice(max(x_index - 2, 0), n)       # [F2] was len(x_index)
    else:
        sl = slice(0, min(x_index + 2, n))

    ar_sl = _ARATIO[sl.start:sl.stop]
    tn_sl = theta_n_deg[sl.start:sl.stop]
    te_sl = theta_e_deg[sl.start:sl.stop]

    tn_val = _interpolate(ar_sl, tn_sl, ar)
    te_val = _interpolate(ar_sl, te_sl, ar)

    return Ln, math.radians(tn_val), math.radians(te_val)


# ──────────────────────────────────────────────────────────────────────────────
# Main contour generator   (all units: mm)
# ──────────────────────────────────────────────────────────────────────────────

def bell_nozzle(aratio, Rt, l_percent, cratio, alpha, Lc):
    """
    Generate the full nozzle / chamber contour.

    Parameters
    ----------
    aratio    : float  — nozzle expansion ratio Ae/At
    Rt        : float  — throat radius [mm]
    l_percent : int    — bell length as % of 15° cone equivalent (60/80/90)
    cratio    : float  — chamber contraction ratio Ac/At
    alpha     : float  — convergent half-angle [degrees]
    Lc        : float  — cylindrical chamber length [mm]

    Returns
    -------
    angles  : (Ln, theta_n, theta_e)
    contour : tuple of 18 lists (x, y, -y) × 6 segments
    R2      : float — convergent-arc radius [mm]

    [F3] fixed: now passes Rt (local parameter) to find_wall_angles,
         not the global variable throat_radius
    """
    # [F3] was: find_wall_angles(aratio, throat_radius, l_percent)
    angles = find_wall_angles(aratio, Rt, l_percent)
    nozzle_length, theta_n, theta_e = angles

    Rc          = Rt * math.sqrt(cratio)
    alpha_rad   = math.radians(alpha)
    entrant_rad = math.radians(-90 - alpha)   # typically −135° → −90°

    # ── 1. Throat entrant arc  (radius = 1.5 Rt, sweep from entrant_rad → −π/2)
    ang = np.linspace(entrant_rad, -math.pi / 2, N_PTS)
    xe  = [1.5 * Rt * math.cos(a)               for a in ang]
    ye  = [1.5 * Rt * math.sin(a) + 2.5 * Rt    for a in ang]

    # ── 2. Convergent diagonal (straight section)
    R2max   = (Rc - ye[0]) / (1 - math.cos(alpha_rad))
    R2      = 0.300 * R2max
    diag_x0 = xe[0];   diag_y0 = ye[0]
    diag_yf = Rc - R2 * (1 - math.cos(alpha_rad))
    diag_xf = (diag_yf - diag_y0) / (-math.tan(alpha_rad)) + diag_x0
    iters   = np.linspace(diag_x0, diag_xf, N_PTS)
    xed     = list(iters)
    yed     = [diag_y0 + abs(i - diag_x0) * math.tan(alpha_rad) for i in iters]

    # ── 3. Combustion-chamber convergent arc
    cca_start = math.radians(90 - alpha)
    cca_end   = math.pi / 2
    ang       = np.linspace(cca_start, cca_end, N_PTS)
    xeca      = [diag_xf + R2 * math.cos(a) - R2 * math.cos(cca_start) for a in ang]
    yeca      = [diag_yf + R2 * math.sin(a) - R2 * math.sin(cca_start) for a in ang]

    # ── 4. Combustion chamber cylinder
    iters = np.linspace(xeca[-1], -Lc, N_PTS)
    xecc  = list(iters)
    yecc  = [yeca[-1]] * N_PTS

    # ── 5. Throat exit arc  (radius = 0.382 Rt)
    ang  = np.linspace(-math.pi / 2, theta_n - math.pi / 2, N_PTS)
    xe2  = [0.382 * Rt * math.cos(a)              for a in ang]
    ye2  = [0.382 * Rt * math.sin(a) + 1.382 * Rt for a in ang]

    # ── 6. Bell section (quadratic Bézier)
    # N point  — end of throat-exit arc  [Eqn. 5]
    Nx = 0.382 * Rt * math.cos(theta_n - math.pi / 2)
    Ny = 0.382 * Rt * math.sin(theta_n - math.pi / 2) + 1.382 * Rt
    # E point  — nozzle exit  [Eqns. 2 & 3]
    Lnp = {60: 0.6, 80: 0.8, 90: 0.9}.get(l_percent, 0.8)
    Ex  = Lnp * (math.sqrt(aratio) - 1) * Rt / math.tan(math.radians(15))
    Ey  = math.sqrt(aratio) * Rt
    # Q point  — Bézier control point  [Eqns. 7–10]
    m1  = math.tan(theta_n);  m2 = math.tan(theta_e)
    C1  = Ny - m1 * Nx;       C2 = Ey - m2 * Ex
    Qx  = (C2 - C1)  / (m1 - m2)
    Qy  = (m1 * C2 - m2 * C1) / (m1 - m2)

    t_vals = np.linspace(0, 1, N_PTS)
    xbell  = [(1-t)**2 * Nx + 2*(1-t)*t * Qx + t**2 * Ex for t in t_vals]
    ybell  = [(1-t)**2 * Ny + 2*(1-t)*t * Qy + t**2 * Ey for t in t_vals]

    # ── Mirror for bottom half
    contour = (
        xe,    ye,    [-y for y in ye],
        xe2,   ye2,   [-y for y in ye2],
        xed,   yed,   [-y for y in yed],
        xeca,  yeca,  [-y for y in yeca],
        xecc,  yecc,  [-y for y in yecc],
        xbell, ybell, [-y for y in ybell],
    )
    return angles, contour, R2


# ──────────────────────────────────────────────────────────────────────────────
# Plotting helpers
# ──────────────────────────────────────────────────────────────────────────────

def _draw_angle_arc(ax, theta_rad, origin, label=r'$\theta$'):
    """Draw a small arc indicating angle theta_rad at origin."""
    length = 50
    sx, sy = origin
    ex = sx + math.cos(-theta_rad) * length * 0.5
    ey = sy + math.sin(-theta_rad) * length * 0.5
    ax.plot([sx, ex], [sy, ey], linewidth=0.5, color='k')
    arc = Arc([sx, sy], 1, 1, angle=0, theta1=0,
              theta2=math.degrees(theta_rad), color='k')
    ax.add_patch(arc)
    ax.text(sx + 0.5, sy + 0.5,
            label + ' = ' + str(round(math.degrees(theta_rad), 1)) + '°')


def _ring(r, h, a=0, n_theta=30, n_height=10):
    """Parametric ring surface for 3-D plot."""
    theta = np.linspace(0, 2 * math.pi, n_theta)
    v     = np.linspace(a, a + h, n_height)
    theta, v = np.meshgrid(theta, v)
    return r * np.cos(theta), r * np.sin(theta), v


def _set_axes_equal_3d(ax):
    limits = np.array([ax.get_xlim3d(), ax.get_ylim3d(), ax.get_zlim3d()])
    origin = np.mean(limits, axis=1)
    radius = 0.5 * np.max(np.abs(limits[:, 1] - limits[:, 0]))
    x, y, z = origin
    ax.set_xlim3d([x - radius, x + radius])
    ax.set_ylim3d([y - radius, y + radius])
    ax.set_zlim3d([z - radius, z + radius])


def plot_nozzle_mm(ax, title, Rt, angles, contour):
    """2-D annotated contour plot in mm."""
    nozzle_length, theta_n, theta_e = angles
    (xe, ye, nye, xe2, ye2, nye2,
     xed, yed, nyed, xeca, yeca, nyeca,
     xecc, yecc, nyecc, xbell, ybell, nybell) = contour

    ax.set_aspect('equal')
    ax.plot(xe,   ye,   linewidth=2.5, color='g', label='throat entrant')
    ax.plot(xe,   nye,  linewidth=2.5, color='g')
    ax.plot(xed,  yed,  linewidth=2.5, color='k', label='convergent diagonal')
    ax.plot(xed,  nyed, linewidth=2.5, color='k')
    ax.plot(xeca, yeca, linewidth=2.5, color='r', label='convergent arc')
    ax.plot(xeca, nyeca,linewidth=2.5, color='r')
    ax.plot(xecc, yecc, linewidth=2.5, color='b', label='chamber cylinder')
    ax.plot(xecc, nyecc,linewidth=2.5, color='b')
    ax.plot(xe2,  ye2,  linewidth=2.5, color='r', label='throat exit arc')
    ax.plot(xe2,  nye2, linewidth=2.5, color='r')
    ax.plot(xbell,ybell,linewidth=2.5, color='b', label='bell')
    ax.plot(xbell,nybell,linewidth=2.5,color='b')

    # annotations
    ax.annotate('', [xe[0], 0], [xe[0], nye[0]],
                arrowprops=dict(lw=0.5, arrowstyle='<-'))
    ax.text(xe[0], nye[0] / 2, f' Ri={round(abs(nye[0]), 1)}', fontsize=9)

    ax.annotate('', [0, 0], [xe[0], 0],
                arrowprops=dict(lw=0.5, arrowstyle='<-'))
    ax.text(xe[0], 0, f' Li={round(abs(xe[0]), 1)}', fontsize=9)

    ax.annotate('', [0, 2.5*Rt], [xe[int(N_PTS/2)], ye[int(N_PTS/2)]],
                arrowprops=dict(lw=0.5, arrowstyle='<-'))
    ax.text(xe[int(N_PTS/2)]/2, ye[int(N_PTS/2)]/2,
            f' 1.5Rt={round(1.5*Rt, 1)}', fontsize=9)

    _draw_angle_arc(ax, theta_n, [xe2[-1], nye2[-1] - 2], r'$\theta_n$')
    _draw_angle_arc(ax, theta_e, [xbell[-1], nybell[-1]], r'$\theta_e$')

    ax.axhline(color='black', lw=0.5, linestyle='dashed')
    ax.axvline(color='black', lw=0.5, linestyle='dashed')
    ax.grid(True, which='major', linestyle='-', linewidth=0.5)
    ax.grid(True, which='minor', linestyle=':', linewidth=0.5)
    ax.minorticks_on()
    ax.set_title(title, fontsize=9)
    ax.legend(fontsize=7, loc='upper left')


def plot3d(ax, contour):
    """3-D surface revolution plot."""
    (xe, ye, _, xe2, ye2, _, xed, yed, _,
     xeca, yeca, _, xecc, yecc, _, xbell, ybell, _) = contour

    x = np.concatenate([xecc, xeca, xed, xe, xe2, xbell])
    y = np.concatenate([yecc, yeca, yed, ye, ye2, ybell])

    thick = 5 * abs(x[1] - x[0])
    for xi, yi in zip(x, y):
        X, Y, Z = _ring(yi, thick, xi)
        ax.plot_surface(X, Y, Z, color='g')

    ax.set_box_aspect([1, 1, 1])
    _set_axes_equal_3d(ax)
    ax.view_init(-170, -15)


def plot_nozzle_inches(contour, angles, Dt_in, Dc_in, De_in, Lc_in, R2_mm, cangle):
    """
    Clean engineering drawing in inches.

    [F6] fixed: dimension label now says 'Cyl. chamber length' not 'Chamber length'
    """
    theta_n, theta_e = angles[1], angles[2]
    R2_in = R2_mm / 25.4

    # convert mm → inches
    def _in(arr):
        return np.array(arr) / 25.4

    xe,    ye,    nye    = _in(contour[0]),  _in(contour[1]),  _in(contour[2])
    xe2,   ye2,   nye2   = _in(contour[3]),  _in(contour[4]),  _in(contour[5])
    xed,   yed,   nyed   = _in(contour[6]),  _in(contour[7]),  _in(contour[8])
    xeca,  yeca,  nyeca  = _in(contour[9]),  _in(contour[10]), _in(contour[11])
    xecc,  yecc,  nyecc  = _in(contour[12]), _in(contour[13]), _in(contour[14])
    xbell, ybell, nybell = _in(contour[15]), _in(contour[16]), _in(contour[17])

    fig, ax = plt.subplots(figsize=(14, 9))

    for xs, ys in [(xe, ye), (xe, nye), (xed, yed), (xed, nyed),
                   (xeca, yeca), (xeca, nyeca), (xecc, yecc), (xecc, nyecc),
                   (xe2, ye2), (xe2, nye2), (xbell, ybell), (xbell, nybell)]:
        ax.plot(xs, ys, linewidth=2.5, color='k')

    FS = 18  # annotation font size

    # throat diameter
    ax.annotate('', [xe[-1], 0.95*nye[-1]], [xe[-1], 0.95*ye[-1]],
                arrowprops=dict(lw=1, arrowstyle='|-|'))
    ax.text(0.1, 0.1, f'{round(Dt_in, 2)} in', fontsize=FS)

    # chamber diameter
    ax.annotate('', [xecc[-1], 0.95*nyecc[-1]], [xecc[-1], 0.95*yecc[-1]],
                arrowprops=dict(lw=1, arrowstyle='|-|'))
    ax.text(xecc[-1] + 0.1, 0.1, f'{round(Dc_in, 2)} in', fontsize=FS)

    # exit diameter
    ax.annotate('', [xbell[-1], 0.95*nybell[-1]], [xbell[-1], 0.95*ybell[-1]],
                arrowprops=dict(lw=1, arrowstyle='|-|'))
    ax.text(4.0, 0.1, f'{round(De_in, 2)} in', fontsize=FS)

    dim_y = 1.25 * ybell[-1]

    # [F6] cylindrical chamber length (not total chamber length)
    cyl_len = Lc_in - abs(xecc[0])
    ax.annotate('', [xecc[-1] - 0.05, 1.2*ybell[-1]],
                    [xecc[0]  + 0.05, 1.2*ybell[-1]],
                arrowprops=dict(lw=1, arrowstyle='|-|'))
    ax.text(xecc[0], dim_y, f'Cyl. Lc={round(cyl_len, 2)} in', fontsize=FS)

    # convergent section length
    conv_len = abs(xecc[0])
    ax.annotate('', [xecc[0] - 0.05, 1.2*ybell[-1]],
                    [0.05,            1.2*ybell[-1]],
                arrowprops=dict(lw=1, arrowstyle='|-|'))
    ax.text(xecc[0] / 2, dim_y, f'Conv.={round(conv_len, 2)} in', fontsize=FS)

    # divergent section length
    ax.annotate('', [-0.05,           1.2*ybell[-1]],
                    [xbell[-1] + 0.05, 1.2*ybell[-1]],
                arrowprops=dict(lw=1, arrowstyle='|-|'))
    ax.text(xbell[-1] / 2, dim_y, f'Div.={round(xbell[-1], 2)} in', fontsize=FS)

    # theta_n arc
    tn_arc_x = [xe2[-1] + 2.0*math.cos(a) for a in np.linspace(0, theta_n, N_PTS)]
    tn_arc_y = [ye2[-1] + 2.0*math.sin(a) for a in np.linspace(0, theta_n, N_PTS)]
    ax.annotate('', [xe2[-1], ye2[-1]], [xe2[-1] + 2.5, ye2[-1]],
                arrowprops=dict(lw=1, arrowstyle='-', color='r'))
    ax.annotate('', [xe2[-1], ye2[-1]],
                    [xe2[-1] + 2.5*math.cos(theta_n), ye2[-1] + 2.5*math.sin(theta_n)],
                arrowprops=dict(lw=1, arrowstyle='-', color='r'))
    ax.plot(tn_arc_x, tn_arc_y, linewidth=1, color='r')
    ax.text(2.5, 1.75,
            r'$\theta_n$ = ' + f'{round(math.degrees(theta_n), 1)}°', fontsize=FS)

    # theta_e arc
    te_arc_x = [xbell[-1] + 1.5*math.cos(a) for a in np.linspace(0, theta_e, N_PTS)]
    te_arc_y = [ybell[-1] + 1.5*math.sin(a) for a in np.linspace(0, theta_e, N_PTS)]
    ax.annotate('', [xbell[-1], ybell[-1]], [xbell[-1] + 2.0, ybell[-1]],
                arrowprops=dict(lw=1, arrowstyle='-', color='r'))
    ax.annotate('', [xbell[-1], ybell[-1]],
                    [xbell[-1] + 2.0*math.cos(theta_e), ybell[-1] + 2.0*math.sin(theta_e)],
                arrowprops=dict(lw=1, arrowstyle='-', color='r'))
    ax.plot(te_arc_x, te_arc_y, linewidth=1, color='r')
    ax.text(6, 5,
            r'$\theta_e$ = ' + f'{round(math.degrees(theta_e), 1)}°', fontsize=FS)

    # R1, R2, Rn labels
    rt = Dt_in / 2.0
    r1_ang = math.radians((270 + 180 + cangle) / 2)
    ax.annotate('', [0, 2.5*rt],
                    [1.5*rt*math.cos(r1_ang), rt*(1.5*math.sin(r1_ang) + 2.5)],
                arrowprops=dict(lw=0.5, arrowstyle='<-', color='b'))
    ax.plot(0, 2.5*rt, '+', color='b')
    ax.text(-1.0, 2.75, r'$R_1$', fontsize=15)

    r2_ang = math.radians((90 + cangle) / 2)
    ax.annotate('', [xecc[0], yecc[0] - R2_in],
                    [xecc[0] + R2_in*math.cos(r2_ang),
                     yecc[0] + R2_in*(math.sin(r2_ang) - 1)],
                arrowprops=dict(lw=0.5, arrowstyle='<-', color='b'))
    ax.plot(xecc[0], yecc[0] - R2_in, '+', color='b')
    ax.text(-3.1, 2.1, r'$R_2$', fontsize=15)

    rn_ang = math.radians((270 + 270 + cangle) / 2)
    ax.annotate('', [0, 1.382*rt],
                    [0.382*rt*math.cos(rn_ang), rt*(0.382*math.sin(rn_ang) + 1.382)],
                arrowprops=dict(lw=0.5, arrowstyle='<-', color='b'))
    ax.plot(0, 1.382*rt, '+', color='b')
    ax.text(0.1, 2.0, r'$R_n$', fontsize=15)

    ax.axhline(color='black', lw=0.5, linestyle='dashed')
    ax.axvline(color='black', lw=0.5, linestyle='dashed')
    ax.grid(True, which='major', linestyle='-', linewidth=0.5)
    ax.grid(True, which='minor', linestyle=':', linewidth=0.5)
    ax.minorticks_on()
    ax.tick_params(axis='both', labelsize=15)
    ax.set_xlabel('Inches', fontsize=18)
    ax.set_ylabel('Inches', fontsize=18)
    ax.set_aspect('equal')
    fig.tight_layout(rect=[0, 0.03, 1, 0.95])

    _here = os.path.dirname(os.path.abspath(__file__))
    out   = os.path.join(_here, 'CSV_DXF_OUTPUTS', 'nozzle_contour_inches_plot.png')
    fig.savefig(out)
    print(f"Saved inches plot → {out}")
    plt.show()


def plot_overview(title, throat_radius, angles, contour):
    """Side-by-side 2-D + 3-D overview figure."""
    fig = plt.figure(1, figsize=(14, 9))
    ax1 = fig.add_subplot(121)
    plot_nozzle_mm(ax1, title, throat_radius, angles, contour)
    ax2 = fig.add_subplot(122, projection='3d')
    plot3d(ax2, contour)
    fig.tight_layout(rect=[0, 0.03, 1, 0.95])

    _here = os.path.dirname(os.path.abspath(__file__))
    out   = os.path.join(_here, 'CSV_DXF_OUTPUTS', 'nozzle_contour_plot.png')
    fig.savefig(out)
    print(f"Saved overview plot → {out}")
    plt.show()


# ──────────────────────────────────────────────────────────────────────────────
# Export helpers
# ──────────────────────────────────────────────────────────────────────────────

def export_nozzle_csv(contour, filename=None):
    """
    Write all contour segments to a single CSV (mm units).

    Segment order: chamber_wall → convergent_arc → convergent_diagonal
                   → throat_entrant → throat_exit → bell

    [F4] fixed: no longer has a hardcoded Windows path as default argument.
    [F5] fixed: segment orientation now matches the DXF (throat→exit).
                Chamber wall, convergent arc, and convergent diagonal are
                reversed here to be consistent.
    """
    if filename is None:
        _here    = os.path.dirname(os.path.abspath(__file__))
        filename = os.path.join(_here, 'CSV_DXF_OUTPUTS', 'nozzle_contour.csv')

    xecc,  yecc  = contour[12], contour[13]
    xeca,  yeca  = contour[9],  contour[10]
    xed,   yed   = contour[6],  contour[7]
    xe,    ye    = contour[0],  contour[1]
    xe2,   ye2   = contour[3],  contour[4]
    xbell, ybell = contour[15], contour[16]

    # [F5] reverse upstream segments so all run throat→exit
    segments = [
        ('chamber_wall',          list(reversed(xecc)),  list(reversed(yecc))),
        ('convergent_arc',        list(reversed(xeca)),  list(reversed(yeca))),
        ('convergent_diagonal',   list(reversed(xed)),   list(reversed(yed))),
        ('throat_entrant_arc',    xe,    ye),
        ('throat_exit_arc',       xe2,   ye2),
        ('bell',                  xbell, ybell),
    ]

    with open(filename, 'w', newline='') as f:
        w = csv.writer(f)
        w.writerow(['segment', 'x_mm', 'y_mm', 'index'])
        for name, xs, ys in segments:
            for i, (x, y) in enumerate(zip(xs, ys)):
                w.writerow([name, float(x), float(y), i])

    print(f"Saved CSV → {filename}")
    return filename


def export_nozzle_dxf(contour):
    """
    Write nozzle contour to DXF (splines, inches).

    [F1] fixed: now divides by 25.4 (mm→in) instead of 1000 (was treating mm as m).
    [F5] fixed: upstream segments reversed so all splines run throat→exit,
                matching the CSV export direction.
    """
    def _to_in(arr):
        return list(np.array(arr) / 25.4)

    xecc,  yecc  = _to_in(contour[12]), _to_in(contour[13])
    xeca,  yeca  = _to_in(contour[9]),  _to_in(contour[10])
    xed,   yed   = _to_in(contour[6]),  _to_in(contour[7])
    xe,    ye    = _to_in(contour[0]),  _to_in(contour[1])
    xe2,   ye2   = _to_in(contour[3]),  _to_in(contour[4])
    xbell, ybell = _to_in(contour[15]), _to_in(contour[16])

    # [F5] reverse upstream segments so all splines run throat→exit
    sections = [
        (list(reversed(xecc)),  list(reversed(yecc))),   # chamber wall
        (list(reversed(xeca)),  list(reversed(yeca))),   # convergent arc
        (list(reversed(xed)),   list(reversed(yed))),    # convergent diagonal
        (xe,    ye),                                      # throat entrant
        (xe2,   ye2),                                     # throat exit
        (xbell, ybell),                                   # bell
    ]

    doc = ezdxf.new('R2010')
    msp = doc.modelspace()
    for xs, ys in sections:
        pts = list(zip(xs, ys))
        if len(pts) >= 3:
            msp.add_spline(pts, degree=3)

    _here = os.path.dirname(os.path.abspath(__file__))
    out   = os.path.join(_here, 'CSV_DXF_OUTPUTS', 'nozzle_contour.dxf')
    doc.saveas(out)
    print(f"Saved DXF → {out}")


# ──────────────────────────────────────────────────────────────────────────────
# Entry point
# ──────────────────────────────────────────────────────────────────────────────

if __name__ == '__main__':
    _here = os.path.dirname(os.path.abspath(__file__))
    yaml_path = os.path.join(_here, '..', 'TCA_params.yaml')

    with open(yaml_path) as f:
        p = yaml.safe_load(f)

    # ── Read inputs from YAML ──────────────────────────────────────────────────
    l_percent     = p['bell_nozzle_l_percent']         # 60 / 80 / 90
    aratio        = p['tca_expansion_ratio']            # Ae / At
    cratio        = p['tca_contraction_ratio']          # Ac / At
    cangle        = p['tca_convergent_half_angle']      # degrees
    Lc_mm         = p['tca_chamber_length'] * 25.4      # in → mm
    throat_radius = p['tca_throat_diameter'] * 25.4 / 2 # in → mm

    # ── Generate contour (all in mm) ──────────────────────────────────────────
    angles, contour, R2 = bell_nozzle(
        aratio, throat_radius, l_percent, cratio, cangle, Lc_mm
    )

    title = (f'Bell Nozzle\n'
             f'[ε = {round(aratio, 1)}, '
             f'Rt = {round(throat_radius, 2)} mm, '
             f'L% = {l_percent}%]')

    # ── Exports ───────────────────────────────────────────────────────────────
    os.makedirs(os.path.join(_here, 'CSV_DXF_OUTPUTS'), exist_ok=True)
    export_nozzle_csv(contour)
    export_nozzle_dxf(contour)

    # ── Plots ─────────────────────────────────────────────────────────────────
    plot_overview(title, throat_radius, angles, contour)

    Dt_in = p['tca_throat_diameter']
    Dc_in = p['tca_chamber_diameter']
    De_in = p['tca_exit_diameter']
    Lc_in = p['tca_chamber_length']
    plot_nozzle_inches(contour, angles, Dt_in, Dc_in, De_in, Lc_in, R2, cangle)
