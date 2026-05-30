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
# Export helpers
# ──────────────────────────────────────────────────────────────────────────────

def export_nozzle_csv(contour, filename):
    """
    Write all contour segments to a single CSV (mm units).
    """
    if not filename:
        _here    = os.path.dirname(os.path.abspath(__file__))
        filename = os.path.join(_here, 'CSV_DXF_OUTPUTS', 'nozzle_contour.csv')
    else:
        os.makedirs(os.path.dirname(os.path.abspath(filename)), exist_ok=True)

    xecc,  yecc  = contour[12], contour[13]
    xeca,  yeca  = contour[9],  contour[10]
    xed,   yed   = contour[6],  contour[7]
    xe,    ye    = contour[0],  contour[1]
    xe2,   ye2   = contour[3],  contour[4]
    xbell, ybell = contour[15], contour[16]

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


def export_nozzle_dxf(contour, filename=None):
    """
    Write nozzle contour to DXF (splines, mm).
    """

    xecc,  yecc  = contour[12], contour[13]
    xeca,  yeca  = contour[9],  contour[10]
    xed,   yed   = contour[6],  contour[7]
    xe,    ye    = contour[0],  contour[1]
    xe2,   ye2   = contour[3],  contour[4]
    xbell, ybell = contour[15], contour[16]

    sections = [
        (list(reversed(xecc)),  list(reversed(yecc))),
        (list(reversed(xeca)),  list(reversed(yeca))),
        (list(reversed(xed)),   list(reversed(yed))),
        (xe,    ye),
        (xe2,   ye2),
        (xbell, ybell),
    ]

    doc = ezdxf.new('R2010')
    msp = doc.modelspace()
    for xs, ys in sections:
        pts = list(zip(xs, ys))
        if len(pts) >= 3:
            msp.add_spline(pts, degree=3)

    if not filename:
        _here    = os.path.dirname(os.path.abspath(__file__))
        filename = os.path.join(_here, 'CSV_DXF_OUTPUTS', 'nozzle_contour.dxf')
    else:
        os.makedirs(os.path.dirname(os.path.abspath(filename)), exist_ok=True)

    doc.saveas(filename)
    print(f"Saved DXF → {filename}")
    return filename


# ──────────────────────────────────────────────────────────────────────────────
# Function Caller Generates the nozzle contour and exports CSV and DXF files.
# ──────────────────────────────────────────────────────────────────────────────

def contour_script(lpercent, aratio, cratio, cangle, Lc_mm, throat_radius,output_filename):
    """
    Generate the full nozzle / chamber contour.

    Parameters
    ----------
    aratio          : float  — nozzle expansion ratio Ae/At
    throat_radius   : float  — throat radius [mm]
    lpercent        : int    — bell length as % of 15° cone equivalent (60/80/90)
    cratio          : float  — chamber contraction ratio Ac/At
    cangle          : float  — convergent half-angle [degrees]
    Lc_mm           : float  — cylindrical chamber length [mm]
    output_filename : str    — path/name for the exported CSV (optional)
    """

    angles, contour, R2 = bell_nozzle(aratio, throat_radius, lpercent, cratio, cangle, Lc_mm)

    if output_filename:
        base     = os.path.splitext(output_filename)[0]  # strip any extension
        csv_file = base + '.csv'
        dxf_file = base + '.dxf'
    else:
        csv_file = None
        dxf_file = None

    export_nozzle_csv(contour, filename=csv_file)
    export_nozzle_dxf(contour, filename=dxf_file)
    return angles, contour