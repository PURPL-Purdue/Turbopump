#!/usr/bin/env python3
"""
========================================================================
2-D TRANSIENT CYLINDRICAL HEAT TRANSFER SOLVER — ROCKET NOZZLE WALL
------------------------------------------------------------------------
Solves the 2-D unsteady conduction equation in cylindrical coordinates
over the rocket nozzle/combustion-chamber wall cross-section.

    rho*cp * dT/dt = 1/r * d/dr(r*k*dT/dr) + d/dz(k*dT/dz)

Author         : Rafael Macia Titos
Date           : 2026

Discretization : FVM, implicit backward-Euler
Linear solver  : Line-by-line TDMA (Thomas algorithm, alternating
                 radial / axial sweeps).  Numba-accelerated if available.

Can be run standalone or called via run() from run.py.

Usage (standalone):
    python HeatTransferSolver_2D.py
    # reads heat_transfer.yaml and TCA_params.yaml from known locations

Usage (from run.py):
    import HeatTransferSolver_2D
    HeatTransferSolver_2D.run(ht_params, tca_params, props_csv_path, script_dir)

------------------------------------------------------------------------
Boundary conditions (explicit, for reference):
    Inner radial face (i=0):
        FIRING  : Robin / Bartz hot-gas convection at h_g(z) and T_aw(z)
        COOLING : Robin / natural convection at h_nat_inner and T_amb
    Outer radial face (i=n-1):
        Robin / natural convection at h_nat_outer and T_amb (always)
    Axial faces (j=0 and j=m-1):
        ADIABATIC.  Encoded by aS[:,0]=0 and aN[:,-1]=0 in the coefficient
        assembly.  Suitable if the inlet/exit ends are well away from the
        hot zone and you don't care about leakage there.
    Step face between steel & aluminum sections:
        Robin / natural convection at h_nat_outer and T_amb on the exposed
        annular area where the outer radius drops.
========================================================================
"""

import os
import sys
import yaml
import numpy as np
import pandas as pd
import matplotlib
import matplotlib.pyplot as plt
import matplotlib.patches as patches
from scipy.interpolate import interp1d, CubicSpline
import time as time_module
import traceback

try:
    import imageio
    IMAGEIO_AVAILABLE = True
except ImportError:
    IMAGEIO_AVAILABLE = False
    print("WARNING: imageio not found. Video recording disabled.")
    print("         Install with: pip install imageio[ffmpeg]")

try:
    from tqdm import tqdm
    TQDM_AVAILABLE = True
except ImportError:
    TQDM_AVAILABLE = False

try:
    import numba
    NUMBA_AVAILABLE = True
except ImportError:
    NUMBA_AVAILABLE = False

LBM_PER_S_TO_KG_PER_S = 0.453592


# =============================================================================
# HELPER FUNCTIONS
# =============================================================================

def uniquetol(x, tol):
    """Return indices of unique sorted values within relative tolerance."""
    if len(x) == 0:
        return np.array([], dtype=int)
    scale   = np.max(np.abs(x))
    abs_tol = tol * scale if scale > 0 else tol
    idx = [0]
    for i in range(1, len(x)):
        if np.abs(x[i] - x[idx[-1]]) > abs_tol:
            idx.append(i)
    return np.array(idx, dtype=int)


def _thomas_r_np(aW, aP, aE, rhs, T_out, c_buf, d_buf):
    """Vectorized Thomas (TDMA) sweep in the radial direction (NumPy)."""
    n = aP.shape[0]
    c_buf[0] = -aE[0] / aP[0]
    d_buf[0] =  rhs[0] / aP[0]
    for i in range(1, n):
        denom    = aP[i] + aW[i] * c_buf[i - 1]
        c_buf[i] = -aE[i] / denom
        d_buf[i] = (rhs[i] + aW[i] * d_buf[i - 1]) / denom
    T_out[-1] = d_buf[-1]
    for i in range(n - 2, -1, -1):
        T_out[i] = d_buf[i] - c_buf[i] * T_out[i + 1]


def _thomas_z_np(aS, aP, aN, rhs, T_out, c_buf, d_buf):
    """Vectorized Thomas (TDMA) sweep in the axial direction (NumPy)."""
    m = aP.shape[1]
    c_buf[:, 0] = -aN[:, 0] / aP[:, 0]
    d_buf[:, 0] =  rhs[:, 0] / aP[:, 0]
    for j in range(1, m):
        denom       = aP[:, j] + aS[:, j] * c_buf[:, j - 1]
        c_buf[:, j] = -aN[:, j] / denom
        d_buf[:, j] = (rhs[:, j] + aS[:, j] * d_buf[:, j - 1]) / denom
    T_out[:, -1] = d_buf[:, -1]
    for j in range(m - 2, -1, -1):
        T_out[:, j] = d_buf[:, j] - c_buf[:, j] * T_out[:, j + 1]


if NUMBA_AVAILABLE:
    @numba.njit(cache=True, fastmath=True, parallel=True)
    def _thomas_r_nb(aW, aP, aE, rhs, T_out, c_buf, d_buf):
        n, m = aP.shape
        for j in numba.prange(m):
            c_buf[0, j] = -aE[0, j] / aP[0, j]
            d_buf[0, j] =  rhs[0, j] / aP[0, j]
            for i in range(1, n):
                denom       = aP[i, j] + aW[i, j] * c_buf[i - 1, j]
                c_buf[i, j] = -aE[i, j] / denom
                d_buf[i, j] = (rhs[i, j] + aW[i, j] * d_buf[i - 1, j]) / denom
            T_out[n - 1, j] = d_buf[n - 1, j]
            for i in range(n - 2, -1, -1):
                T_out[i, j] = d_buf[i, j] - c_buf[i, j] * T_out[i + 1, j]

    @numba.njit(cache=True, fastmath=True, parallel=True)
    def _thomas_z_nb(aS, aP, aN, rhs, T_out, c_buf, d_buf):
        n, m = aP.shape
        for i in numba.prange(n):
            c_buf[i, 0] = -aN[i, 0] / aP[i, 0]
            d_buf[i, 0] =  rhs[i, 0] / aP[i, 0]
            for j in range(1, m):
                denom       = aP[i, j] + aS[i, j] * c_buf[i, j - 1]
                c_buf[i, j] = -aN[i, j] / denom
                d_buf[i, j] = (rhs[i, j] + aS[i, j] * d_buf[i, j - 1]) / denom
            T_out[i, m - 1] = d_buf[i, m - 1]
            for j in range(m - 2, -1, -1):
                T_out[i, j] = d_buf[i, j] - c_buf[i, j] * T_out[i, j + 1]

    _thomas_r = _thomas_r_nb
    _thomas_z = _thomas_z_nb
else:
    _thomas_r = _thomas_r_np
    _thomas_z = _thomas_z_np


def _eval_poly(coeffs_per_mat, mat_id, T, T_clamp_per_mat, out):
    """Evaluate cubic poly value(T) = c0 + c1·T + c2·T² + c3·T³ per cell.

    coeffs_per_mat   : (nmat, 4) array
    mat_id           : (n, m) int array
    T                : (n, m) float array (will be clamped per material)
    T_clamp_per_mat  : (nmat,) array; per-cell T is min(T, T_clamp[mat_id])
    out              : (n, m) preallocated output array
    """
    Tc = np.minimum(T, T_clamp_per_mat[mat_id])
    c0 = coeffs_per_mat[mat_id, 0]
    c1 = coeffs_per_mat[mat_id, 1]
    c2 = coeffs_per_mat[mat_id, 2]
    c3 = coeffs_per_mat[mat_id, 3]
    np.multiply(c3, Tc, out=out); out += c2; out *= Tc; out += c1; out *= Tc; out += c0


# =============================================================================
# MAIN SOLVER FUNCTION
# =============================================================================

def run(ht_params, tca_params, props_csv_path, script_dir):
    """Run the 2-D transient heat transfer simulation."""

    # =========================================================================
    # 1. USER SETTINGS  (from ht_params / tca_params) + validation
    # =========================================================================

    dt      = float(ht_params.get('dt',      1e-4))
    runtime = float(ht_params.get('runtime', 2.0))
    tf      = float(ht_params.get('tf',      10.0))
    n       = int(ht_params.get('n_radial',  100))

    mtype1    = int(ht_params.get('wall_material',   1))
    mtype_ins = int(ht_params.get('insert_material', 2))
    mtype_al  = int(ht_params.get('al_material',     3))

    ins_x_start = float(ht_params.get('insert_x_start',   0.28685))
    ins_length  = float(ht_params.get('insert_length',    0.1016))
    ins_thick   = float(ht_params.get('insert_thickness', 0.03992))

    al_x_start = ins_x_start + ins_length
    al_thick   = float(ht_params.get('al_thickness', 0.020))

    h_nat       = float(ht_params.get('h_nat_outer', 5.0))
    h_nat_inner = float(ht_params.get('h_nat_inner', 10.0))
    T_amb       = float(ht_params.get('T_ambient',   300.0))

    mdot_lbm = float(tca_params.get('turbopump_mdot', 20.0))
    mdot     = mdot_lbm * LBM_PER_S_TO_KG_PER_S

    steel_thick = float(ht_params.get('steel_thickness', 0.03))

    plot_interval_cfg = int(ht_params.get('plot_interval',  0))
    record_video      = bool(ht_params.get('record_video',  True))
    video_filename    = str(ht_params.get('video_filename', 'nozzle_heat_transfer.mp4'))
    video_fps         = int(ht_params.get('video_fps',      15))
    video_quality     = int(ht_params.get('video_quality',  8))
    show_live_plot    = bool(ht_params.get('show_live_plot', True))

    results_filename  = str(ht_params.get('results_filename',  'results.npz'))
    results_snapshots = int(ht_params.get('results_snapshots', 200))

    solver_tol      = float(ht_params.get('solver_tolerance',    1e-6))
    solver_max_iter = int(ht_params.get('solver_max_iterations', 10))

    # ---- Input validation -------------------------------------------------
    _problems = []
    if dt <= 0:           _problems.append(f'dt must be > 0 (got {dt})')
    if runtime <= 0:      _problems.append(f'runtime must be > 0 (got {runtime})')
    if tf < runtime:      _problems.append(f'tf ({tf}) must be >= runtime ({runtime})')
    if n < 5:             _problems.append(f'n_radial must be >= 5 (got {n})')
    if ins_length <= 0:   _problems.append(f'insert_length must be > 0 (got {ins_length})')
    if ins_thick <= 0:    _problems.append(f'insert_thickness must be > 0 (got {ins_thick})')
    if al_thick <= 0:     _problems.append(f'al_thickness must be > 0 (got {al_thick})')
    if steel_thick <= 0:  _problems.append(f'steel_thickness must be > 0 (got {steel_thick})')
    if mdot <= 0:         _problems.append(f'mdot must be > 0 (got {mdot})')
    if _problems:
        raise ValueError('Invalid heat_transfer.yaml settings:\n  - ' + '\n  - '.join(_problems))

    contour_rel  = ht_params.get('contour_csv', 'contour.csv')
    contour_path = os.path.normpath(os.path.join(script_dir, contour_rel))
    video_path   = os.path.join(script_dir, video_filename)
    results_path = os.path.join(script_dir, results_filename)

    # =========================================================================
    # 2. LOAD GEOMETRY
    # =========================================================================

    contour_tbl    = pd.read_csv(contour_path)
    x_raw_unsorted = (contour_tbl['x_mm'].values - contour_tbl['x_mm'].min()) / 1000.0
    I              = np.argsort(x_raw_unsorted)
    x_raw          = x_raw_unsorted[I]
    r_raw          = contour_tbl['y_mm'].values[I] / 1000.0

    idx_unique = uniquetol(x_raw, 1e-12)
    x_contour  = x_raw[idx_unique]
    r_contour  = r_raw[idx_unique]
    m          = len(x_contour)

    if ins_x_start + ins_length > x_contour[-1]:
        raise ValueError(f'insert (x_start={ins_x_start}, length={ins_length}) extends '
                         f'beyond nozzle exit at x={x_contour[-1]:.4f} m')

    print(f'Geometry loaded   : {m} axial stations.')
    print(f'Duplicates removed: {len(x_raw) - m} points.')

    # =========================================================================
    # 3. KEY GEOMETRIC FEATURES
    # =========================================================================

    throat_idx = int(np.argmin(r_contour))
    r_throat   = r_contour[throat_idx]
    Throat_D   = 2.0 * r_throat
    Throat_A   = np.pi * r_throat**2

    slope  = np.diff(r_contour) / (np.diff(x_contour) + 1e-9)
    cc_arr = np.where(slope < -1e-4)[0]
    cc_idx = int(cc_arr[0]) if len(cc_arr) > 0 else 0

    print(f'CC end  : idx={cc_idx}  x={x_contour[cc_idx]:.3f} m')
    print(f'Throat  : idx={throat_idx}  x={x_contour[throat_idx]:.3f} m')

    # =========================================================================
    # 4. LOAD THERMODYNAMIC PROPERTIES
    # =========================================================================

    param = pd.read_csv(props_csv_path)

    pos = np.array([0.0,
                    x_contour[cc_idx],
                    x_contour[throat_idx],
                    x_contour[-1]])

    def iprop(p):
        p  = np.asarray(p, dtype=float)
        f1 = interp1d(pos[0:2], p[0:2], kind='linear', fill_value='extrapolate')
        part1 = f1(x_contour[:cc_idx])
        cs    = CubicSpline(pos[1:], p[1:])
        part2 = cs(x_contour[cc_idx:])
        return np.concatenate([part1, part2])

    gamma_g = iprop(param['Gamma'].values)
    Pr_g    = iprop(param['Prandtl'].values)
    mu_g    = iprop(param['Viscosity_Pa_s'].values)
    cp_g    = iprop(param['Cp_J_kgK'].values)

    # =========================================================================
    # 5. MACH-NUMBER DISTRIBUTION
    # =========================================================================

    M              = np.ones(m)
    M[:throat_idx] = 0.3
    M[throat_idx:] = 2.0

    for it in range(1000):
        term   = 1.0 + (gamma_g - 1.0) / 2.0 * M**2
        AoA    = (((gamma_g + 1.0) / 2.0) ** (-(gamma_g + 1.0) / (2.0*(gamma_g - 1.0)))
                  * term ** ((gamma_g + 1.0) / (2.0*(gamma_g - 1.0))) / M)
        res_M  = np.pi * r_contour**2 / Throat_A - AoA
        dM_eps = 1e-6
        term_p = 1.0 + (gamma_g - 1.0) / 2.0 * (M + dM_eps)**2
        AoA_dM = (((gamma_g + 1.0) / 2.0) ** (-(gamma_g + 1.0) / (2.0*(gamma_g - 1.0)))
                  * term_p ** ((gamma_g + 1.0) / (2.0*(gamma_g - 1.0))) / (M + dM_eps))
        M += 0.5 * res_M / ((AoA_dM - AoA) / dM_eps)
        if np.all(np.abs(res_M) < 1e-6):
            break

    print(f'Mach converged: {it + 1} iters, max res = {np.max(np.abs(res_M)):.2e}')

    # =========================================================================
    # 6. BARTZ CORRELATION
    # =========================================================================

    inj     = param['Station'] == 'Injector'
    Pc      = float(param.loc[inj, 'Pressure_Pa'].values[0])
    Tc      = float(param.loc[inj, 'Temperature_K'].values[0])
    if Pc <= 0:
        raise ValueError(f'Chamber pressure must be > 0 (got {Pc} Pa)')
    cstar   = Pc * Throat_A / mdot
    T_g     = Tc / (1.0 + (gamma_g - 1.0) / 2.0 * M**2)
    hG_base = ((0.026 / Throat_D**0.2)
               * (mu_g**0.2 * cp_g / Pr_g**0.6)
               * (Pc / cstar)**0.8
               * (Throat_A / (np.pi * r_contour**2))**0.9)
    Taw     = T_g * (1.0 + (gamma_g - 1.0) / 2.0 * M**2 * Pr_g**(1.0/3.0))

    # =========================================================================
    # 7. FVM GRID
    # =========================================================================

    al_start_j = int(np.searchsorted(x_contour, al_x_start))
    if al_start_j >= m:
        al_start_j = m - 1
    r_outer_steel = steel_thick
    r_outer_al    = al_thick

    thickness = np.where(x_contour >= al_x_start,
                         r_outer_al   - r_contour,
                         r_outer_steel - r_contour)
    thickness = np.maximum(thickness, 1e-4)

    print(f'Steel outer radius   : {r_outer_steel:.4f} m')
    print(f'Aluminum outer radius: {r_outer_al:.4f} m')
    dr = thickness / n

    r_f = np.zeros((n + 1, m))
    for j in range(m):
        r_f[:, j] = np.linspace(r_contour[j], r_contour[j] + thickness[j], n + 1)

    rw = r_f[:n,    :]
    re = r_f[1:n+1, :]

    dx_vec = np.diff(x_contour)
    dx_vec = np.append(dx_vec, dx_vec[-1])
    Vp = np.pi * (re**2 - rw**2) * dx_vec
    Aw = 2.0 * np.pi * rw * dx_vec
    Ae = 2.0 * np.pi * re * dx_vec
    Ar = np.pi * (re**2 - rw**2)

    dx_S       = np.zeros(m)
    dx_N       = np.zeros(m)
    dx_S[1:m]  = (dx_vec[:m-1] + dx_vec[1:m]) / 2.0
    dx_N[:m-1] = (dx_vec[:m-1] + dx_vec[1:m]) / 2.0
    dx_S[0]    = dx_vec[0]
    dx_N[m-1]  = dx_vec[m-1]

    # =========================================================================
    # 8. MATERIAL PROPERTY POLYNOMIALS  (loaded from materials.yaml)
    # =========================================================================

    mat_yaml_path = os.path.join(script_dir, 'materials.yaml')
    with open(mat_yaml_path) as _f:
        _mat_data = yaml.safe_load(_f)
    mat_names  = []
    mat_T_melt = []
    k_polys, rho_polys, cp_polys, T_maxes = [], [], [], []
    for _m in _mat_data['materials']:
        k_polys.append(_m['k_poly'])
        rho_polys.append(_m['rho_poly'])
        cp_polys.append(_m['cp_poly'])
        T_maxes.append(float(_m['T_max']))
        mat_names.append(_m['name'])
        mat_T_melt.append(float(_m.get('T_melt', float('nan'))))
    nmat = len(mat_names)

    for _idx, _name in [(mtype1, 'wall_material'),
                        (mtype_ins, 'insert_material'),
                        (mtype_al, 'al_material')]:
        if not (0 <= _idx < nmat):
            raise ValueError(f'{_name} index {_idx} is out of range '
                             f'(materials.yaml has {nmat} entries: 0..{nmat - 1})')

    k_coeffs   = np.asarray(k_polys,   dtype=np.float64)   # (nmat, 4)
    rho_coeffs = np.asarray(rho_polys, dtype=np.float64)
    cp_coeffs  = np.asarray(cp_polys,  dtype=np.float64)
    T_clamp    = np.asarray(T_maxes,   dtype=np.float64)   # (nmat,)
    T_melt_arr = np.asarray(mat_T_melt, dtype=np.float64)

    print(f'Materials loaded  : {nmat} entries from materials.yaml')
    print('=' * 70)
    print(f'  Wall   : [{mtype1}]  {mat_names[mtype1]}')
    print(f'  Insert : [{mtype_ins}]  {mat_names[mtype_ins]}')
    print(f'  Exit   : [{mtype_al}]  {mat_names[mtype_al]}')
    print('=' * 70)

    # =========================================================================
    # 9. STATIC PRECOMPUTATIONS  (region masks, mat_id, area arrays)
    # =========================================================================

    x_2d  = np.tile(x_contour, (n, 1))
    rc_2d = np.tile(r_contour, (n, 1))

    ins_mask  = ((x_2d >= ins_x_start) &
                 (x_2d <= ins_x_start + ins_length) &
                 (r_f[:n, :] - r_throat <= ins_thick) &
                 (r_f[:n, :] >= rc_2d))
    al_mask   = ((x_2d >= al_x_start) &
                 (r_f[:n, :] - r_throat <= al_thick) &
                 (r_f[:n, :] >= rc_2d))
    base_mask = ~ins_mask & ~al_mask

    mat_id = np.full((n, m), mtype1, dtype=np.intp)
    mat_id[ins_mask] = mtype_ins
    mat_id[al_mask]  = mtype_al

    dr_2d     = np.tile(dr,   (n, 1))
    dx_S_2d   = np.tile(dx_S, (n, 1))
    dx_N_2d   = np.tile(dx_N, (n, 1))
    Aw_ov_dr  = Aw / dr_2d
    Ae_ov_dr  = Ae / dr_2d
    Ar_ov_dxS = Ar / dx_S_2d
    Ar_ov_dxN = Ar / dx_N_2d

    firing_steps = round(runtime / dt)
    total_steps  = max(1, round(tf / dt))

    if plot_interval_cfg > 0:
        plot_interval = plot_interval_cfg
    else:
        plot_interval = max(1, total_steps // 200)

    snapshot_stride = max(1, total_steps // max(1, results_snapshots))

    print(f'Time steps        : {total_steps} total ({firing_steps} firing)')
    print(f'Plot interval     : every {plot_interval} steps '
          f'(~{max(1, total_steps // plot_interval)} frames)')
    print(f'Numba acceleration: {"ON" if NUMBA_AVAILABLE else "OFF (install numba for ~5x speedup)"}')

    # =========================================================================
    # 10. INITIALISATION  (preallocate everything used in the hot loop)
    # =========================================================================

    T_2d     = T_amb * np.ones((n, m))
    T_2d_hot = None
    sim_time = 0.0
    step_count = 0

    k_field = np.empty((n, m))
    rho_f   = np.empty((n, m))
    cp_f    = np.empty((n, m))

    kW_nb = np.empty((n, m))
    kE_nb = np.empty((n, m))
    kS_nb = np.empty((n, m))
    kN_nb = np.empty((n, m))

    T_W = np.empty((n, m))
    T_E = np.empty((n, m))
    T_S = np.empty((n, m))
    T_N = np.empty((n, m))

    aP    = np.empty((n, m))
    src   = np.empty((n, m))
    Bm    = np.empty((n, m))
    R_mat = np.empty((n, m))

    rhs_r = np.empty((n, m))
    rhs_z = np.empty((n, m))
    T_new = np.empty((n, m))
    c_buf = np.empty((n, m))
    d_buf = np.empty((n, m))

    # Trigger Numba compile now so it doesn't pollute the first iteration timing.
    if NUMBA_AVAILABLE:
        _ones = np.ones((n, m))
        _thomas_r(_ones, _ones, _ones, _ones, T_new, c_buf, d_buf)
        _thomas_z(_ones, _ones, _ones, _ones, T_new, c_buf, d_buf)

    # ---- Time-history storage --------------------------------------------
    n_snaps = total_steps // snapshot_stride + 2
    T_history    = np.empty((n_snaps, n, m), dtype=np.float32)
    t_array      = np.empty(n_snaps,         dtype=np.float64)
    inner_hist   = np.empty((n_snaps, m),    dtype=np.float32)
    outer_hist   = np.empty((n_snaps, m),    dtype=np.float32)
    snap_count   = 0
    melt_events  = []   # list of {'t':…, 'material':…, 'max_T':…}
    _melt_warned = set()

    def _record_snapshot(t):
        nonlocal snap_count
        T_history[snap_count]  = T_2d.astype(np.float32, copy=False)
        inner_hist[snap_count] = T_2d[0,  :].astype(np.float32, copy=False)
        outer_hist[snap_count] = T_2d[-1, :].astype(np.float32, copy=False)
        t_array[snap_count]    = t
        snap_count += 1

    _record_snapshot(0.0)

    # =========================================================================
    # 11. FIGURE & VIDEO SETUP
    # =========================================================================

    draw_frames = show_live_plot or record_video

    Z_grid, R_norm = np.meshgrid(x_contour, np.linspace(0.0, 1.0, n))
    R_grid = np.zeros_like(R_norm)
    for j in range(m):
        R_grid[:, j] = r_contour[j] + R_norm[:, j] * thickness[j]

    fig = None
    h_plot = h_wall_line = h_taw_line = None
    ax1 = ax2 = None
    if draw_frames:
        if show_live_plot:
            plt.ion()
        else:
            matplotlib.use('Agg', force=True)
        fig = plt.figure(figsize=(12, 8), facecolor='w')
        gs  = fig.add_gridspec(2, 1, hspace=0.4)
        ax1 = fig.add_subplot(gs[0])
        ax2 = fig.add_subplot(gs[1])

        h_plot = ax1.pcolormesh(Z_grid, R_grid, T_2d,
                                shading='gouraud', cmap='hot', vmin=300, vmax=3000)
        cb = fig.colorbar(h_plot, ax=ax1)
        cb.set_label('Temperature [K]')
        ax1.set_xlabel('Axial Position [m]')
        ax1.set_ylabel('Radius [m]')
        ax1.set_title('2D Wall Temperature Distribution')
        ax1.set_aspect('equal')
        ax1.autoscale_view()
        ins_rect = patches.Rectangle((ins_x_start, r_throat), ins_length, ins_thick,
                                       edgecolor='w', facecolor='none',
                                       linewidth=0.5, linestyle=':')
        ax1.add_patch(ins_rect)
        ax1.axvline(al_x_start, color='cyan', linewidth=0.6, linestyle=':', zorder=4)

        h_wall_line, = ax2.plot(x_contour, T_2d[0, :], 'r-', linewidth=2,
                                label='Inner Wall')
        ax2.axvline(ins_x_start,              color='b',    linestyle='--', label=f'{mat_names[mtype_ins]} start')
        ax2.axvline(ins_x_start + ins_length, color='b',    linestyle='--', label=f'{mat_names[mtype_ins]} end')
        ax2.axvline(al_x_start,               color='cyan', linestyle='--', label=f'{mat_names[mtype_al]} start')
        h_taw_line, = ax2.plot(x_contour, Taw, 'k:', linewidth=1.0, label='Gas Recovery (Taw)')
        ax2.grid(True)
        ax2.set_xlabel('Axial Position [m]')
        ax2.set_ylabel('Inner Wall Temperature [K]')
        ax2.set_title('Wall Temperature Profile')
        ax2.set_xlim([x_contour[0], x_contour[-1]])
        ax2.set_ylim([300, 3500])
        ax2.legend(loc='upper right', fontsize=8)
        plt.tight_layout()
        if show_live_plot:
            plt.pause(0.05)

    vid_writer = None
    if record_video:
        if IMAGEIO_AVAILABLE and draw_frames:
            vid_writer = imageio.get_writer(video_path, fps=video_fps,
                                            quality=video_quality,
                                            macro_block_size=None)
            print(f'Video: recording to {os.path.basename(video_path)}')
        elif not IMAGEIO_AVAILABLE:
            print('Video: DISABLED (imageio[ffmpeg] not installed).')
    else:
        print('Video: DISABLED (record_video=false in heat_transfer.yaml)')
        if os.path.exists(video_path):
            os.remove(video_path)
            print(f'       Removed stale file: {os.path.basename(video_path)}')

    # =========================================================================
    # 12. MAIN TRANSIENT LOOP
    # =========================================================================

    print(f'\nStarting transient simulation '
          f'(line-by-line TDMA, mdot={mdot:.3f} kg/s)...')
    t_wall_start = time_module.time()

    j_exp = al_start_j - 1
    exp_face_active = (0 <= j_exp < m - 1)

    progress = None
    if TQDM_AVAILABLE:
        progress = tqdm(total=total_steps, unit='step', dynamic_ncols=True)

    try:
        while sim_time < tf:

            # (a) Material properties (fused single eval per property)
            _eval_poly(k_coeffs,   mat_id, T_2d, T_clamp, k_field)
            _eval_poly(rho_coeffs, mat_id, T_2d, T_clamp, rho_f)
            _eval_poly(cp_coeffs,  mat_id, T_2d, T_clamp, cp_f)

            # (b) Neighbor-conductivity arrays — in-place slice copies
            kW_nb[0, :]    = k_field[0, :]
            kW_nb[1:, :]   = k_field[:-1, :]
            kE_nb[:-1, :]  = k_field[1:, :]
            kE_nb[-1, :]   = k_field[-1, :]
            kS_nb[:, 0]    = k_field[:, 0]
            kS_nb[:, 1:]   = k_field[:, :-1]
            kN_nb[:, :-1]  = k_field[:, 1:]
            kN_nb[:, -1]   = k_field[:, -1]

            # Harmonic-mean face conductivities
            denom_w = k_field + kW_nb + 1e-30
            denom_e = k_field + kE_nb + 1e-30
            denom_s = k_field + kS_nb + 1e-30
            denom_n = k_field + kN_nb + 1e-30
            kw_h = (2.0 * k_field * kW_nb) / denom_w
            ke_h = (2.0 * k_field * kE_nb) / denom_e
            ks_h = (2.0 * k_field * kS_nb) / denom_s
            kn_h = (2.0 * k_field * kN_nb) / denom_n

            # Coefficients (Robin BCs handled below)
            aW = kw_h * Aw_ov_dr
            aE = ke_h * Ae_ov_dr
            aS = ks_h * Ar_ov_dxS
            aN = kn_h * Ar_ov_dxN

            # Adiabatic axial ends + zero "neighbor" beyond radial bounds
            aW[0,   :] = 0.0
            aE[-1,  :] = 0.0
            aS[:,   0] = 0.0
            aN[:,  -1] = 0.0

            # Aluminum-step exposed face: cut internal axial link
            if exp_face_active:
                aN[:, j_exp] = 0.0

            # Inner-wall Robin: Bartz during firing, natural conv during cooldown
            if step_count < firing_steps:
                term_Ma  = 1.0 + (gamma_g - 1.0) / 2.0 * M**2
                sigma_v  = 1.0 / ((0.5 * (T_2d[0, :] / T_g) * term_Ma
                                   + 0.5)**0.68 * term_Ma**0.12)
                hg_vec   = hG_base * sigma_v
                Tref_vec = Taw
            else:
                hg_vec   = h_nat_inner * np.ones(m)
                Tref_vec = T_amb       * np.ones(m)

            bin_row  = hg_vec * Aw[0,   :]
            bout_row = h_nat  * Ae[-1,  :]

            src.fill(0.0)
            src[0,  :]  = bin_row * Tref_vec
            src[-1, :] += bout_row * T_amb

            AP0 = rho_f * cp_f * Vp / dt
            aP_bc = np.zeros((n, m))
            aP_bc[0,  :]  = bin_row
            aP_bc[-1, :] += bout_row

            if exp_face_active:
                h_face            = h_nat * Ar[:, j_exp]
                aP_bc[:, j_exp]  += h_face
                src[:,   j_exp]  += h_face * T_amb

            aP[:] = aW + aE + aS + aN + AP0 + aP_bc
            Bm[:] = AP0 * T_2d + src
            B_norm = np.max(np.abs(Bm))

            # (c) Line-by-line TDMA
            converged = False
            for tdma_iter in range(1, solver_max_iter + 1):

                # Radial sweep
                T_S[:, 0]   = T_2d[:, 0]
                T_S[:, 1:]  = T_2d[:, :-1]
                T_N[:, :-1] = T_2d[:, 1:]
                T_N[:, -1]  = T_2d[:, -1]
                np.add(Bm, aS * T_S, out=rhs_r); rhs_r += aN * T_N
                _thomas_r(aW, aP, aE, rhs_r, T_new, c_buf, d_buf)
                T_2d, T_new = T_new, T_2d

                # Axial sweep
                T_W[0, :]   = T_2d[0, :]
                T_W[1:, :]  = T_2d[:-1, :]
                T_E[:-1, :] = T_2d[1:, :]
                T_E[-1, :]  = T_2d[-1, :]
                np.add(Bm, aW * T_W, out=rhs_z); rhs_z += aE * T_E
                _thomas_z(aS, aP, aN, rhs_z, T_new, c_buf, d_buf)
                T_2d, T_new = T_new, T_2d

                # Residual
                T_W[0, :]   = T_2d[0, :]
                T_W[1:, :]  = T_2d[:-1, :]
                T_E[:-1, :] = T_2d[1:, :]
                T_E[-1, :]  = T_2d[-1, :]
                T_S[:, 0]   = T_2d[:, 0]
                T_S[:, 1:]  = T_2d[:, :-1]
                T_N[:, :-1] = T_2d[:, 1:]
                T_N[:, -1]  = T_2d[:, -1]
                R_mat[:] = aP*T_2d - aW*T_W - aE*T_E - aS*T_S - aN*T_N - Bm
                if np.max(np.abs(R_mat)) < solver_tol * B_norm:
                    converged = True
                    break

            if not converged:
                ratio = np.max(np.abs(R_mat)) / max(B_norm, 1e-30)
                print(f'WARNING: TDMA did not converge at t={sim_time:.4f} s  '
                      f'(max|R|/|B| = {ratio:.2e}, iters={tdma_iter})')

            sim_time   += dt
            step_count += 1

            if step_count == firing_steps:
                T_2d_hot = T_2d.copy()

            # ---- Melt-exceedance check (one warning per material) ---------
            for _midx in (mtype1, mtype_ins, mtype_al):
                if _midx in _melt_warned:
                    continue
                _Tm = T_melt_arr[_midx]
                if np.isnan(_Tm):
                    continue
                _Tmax_mat = T_2d[mat_id == _midx].max()
                if _Tmax_mat > _Tm:
                    msg = (f'\nWARNING: {mat_names[_midx]} exceeded melt T={_Tm:.0f} K '
                           f'at t={sim_time:.3f} s (peak {_Tmax_mat:.0f} K)')
                    if progress is not None:
                        progress.write(msg)
                    else:
                        print(msg)
                    melt_events.append({'t': sim_time, 'material': mat_names[_midx],
                                        'max_T': float(_Tmax_mat), 'T_melt': float(_Tm)})
                    _melt_warned.add(_midx)

            # ---- Snapshot recording ---------------------------------------
            if step_count % snapshot_stride == 0 and snap_count < n_snaps:
                _record_snapshot(sim_time)

            # ---- Live plot / video frame ----------------------------------
            if draw_frames and (step_count % plot_interval == 0):
                h_plot.set_array(T_2d)
                state_str = 'FIRING' if step_count < firing_steps else 'COOLING'
                ax1.set_title(f'{state_str} | t = {sim_time:.3f} s')
                h_wall_line.set_ydata(T_2d[0, :])
                if step_count >= firing_steps:
                    h_taw_line.set_ydata(T_amb * np.ones_like(x_contour))
                fig.canvas.draw()
                if show_live_plot:
                    fig.canvas.flush_events()
                if vid_writer is not None:
                    buf = np.asarray(fig.canvas.buffer_rgba())
                    img = buf[..., :3]
                    vid_writer.append_data(img)

            # ---- Progress bar update --------------------------------------
            if progress is not None:
                progress.update(1)
                if step_count % 25 == 0:
                    progress.set_postfix(
                        t=f'{sim_time:6.3f}s',
                        Tmax=f'{T_2d.max():6.0f}K',
                        iters=tdma_iter,
                        snaps=snap_count,
                    )

    except Exception:
        print(f'\nSimulation error at t = {sim_time:.4f} s:')
        traceback.print_exc()
        if vid_writer is not None:
            try:
                vid_writer.close()
            except Exception:
                pass
        if progress is not None:
            progress.close()
        raise

    if progress is not None:
        progress.close()

    elapsed = time_module.time() - t_wall_start
    print(f'\nSimulation complete.  Elapsed: {elapsed:.1f} s '
          f'({total_steps/elapsed:.0f} steps/s)')
    if vid_writer is not None:
        vid_writer.close()

    # ---- Final snapshot ---------------------------------------------------
    if snap_count < n_snaps:
        _record_snapshot(sim_time)

    # =========================================================================
    # 13. SAVE RESULTS
    # =========================================================================

    T_hot_save = T_2d_hot if T_2d_hot is not None else T_2d
    config_summary = {
        'dt': dt, 'runtime': runtime, 'tf': tf, 'n_radial': n,
        'wall_material': mtype1, 'insert_material': mtype_ins, 'al_material': mtype_al,
        'insert_x_start': ins_x_start, 'insert_length': ins_length,
        'insert_thickness': ins_thick, 'al_thickness': al_thick,
        'steel_thickness': steel_thick,
        'h_nat_outer': h_nat, 'h_nat_inner': h_nat_inner, 'T_ambient': T_amb,
        'mdot_kg_s': mdot, 'Pc_Pa': Pc, 'Tc_K': Tc,
    }
    np.savez_compressed(
        results_path,
        T_2d_hot=T_hot_save.astype(np.float32),
        T_2d_final=T_2d.astype(np.float32),
        T_history=T_history[:snap_count],
        t_array=t_array[:snap_count],
        inner_hist=inner_hist[:snap_count],
        outer_hist=outer_hist[:snap_count],
        R_grid=R_grid.astype(np.float32),
        Z_grid=Z_grid.astype(np.float32),
        x_contour=x_contour,
        r_contour=r_contour,
        mat_id=mat_id.astype(np.int8),
        mat_names=np.array(mat_names),
        mat_T_melt=T_melt_arr,
        Taw=Taw,
        runtime=runtime,
        melt_events_t=np.array([e['t'] for e in melt_events]),
        melt_events_mat=np.array([e['material'] for e in melt_events]),
        melt_events_Tmax=np.array([e['max_T'] for e in melt_events]),
        config=np.array(yaml.safe_dump(config_summary)),
    )
    print(f'Results: saved {snap_count} snapshots to {os.path.basename(results_path)} '
          f'({os.path.getsize(results_path)/1e6:.1f} MB)')
    print(f'         View with:  python view_results.py "{results_path}"')

    # =========================================================================
    # 14. FINAL PLOTS
    # =========================================================================

    T_plot = T_hot_save
    t_plot = runtime if T_2d_hot is not None else tf

    plt.ioff()
    fig2, ax_f = plt.subplots(figsize=(12, 5), facecolor='w')
    pcm = ax_f.pcolormesh(Z_grid, R_grid, T_plot, shading='gouraud', cmap='hot')
    fig2.colorbar(pcm, ax=ax_f)
    ins_rect2 = patches.Rectangle((ins_x_start, r_throat), ins_length, ins_thick,
                                    edgecolor='g', facecolor='none',
                                    linewidth=1.5, linestyle='--')
    ax_f.add_patch(ins_rect2)
    ax_f.set_aspect('equal')
    ax_f.autoscale_view()
    ax_f.set_title(f'Peak Temperature Field  (t = {t_plot:.2f} s — end of firing)')
    ax_f.set_xlabel('Axial Position [m]')
    ax_f.set_ylabel('Radius [m]')
    plt.tight_layout()

    fig3, ax3 = plt.subplots(figsize=(10, 5), facecolor='w')
    ax3.plot(x_contour, T_plot[0,   :], 'r-',  linewidth=2, label='Inner wall (gas side)')
    ax3.plot(x_contour, T_plot[-1,  :], 'b--', linewidth=2, label='Outer wall (ambient side)')
    ax3.axvline(ins_x_start,              color='g',    linestyle='--', linewidth=1.2, label=f'{mat_names[mtype_ins]} start')
    ax3.axvline(ins_x_start + ins_length, color='g',    linestyle='--', linewidth=1.2, label=f'{mat_names[mtype_ins]} end')
    ax3.axvline(al_x_start,               color='cyan', linestyle='--', linewidth=1.2, label=f'{mat_names[mtype_al]} start')

    _x0 = x_contour[0]
    _x1 = x_contour[-1]
    _melt_styles = [
        (mtype1,    'orange', mat_names[mtype1],    (_x0,          ins_x_start)),
        (mtype_ins, 'purple', mat_names[mtype_ins], (ins_x_start,  ins_x_start + ins_length)),
        (mtype_al,  'teal',   mat_names[mtype_al],  (al_x_start,   _x1)),
    ]
    _seen = set()
    for _idx, _color, _name, (_xa, _xb) in _melt_styles:
        _Tm = T_melt_arr[_idx]
        if np.isnan(_Tm) or _idx in _seen or _xa >= _xb:
            continue
        _seen.add(_idx)
        ax3.plot([_xa, _xb], [_Tm, _Tm], color=_color, linestyle=':', linewidth=1.2,
                 label=f'{_name} melt  {_Tm:.0f} K')

    ax3.set_xlabel('Axial Position [m]', fontsize=12)
    ax3.set_ylabel('Temperature [K]',    fontsize=12)
    ax3.set_title(f'Peak Wall Temperature Profile  (t = {t_plot:.2f} s — end of firing)', fontsize=13)
    ax3.legend(loc='best', fontsize=8)
    ax3.grid(True)
    plt.tight_layout()
    plt.show()


# =============================================================================
# Standalone entry point
# =============================================================================

if __name__ == '__main__':
    _script_dir = os.path.dirname(os.path.abspath(__file__))

    _ht_yaml = os.path.join(_script_dir, 'heat_transfer.yaml')
    with open(_ht_yaml) as _f:
        _ht_params = yaml.safe_load(_f) or {}

    _tca_rel = _ht_params.get('tca_params_path', '../../TCA_params.yaml')
    _tca_yaml = _tca_rel if os.path.isabs(_tca_rel) else os.path.normpath(os.path.join(_script_dir, _tca_rel))
    if not os.path.exists(_tca_yaml):
        _alt = os.path.join(os.getcwd(), 'TCA_params.yaml')
        if os.path.exists(_alt):
            _tca_yaml = _alt
    if not os.path.exists(_tca_yaml):
        raise FileNotFoundError(
            f'TCA_params.yaml not found.\n'
            f'  Tried: {_tca_yaml}\n'
            f'  Set tca_params_path in heat_transfer.yaml to point to it.'
        )
    with open(_tca_yaml) as _f:
        _tca_params = yaml.safe_load(_f) or {}

    _props_csv = os.path.join(_script_dir, 'turbopump_properties.csv')
    if not os.path.exists(_props_csv):
        raise FileNotFoundError(
            f'turbopump_properties.csv not found at {_props_csv}\n'
            'Run Bartz_Values.py first, or use run.py to run the full pipeline.'
        )

    run(_ht_params, _tca_params, _props_csv, _script_dir)
