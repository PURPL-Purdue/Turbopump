#!/usr/bin/env python3
"""
========================================================================
2-D TRANSIENT CYLINDRICAL HEAT TRANSFER SOLVER — ROCKET NOZZLE WALL
------------------------------------------------------------------------
Solves the 2-D unsteady conduction equation in cylindrical coordinates
over the rocket nozzle/combustion-chamber wall cross-section.

    rho*cp * dT/dt = 1/r * d/dr(r*k*dT/dr) + d/dz(k*dT/dz)

Discretization : FVM, implicit backward-Euler
Linear solver  : Red-Black SOR

Can be run standalone or called via run() from run.py.

Usage (standalone):
    python HeatTranferSolver_2D.py
    # reads heat_transfer.yaml and TCA_params.yaml from known locations

Usage (from run.py):
    import HeatTranferSolver_2D
    HeatTranferSolver_2D.run(ht_params, tca_params, props_csv_path, script_dir)
========================================================================
"""

import os
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


def poly4(row, T):
    """Cubic polynomial: row[0] + row[1]*T + row[2]*T^2 + row[3]*T^3"""
    return row[0] + row[1]*T + row[2]*T**2 + row[3]*T**3


def _thomas_r(aW, aP, aE, rhs):
    """Vectorized Thomas (TDMA) sweep in the radial direction.

    Solves m independent tridiagonal systems of size n simultaneously
    (one per axial column j).

    System per column j, row i:
        aP[i]*T[i] - aW[i]*T[i-1] - aE[i]*T[i+1] = rhs[i]

    Boundary conditions are already encoded:
        aW[0,:] = 0  (inner wall),  aE[-1,:] = 0  (outer wall)
    """
    n = aP.shape[0]
    c = np.empty_like(aP)    # modified upper / pivot ratio
    d = np.empty_like(rhs)   # modified RHS

    c[0] = -aE[0] / aP[0]
    d[0] =  rhs[0] / aP[0]
    for i in range(1, n):
        denom = aP[i] + aW[i] * c[i - 1]   # aP[i] - (-aW[i])*c'[i-1]
        c[i]  = -aE[i] / denom
        d[i]  = (rhs[i] + aW[i] * d[i - 1]) / denom

    T = np.empty_like(aP)
    T[-1] = d[-1]
    for i in range(n - 2, -1, -1):
        T[i] = d[i] - c[i] * T[i + 1]
    return T


def _thomas_z(aS, aP, aN, rhs):
    """Vectorized Thomas (TDMA) sweep in the axial direction.

    Solves n independent tridiagonal systems of size m simultaneously
    (one per radial row i).

    System per row i, column j:
        aP[j]*T[j] - aS[j]*T[j-1] - aN[j]*T[j+1] = rhs[j]

    Boundary conditions are already encoded:
        aS[:,0] = 0  (axial left),  aN[:,-1] = 0  (axial right)
    """
    m = aP.shape[1]
    c = np.empty_like(aP)
    d = np.empty_like(rhs)

    c[:, 0] = -aN[:, 0] / aP[:, 0]
    d[:, 0] =  rhs[:, 0] / aP[:, 0]
    for j in range(1, m):
        denom   = aP[:, j] + aS[:, j] * c[:, j - 1]
        c[:, j] = -aN[:, j] / denom
        d[:, j] = (rhs[:, j] + aS[:, j] * d[:, j - 1]) / denom

    T = np.empty_like(aP)
    T[:, -1] = d[:, -1]
    for j in range(m - 2, -1, -1):
        T[:, j] = d[:, j] - c[:, j] * T[:, j + 1]
    return T


# =============================================================================
# MAIN SOLVER FUNCTION
# =============================================================================

def run(ht_params, tca_params, props_csv_path, script_dir):
    """
    Run the 2-D transient heat transfer simulation.

    Parameters
    ----------
    ht_params : dict
        Contents of heat_transfer.yaml.
    tca_params : dict
        Contents of TCA_params.yaml.  Uses:
            turbopump_mdot  [lbm/s]  -> converted to kg/s internally
    props_csv_path : str
        Absolute path to turbopump_properties.csv (written by Bartz_Values).
    script_dir : str
        Directory used for resolving relative paths (contour.csv, video output).
    """

    # =========================================================================
    # 1. USER SETTINGS  (from ht_params / tca_params)
    # =========================================================================

    # Time integration
    dt      = float(ht_params.get('dt',      1e-4))
    runtime = float(ht_params.get('runtime', 2.0))
    tf      = float(ht_params.get('tf',      10.0))

    # Spatial
    n = int(ht_params.get('n_radial', 100))

    # Materials (0-based index)
    mtype1    = int(ht_params.get('wall_material',   1))
    mtype_ins = int(ht_params.get('insert_material', 2))
    mtype_al  = int(ht_params.get('al_material',     3))

    # Insert geometry (graphite)
    ins_x_start = float(ht_params.get('insert_x_start',   0.28685))
    ins_length  = float(ht_params.get('insert_length',    0.1016))
    ins_thick   = float(ht_params.get('insert_thickness', 0.03992))

    # Aluminum exit section geometry (axially connected to graphite end)
    al_x_start = ins_x_start + ins_length
    al_thick   = float(ht_params.get('al_thickness', 0.020))

    # Boundary conditions
    h_nat       = float(ht_params.get('h_nat_outer', 5.0))
    h_nat_inner = float(ht_params.get('h_nat_inner', 10.0))
    T_amb       = float(ht_params.get('T_ambient',   300.0))

    # Mass flow: read from TCA_params in lbm/s, convert to kg/s
    mdot_lbm = float(tca_params.get('turbopump_mdot', 20.0))
    mdot     = mdot_lbm * LBM_PER_S_TO_KG_PER_S   # [kg/s]

    # Wall section thicknesses
    steel_thick = float(ht_params.get('steel_thickness', 0.03))

    # Visualization
    plot_interval  = int(ht_params.get('plot_interval',  20))
    record_video   = bool(ht_params.get('record_video',  True))
    video_filename = str(ht_params.get('video_filename', 'nozzle_heat_transfer.mp4'))
    video_fps      = int(ht_params.get('video_fps',      15))
    video_quality  = int(ht_params.get('video_quality',  8))

    # Solver
    solver_tol      = float(ht_params.get('solver_tolerance',    1e-6))
    solver_max_iter = int(ht_params.get('solver_max_iterations', 10))

    # Resolve file paths relative to script_dir
    contour_rel  = ht_params.get('contour_csv', '../contour.csv')
    contour_path = os.path.normpath(os.path.join(script_dir, contour_rel))
    video_path   = os.path.join(script_dir, video_filename)

    # =========================================================================
    # 2. LOAD GEOMETRY
    # =========================================================================

    contour_tbl    = pd.read_csv(contour_path)
    x_raw_unsorted = (contour_tbl['x'].values - contour_tbl['x'].min()) / 1000.0
    I              = np.argsort(x_raw_unsorted)
    x_raw          = x_raw_unsorted[I]
    r_raw          = contour_tbl['y'].values[I] / 1000.0

    idx_unique = uniquetol(x_raw, 1e-12)
    x_contour  = x_raw[idx_unique]
    r_contour  = r_raw[idx_unique]
    m          = len(x_contour)

    print(f'Geometry loaded   : {m} axial stations.')
    print(f'Duplicates removed: {len(x_raw) - m} points.')

    # =========================================================================
    # 3. KEY GEOMETRIC FEATURES
    # =========================================================================

    throat_idx = int(np.argmin(r_contour))
    r_throat   = r_contour[throat_idx]
    Throat_D   = 2.0 * r_throat
    Throat_A   = np.pi * r_throat**2

    slope    = np.diff(r_contour) / (np.diff(x_contour) + 1e-9)
    cc_arr   = np.where(slope < -1e-4)[0]
    cc_idx   = int(cc_arr[0]) if len(cc_arr) > 0 else 0

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
        """Interpolate 4-point property array over all m axial stations."""
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

    # Each section has a constant (cylindrical) outer radius — stepped at the
    # steel/aluminum boundary.
    #
    #   Steel body   : r_outer_steel = r_CC + steel_thick  (constant along entire steel region)
    #   Aluminum exit: r_outer_al    = r(al_x_start) + al_thick (constant along aluminum region)
    #
    # wall thickness at station j = outer_radius(section) - r_contour[j]
    al_start_j    = int(np.searchsorted(x_contour, al_x_start))
    if al_start_j >= m:
        al_start_j = m - 1
    r_outer_steel = r_contour[cc_idx] + steel_thick
    r_outer_al    = r_contour[al_start_j] + al_thick

    thickness = np.where(x_contour >= al_x_start,
                         r_outer_al   - r_contour,
                         r_outer_steel - r_contour)
    thickness = np.maximum(thickness, 1e-4)   # guard against zero/negative cells

    print(f'Steel outer radius : {r_outer_steel:.4f} m')
    print(f'Aluminum outer radius: {r_outer_al:.4f} m')
    dr        = thickness / n

    r_f = np.zeros((n + 1, m))
    for j in range(m):
        r_f[:, j] = np.linspace(r_contour[j], r_contour[j] + thickness[j], n + 1)

    rw = r_f[:n,    :]
    re = r_f[1:n+1, :]

    dx_vec        = np.diff(x_contour)
    dx_vec        = np.append(dx_vec, dx_vec[-1])
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
    # Each material row: [k_poly, rho_poly, cp_poly, [T_max, 0, 0, 0]]
    # Polynomial:  value(T) = c0 + c1*T + c2*T^2 + c3*T^3

    mat_yaml_path = os.path.join(script_dir, 'materials.yaml')
    with open(mat_yaml_path) as _f:
        _mat_data = yaml.safe_load(_f)
    materials    = []
    mat_names    = []
    mat_T_melt   = []
    for _m in _mat_data['materials']:
        materials.append(np.array([
            _m['k_poly'],
            _m['rho_poly'],
            _m['cp_poly'],
            [_m['T_max'], 0.0, 0.0, 0.0],
        ]))
        mat_names.append(_m['name'])
        mat_T_melt.append(float(_m.get('T_melt', float('nan'))))
    print(f'Materials loaded  : {len(materials)} entries from materials.yaml')

    Mat1   = materials[mtype1]
    MatIns = materials[mtype_ins]
    MatAl  = materials[mtype_al]

    # =========================================================================
    # 9. STATIC PRECOMPUTATIONS
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

    dr_2d     = np.tile(dr,   (n, 1))
    dx_S_2d   = np.tile(dx_S, (n, 1))
    dx_N_2d   = np.tile(dx_N, (n, 1))
    Aw_ov_dr  = Aw / dr_2d
    Ae_ov_dr  = Ae / dr_2d
    Ar_ov_dxS = Ar / dx_S_2d
    Ar_ov_dxN = Ar / dx_N_2d

    firing_steps = round(runtime / dt)

    print('Static precomputation done.')

    # =========================================================================
    # 10. INITIALISATION
    # =========================================================================

    T_2d       = T_amb * np.ones((n, m))
    T_2d_hot   = None   # snapshot at end of firing (t = runtime), used in final plots
    sim_time   = 0.0
    step_count = 0

    # =========================================================================
    # 11. FIGURE & VIDEO SETUP
    # =========================================================================

    plt.ion()
    fig = plt.figure(figsize=(12, 8), facecolor='w')
    gs  = fig.add_gridspec(2, 1, hspace=0.4)
    ax1 = fig.add_subplot(gs[0])
    ax2 = fig.add_subplot(gs[1])

    Z_grid, R_norm = np.meshgrid(x_contour, np.linspace(0.0, 1.0, n))
    R_grid = np.zeros_like(R_norm)
    for j in range(m):
        R_grid[:, j] = r_contour[j] + R_norm[:, j] * thickness[j]

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
    ax2.axvline(ins_x_start,              color='b',    linestyle='--', label='Graphite start')
    ax2.axvline(ins_x_start + ins_length, color='b',    linestyle='--', label='Graphite end')
    ax2.axvline(al_x_start,               color='cyan', linestyle='--', label='Al start')
    h_taw_line, = ax2.plot(x_contour, Taw, 'k:', linewidth=1.0,
                           label='Gas Recovery (Taw)')
    ax2.grid(True)
    ax2.set_xlabel('Axial Position [m]')
    ax2.set_ylabel('Inner Wall Temperature [K]')
    ax2.set_title('Wall Temperature Profile')
    ax2.set_xlim([x_contour[0], x_contour[-1]])
    ax2.set_ylim([300, 3500])
    ax2.legend(loc='upper right')
    plt.tight_layout()
    plt.pause(0.05)

    vid_writer = None
    if record_video:
        if IMAGEIO_AVAILABLE:
            vid_writer = imageio.get_writer(video_path, fps=video_fps,
                                            quality=video_quality,
                                            macro_block_size=None)
            print(f'Video: recording to {os.path.basename(video_path)}')
        else:
            print('Video: DISABLED (imageio[ffmpeg] not installed).')
    else:
        print('Video: DISABLED (record_video=false in heat_transfer.yaml)')
        if os.path.exists(video_path):
            os.remove(video_path)
            print(f'       Removed stale file: {os.path.basename(video_path)}')

    # =========================================================================
    # 12. MAIN TRANSIENT LOOP
    # =========================================================================

    print(f'\nStarting transient simulation (line-by-line TDMA, '
          f'mdot={mdot:.3f} kg/s)...')
    t_wall_start = time_module.time()

    try:
        while sim_time < tf:

            # (a) Conductivity field
            T_c_base = np.minimum(T_2d, Mat1[3, 0])
            T_c_ins  = np.minimum(T_2d, MatIns[3, 0])
            T_c_al   = np.minimum(T_2d, MatAl[3, 0])
            k_field  = (base_mask * poly4(Mat1[0,   :], T_c_base) +
                        ins_mask  * poly4(MatIns[0, :], T_c_ins)  +
                        al_mask   * poly4(MatAl[0,  :], T_c_al))

            # (b) Coefficient assembly
            rho_f = (base_mask * poly4(Mat1[1,   :], T_c_base) +
                     ins_mask  * poly4(MatIns[1, :], T_c_ins)  +
                     al_mask   * poly4(MatAl[1,  :], T_c_al))
            cp_f  = (base_mask * poly4(Mat1[2,   :], T_c_base) +
                     ins_mask  * poly4(MatIns[2, :], T_c_ins)  +
                     al_mask   * poly4(MatAl[2,  :], T_c_al))

            kW_nb = np.vstack([k_field[0:1,   :], k_field[:n-1, :]])
            kE_nb = np.vstack([k_field[1:n,   :], k_field[n-1:n, :]])
            kS_nb = np.hstack([k_field[:, 0:1],   k_field[:, :m-1]])
            kN_nb = np.hstack([k_field[:, 1:m],   k_field[:, m-1:m]])
            kw_h  = 2.0 * k_field * kW_nb / (k_field + kW_nb + 1e-30)
            ke_h  = 2.0 * k_field * kE_nb / (k_field + kE_nb + 1e-30)
            ks_h  = 2.0 * k_field * kS_nb / (k_field + kS_nb + 1e-30)
            kn_h  = 2.0 * k_field * kN_nb / (k_field + kN_nb + 1e-30)

            aW = kw_h * Aw_ov_dr
            aE = ke_h * Ae_ov_dr
            aS = ks_h * Ar_ov_dxS
            aN = kn_h * Ar_ov_dxN

            aW[0,   :] = 0.0
            aE[n-1, :] = 0.0
            aS[:,   0] = 0.0
            aN[:, m-1] = 0.0

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
            bout_row = h_nat  * Ae[n-1, :]

            src         = np.zeros((n, m))
            src[0,  :]  = bin_row * Tref_vec
            src[n-1,:] += bout_row * T_amb

            AP0          = rho_f * cp_f * Vp / dt
            aP_bc        = np.zeros((n, m))
            aP_bc[0,  :] = bin_row
            aP_bc[n-1,:] += bout_row
            aP = aW + aE + aS + aN + AP0 + aP_bc

            Bm     = AP0 * T_2d + src
            B_norm = np.max(np.abs(Bm))

            # (c) Line-by-line TDMA (Thomas algorithm)
            #
            # Each outer iteration:
            #   Radial sweep  — for every axial column j, solve the n×1
            #                   tridiagonal system in i exactly (Thomas),
            #                   treating the axial neighbours as explicit.
            #   Axial sweep   — for every radial row i, solve the m×1
            #                   tridiagonal system in j exactly (Thomas),
            #                   treating the radial neighbours as explicit
            #                   (already updated by the radial sweep).
            #
            # Typically converges in 3–5 iterations vs. up to 500 SOR sweeps.
            converged = False
            for tdma_iter in range(1, solver_max_iter + 1):

                # -- Radial sweep (i-direction, explicit axial neighbours) --
                T_S = np.hstack([T_2d[:, 0:1],  T_2d[:, :m-1]])
                T_N = np.hstack([T_2d[:, 1:m],  T_2d[:, m-1:m]])
                rhs_r = Bm + aS * T_S + aN * T_N
                T_2d  = _thomas_r(aW, aP, aE, rhs_r)

                # -- Axial sweep (j-direction, explicit radial neighbours) --
                T_W = np.vstack([T_2d[0:1, :],  T_2d[:n-1, :]])
                T_E = np.vstack([T_2d[1:n, :],  T_2d[n-1:n, :]])
                rhs_z = Bm + aW * T_W + aE * T_E
                T_2d  = _thomas_z(aS, aP, aN, rhs_z)

                # Residual check (same criterion as before)
                T_W_r = np.vstack([T_2d[0:1, :],  T_2d[:n-1, :]])
                T_E_r = np.vstack([T_2d[1:n, :],  T_2d[n-1:n, :]])
                T_S_r = np.hstack([T_2d[:, 0:1],  T_2d[:, :m-1]])
                T_N_r = np.hstack([T_2d[:, 1:m],  T_2d[:, m-1:m]])
                R_mat = (aP*T_2d - aW*T_W_r - aE*T_E_r
                         - aS*T_S_r - aN*T_N_r - Bm)
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

            if step_count % plot_interval == 0:
                h_plot.set_array(T_2d)
                state_str = 'FIRING' if step_count < firing_steps else 'COOLING'
                ax1.set_title(f'{state_str} | t = {sim_time:.3f} s')
                h_wall_line.set_ydata(T_2d[0, :])
                if step_count >= firing_steps:
                    h_taw_line.set_ydata(T_amb * np.ones_like(x_contour))
                fig.canvas.draw()
                fig.canvas.flush_events()
                if vid_writer is not None:
                    buf        = np.frombuffer(fig.canvas.buffer_rgba(), dtype=np.uint8)
                    w_px, h_px = fig.canvas.get_width_height()
                    actual     = len(buf) // 4
                    if actual != w_px * h_px:
                        scale = int(round((actual / (w_px * h_px)) ** 0.5))
                        w_px *= scale
                        h_px *= scale
                    img = buf.reshape(h_px, w_px, 4)[:, :, :3]
                    vid_writer.append_data(img)

    except Exception:
        print(f'\nSimulation error at t = {sim_time:.4f} s:')
        traceback.print_exc()
        if vid_writer is not None:
            try:
                vid_writer.close()
            except Exception:
                pass
        raise

    elapsed = time_module.time() - t_wall_start
    print(f'\nSimulation complete.  Elapsed: {elapsed:.1f} s')
    if vid_writer is not None:
        vid_writer.close()

    # =========================================================================
    # 13. FINAL TEMPERATURE FIELD
    # =========================================================================

    T_plot = T_2d_hot if T_2d_hot is not None else T_2d
    t_plot = runtime  if T_2d_hot is not None else tf

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

    # =========================================================================
    # 14. INNER / OUTER WALL TEMPERATURE PROFILE
    # =========================================================================

    fig3, ax3 = plt.subplots(figsize=(10, 5), facecolor='w')
    ax3.plot(x_contour, T_plot[0,   :], 'r-',  linewidth=2, label='Inner wall (gas side)')
    ax3.plot(x_contour, T_plot[n-1, :], 'b--', linewidth=2, label='Outer wall (ambient side)')
    ax3.axvline(ins_x_start,              color='g',    linestyle='--', linewidth=1.2, label='Graphite start')
    ax3.axvline(ins_x_start + ins_length, color='g',    linestyle='--', linewidth=1.2, label='Graphite end')
    ax3.axvline(al_x_start,               color='cyan', linestyle='--', linewidth=1.2, label='Al start')

    # Melting / sublimation temperature reference lines — each drawn only over
    # the x-range where that material is actually present.
    #   Steel   : chamber inlet  → graphite insert start
    #   Graphite: insert start   → insert end  (= al_x_start)
    #   Aluminum: al_x_start     → nozzle exit
    _x0 = x_contour[0]
    _x1 = x_contour[-1]
    _melt_styles = [
        (mtype1,    'orange', mat_names[mtype1],    (_x0,          ins_x_start)),
        (mtype_ins, 'purple', mat_names[mtype_ins], (ins_x_start,  ins_x_start + ins_length)),
        (mtype_al,  'teal',  mat_names[mtype_al],  (al_x_start,   _x1)),
    ]
    _seen = set()
    for _idx, _color, _name, (_xa, _xb) in _melt_styles:
        _Tm = mat_T_melt[_idx]
        if np.isnan(_Tm) or _idx in _seen or _xa >= _xb:
            continue
        _seen.add(_idx)
        ax3.plot([_xa, _xb], [_Tm, _Tm], color=_color, linestyle=':', linewidth=1.2,
                 label=f'{_name} melt  {_Tm:.0f} K')

    ax3.set_xlabel('Axial Position [m]', fontsize=12)
    ax3.set_ylabel('Temperature [K]',    fontsize=12)
    ax3.set_title(f'Peak Wall Temperature Profile  (t = {t_plot:.2f} s — end of firing)', fontsize=13)
    ax3.legend(loc='best')
    ax3.grid(True)
    plt.tight_layout()
    plt.show()


# =============================================================================
# Standalone entry point
# =============================================================================

if __name__ == '__main__':
    import yaml

    _script_dir = os.path.dirname(os.path.abspath(__file__))

    # heat_transfer.yaml lives next to this script
    _ht_yaml = os.path.join(_script_dir, 'heat_transfer.yaml')
    with open(_ht_yaml) as _f:
        _ht_params = yaml.safe_load(_f) or {}

    # TCA_params.yaml is two levels up: .../TCA/TCA_params.yaml
    _tca_yaml = os.path.normpath(os.path.join(_script_dir, '..', '..', 'TCA_params.yaml'))
    if not os.path.exists(_tca_yaml):
        _tca_yaml = os.path.join(os.getcwd(), 'TCA_params.yaml')
    with open(_tca_yaml) as _f:
        _tca_params = yaml.safe_load(_f) or {}

    # turbopump_properties.csv is expected in the same directory
    _props_csv = os.path.join(_script_dir, 'turbopump_properties.csv')
    if not os.path.exists(_props_csv):
        raise FileNotFoundError(
            f'turbopump_properties.csv not found at {_props_csv}\n'
            'Run Bartz_Values.py first, or use run.py to run the full pipeline.'
        )

    run(_ht_params, _tca_params, _props_csv, _script_dir)
