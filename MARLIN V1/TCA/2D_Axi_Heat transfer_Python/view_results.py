#!/usr/bin/env python3
"""
view_results.py  —  Interactive viewer for HeatTransferSolver_2D results
-------------------------------------------------------------------------
Loads a results.npz file written by HeatTransferSolver_2D and shows:

    Top panel    : r-z temperature heatmap (animated via time slider)
    Middle panel : inner/outer wall temperature vs axial position
    Bottom panel : peak inner-wall temperature vs time + firing-end marker
                   + melt-exceedance event markers

Usage:
    python view_results.py [results.npz]

If no path is given, defaults to ./results.npz next to this script.
"""

import os
import sys
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.widgets import Slider, Button


def _default_results_path():
    here = os.path.dirname(os.path.abspath(__file__))
    return os.path.join(here, 'results.npz')


def main(results_path):
    if not os.path.exists(results_path):
        print(f'ERROR: file not found: {results_path}')
        print('       Run `python run.py` first to generate results.npz.')
        sys.exit(1)

    data = np.load(results_path, allow_pickle=False)
    T_history  = data['T_history']         # (n_snaps, n, m)
    t_array    = data['t_array']           # (n_snaps,)
    inner_hist = data['inner_hist']        # (n_snaps, m)
    outer_hist = data['outer_hist']        # (n_snaps, m)
    R_grid     = data['R_grid']
    Z_grid     = data['Z_grid']
    x_contour  = data['x_contour']
    r_contour  = data['r_contour']
    mat_id     = data['mat_id']
    mat_names  = [str(s) for s in data['mat_names']]
    T_melt     = data['mat_T_melt']
    Taw        = data['Taw']
    runtime    = float(data['runtime'])
    melt_ts    = data['melt_events_t']     if 'melt_events_t' in data.files else np.array([])
    melt_mats  = data['melt_events_mat']   if 'melt_events_mat' in data.files else np.array([])
    melt_Tmax  = data['melt_events_Tmax']  if 'melt_events_Tmax' in data.files else np.array([])

    n_snaps = T_history.shape[0]
    print(f'Loaded {n_snaps} snapshots from {os.path.basename(results_path)}')
    print(f'  Time range : {t_array[0]:.3f} s  ->  {t_array[-1]:.3f} s')
    print(f'  Grid       : {T_history.shape[1]} radial x {T_history.shape[2]} axial')
    print(f'  Materials  : {", ".join(mat_names)}')
    print(f'  Firing ends: t = {runtime:.3f} s')
    if len(melt_ts) > 0:
        print(f'  Melt events: {len(melt_ts)}')
        for _t, _mat, _Tm in zip(melt_ts, melt_mats, melt_Tmax):
            print(f'     - {_mat} at t={float(_t):.3f}s  peak={float(_Tm):.0f}K')

    # ------------------------------------------------------------------
    # Figure layout
    # ------------------------------------------------------------------
    fig = plt.figure(figsize=(13, 9), facecolor='w')
    gs  = fig.add_gridspec(3, 1, hspace=0.45, top=0.95, bottom=0.12)
    ax_heat    = fig.add_subplot(gs[0])
    ax_profile = fig.add_subplot(gs[1])
    ax_peak    = fig.add_subplot(gs[2])

    # Heatmap
    T_init = T_history[0]
    vmin   = 300.0
    vmax   = max(float(T_history.max()), vmin + 1.0)
    pcm = ax_heat.pcolormesh(Z_grid, R_grid, T_init,
                             shading='gouraud', cmap='hot', vmin=vmin, vmax=vmax)
    cb = fig.colorbar(pcm, ax=ax_heat, pad=0.01)
    cb.set_label('Temperature [K]')
    ax_heat.set_aspect('equal')
    ax_heat.set_xlabel('Axial Position [m]')
    ax_heat.set_ylabel('Radius [m]')
    ax_heat.autoscale_view()
    title_heat = ax_heat.set_title('')

    # Profile
    line_inner, = ax_profile.plot(x_contour, inner_hist[0], 'r-', lw=2, label='Inner wall')
    line_outer, = ax_profile.plot(x_contour, outer_hist[0], 'b--', lw=2, label='Outer wall')
    ax_profile.plot(x_contour, Taw, 'k:', lw=1.0, label='Gas T_aw')
    ax_profile.set_xlabel('Axial Position [m]')
    ax_profile.set_ylabel('Temperature [K]')
    ax_profile.set_xlim(x_contour[0], x_contour[-1])
    y_max = max(float(inner_hist.max()), float(outer_hist.max()), float(Taw.max())) * 1.05
    ax_profile.set_ylim(290, max(3500.0, y_max))
    ax_profile.grid(True)
    ax_profile.legend(loc='upper right', fontsize=8)
    title_profile = ax_profile.set_title('')

    # Material melting reference lines
    seen = set()
    color_cycle = ['orange', 'purple', 'teal', 'olive', 'brown', 'magenta', 'gray']
    for k, _name in enumerate(mat_names):
        if k in seen:
            continue
        seen.add(k)
        if np.isnan(T_melt[k]):
            continue
        if not np.any(mat_id == k):
            continue
        ax_profile.axhline(float(T_melt[k]),
                           color=color_cycle[k % len(color_cycle)],
                           linestyle=':', lw=1.0,
                           label=f'{_name} melt {float(T_melt[k]):.0f}K')

    # Peak history (inner-wall maximum over axial stations vs time)
    peak_inner_t = inner_hist.max(axis=1)
    ax_peak.plot(t_array, peak_inner_t, 'r-', lw=1.5, label='max inner-wall T')
    ax_peak.axvline(runtime, color='k', linestyle='--', lw=1, label='firing end')
    for _t, _Tm in zip(melt_ts, melt_Tmax):
        ax_peak.plot(float(_t), float(_Tm), 'rx', markersize=10, markeredgewidth=2)
    if len(melt_ts) > 0:
        ax_peak.plot([], [], 'rx', markersize=10, markeredgewidth=2, label='melt event')
    ax_peak.set_xlabel('Time [s]')
    ax_peak.set_ylabel('Peak Temperature [K]')
    ax_peak.grid(True)
    ax_peak.legend(loc='best', fontsize=8)
    ax_peak.set_title('Peak inner-wall temperature history')
    cursor = ax_peak.axvline(t_array[0], color='b', lw=1.0, alpha=0.7)

    # ------------------------------------------------------------------
    # Slider + buttons
    # ------------------------------------------------------------------
    ax_slider = fig.add_axes([0.13, 0.04, 0.62, 0.03])
    slider = Slider(ax_slider, 'Time', 0, n_snaps - 1, valinit=0, valstep=1)

    ax_play  = fig.add_axes([0.78, 0.035, 0.06, 0.04])
    btn_play = Button(ax_play, 'Play')
    ax_reset = fig.add_axes([0.86, 0.035, 0.07, 0.04])
    btn_reset = Button(ax_reset, 'Reset')

    state = {'playing': False, 'timer': None}

    def update(val):
        idx = int(slider.val)
        idx = max(0, min(idx, n_snaps - 1))
        T = T_history[idx]
        pcm.set_array(T)
        line_inner.set_ydata(inner_hist[idx])
        line_outer.set_ydata(outer_hist[idx])
        t = float(t_array[idx])
        phase = 'FIRING' if t <= runtime else 'COOLING'
        title_heat.set_text(f'{phase}  |  t = {t:.3f} s   (snapshot {idx+1}/{n_snaps})')
        title_profile.set_text(f'Wall profile at t = {t:.3f} s')
        cursor.set_xdata([t, t])
        fig.canvas.draw_idle()

    slider.on_changed(update)

    def _tick(event=None):
        if not state['playing']:
            return
        nxt = int(slider.val) + 1
        if nxt >= n_snaps:
            nxt = 0
        slider.set_val(nxt)

    def _toggle_play(event):
        state['playing'] = not state['playing']
        btn_play.label.set_text('Pause' if state['playing'] else 'Play')
        if state['playing']:
            if state['timer'] is None:
                state['timer'] = fig.canvas.new_timer(interval=80)
                state['timer'].add_callback(_tick)
            state['timer'].start()
        else:
            if state['timer'] is not None:
                state['timer'].stop()
        fig.canvas.draw_idle()

    def _reset(event):
        if state['playing']:
            _toggle_play(event)
        slider.set_val(0)

    btn_play.on_clicked(_toggle_play)
    btn_reset.on_clicked(_reset)

    update(0)
    plt.show()


if __name__ == '__main__':
    if len(sys.argv) > 1:
        path = sys.argv[1]
    else:
        path = _default_results_path()
    main(path)
