"""
Printing and plotting of run results. No physics here; the run scripts in
scripts/ call these.
"""
import numpy as np
import matplotlib.pyplot as plt

from heat_transfer import coolant_side
from util import provenance as prov

MM, PA_TO_PSI = 1e3, 1 / 6894.76


def print_contour_dims(d):
    mm = lambda v: v * 1e3
    print("\n--- Contour (sketch dimensions, mm) ---")
    print(f"  Chamber radius Rc              {mm(d['Rc']):8.2f}")
    print(f"  Throat radius Rt               {mm(d['Rt']):8.2f}")
    print(f"  Exit radius Re                 {mm(d['Re']):8.2f}   (eps {d['eps']:.3f})")
    print(f"  Cylinder length                {mm(d['L_cyl']):8.2f}")
    print(f"  Convergent cone half-angle     {d['theta_conv_deg']:8.1f} deg")
    print(f"  Upstream throat arc radius     {mm(d['R_up']):8.2f}  (centre {mm(d['Rt']+d['R_up']):.2f} above axis at throat)")
    print(f"  Downstream throat arc radius   {mm(d['R_dn']):8.2f}  (centre {mm(d['Rt']+d['R_dn']):.2f} above axis at throat)")
    print(f"  Bell                           {d['bell_percent']:8d} % of 15 deg cone length")
    print(f"  Bell initial angle theta_n     {d['theta_n_deg']:8.2f} deg")
    print(f"  Bell exit angle theta_e        {d['theta_e_deg']:8.2f} deg")
    print(f"  Bell length, throat -> exit    {mm(d['L_bell']):8.2f}")
    for k in ("N", "Q", "E"):
        print(f"  Bell point {k} (x, r)            {mm(d[k][0]):8.2f}, {mm(d[k][1]):.2f}")
    print(f"  Injector face -> throat        {mm(d['x_t']):8.2f}")
    print(f"  Injector face -> exit          {mm(d['L_total']):8.2f}")


def throat_index(x, dims):
    return int(np.argmin(np.abs(x - dims["x_t"])))


def print_summary(x, dims, out, jk, T_limit, limit_name):
    """Throat state, hottest wall vs the limit, coolant in -> out."""
    i_t = throat_index(x, dims)
    i_max = int(np.nanargmax(out["T_wg"]))
    psi = PA_TO_PSI
    print(f"  Throat (x = {x[i_t]*MM:.1f} mm):")
    print(f"    T_wg {out['T_wg'][i_t]:6.0f} K   T_wl {out['T_wl'][i_t]:6.0f} K   T_l {out['T_l'][i_t]:6.0f} K")
    print(f"    q    {out['q'][i_t]/1e6:6.1f} MW/m2   h_g {out['h_g'][i_t]/1e3:5.1f}   "
          f"h_l {out['h_l'][i_t]/1e3:5.1f} kW/m2K")
    print(f"    v_l  {out['v_l'][i_t]:6.1f} m/s   Re {out['Re_l'][i_t]:8.0f}")
    if jk.fins:
        print(f"    fins phi {out['phi'][i_t]:5.2f}   eta {out['eta_fin'][i_t]:5.3f}   "
              f"h_l phi {out['h_l'][i_t]*out['phi'][i_t]/1e3:5.1f} kW/m2K")
    # 0.5 K slack so a channel-solver run that sits exactly on its target (to the
    # root-finder's tolerance) doesn't read as EXCEEDED. Display only.
    print(f"  Max T_wg {out['T_wg'][i_max]:.0f} K at x = {x[i_max]*MM:.1f} mm, {limit_name} "
          f"{T_limit:.0f} K -> " + ("OK" if out["T_wg"][i_max] <= T_limit + 0.5 else "EXCEEDED"))
    print(f"  Coolant in {jk.T_l_in:.0f} K / {jk.p_l_in*psi:.0f} psi -> out {out['T_l_out']:.0f} K / "
          f"{out['p_l_out']*psi:.0f} psi   (dT {out['T_l_out'] - jk.T_l_in:.0f} K, "
          f"dp {(jk.p_l_in - out['p_l_out'])*psi:.0f} psi)")


def print_checks(out, jk, cfg):
    """Validity warnings: correlation Re range, property table range, injection pressure."""
    fluid = jk.fluid
    Re_lo, Re_min = np.nanmin(out["Re_l"]), coolant_side.Re_min(jk.hl_name)
    print(f"  Re range {Re_lo:.0f} - {np.nanmax(out['Re_l']):.0f}"
          + (f"   ** below {Re_min:,.0f}: {jk.hl_name} out of range **" if Re_lo < Re_min else ""))
    n_lam = int(np.sum(out["Re_l"] <= coolant_side.RE_LAM))
    n_tr = int(np.sum((out["Re_l"] > coolant_side.RE_LAM) & (out["Re_l"] < coolant_side.RE_TURB)))
    if n_lam or n_tr:
        print(f"  {n_tr} stations in transition (Re {coolant_side.RE_LAM:.0f}-{coolant_side.RE_TURB:.0f}), "
              f"{n_lam} laminar; entrance effects {'on' if jk.entrance else 'off'}")
    T_wl_max = np.nanmax(out["T_wl"])
    if T_wl_max > fluid.T_MAX:
        print(f"  ** T_wl reaches {T_wl_max:.0f} K, above the {fluid.NAME} property table "
              f"({fluid.T_MAX:.0f} K; {fluid.LIMIT_NOTE}) - wall-side properties clamped **")
    if out["T_l_out"] > fluid.T_MAX:
        print(f"  ** coolant bulk reaches {out['T_l_out']:.0f} K, above the {fluid.NAME} "
              f"property table ({fluid.LIMIT_NOTE}) **")
    if out["p_l_out"] < cfg.p0_Pa:
        print(f"  ** coolant leaves at {out['p_l_out']*PA_TO_PSI:.0f} psi, below chamber "
              f"pressure ({cfg.p0_psi:.0f} psi): it can't be injected **")


def print_station_table(x, r, dims, chan, gas, out, status=None, spacing=0.010):
    """Every ~10 mm plus the throat and the exit. status marks clamped heights (channel solver)."""
    i_t = throat_index(x, dims)
    fins = not np.all(out["phi"] == 1)              # fin columns only when fins are on
    print("\n   x[mm]   r[mm]   w[mm]   h[mm]  land[mm]  T_r[K]  T_wg[K]  T_wl[K]  T_l[K]  q[MW/m2]"
          "  h_g[kW/m2K]  h_l[kW/m2K]" + ("   phi    eta" if fins else "") + "  v_l[m/s]  p_l[psi]")
    step = max(1, int(round(spacing / np.mean(np.diff(x)))))
    for i in sorted(set(range(0, len(x), step)) | {i_t, len(x) - 1}):
        tag = " <- throat" if i == i_t else (
            f"  ({status[i]})" if status is not None and status[i] != "ok" else "")
        print(f"  {x[i]*MM:6.1f}  {r[i]*MM:6.2f}  {chan['w'][i]*MM:6.3f}  {chan['h'][i]*MM:6.3f}"
              f"   {chan['land'][i]*MM:6.2f}   {gas['T_r'][i]:6.0f}   {out['T_wg'][i]:6.0f}"
              f"   {out['T_wl'][i]:6.0f}  {out['T_l'][i]:6.0f}    {out['q'][i]/1e6:5.1f}"
              f"       {out['h_g'][i]/1e3:6.2f}       {out['h_l'][i]/1e3:6.1f}"
              + (f"  {out['phi'][i]:5.2f}  {out['eta_fin'][i]:5.3f}" if fins else "")
              + f"     {out['v_l'][i]:5.1f}    {out['p_l'][i]*PA_TO_PSI:6.0f}{tag}")


def plot_run(x, r, dims, chan, out, title, T_limit, limit_name):
    """Contour, channel, temperatures, heat flux, h_g / h_l, coolant pressure vs x."""
    fig, ax = plt.subplots(6, 1, sharex=True, figsize=(11, 14))
    xm = x * MM
    ax[0].plot(xm, r * MM, "k"); ax[0].set_ylabel("r [mm]"); ax[0].set_aspect("equal")
    ax[1].plot(xm, chan["h"] * MM, label="h (height)"); ax[1].plot(xm, chan["w"] * MM, label="w (width)")
    ax[1].set_ylabel("channel [mm]")
    ax[2].plot(xm, out["T_wg"], label="T_wg (gas-side wall)")
    ax[2].plot(xm, out["T_wl"], label="T_wl (coolant-side wall)")
    ax[2].plot(xm, out["T_l"], label="T_l (coolant bulk)")
    ax[2].axhline(T_limit, color="r", ls="--", lw=0.8, label=f"{limit_name} {T_limit:.0f} K")
    ax[2].set_ylabel("T [K]")
    ax[3].plot(xm, out["q"] / 1e6); ax[3].set_ylabel("q [MW/m²]")
    ax[4].plot(xm, out["h_g"] / 1e3, label="h_g"); ax[4].plot(xm, out["h_l"] / 1e3, label="h_l")
    if not np.all(out["phi"] == 1):
        ax[4].plot(xm, out["h_l"] * out["phi"] / 1e3, label="h_l·phi (with fins)")
    ax[4].set_ylabel("h [kW/m²K]")
    ax[5].plot(xm, out["p_l"] * PA_TO_PSI); ax[5].set_ylabel("p_l [psi]")
    ax[5].set_xlabel("x from injector face [mm]")
    for a in ax:
        a.axvline(dims["x_t"] * MM, color="k", ls=":", lw=0.8); a.grid(alpha=0.3)
        if a.get_legend_handles_labels()[0]:                  # legend outside, right of the plot
            a.legend(fontsize=8, loc="upper left", bbox_to_anchor=(1.01, 1), borderaxespad=0)
    fig.suptitle(title, fontsize=10)
    plt.tight_layout(rect=(0, 0.015, 1, 1)); prov.stamp_figure(fig); plt.show()
