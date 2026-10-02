"""NPSH charts for the inducer, built on bladeParams.

Chart 1: NPSH required vs. tip diameter, one curve per blade cavitation number
         sigma_b (chosen value plus SIGMA_B_STEPS above and below), against the
         flat NPSH available.
Chart 2: NPSH vs. shaft speed with lines of constant suction specific speed.

Units follow bladeParams (gpm, rpm, ft); diameters are plotted in inches.
PNGs are saved next to this script.

    python npshCharts.py
"""

from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np

import BladeParams as bp

# --- Chart settings ---------------------------------------------------------
SIGMA_B = 0.05          # assumed blade cavitation number at the tip
SIGMA_B_STEP = 0.01     # spacing of the comparison curves
SIGMA_B_STEPS = 4       # curves on each side of SIGMA_B
LAMBDA_C = 1.2          # C_m^2/2g coefficient in NPSH_r (see bp.npsh_required)
D_TIP_RANGE_IN = (0.4, 2.0)

SPEED_RANGE_RPM = (5_000, 60_000)
N_SS_LINES = [5_000, 10_000, 15_000, 20_000, 30_000, 40_000]

OUT_DIR = Path(__file__).resolve().parent

# --- Colors (sequential blue ramp for magnitude, orange for the chosen curve)
RAMP = ["#86b6ef", "#6da7ec", "#5598e7", "#3987e5", "#2a78d6",
        "#256abf", "#1c5cab", "#184f95", "#104281", "#0d366b"]
CHOSEN = "#eb6834"
INK = "#0b0b0b"
INK_2 = "#52514e"
MUTED = "#898781"
GRID = "#e1e0d9"
AXIS = "#c3c2b7"


def ramp(k: int) -> list[str]:
    """k evenly spaced steps from the blue ramp, light to dark."""
    idx = np.linspace(0, len(RAMP) - 1, k).round().astype(int)
    return [RAMP[i] for i in idx]


def style_axes(ax) -> None:
    ax.grid(True, color=GRID, linewidth=0.8)
    ax.set_axisbelow(True)
    for side in ("top", "right"):
        ax.spines[side].set_visible(False)
    for side in ("left", "bottom"):
        ax.spines[side].set_color(AXIS)
    ax.tick_params(colors=INK_2)
    ax.xaxis.label.set_color(INK_2)
    ax.yaxis.label.set_color(INK_2)
    ax.title.set_color(INK)


def label_line_end(ax, x, y, y_max, text, bold=False) -> None:
    """Direct label where a line leaves the plot: right end or top edge."""
    x, y = np.asarray(x), np.asarray(y)
    j = np.nonzero(y <= y_max)[0][-1]
    at_right = j == len(x) - 1
    ax.annotate(
        text, (x[j], y[j]),
        xytext=(4, 0) if at_right else (-3, -2),
        textcoords="offset points", fontsize=8,
        va="center" if at_right else "top",
        ha="left" if at_right else "right",
        color=INK if bold else INK_2,
        fontweight="bold" if bold else "normal",
    )


def add_legend(ax) -> None:
    ax.legend(loc="upper left", bbox_to_anchor=(1.02, 1.0), fontsize=8,
              frameon=False, labelcolor=INK_2)


def plot_npsh_vs_tip_diameter(inp: bp.InducerInputs) -> plt.Figure:
    q, n, rho = inp.flow_gpm, inp.speed_rpm, inp.hub_tip_ratio
    d_in = np.linspace(*D_TIP_RANGE_IN, 400)
    d_ft = d_in / 12.0

    offsets = np.arange(-SIGMA_B_STEPS, SIGMA_B_STEPS + 1)
    sigmas = SIGMA_B + offsets * SIGMA_B_STEP
    if sigmas.min() <= 0:
        raise ValueError("SIGMA_B_STEP too large: sweep reaches sigma_b <= 0")
    colors = ramp(len(sigmas))

    y_max = 3.0 * inp.npsh_ft
    fig, ax = plt.subplots(figsize=(12, 6.5))
    for sigma, color, off in zip(sigmas, colors, offsets):
        npsh = [bp.npsh_required(d, q, n, rho, sigma, LAMBDA_C) for d in d_ft]
        chosen = off == 0
        ax.plot(
            d_in, npsh,
            color=CHOSEN if chosen else color,
            linewidth=2.5 if chosen else 1.5,
            zorder=3 if chosen else 2,
            label=rf"$\sigma_b$ = {sigma:.3f}" + (" (assumed)" if chosen else ""),
        )
        label_line_end(ax, d_in, npsh, y_max, f"{sigma:.3f}", bold=chosen)

        # Minimum-NPSH tip diameter for the assumed sigma_b.
        if chosen:
            i = int(np.argmin(npsh))
            ax.plot(d_in[i], npsh[i], "o", ms=8, color=CHOSEN,
                    mec="white", mew=2, zorder=4)
            ax.annotate(
                f"min NPSH$_r$ = {npsh[i]:.1f} ft\nat D$_t$ = {d_in[i]:.3f} in",
                (d_in[i], npsh[i]), xytext=(60, -50),
                textcoords="offset points", ha="left", fontsize=9, color=INK,
                arrowprops=dict(arrowstyle="-", color=INK_2, lw=0.8),
            )

    ax.axhline(inp.npsh_ft, color=INK, linestyle="--", linewidth=1.5,
               label=f"NPSH available = {inp.npsh_ft:.1f} ft", zorder=2.5)

    # Tip diameter bladeParams sizes from N_ss (for reference).
    d_design_in = bp.size_inducer(inp).tip_diameter_ft * 12.0
    ax.axvline(d_design_in, color=MUTED, linestyle=":", linewidth=1.5,
               label=f"bladeParams D$_t$ = {d_design_in:.3f} in")

    ax.set_xlim(D_TIP_RANGE_IN[0], D_TIP_RANGE_IN[1] * 1.06)
    ax.set_ylim(0, y_max)
    ax.set_xlabel("Inducer tip diameter D$_t$ (in)")
    ax.set_ylabel("NPSH required (ft)")
    ax.set_title(
        f"NPSH required vs. tip diameter  "
        f"(Q = {q:.1f} gpm, n = {n:,.0f} rpm, "
        rf"hub/tip = {rho}, $\lambda_c$ = {LAMBDA_C})",
        loc="left",
    )
    style_axes(ax)
    add_legend(ax)
    fig.tight_layout()
    return fig


def plot_npsh_vs_speed(inp: bp.InducerInputs) -> plt.Figure:
    q = inp.flow_gpm
    speeds = np.linspace(*SPEED_RANGE_RPM, 400)
    colors = ramp(len(N_SS_LINES))

    fig, ax = plt.subplots(figsize=(12, 6.5))
    y_max = 3.0 * inp.npsh_ft
    for n_ss, color in zip(N_SS_LINES, colors):
        npsh = bp.npsh_for_suction_specific_speed(speeds, q, n_ss)
        ax.plot(speeds, npsh, color=color, linewidth=1.5,
                label=f"N$_{{ss}}$ = {n_ss:,}")
        label_line_end(ax, speeds, npsh, y_max, f"{n_ss:,}")

    ax.axhline(inp.npsh_ft, color=INK, linestyle="--", linewidth=1.5,
               label=f"NPSH available = {inp.npsh_ft:.1f} ft")

    n_ss_design = bp.suction_specific_speed(inp.speed_rpm, q, inp.npsh_ft)
    ax.plot(inp.speed_rpm, inp.npsh_ft, "o", ms=8, color=CHOSEN,
            mec="white", mew=2, zorder=4,
            label=f"Design point (N$_{{ss}}$ = {n_ss_design:,.0f})")
    ax.annotate(
        f"{inp.speed_rpm:,.0f} rpm\nN$_{{ss}}$ = {n_ss_design:,.0f}",
        (inp.speed_rpm, inp.npsh_ft), xytext=(8, -26),
        textcoords="offset points", fontsize=9, color=INK,
    )

    ax.set_xlim(SPEED_RANGE_RPM[0], SPEED_RANGE_RPM[1] * 1.06)
    ax.set_ylim(0, y_max)
    ax.xaxis.set_major_formatter(plt.FuncFormatter(lambda x, _: f"{x:,.0f}"))
    ax.set_xlabel("Shaft speed n (rpm)")
    ax.set_ylabel("NPSH (ft)")
    ax.set_title(
        f"NPSH vs. shaft speed at constant N$_{{ss}}$  (Q = {q:.1f} gpm)",
        loc="left",
    )
    style_axes(ax)
    add_legend(ax)
    fig.tight_layout()
    return fig


if __name__ == "__main__":
    inp = bp.DESIGN_INPUTS
    charts = {
        "npsh_vs_tip_diameter.png": plot_npsh_vs_tip_diameter(inp),
        "npsh_vs_shaft_speed.png": plot_npsh_vs_speed(inp),
    }
    for name, fig in charts.items():
        fig.savefig(OUT_DIR / name, dpi=200)
        print(f"saved {OUT_DIR / name}")
    plt.show()
