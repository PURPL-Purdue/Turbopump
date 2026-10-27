from dataclasses import dataclass
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.patches import Arc, Polygon, Rectangle
from pint import Quantity as Q_

from Inducer import Inducer, constraints, design_points

@dataclass
class MeridionalViewSettings:
    solidity: float = 2.75                  #C_i / P_i at the tip (input)
    n_blades: int = 3                       #Z
    solidity_limits: tuple = (1.5, 3.0)     #allowable C_i / P_i at the tip
    attack_limits: tuple = (Q_(3.0, 'deg'), Q_(8.0, 'deg'))   #allowable incidence at the tip
    thickness_ratio: float = 0.06           #drawn vane thickness as a fraction of P_i


@dataclass
class MeridionalViewGeometry:
    circumference: Q_[float]
    pitch: Q_[float]
    chord: Q_[float]
    axial_length: Q_[float]
    tangential_extent: Q_[float]
    flow_angle: Q_[float]
    attack_angle: Q_[float]
    blade_angle: Q_[float]
    n_blades: int
    solidity: float


def meridional_view_geometry(inducer: Inducer, settings: MeridionalViewSettings) -> MeridionalViewGeometry:
    gamma = inducer.gamma_1.to('deg')                # from the velocity triangle
    circumference = np.pi * inducer.D_tip.to('m')
    pitch = circumference / settings.n_blades
    chord = settings.solidity * pitch
    axial_length = inducer.C.axial_length_ratio * inducer.D_tip.to('m')

    sin_theta = (axial_length / chord).to('').magnitude
    if sin_theta > 1.0:
        raise ValueError(f"Solidity {settings.solidity} is too low: the chord is shorter than "
                         f"the target axial length of {axial_length.to('mm'):.2f}.")
    theta = Q_(float(np.degrees(np.arcsin(sin_theta))), 'deg')
    alpha = theta - gamma

    lo, hi = settings.attack_limits
    sol_lo, sol_hi = settings.solidity_limits
    if not sol_lo <= settings.solidity <= sol_hi:
        raise ValueError(f"Solidity {settings.solidity} is outside the limits {sol_lo} to {sol_hi}.")
    if not lo <= alpha <= hi:
        raise ValueError(f"Solidity {settings.solidity} gives incidence {alpha:.2f}, outside the limits "
                         f"{lo:~.1f} to {hi:~.1f} (blade angle {theta:.2f}, flow angle {gamma:.2f}).")

    return MeridionalViewGeometry(
        circumference=circumference,
        pitch=pitch,
        chord=chord,
        axial_length=axial_length,
        tangential_extent=chord * np.cos(theta.to('rad').magnitude),
        flow_angle=gamma,
        attack_angle=alpha,
        blade_angle=theta,
        n_blades=settings.n_blades,
        solidity=settings.solidity,
    )


def print_meridional_view(geo: MeridionalViewGeometry, name: str = 'Inducer'):
    print(f"\n===== {name} Meridional View (tip diameter) =====")
    print(f"Number of Vanes (Z)      = {geo.n_blades}")
    print(f"Solidity (C_i / P_i)     = {geo.solidity:.2f}")
    print(f"Circumference (pi d_t)   = {geo.circumference.to('mm'):.3f}")
    print(f"Pitch (P_i)              = {geo.pitch.to('mm'):.3f}")
    print(f"Chord (C_i)              = {geo.chord.to('mm'):.3f}")
    print(f"Flow Angle (gamma)       = {geo.flow_angle.to('deg'):.2f}")
    print(f"Angle of Attack (alpha)  = {geo.attack_angle.to('deg'):.2f}")
    print(f"Blade Angle (theta)      = {geo.blade_angle.to('deg'):.2f}")
    print(f"Axial Length (L_a)       = {geo.axial_length.to('mm'):.3f}")
    print(f"Axial Length Ratio L/D_t = {(geo.axial_length / (geo.circumference / np.pi)).to('').magnitude:.3f}")
    print(f"Vane Wrap (C_i cos/pi d) = {(geo.tangential_extent / geo.circumference).magnitude * 360.0:.1f} deg")


def plot_meridional_view(inducer: Inducer, name: str = 'Inducer',
                        settings: MeridionalViewSettings = MeridionalViewSettings()):
    geo = meridional_view_geometry(inducer, settings)
    w = geo.circumference.to('mm').magnitude
    p = geo.pitch.to('mm').magnitude
    c = geo.chord.to('mm').magnitude
    h = geo.axial_length.to('mm').magnitude          # strip height = vane axial length
    theta = geo.blade_angle.to('rad').magnitude
    gamma = geo.flow_angle.to('rad').magnitude
    th_deg = geo.blade_angle.to('deg').magnitude
    gm_deg = geo.flow_angle.to('deg').magnitude
    cx = c * np.cos(theta)                           # tangential projection of one vane
    tb = settings.thickness_ratio * p                # drawn vane thickness
    u = 0.04 * w                                     # base spacing unit for annotations

    fig, ax = plt.subplots(figsize=(11, 5.8))
    fig.subplots_adjust(left=0.02, right=0.98, top=0.98, bottom=0.16)
    strip = Rectangle((0, 0), w, h, fill=False, lw=1.6, ec='black', zorder=4)
    ax.add_patch(strip)

    # Vanes: leading edge on the bottom (inlet) edge, trailing edge at the top.
    # Copies shifted by one circumference show the wrap through the left edge.
    le_x = [w - tb - j * p for j in range(geo.n_blades)]
    for x0 in le_x:
        for k in (-1, 0, 1):
            xs = x0 + k * w
            vane = Polygon([(xs, 0), (xs + tb, 0), (xs + tb - cx, h), (xs - cx, h)],
                           closed=True, fc='lightgray', ec='black', lw=1.0, zorder=3)
            ax.add_patch(vane)
            vane.set_clip_path(Rectangle((0, 0), w, h, transform=ax.transData))

    def dim(p1, p2, color='black'):
        ax.annotate('', xy=p2, xytext=p1,
                    arrowprops=dict(arrowstyle='<|-|>', color=color, lw=1.0, mutation_scale=9,
                                    shrinkA=0, shrinkB=0))

    def label(text, xy, ha='center', va='center', size=9, bbox=True, **kw):
        ax.text(*xy, text, ha=ha, va=va, fontsize=size, fontweight='bold',
                bbox=dict(facecolor='white', edgecolor='none', pad=1.5) if bbox else None,
                zorder=6, **kw)

    # pi d_t: full circumference, below the strip
    y_pd = -2.8 * u
    for x in (0, w):
        ax.plot([x, x], [0, y_pd - 0.4 * u], color='black', lw=0.8)
    dim((0, y_pd), (w, y_pd))
    label(r'$\pi\,d_t$', (w / 2, y_pd), size=11)

    # P_i: tangential pitch between the leading edges of two adjacent vanes
    y_p = -1.3 * u
    x_a, x_b = le_x[0] + tb / 2, le_x[1] + tb / 2
    for x in (x_a, x_b):
        ax.plot([x, x], [0, y_p - 0.4 * u], color='black', lw=0.8)
    dim((x_b, y_p), (x_a, y_p))
    label(r'$P_i$', ((x_a + x_b) / 2, y_p), size=11)

    # C_i: chord dimension, parallel to the vanes, offset above the strip
    nx, ny = np.sin(theta), np.cos(theta)
    d_off = h / ny + 1.6 * u
    le = np.array([le_x[0] + tb / 2, 0.0])
    te = np.array([le_x[0] + tb / 2 - cx, h])
    for pt in (le, te):
        ax.plot([pt[0], pt[0] + d_off * nx], [pt[1], pt[1] + d_off * ny],
                color='black', lw=0.8, zorder=2)
    dim(le + d_off * np.array([nx, ny]), te + d_off * np.array([nx, ny]))
    mid = (le + te) / 2 + d_off * np.array([nx, ny])
    label(r'$C_i$', mid, size=11)

    # theta: blade angle measured from the tangential line at a leading edge
    x_t = le_x[1] + tb / 2
    r_arc = 0.9 * p
    ax.plot([x_t - r_arc * 1.15, x_t], [0, 0], color='black', lw=0.8, zorder=5)
    ax.add_patch(Arc((x_t, 0), 2 * r_arc, 2 * r_arc, theta1=180 - th_deg, theta2=180,
                     color='black', lw=1.0, zorder=5))
    ax.text(x_t - r_arc * 1.05, 0.55 * u, r'$\theta$', fontsize=13, ha='center', va='bottom', zorder=6)

    # Relative inlet flow and angle of attack at another leading edge
    x_f = le_x[2] + tb / 2
    r_f = 0.55 * p
    ax.annotate('', xy=(x_f - r_f * np.cos(gamma), r_f * np.sin(gamma)), xytext=(x_f, 0),
                arrowprops=dict(arrowstyle='-|>', color='black', lw=1.3, mutation_scale=10,
                                shrinkA=0, shrinkB=0), zorder=5)
    r_alpha = 0.40 * p
    ax.add_patch(Arc((x_f, 0), 2 * r_alpha, 2 * r_alpha, theta1=180 - th_deg, theta2=180 - gm_deg,
                     color='black', lw=1.0, zorder=5))
    ang_x = x_f - 0.40 * p
    mid_g = (gamma + theta) / 2
    label('ANGLE\nOF ATTACK', (ang_x, -0.9 * u), va='top', size=8, bbox=False)
    ax.annotate('', xy=(x_f - r_alpha * np.cos(mid_g), r_alpha * np.sin(mid_g)),
                xytext=(ang_x, -0.7 * u),
                arrowprops=dict(arrowstyle='-', color='black', lw=0.7))
    rel_x = x_f + 0.20 * p
    label('RELATIVE\nINLET FLOW', (rel_x, -0.9 * u), va='top', size=8, bbox=False)
    ax.annotate('', xy=(x_f - 0.5 * r_f * np.cos(gamma), 0.5 * r_f * np.sin(gamma)),
                xytext=(rel_x - 0.5 * u, -0.7 * u),
                arrowprops=dict(arrowstyle='-', color='black', lw=0.7))

    # Rotation and axial flow arrows
    ax.annotate('', xy=(0.12 * w, -0.9 * u), xytext=(0.0, -0.9 * u),
                arrowprops=dict(arrowstyle='-|>', color='black', lw=1.5, mutation_scale=12))
    label('ROTATION', (0.06 * w, -1.7 * u), size=8, bbox=False)
    y_top = h + d_off * ny + 0.6 * u
    ax.annotate('', xy=(0.02 * w, h + 6.0 * u), xytext=(0.02 * w, h + 1.0 * u),
                arrowprops=dict(arrowstyle='-|>', color='black', lw=1.5, mutation_scale=12))
    label('AXIAL FLOW', (0.02 * w, h + 7.0 * u), size=9, bbox=False)

    ax.set_xlim(-0.06 * w, 1.06 * w)
    ax.set_ylim(-4.2 * u, max(y_top, h + 7.0 * u) + 1.5 * u)
    ax.set_aspect('equal', adjustable='box')
    ax.axis('off')

    fig.suptitle(f'{name}: MERIDIONAL VIEW OF BLADES AT THE TIP DIAMETER',
                 fontsize=11, fontweight='bold', y=0.115)
    fig.text(0.5, 0.05,
             f'Z = {geo.n_blades},  solidity = {geo.solidity:.2f},  '
             rf'$\gamma$ = {gm_deg:.2f}°,  $\alpha$ = {geo.attack_angle.magnitude:.2f}°,  '
             rf'$\theta$ = {th_deg:.2f}°,  $C_i$ = {c:.1f} mm,  $P_i$ = {p:.1f} mm,  $L_a$ = {h:.1f} mm',
             ha='center', fontsize=9)
    return fig


if __name__ == '__main__':
    settings = MeridionalViewSettings()
    for name, design_point in design_points.items():
        inducer = Inducer(design_point, constraints)
        print_meridional_view(meridional_view_geometry(inducer, settings), name)
        plot_meridional_view(inducer, name, settings)
    plt.show()