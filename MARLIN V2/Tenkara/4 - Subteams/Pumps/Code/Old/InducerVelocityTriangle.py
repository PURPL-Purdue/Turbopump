"""Plot inducer tip velocity triangles calculated by Inducer.py."""

from math import hypot

import matplotlib.pyplot as plt
from matplotlib.patches import Arc

from Inducer import Inducer, constraints, design_points, print_parameters


def plot_velocity_triangles(inducer: Inducer, name: str = 'Inducer'):
    """Return inlet/outlet figures; all flow calculations live in Inducer.py.

    Assumes uniform annulus flow, zero inlet prewhirl, constant tip radius,
    and ideal Euler head rise for outlet swirl. Gamma is a flow angle
    measured from the negative tangent, not a blade metal angle.
    """
    figures = []
    stations = (
        ('Inlet', inducer.C_m1, inducer.C_u1, inducer.W_1, inducer.gamma_1),
        ('Outlet', inducer.C_m2, inducer.C_u2, inducer.W_2, inducer.gamma_2),
    )
    for station, (label, cm, cu, w, gamma) in enumerate(stations, start=1):
        u = inducer.U_tip.to('m/s').magnitude
        cm, cu, w = (v.to('m/s').magnitude for v in (cm, cu, w))
        angle = gamma.to('deg').magnitude
        fig, ax = plt.subplots(figsize=(9, 5))
        vectors = (
            ((0, 0), (u, 0), 'U', u, 'green', (0, -16)),
            ((0, 0), (cu, cm), 'C', hypot(cu, cm), 'blue', (-65, 12)),
            ((u, 0), (cu, cm), 'W', w, 'red', (0, 12)),
        )
        def component_label(symbol, value, position, offset, color):
            text = f'{symbol}{station} = {value:.2f} m/s'
            ax.annotate(text, xy=position, xytext=offset,
                        textcoords='offset points', ha='center', va='center',
                        color=color, fontsize=9,
                        bbox=dict(facecolor='white', edgecolor='none', alpha=0.85, pad=2))
            return text

        for start, end, symbol, value, color, offset in vectors:
            ax.annotate('', xy=end, xytext=start,
                        arrowprops=dict(arrowstyle='->', color=color, lw=2))
            midpoint = ((start[0] + end[0]) / 2, (start[1] + end[1]) / 2)
            text = component_label(symbol, value, midpoint, offset, color)
            ax.plot([], [], color=color, label=text)
        ax.plot([cu, cu], [0, cm], '--', color='orange', label=f'Cm{station} = {cm:.2f} m/s')
        component_label('Cm', cm, (cu, cm / 2), (65, 0), 'darkorange')
        ax.plot([0, cu], [cm, cm], '--', color='purple', label=f'Cu{station} = {cu:.2f} m/s')
        component_label('Cu', cu, (cu / 2, cm), (0, 16), 'purple')
        ax.set_title(f'{name} {label} Velocity Triangle\n'
                     f'U = {u:.2f}, Cm = {cm:.2f}, Cu = {cu:.2f}, W = {w:.2f} m/s')
        radius = 0.25 * w
        ax.add_patch(Arc((u, 0), 2 * radius, 2 * radius,
                         theta1=180-angle, theta2=180, color='black'))
        ax.text(u - radius, -0.08 * max(u, cm),
                rf'$\gamma_{station} = {angle:.2f}^\circ$', ha='center')
        ax.update_datalim([(0, 0), (u, 0), (cu, cm)])
        ax.autoscale_view()
        ax.margins(0.2)
        ax.set_aspect('equal', adjustable='box')
        ax.set_xlabel('Tangential velocity [m/s]')
        ax.set_ylabel('Meridional velocity [m/s]')
        ax.grid(True, alpha=0.3)
        ax.legend()
        fig.tight_layout()
        figures.append(fig)
    return tuple(figures)


if __name__ == '__main__':
    for name, design_point in design_points.items():
        inducer = Inducer(design_point, constraints)
        print_parameters(inducer, name)
        plot_velocity_triangles(inducer, name)
    plt.show()
