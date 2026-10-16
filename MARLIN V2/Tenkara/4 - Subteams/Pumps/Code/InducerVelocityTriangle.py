"""Plot inducer tip velocity triangles using the sizing results from Inducer.py."""

import math

import matplotlib.pyplot as plt
from matplotlib.patches import Arc

from Inducer import DESIGN_INPUTS, InducerResults, get_velocity_triangles, size_inducer


def plot_velocity_triangles(results: InducerResults) -> None:
    """Draw separate inlet/outlet graphs; all flow calculations live in Inducer.py.

    Gamma is the relative flow angle from the tangent, not the blade metal
    angle. Velocities from the US-unit sizing model are displayed in m/s.
    """
    inlet, outlet = get_velocity_triangles(results)
    stations = (
        (1, "Inlet", inlet, results.inlet_flow_angle_rad),
        (2, "Outlet", outlet, results.outlet_flow_angle_rad),
    )
    for station, name, triangle, flow_angle in stations:
        angle_deg = math.degrees(flow_angle)
        triangle.Plot(
            title=f"Inducer {name} Velocity Triangle (flow angle = {angle_deg:.2f} deg)",
            station=station,
        )
        figure = plt.gcf()
        figure.canvas.manager.set_window_title(f"Inducer {name} Velocity Triangle")
        ax = plt.gca()

        # Mark the calculated relative flow angle between the tangent and W.
        u_plot = triangle.u.to("m/s").magnitude
        radius = 0.25 * triangle.w.to("m/s").magnitude
        ax.add_patch(Arc(
            (u_plot, 0), 2 * radius, 2 * radius,
            theta1=180 - angle_deg, theta2=180,
            color="black", linewidth=1.2,
        ))
        half_angle = flow_angle / 2
        label_radius = 1.6 * radius
        ax.text(
            u_plot - label_radius * math.cos(half_angle),
            label_radius * math.sin(half_angle),
            rf"$\gamma_{station} = {angle_deg:.2f}^\circ$",
            ha="center", va="center", fontsize=10,
            bbox=dict(facecolor="white", edgecolor="none", alpha=0.8, pad=1),
        )
        figure.tight_layout()


if __name__ == "__main__":
    plot_velocity_triangles(size_inducer(DESIGN_INPUTS))
    plt.show()
