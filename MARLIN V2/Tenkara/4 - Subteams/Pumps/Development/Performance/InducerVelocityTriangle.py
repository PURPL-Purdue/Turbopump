
"""Inducer inlet and outlet velocity triangles with no inlet preswirl."""

import numpy as np
from pint import Quantity as Q_

from VelocityTriangle import VelocityTriangle


def GetInletVelocities(
    flow_Q: Q_, speed_N: Q_, d_inlet: Q_, d_hub: Q_
) -> VelocityTriangle:
    """Return the inlet velocity triangle at the inducer tip radius.

    Supply Pint quantities: volumetric flow rate in L/s, shaft speed in rpm,
    inlet outer diameter and hub diameter in meters. Assumes uniform axial
    flow, no inlet preswirl, and no blade blockage.
    """
    flow_Q = flow_Q.to("m^3/s")
    speed_N = speed_N.to("rad/s")
    d_inlet = d_inlet.to("m")
    d_hub = d_hub.to("m")

    if not all(np.isfinite(value.magnitude) for value in
               (flow_Q, speed_N, d_inlet, d_hub)):
        raise ValueError("Inputs must be finite scalar quantities.")
    if d_inlet.magnitude <= 0 or not 0 <= d_hub.magnitude < d_inlet.magnitude:
        raise ValueError("Diameters must satisfy 0 <= d_hub < d_inlet.")
    if flow_Q.magnitude < 0 or speed_N.magnitude < 0:
        raise ValueError("Flow rate and rotational speed must be nonnegative.")

    inlet_area = np.pi / 4 * (d_inlet**2 - d_hub**2)
    c_m = (flow_Q / inlet_area).to("m/s")  #meridional velocity 
    c_u = Q_(0, "m/s") #no preswirl#
    u = ((speed_N * d_inlet)/ 2).to("m/s") #tip velocity
    w = np.sqrt(c_m**2 + (u - c_u)**2)

    return VelocityTriangle(u=u, c_m=c_m, c_u=c_u, w=w)


def GetOutletVelocities(
    flow_Q: Q_, speed_N: Q_, d_outlet: Q_, d_hub: Q_, head_rise: Q_
) -> VelocityTriangle:
    """Return the outlet velocity triangle at the inducer tip radius.

    Supply Pint quantities for volumetric flow rate, shaft speed, outlet tip
    and hub diameters, and ideal inducer head rise (length, e.g. meters).
    Assumes uniform axial flow, no blockage, and no inlet preswirl. Euler's
    pump equation gives c_u = g * head_rise / u. Hydraulic losses are neglected;
    measured head must be converted to ideal head before calling this function.
    """
    head_rise = head_rise.to("m")
    if not np.isfinite(head_rise.magnitude) or head_rise.magnitude < 0:
        raise ValueError("Head rise must be finite and nonnegative.")

    # The annular-area and tip-speed calculations also apply at the outlet.
    axial = GetInletVelocities(flow_Q, speed_N, d_outlet, d_hub)
    if axial.u.magnitude == 0:
        raise ValueError("Outlet calculation requires positive rotational speed.")

    gravity = Q_(9.80665, "m/s^2")
    c_u = (gravity * head_rise / axial.u).to("m/s")
    w = np.sqrt(axial.c_m**2 + (axial.u - c_u)**2)
    return VelocityTriangle(u=axial.u, c_m=axial.c_m, c_u=c_u, w=w)


def GetInletBladeAngle(c_m: Q_, u: Q_) -> Q_:
    """Return tip angle from the tangent, assuming no preswirl or incidence.

    This relative flow angle equals the blade metal angle at zero incidence.
    """
    c_m_value = c_m.to("m/s").magnitude
    u_value = u.to("m/s").magnitude
    if c_m_value == 0 and u_value == 0:
        raise ValueError("Blade angle is undefined when both velocities are zero.")
    return Q_(np.arctan2(c_m_value, u_value), "rad").to("degree")


if __name__ == "__main__":
    import matplotlib.pyplot as plt
    from matplotlib.patches import Arc

    # Illustrative inputs; replace with your inducer design values.
    inlet = GetInletVelocities(
        flow_Q=Q_(1.6, "L/s"),
        speed_N=Q_(40000, "rpm"),
        d_inlet=Q_(0.01750748707, "m"),
        d_hub=Q_(0.005252246121, "m"),
    )
    beta_1 = GetInletBladeAngle(inlet.c_m, inlet.u)
    print(f"Inlet blade angle from tangent (zero incidence): {beta_1.magnitude:.2f} deg")
    inlet.Plot(
        title=f"Inducer Inlet Velocity Triangle (beta = {beta_1.magnitude:.2f} deg)",
        station=1,
        show_c=False,
    )
    inlet_figure = plt.gcf()
    inlet_figure.canvas.manager.set_window_title("Inducer Inlet Velocity Triangle")
    # Mark beta at the tip of U, between the tangent and the W leg.
    ax = plt.gca()
    u_plot = inlet.u.to("m/s").magnitude
    radius = 0.25 * inlet.w.to("m/s").magnitude
    ax.add_patch(Arc(
        (u_plot, 0), 2 * radius, 2 * radius,
        theta1=180 - beta_1.magnitude, theta2=180,
        color="black", linewidth=1.2,
    ))
    half_beta = beta_1.to("rad").magnitude / 2
    label_radius = 1.6 * radius
    ax.text(
        u_plot - label_radius * np.cos(half_beta),
        label_radius * np.sin(half_beta),
        rf"$\beta_1 = {beta_1.magnitude:.2f}^\circ$",
        ha="center", va="center", fontsize=10,
        bbox=dict(facecolor="white", edgecolor="none", alpha=0.8, pad=1),
    )
    # Illustrative outlet geometry and ideal head; replace with design values.
    outlet = GetOutletVelocities(
        flow_Q=Q_(1.6, "L/s"),
        speed_N=Q_(40000, "rpm"),
        d_outlet=Q_(0.01750748707, "m"),
        d_hub=Q_(0.005252246121, "m"),
        head_rise=Q_(10, "m"),
    )
    outlet.Plot(title="Inducer Outlet Velocity Triangle", station=2)
    outlet_figure = plt.gcf()
    outlet_figure.canvas.manager.set_window_title("Inducer Outlet Velocity Triangle")
    print(f"Outlet meridional velocity: {outlet.c_m.to('m/s'):.2f}")
    print(f"Outlet tangential velocity: {outlet.c_u.to('m/s'):.2f}")
    # Each Plot call creates its own figure; display both separate windows.
    plt.show()
