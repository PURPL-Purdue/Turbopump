# Create hub geometry for inducer

import json
from pathlib import Path

import numpy as np
import plotly.graph_objects as go
from pint import Quantity as Q_

g = Q_(9.81, 'm/s^2') # pint Quantity, m/s^2
MIN_INLET_HUB_TO_TIP_RATIO = 0.15 # Gülich recommended. Unit: dimensionless
MAX_INLET_HUB_TO_TIP_RATIO = 0.6 # Idk where ask Aarya. Unit: dimensionless

with open(Path(__file__).with_name("design_parameters.json"), "r", encoding="utf-8") as file:
    design_parameters = json.load(file)

lox_parameters = design_parameters["propellants"]["LOX"]

# Values below are plain floats read from the JSON, in the units noted (the JSON stores each unit next to its value)
impeller_tip_dia = lox_parameters["impeller_inlet_geometry"]["tip_diameter"]["value"] # Unit: m
impeller_hub_dia = lox_parameters["impeller_inlet_geometry"]["hub_diameter"]["value"] # Unit: m
volumetric_flow_rate = lox_parameters["operating_point"]["volumetric_flow_rate"]["value"] # Unit: m^3/s
shaft_speed = lox_parameters["operating_point"]["shaft_speed"]["value"] # Unit: rpm
lambda_c = lox_parameters["inducer_design"]["sizing_constants"]["lambda_c"] # Unit: dimensionless
lambda_w = lox_parameters["inducer_design"]["sizing_constants"]["lambda_w"] # Unit: dimensionless
axial_length_ratio = lox_parameters["inducer_design"]["sizing_constants"]["axial_length_ratio"] # L_ax / D_tip, unit: dimensionless
# psi_max = lox_parameters["inducer_design"]["sizing_constants"]["psi_max"] # head coefficient, unit: dimensionless
# req_head_rise = lox_parameters["inducer_design"]["required_head_rise"]["value"] # Unit: m

inducer_tip_dia = impeller_tip_dia # Assuming constant inducer tip diameter, constrained by impeller tip. Unit: m
inducer_axial_length = axial_length_ratio * inducer_tip_dia # Unit: m

hub_to_tip_ratios = np.linspace(MIN_INLET_HUB_TO_TIP_RATIO, MAX_INLET_HUB_TO_TIP_RATIO, 200) # Unit: dimensionless
k_n = 1 - hub_to_tip_ratios**2 # Unit: dimensionless
# The 3.25 coefficient is empirical: it assumes Q in m^3/s and n in rpm, and gives a diameter in m
optimal_tip_diameters = ( # Unit: m
    3.25
    * (volumetric_flow_rate / (shaft_speed * k_n)) ** (1 / 3)
    * ((lambda_c + lambda_w) / lambda_w) ** (1 / 6)
)


def GetHubToTipRatioAtTipDiameter(tip_diameter: float) -> float | None:
    """Interpolate the plotted curve's hub-to-tip ratio (dimensionless) at a given tip diameter (m)."""
    if not optimal_tip_diameters[0] <= tip_diameter <= optimal_tip_diameters[-1]:
        return None
    return float(np.interp(tip_diameter, optimal_tip_diameters, hub_to_tip_ratios))


def PlotOptimalTipDiameter() -> go.Figure:
    """Plot Gülich optimum tip diameter against hub-to-tip ratio."""
    fig = go.Figure()
    fig.add_trace(go.Scatter(
        x=hub_to_tip_ratios, # dimensionless
        y=optimal_tip_diameters * 1000, # m -> mm
        mode="lines",
        name="Gülich optimal tip diameter",
    ))
    fig.add_trace(go.Scatter(
        x=[MIN_INLET_HUB_TO_TIP_RATIO, MAX_INLET_HUB_TO_TIP_RATIO], # dimensionless
        y=[inducer_tip_dia * 1000, inducer_tip_dia * 1000], # m -> mm
        mode="lines",
        name="Impeller-constrained inducer tip diameter",
        line=dict(dash="dash"),
    ))
    fig.update_layout(
        title="LOX Inducer Optimal Tip Diameter vs Hub-to-Tip Ratio",
        xaxis_title="Hub-to-tip ratio",
        yaxis_title="Tip diameter (mm)",
        template="plotly_white",
    )
    return fig


def PlotMeridionalFlowArea(hub_tip_ratio: float) -> go.Figure:
    """Plot the inducer's meridional flow passage for the selected inlet hub-to-tip ratio (dimensionless). All plot coordinates are in mm."""
    if not 0 < hub_tip_ratio < 1:
        raise ValueError("Hub-to-tip ratio must be between 0 and 1.")

    length_mm = inducer_axial_length * 1000 # m -> mm
    tip_radius_mm = inducer_tip_dia * 500 # diameter in m -> radius in mm (x1000 / 2)
    inlet_hub_radius_mm = inducer_tip_dia * hub_tip_ratio * 500 # m -> mm radius
    outlet_hub_radius_mm = impeller_hub_dia * 500 # m -> mm radius

    fig = go.Figure()
    fig.add_trace(go.Scatter(
        x=[0, length_mm],
        y=[tip_radius_mm, tip_radius_mm],
        mode="lines",
        name="Shroud",
        line=dict(color="royalblue"),
    ))
    fig.add_trace(go.Scatter(
        x=[0, length_mm],
        y=[inlet_hub_radius_mm, outlet_hub_radius_mm],
        mode="lines",
        name="Hub",
        line=dict(color="darkorange"),
        fill="tonexty",
        fillcolor="rgba(100, 149, 237, 0.2)",
    ))
    fig.add_trace(go.Scatter(
        x=[0, 0],
        y=[inlet_hub_radius_mm, tip_radius_mm],
        mode="lines",
        name="Inlet",
        line=dict(color="grey", dash="dot"),
    ))
    fig.add_trace(go.Scatter(
        x=[length_mm, length_mm],
        y=[outlet_hub_radius_mm, tip_radius_mm],
        mode="lines",
        name="Outlet",
        line=dict(color="grey", dash="dot"),
    ))
    fig.update_layout(
        title=f"LOX Inducer Meridional Flow Area (hub-to-tip ratio = {hub_tip_ratio:.4f})",
        xaxis_title="Axial position (mm)",
        yaxis_title="Radius (mm)",
        template="plotly_white",
    )
    fig.update_yaxes(scaleanchor="x", scaleratio=1)
    
    # The angle is between a horizontal line (the axis direction) through the inlet hub point and the hub line.
    hub_taper_angle = GetHubTaperAngle(hub_tip_ratio) # Unit: degrees
    angle_rad = np.radians(hub_taper_angle) # Unit: radians
    arc_radius_mm = 0.4 * length_mm # Unit: mm
    arc_angles = np.linspace(0, angle_rad, 30) # Unit: radians
    fig.add_trace(go.Scatter(
        x=[0, length_mm],
        y=[inlet_hub_radius_mm, inlet_hub_radius_mm],
        mode="lines",
        name="Axis-parallel reference",
        line=dict(color="black", dash="dash", width=1),
    ))
    fig.add_trace(go.Scatter(
        x=arc_radius_mm * np.cos(arc_angles), # mm
        y=inlet_hub_radius_mm + arc_radius_mm * np.sin(arc_angles), # mm
        mode="lines",
        name="Hub taper angle",
        line=dict(color="black", width=2),
    ))
    fig.add_annotation(
        x=arc_radius_mm * np.cos(angle_rad / 2), # mm
        y=inlet_hub_radius_mm + arc_radius_mm * np.sin(angle_rad / 2), # mm
        text=f"Hub taper angle = {hub_taper_angle:.2f}°",
        showarrow=True,
        arrowhead=2,
        ax=60,
        ay=-30,
    )
    return fig

# def GetInletHubDiameter(tip_diameter: float) -> float:
#     """Calculate the inducer inlet hub diameter; tip_diameter and result are in meters."""
#     tip_diameter_q = Q_(tip_diameter, "m")
#     outlet_hub_diameter_q = Q_(impeller_hub_dia, "m")
#     if not 0 < outlet_hub_diameter_q < tip_diameter_q:
#         raise ValueError("The impeller hub diameter must be positive and smaller than the inducer tip diameter.")

#     angular_speed = Q_(shaft_speed, "rpm").to("rad/s")
#     tip_speed = angular_speed * tip_diameter_q / 2
#     psi_target = ((g * Q_(req_head_rise, "m")) / tip_speed**2).to("").m
#     area_ratio = np.sqrt(psi_max / psi_target)
#     outlet_area = (np.pi * (tip_diameter_q**2 - outlet_hub_diameter_q**2) / 4).to("m^2")
#     inlet_area = outlet_area / area_ratio
#     inlet_hub_diameter_squared = (tip_diameter_q**2 - 4 * inlet_area / np.pi).to("m^2")
#     if inlet_hub_diameter_squared.magnitude < 0:
#         raise ValueError("The required inlet flow area exceeds the inducer tip area.")
#     return np.sqrt(inlet_hub_diameter_squared).to("m").magnitude

def GetHubTaperAngle(hub_tip_ratio: float) -> float:
    """Hub taper angle in degrees; positive when the hub radius grows toward the outlet. hub_tip_ratio is dimensionless."""
    inlet_hub_radius = inducer_tip_dia * hub_tip_ratio / 2 # Unit: m
    outlet_hub_radius = impeller_hub_dia / 2 # Unit: m (same unit as the axial length, so the ratio is dimensionless)
    return float(np.degrees(np.arctan((outlet_hub_radius - inlet_hub_radius) / inducer_axial_length)))

HUB_GEOMETRY_PATH = Path(__file__).with_name("hub_geometry.json")


def SaveHubGeometry(axial_length: float, tip_diameter: float, inlet_hub_diameter: float, path: Path = HUB_GEOMETRY_PATH) -> None:
    """Write the inducer axial length, tip diameter and inlet hub diameter (all in m) to a JSON file, each stored with its unit."""
    output = {
        "inducer_axial_length": {"value": float(axial_length), "unit": "m"},
        "inducer_tip_diameter": {"value": float(tip_diameter), "unit": "m"},
        "inducer_inlet_hub_diameter": {"value": float(inlet_hub_diameter), "unit": "m"},
    }
    with open(path, "w", encoding="utf-8") as file:
        json.dump(output, file, indent=4)


if __name__ == "__main__":
    print(f"Inducer Axial Length: {inducer_axial_length * 1000:.4f} mm")
    
    inlet_hub_tip_ratio = GetHubToTipRatioAtTipDiameter(inducer_tip_dia) # Unit: dimensionless

    if inlet_hub_tip_ratio <= 0 or inlet_hub_tip_ratio >= 1:
        raise ValueError("The inlet hub-to-tip ratio must be between 0 and 1.")

    print(f"\n--- Initial Straight Hub Geometry ---")
    print(f"Hub-to-tip ratio at the impeller-constrained tip diameter: {inlet_hub_tip_ratio:.4f}")
    PlotOptimalTipDiameter().show()

    inducer_hub_dia = inducer_tip_dia * inlet_hub_tip_ratio # Unit: m
    print(f"Inducer tip diameter (constant): {inducer_tip_dia * 1000:.4f} mm")
    print(f"Inducer inlet hub diameter: {inducer_hub_dia * 1000:.4f} mm")
    SaveHubGeometry(inducer_axial_length, inducer_tip_dia, inducer_hub_dia)
    PlotMeridionalFlowArea(inlet_hub_tip_ratio).show()

    # hub_taper_angle = GetHubTaperAngle(inlet_hub_tip_ratio) # Unit: degrees
    # print(f"Hub taper angle: {hub_taper_angle:.2f} deg")