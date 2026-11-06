# Zero-blockage inlet meanline velocity triangle of an inducer (NASA SP-8052 secs. 2.1.9, 2.2.1, 3.1.9).
import json
import sys
from pathlib import Path

import numpy as np
import plotly.graph_objects as go
from pint import Quantity as Q_

if not __package__:
    # Run directly as a script: put Code/ first and drop this folder so the Inducer package (not Inducer.py) is imported
    _here = Path(__file__).resolve().parent
    sys.path[:] = [p for p in sys.path if Path(p or '.').resolve() != _here]
    sys.path.insert(0, str(_here.parent))

from Inducer.Inducer import Inducer, InducerDesignPoint, InducerInputGeometry
from VelocityTriangle import PlotTriangles, VelocityTriangle

INCIDENCE = 2 # deg
MEANLINE_TRIANGLE_PATH = Path(__file__).with_name("meanline_triangle.json")

class MeanlineTriangle:
    """Inlet triangle on the rms-radius meanline with no blade blockage and zero prewhirl.

    c_m = Q / A_1 is uniform over the annulus (eq. 23), so the flow angle is gamma_1 = arctan(c_m / u_ms).
    With an incidence ratio alpha/beta the blade angle is beta_1 = gamma_1 + alpha_1 = gamma_1 / (1 - alpha/beta)
    (secs. 2.1.9 / 3.1.9, ~0.425); None means zero incidence.
    Alternatively a fixed incidence angle can be given, so beta_1 = gamma_1 + incidence; incidence_ratio is then
    derived as incidence / beta_1. Give at most one of the two.
    """

    def __init__(self, inducer: Inducer, incidence_ratio: float | None = None, incidence: Q_ | None = None):
        if incidence_ratio is not None and incidence is not None:
            raise ValueError("Give either incidence_ratio or incidence, not both.")
        self.inducer = inducer

        self.R_t = inducer.R_tip.to('m').magnitude
        self.R_h1 = inducer.R_hub_1.to('m').magnitude
        self.omega = inducer.DP.N_shaft.to('rad/s').magnitude
        Q = inducer.DP.Q.to('m^3/s').magnitude

        self.r_ms = np.sqrt((self.R_t**2 + self.R_h1**2) / 2.0)  # rms radius, splits the annulus into equal areas
        self.c_m = Q / (np.pi * (self.R_t**2 - self.R_h1**2))
        self.gamma1 = np.arctan(self.c_m / (self.omega * self.r_ms))
        if incidence is not None:
            self.alpha1 = incidence.to('rad').magnitude
            self.beta1 = self.gamma1 + self.alpha1
            incidence_ratio = self.alpha1 / self.beta1
        elif incidence_ratio is None:
            self.beta1, self.alpha1 = self.gamma1, 0.0
        else:
            self.beta1 = self.gamma1 / (1.0 - incidence_ratio)
            self.alpha1 = incidence_ratio * self.beta1
        self.incidence_ratio = incidence_ratio

    @property
    def lead(self) -> Q_:
        """Zero-swirl flow lead lambda = c_m / omega (eq. 39 with psi = 0)."""
        return Q_(self.c_m / self.omega, 'm')

    def Triangle(self, notation: str = '<sub>1</sub>') -> VelocityTriangle:
        return VelocityTriangle(
            u=Q_(self.omega * self.r_ms, 'm/s'),
            c_m=Q_(self.c_m, 'm/s'),
            c_u=Q_(0.0, 'm/s'),
            notation=notation,
            incidence=Q_(self.alpha1, 'rad') if self.alpha1 else None,
        )

    def Plot(self, title: str = 'Inducer inlet meanline velocity triangle, zero blade blockage') -> go.Figure:
        return PlotTriangles([self.Triangle()], names=['Inlet (no blockage)'], title=title)

    def Save(self, path: Path = MEANLINE_TRIANGLE_PATH) -> None:
        """Write V, W, U, the flow angle beta1 (gamma_1 in SP-8052), the incidence and the blade angle to a JSON file, each with its unit."""
        triangle = self.Triangle()
        output = {
            "absolute_velocity": {"value": float(triangle.c.to('m/s').m), "unit": "m/s"},
            "relative_velocity": {"value": float(triangle.w.to('m/s').m), "unit": "m/s"},
            "rotational_velocity": {"value": float(triangle.u.to('m/s').m), "unit": "m/s"},
            "beta1": {"value": float(np.degrees(self.gamma1)), "unit": "deg"},
            "incidence": {"value": float(np.degrees(self.alpha1)), "unit": "deg"},
            "beta1_blade": {"value": float(np.degrees(self.beta1)), "unit": "deg"},
        }
        with open(path, "w", encoding="utf-8") as file:
            json.dump(output, file, indent=4)


def _quantity(data: dict) -> Q_:
    return Q_(data['value'], data['unit'])


def InducerFromJson(propellant: str = 'LOX') -> Inducer:
    """Inducer from design_parameters.json (flow, speed, head, blade count, impeller hub) and hub_geometry.json (tip, inlet hub, length).

    hub_geometry.json holds a single (LOX) geometry, so it only fits propellant 'LOX'.
    """
    here = Path(__file__).resolve().parent
    with (here / 'design_parameters.json').open(encoding='utf-8') as file:
        design_parameters = json.load(file)
    with (here / 'hub_geometry.json').open(encoding='utf-8') as file:
        hub_geometry = json.load(file)

    record = design_parameters['propellants'][propellant]
    operating_point = record['operating_point']
    inducer_design = record['inducer_design']
    sizing_constants = inducer_design['sizing_constants']
    tip_diameter = _quantity(hub_geometry['inducer_tip_diameter'])
    axial_length = _quantity(hub_geometry['inducer_axial_length'])

    return Inducer(
        InducerDesignPoint(
            Q=_quantity(operating_point['volumetric_flow_rate']),
            N_shaft=_quantity(operating_point['shaft_speed']),
            H_rise=_quantity(inducer_design['required_head_rise']),
        ),
        InducerInputGeometry(
            d_tip=tip_diameter,
            d_hub_1=_quantity(hub_geometry['inducer_inlet_hub_diameter']),
            d_hub_2=_quantity(record['impeller_inlet_geometry']['hub_diameter']),
            Z_blade=inducer_design['blade_parameters']['blade_count'],
            axial_length_ratio=(axial_length / tip_diameter).to('').magnitude,
        ),
        lambda_c=sizing_constants['lambda_c'],
        lambda_w=sizing_constants['lambda_w'],
    )


if __name__ == "__main__":
    inducer = InducerFromJson()
    ml =  MeanlineTriangle(inducer, incidence=Q_(INCIDENCE, 'deg'))
    print(f"--- LOX meanline values at incidence {INCIDENCE:g} deg ---")
    print(f"c_m = {ml.c_m:.2f} m/s, gamma1 = {np.degrees(ml.gamma1):.2f} deg, "
          f"beta1B = {np.degrees(ml.beta1):.2f} deg, "
          f"incidence = {np.degrees(ml.alpha1):.2f} deg")
    ml.Save()
    ml.Plot("LOX inducer inlet meanline triangle, zero blade blockage").show()
    