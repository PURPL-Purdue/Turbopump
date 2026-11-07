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

from Inducer.Inducer import Inducer
from VelocityTriangle import PlotTriangles, VelocityTriangle

INCIDENCE = 2 # deg

MEANLINE_TRIANGLE_PATH = Path(__file__).with_name("meanline_triangle.json")

class MeanlineTriangle:
    """Inlet triangle on the rms-radius meanline with no blade blockage and zero prewhirl.

    c_m = Q / A_1 is uniform over the annulus (eq. 23), so the flow angle is gamma_1 = arctan(c_m / u_ms).
    The blade angle is beta_1b = gamma_1 + incidence (secs. 2.1.9 / 3.1.9). None means zero incidence.
    """

    def __init__(self, inducer: Inducer, incidence: Q_ | None = None):
        self.inducer = inducer

        self.R_t = inducer.R_tip.to('m').magnitude
        self.R_h1 = inducer.R_hub_1.to('m').magnitude
        self.omega = inducer.DP.N_shaft.to('rad/s').magnitude
        Q = inducer.DP.Q.to('m^3/s').magnitude

        self.r_ms = np.sqrt((self.R_t**2 + self.R_h1**2) / 2.0)  # rms radius, splits the annulus into equal areas
        self.c_m = Q / (np.pi * (self.R_t**2 - self.R_h1**2))
        self.gamma1 = np.arctan(self.c_m / (self.omega * self.r_ms))
        if incidence is not None:
            self.incidence = incidence.to('rad').magnitude
        else:
            self.incidence = 0.0
        self.beta1b = self.gamma1 + self.incidence

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
            incidence=Q_(self.incidence, 'rad') if self.incidence else None,
        )

    def Plot(self, title: str = 'Inducer inlet meanline velocity triangle, zero blade blockage') -> go.Figure:
        return PlotTriangles([self.Triangle()], names=['Inlet (no blockage)'], title=title)

    def Save(self, path: Path = MEANLINE_TRIANGLE_PATH) -> None:
        """Write V, W, U, the flow angle gamma1, incidence and blade angle beta1b to JSON, each with its unit."""
        triangle = self.Triangle()
        output = {
            "absolute_velocity": {"value": float(triangle.c.to('m/s').m), "unit": "m/s"},
            "relative_velocity": {"value": float(triangle.w.to('m/s').m), "unit": "m/s"},
            "rotational_velocity": {"value": float(triangle.u.to('m/s').m), "unit": "m/s"},
            "gamma1": {"value": float(np.degrees(self.gamma1)), "unit": "deg"},
            "incidence": {"value": float(np.degrees(self.incidence)), "unit": "deg"},
            "beta1b": {"value": float(np.degrees(self.beta1b)), "unit": "deg"},
        }
        with open(path, "w", encoding="utf-8") as file:
            json.dump(output, file, indent=4)


if __name__ == "__main__":
    import Main
    import Pump1D
    propellant = Main.Propellant.LOX
    pump = Main.CreatePump(propellant)
    pump.AddInducer(Pump1D.GetInducerlessNPSH_r(pump), geometry=Main.basic_inducer_geometry[propellant])
    ml = MeanlineTriangle(pump.Inducer, incidence=Q_(INCIDENCE, 'deg'))
    print(f"--- {propellant.name} meanline values at incidence {INCIDENCE:g} deg ---")
    print(f"c_m = {ml.c_m:.2f} m/s, gamma1 = {np.degrees(ml.gamma1):.2f} deg, "
          f"beta1b = {np.degrees(ml.beta1b):.2f} deg, "
          f"incidence = {np.degrees(ml.incidence):.2f} deg")
    # ml.Save()
    ml.Plot(f"{propellant.name} inducer inlet meanline triangle, zero blade blockage").show()
