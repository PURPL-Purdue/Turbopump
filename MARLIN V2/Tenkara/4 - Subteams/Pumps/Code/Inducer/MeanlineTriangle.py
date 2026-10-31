# Zero-blockage inlet meanline velocity triangle of an inducer (NASA SP-8052 secs. 2.1.9, 2.2.1, 3.1.9).
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


class MeanlineTriangle:
    """Inlet triangle on the rms-radius meanline with no blade blockage and zero prewhirl.

    c_m = Q / A_1 is uniform over the annulus (eq. 23), so the flow angle is gamma_1 = arctan(c_m / u_ms).
    With an incidence ratio alpha/beta the blade angle is beta_1 = gamma_1 + alpha_1 = gamma_1 / (1 - alpha/beta)
    (secs. 2.1.9 / 3.1.9, ~0.425); None means zero incidence.
    """

    def __init__(self, inducer: Inducer, incidence_ratio: float | None = None):
        self.inducer = inducer
        self.incidence_ratio = incidence_ratio

        self.R_t = inducer.D_tip.to('m').magnitude / 2.0
        self.R_h1 = inducer.R_hub_1.to('m').magnitude
        self.omega = inducer.DP.N_shaft.to('rad/s').magnitude
        Q = inducer.DP.Q.to('m^3/s').magnitude

        self.r_ms = np.sqrt((self.R_t**2 + self.R_h1**2) / 2.0)  # rms radius, splits the annulus into equal areas
        self.c_m = Q / (np.pi * (self.R_t**2 - self.R_h1**2))
        self.gamma1 = np.arctan(self.c_m / (self.omega * self.r_ms))
        if incidence_ratio is None:
            self.beta1, self.alpha1 = self.gamma1, 0.0
        else:
            self.beta1 = self.gamma1 / (1.0 - incidence_ratio)
            self.alpha1 = incidence_ratio * self.beta1

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
        )

    def Plot(self, title: str = 'Inducer inlet meanline velocity triangle, zero blade blockage') -> go.Figure:
        return PlotTriangles([self.Triangle()], names=['Inlet (no blockage)'], title=title)


if __name__ == "__main__":
    import Main
    print("--- Meanline values at incidence ratio of 0.425 ---")
    for propellant, H in ((Main.Propellant.LOX, 65), (Main.Propellant.IPA, 80)):
        pump = Main.CreatePump(propellant)
        pump.AddInducer(Q_(H, 'm'))
        ml = MeanlineTriangle(pump.Inducer, incidence_ratio=0.425)
        print(f"{propellant}: c_m = {ml.c_m:.2f} m/s, gamma1 = {np.degrees(ml.gamma1):.2f} deg, "
              f"beta1 = {np.degrees(ml.beta1):.2f} deg")
        ml.Plot(f"{propellant} inducer inlet meanline triangle, zero blade blockage").show()
