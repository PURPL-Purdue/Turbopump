import numpy as np
from pint import Quantity as Q_
from dataclasses import dataclass

@dataclass
class InducerDesignPoint:
    Q: Q_[float]                        # volumetric flow rate
    N_shaft: Q_[float]                  # shaft speed
    H_rise: Q_[float]                   # head the inducer has to supply

@dataclass
class InducerInputGeometry:
    d_tip: Q_[float]                    # tip diameter, constant from inlet to outlet
    d_hub_2: Q_[float]                  # outlet hub diameter
    Z_blade: int = 2                    # number of blades
    hub_tip_ratio: float = 0.5          # inlet d_hub_1 / d_tip
    axial_length_ratio: float = 0.6     # L_ax / d_tip

class Inducer:
    """Constant-tip inducer whose hub tapers linearly from R_hub_1 at the inlet to R_hub_2 at the outlet over L_ax.

    Suction metrics from the sheet relations, with the design cavitation parameter lambda_w (the sheet uses 0.03 to 0.06)
    and the empirical constant lambda_c.
    """

    def __init__(self, design_point: InducerDesignPoint, geometry: InducerInputGeometry,
                 lambda_c: float = 1.1, lambda_w: float = 0.03):
        self.DP = design_point
        self.Z_blade = geometry.Z_blade
        self.R_tip: Q_[float] = geometry.d_tip / 2
        self.R_hub_1: Q_[float] = geometry.hub_tip_ratio * self.R_tip
        self.R_hub_2: Q_[float] = geometry.d_hub_2 / 2
        self.L_ax: Q_[float] = geometry.axial_length_ratio * geometry.d_tip

        rho = geometry.hub_tip_ratio
        # phi = sqrt(lambda_w / (2 (lambda_c + lambda_w)))
        self.FlowCoefficient: float = float(np.sqrt(lambda_w / (2 * (lambda_c + lambda_w))))
        # N_ss = 98 / (lambda_c + lambda_w)^(1/4) sqrt((1 - rho^2) / lambda_w), with N in rpm, Q in m^3/s, NPSH in m
        self.N_ss: float = float(98 / (lambda_c + lambda_w)**0.25 * np.sqrt((1 - rho**2) / lambda_w))
        # backed out of N_ss = N sqrt(Q) / NPSH_3^(3/4)
        self.NPSH_3: Q_[float] = Q_((self.DP.N_shaft.to('rpm').m * np.sqrt(self.DP.Q.to('m^3/s').m) / self.N_ss)**(4 / 3), 'm')
