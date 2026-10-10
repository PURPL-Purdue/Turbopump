import numpy as np
from pint import Quantity as Q_
from dataclasses import dataclass

g = Q_(9.8, 'm/s^2')

@dataclass
class InducerDesignPoint:
    Q: Q_[float]
    N: Q_[float]
    NPSH_a: Q_[float]
    NPSH_r: Q_[float]

@dataclass
class InducerInletDesignConstraints:
    hub_tip_ratio: float = 0.3   #rho = D_h/D_t

class Inducer:
    def __init__(self, design_point: InducerDesignPoint, constraints: InducerInletDesignConstraints, g: Q_[float] = g):
        self.DP = design_point
        self.C = constraints
        self.g = g

        # Head rise across inducer (Delta H = NPSHr - NPSHa)
        self.Delta_H = (self.DP.NPSH_r - self.DP.NPSH_a).to('m')

        # ---------------------------------------------------------
        # Suction Specific Speed
        # ---------------------------------------------------------
        q_gpm = self.DP.Q.to('gallon/min').magnitude
        n_rpm = self.DP.N.to('rpm').magnitude
        npsh_ft = self.DP.NPSH_a.to('ft').magnitude
 
        # Standard Suction Specific Speed
        self.N_ss: float = n_rpm * np.sqrt(q_gpm) / (npsh_ft ** 0.75)

        # Corrected Suction Specific Speed
        self.N_ss_corrected: float = self.N_ss / np.sqrt(1.0 - self.C.hub_tip_ratio**2)

        # ---------------------------------------------------------
        # Optimum Flow Coefficient (Brumfield criterion)
        # ---------------------------------------------------------
        # 3574 is the empirical US-unit constant for optimum cavitation inception
        b_ratio = 3574.0 / self.N_ss_corrected
        self.phi: float = float(b_ratio / ((1.0 + np.sqrt(1.0 + 6.0 * (b_ratio**2))) / 2.0))

        # ---------------------------------------------------------
        # Inlet Diameter Sizing
        # ---------------------------------------------------------
        # D_t [ft] = 0.37843 * ( Q[gpm] / ( (1 - rho^2) * N[rpm] * phi ) )^(1/3)
        d_t_ft = 0.37843 * (q_gpm / ((1.0 - self.C.hub_tip_ratio**2) * n_rpm * self.phi)) ** (1.0 / 3.0)

        self.D_tip: Q_[float] = Q_(d_t_ft, 'ft').to('m')
        self.D_hub: Q_[float] = self.D_tip * self.C.hub_tip_ratio