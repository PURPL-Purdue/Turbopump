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
    sigma: float = 0.06          #design cavitation parameter (sheet uses 0.03 to 0.06)
    K: float = 1.1               #empirical constant in phi = sqrt(sigma / (2 (K + sigma)))
    psi_1: float = 0.08          #assumed baseline head coefficient of the untapered inducer
    axial_length_ratio: float = 1.0   #L_a / D_t (axial length as a multiple of tip diameter)

class Inducer:
    def __init__(self, design_point: InducerDesignPoint, constraints: InducerInletDesignConstraints, g: Q_[float] = g):
        self.DP = design_point
        self.C = constraints
        self.g = g

        # Head rise across inducer (Delta H = NPSHr - NPSHa)
        self.Delta_H = (self.DP.NPSH_r - self.DP.NPSH_a).to('m')

        # ---------------------------------------------------------
        # Unit conversions (sheet methodology is SI with N in rpm)
        # ---------------------------------------------------------
        q_m3s = self.DP.Q.to('m^3/s').magnitude
        n_rpm = self.DP.N.to('rpm').magnitude
        sigma = self.C.sigma
        K = self.C.K
        rho = self.C.hub_tip_ratio

        # ---------------------------------------------------------
        # Flow Coefficient
        # ---------------------------------------------------------
        # phi = sqrt( sigma / (2 (K + sigma)) )
        self.phi: float = float(np.sqrt(sigma / (2.0 * (K + sigma))))

        # ---------------------------------------------------------
        # Inlet Diameter Sizing
        # ---------------------------------------------------------
        # D_t [m] = 3.25 * ( Q[m^3/s] / ( N[rpm] (1 - rho^2) ) )^(1/3) * ( (K + sigma) / sigma )^(1/6)
        # (3.25 is dimensional: it only holds for Q in m^3/s, N in rpm, D_t in m)
        d_t_m = 3.25 * (q_m3s / (n_rpm * (1.0 - rho**2))) ** (1.0 / 3.0) * ((K + sigma) / sigma) ** (1.0 / 6.0)

        self.D_tip: Q_[float] = Q_(d_t_m, 'm')
        self.D_hub: Q_[float] = self.D_tip * rho

        # ---------------------------------------------------------
        # Suction Specific Speed
        # ---------------------------------------------------------
        # Nss = ( 98 / (K + sigma)^(1/4) ) * sqrt( (1 - rho^2) / sigma )
        # Units: N in rpm, Q in m^3/s, NPSH in m (not the US-unit convention)
        self.N_ss: float = (98.0 / (K + sigma) ** 0.25) * np.sqrt((1.0 - rho**2) / sigma)

        # ---------------------------------------------------------
        # Estimated NPSH3
        # ---------------------------------------------------------
        # Invert Nss = N sqrt(Q) / NPSH^(3/4)  ->  NPSH3 = ( N sqrt(Q) / Nss )^(4/3)
        npsh3_m = (n_rpm * np.sqrt(q_m3s) / self.N_ss) ** (4.0 / 3.0)
        self.NPSH_3: Q_[float] = Q_(npsh3_m, 'm')

        # ---------------------------------------------------------
        # Target Head Coefficient
        # ---------------------------------------------------------
        # U_t = omega * R_t,  psi_2 = g * H / U_t^2
        r_t_m = self.D_tip.to('m').magnitude / 2.0
        omega = 2.0 * np.pi * n_rpm / 60.0
        self.U_tip: Q_[float] = Q_(omega * r_t_m, 'm/s')
        self.psi_2: float = float(
            self.g.to('m/s^2').magnitude * self.Delta_H.to('m').magnitude / self.U_tip.to('m/s').magnitude ** 2
        )

        # ---------------------------------------------------------
        # Hub Taper Angle (1D area continuity approximation)
        # ---------------------------------------------------------
        # R_h2 = sqrt( R_t^2 - (R_t^2 - R_h1^2) * sqrt(psi_1 / psi_2) )
        # alpha = atan( (R_h2 - R_h1) / L_a ),  L_a = axial_length_ratio * D_t
        r_h1_m = self.D_hub.to('m').magnitude / 2.0
        r_h2_m = np.sqrt(r_t_m**2 - (r_t_m**2 - r_h1_m**2) * np.sqrt(self.C.psi_1 / self.psi_2))
        l_ax_m = self.C.axial_length_ratio * self.D_tip.to('m').magnitude

        self.R_hub_1: Q_[float] = Q_(r_h1_m, 'm')
        self.R_hub_2: Q_[float] = Q_(float(r_h2_m), 'm')
        self.L_ax: Q_[float] = Q_(l_ax_m, 'm')
        self.alpha_taper: Q_[float] = Q_(float(np.degrees(np.arctan((r_h2_m - r_h1_m) / l_ax_m))), 'deg')