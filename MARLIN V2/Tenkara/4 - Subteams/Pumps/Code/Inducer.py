import numpy as np
from pint import Quantity as Q_
from dataclasses import dataclass

g = Q_(9.81, 'm/s^2')

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

        # Uniform annulus flow, constant tip radius, zero inlet prewhirl,
        # and ideal Euler head: g * Delta_H = U_tip * (C_u2 - C_u1).
        self.A_inlet = (np.pi * (self.D_tip**2 - self.D_hub**2) / 4).to('m^2')
        self.A_outlet = (np.pi * ((self.D_tip / 2)**2 - self.R_hub_2**2)).to('m^2')
        self.C_m1 = (self.DP.Q / self.A_inlet).to('m/s')
        self.C_m2 = (self.DP.Q / self.A_outlet).to('m/s')
        self.C_u1 = Q_(0.0, 'm/s')
        self.C_u2 = (self.g * self.Delta_H / self.U_tip).to('m/s')
        self.W_1 = np.sqrt(self.C_m1**2 + (self.U_tip - self.C_u1)**2)
        self.W_2 = np.sqrt(self.C_m2**2 + (self.U_tip - self.C_u2)**2)
        # Relative flow angles measured from the negative tangential direction.
        self.gamma_1 = np.arctan2(self.C_m1, self.U_tip - self.C_u1).to('deg')
        self.gamma_2 = np.arctan2(self.C_m2, self.U_tip - self.C_u2).to('deg')


# Shared nominal inputs for sizing and velocity-triangle scripts.
design_points = {
    'LOx': InducerDesignPoint(
        Q=Q_(1.69, 'L/s'), N=Q_(35000, 'rpm'),
        NPSH_a=Q_(105.46, 'm'), NPSH_r=Q_(150.0, 'm'),
    ),
    'IPA': InducerDesignPoint(
        Q=Q_(2.739, 'L/s'), N=Q_(35000, 'rpm'),
        NPSH_a=Q_(105.46, 'm'), NPSH_r=Q_(150.0, 'm'),
    ),
}
constraints = InducerInletDesignConstraints(
    hub_tip_ratio=0.5, sigma=0.03, K=1.1, psi_1=0.08, axial_length_ratio=1.0,
)


def print_parameters(inducer: Inducer, name: str = 'Inducer') -> None:
    """Print design inputs, sizing results, and inlet/outlet triangle components."""
    print(f'\n===== {name} Inducer =====')

    def print_group(title, parameters):
        print(f'\n{title}')
        for label, value in parameters:
            print(f'{label:<36} = {value:.4f}')

    print_group('Design inputs', vars(inducer.DP).items())
    print_group('Design constraints (dimensionless)', vars(inducer.C).items())
    print_group('Sizing results', (
        (label, getattr(inducer, label)) for label in (
            'g', 'Delta_H', 'phi', 'D_tip', 'D_hub', 'N_ss', 'NPSH_3',
            'psi_2', 'R_hub_1', 'R_hub_2', 'L_ax', 'alpha_taper',
        )
    ))
    print_group('Inlet velocity triangle', (
        ('A_inlet', inducer.A_inlet),
        ('U_tip (tip velocity)', inducer.U_tip),
        ('C_m1 (meridional velocity)', inducer.C_m1),
        ('C_u1 (tangential velocity)', inducer.C_u1),
        ('C_1 (absolute velocity)', np.sqrt(inducer.C_m1**2 + inducer.C_u1**2)),
        ('W_1 (relative velocity)', inducer.W_1),
        ('gamma_1 (relative flow angle)', inducer.gamma_1),
    ))
    print_group('Outlet velocity triangle', (
        ('A_outlet', inducer.A_outlet),
        ('U_tip (tip velocity)', inducer.U_tip),
        ('C_m2 (meridional velocity)', inducer.C_m2),
        ('C_u2 (tangential velocity)', inducer.C_u2),
        ('C_2 (absolute velocity)', np.sqrt(inducer.C_m2**2 + inducer.C_u2**2)),
        ('W_2 (relative velocity)', inducer.W_2),
        ('gamma_2 (relative flow angle)', inducer.gamma_2),
    ))


if __name__ == '__main__':
    for name, design_point in design_points.items():
        print_parameters(Inducer(design_point, constraints), name)
