import numpy as np
import plotly.graph_objects as go
from pint import Quantity as Q_
from dataclasses import dataclass
from functools import cache
from scipy.optimize import brentq
from MeridionalContour import Polyline

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

@dataclass
class InducerGeometry:
    """Blades for the radial equilibrium that sets the blade angles (NASA SP-8052 sec. 2.2.1)"""
    Z_blade: int = 3                    # number of blades
    thk: Q_[float] = Q_(0, 'mm')        # constant normal blade thickness, blocks part of the annulus; 0: no blockage
    eta_bl: float = 0.85                # blade efficiency, g * Delta_H = eta_bl * Euler head (SP-8052 eq. 47, sec. 2.2.6)
    incidence_ratio: float = 0.0        # alpha / beta (~0.425, SP-8052 secs. 2.1.9, 3.1.9): offsets the flow angle from the
                                        # blade angle at the inlet, fading linearly to zero at the outlet; 0: no incidence

class Inducer:
    def __init__(self, design_point: InducerDesignPoint, constraints: InducerInletDesignConstraints,
                 geometry: InducerGeometry = InducerGeometry()):
        self.DP = design_point
        self.C = constraints
        self.Geometry = geometry

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
            g.to('m/s^2').magnitude * self.Delta_H.to('m').magnitude / self.U_tip.to('m/s').magnitude ** 2
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
        self.C_u2 = (g * self.Delta_H / self.U_tip).to('m/s')
        self.W_1 = np.sqrt(self.C_m1**2 + (self.U_tip - self.C_u1)**2)
        self.W_2 = np.sqrt(self.C_m2**2 + (self.U_tip - self.C_u2)**2)
        # Relative flow angles measured from the negative tangential direction.
        self.gamma_1 = np.arctan2(self.C_m1, self.U_tip - self.C_u1).to('deg')
        self.gamma_2 = np.arctan2(self.C_m2, self.U_tip - self.C_u2).to('deg')

        # ---------------------------------------------------------
        # Meridional Geometry
        # ---------------------------------------------------------
        self.BuildMeridian()

        # ---------------------------------------------------------
        # Blade Angles (radial equilibrium, NASA SP-8052 sec. 2.2.1)
        # ---------------------------------------------------------
        # At the rms radius (sec. 3.1.12), from the tangential direction.
        # Zero-swirl inlet (c_u = 0), uniform c_a = Q/A1 (eq. 23); flow lead lambda = c_a/omega (eq. 39 with psi = 0)
        gamma_1_rms = np.arctan(self._lam_flow(0.0) / self._r_ms(0.0))
        # gamma = beta - alpha, alpha/beta from sec. 2.1.9 / 3.1.9
        self.beta_1: Q_[float] = Q_(gamma_1_rms / (1.0 - self.Geometry.incidence_ratio), 'rad').to('deg')
        self.alpha_1: Q_[float] = self.Geometry.incidence_ratio * self.beta_1
        # The blade angle follows beta(z) = beta_1 + (beta_2 - beta_1) (z/L_ax)^2 (eq. 74), with beta_2 solved so the
        # mass-averaged outlet Euler head equals g * Delta_H / eta_bl (eq. 47)
        self.beta_2: Q_[float] = Q_(self._solve_beta2(), 'rad').to('deg')

    def BuildMeridian(self):
        """
        Meridional curves, (z, r) in m with z = 0 at the inducer inlet:
            shroud:         casing at the constant tip radius D_tip / 2
            hub:            linear taper from R_hub_1 to R_hub_2 over L_ax
            inlet, outlet:  stations 1 and 2, straight from hub to tip at z = 0 and z = L_ax
            mean:           rms radius sqrt((r_tip^2 + r_hub^2) / 2), which splits the annulus into equal areas
                            (NASA SP-8052, Glossary). With the uniform c_m assumed above, half the flow passes on each
                            side of it, and it is the radius at which the blade angles beta_1 and beta_2 are set
                            (SP-8052 sec. 3.1.12). MeanCurve(HubContour, ShroudContour) would give the geometric
                            meanline instead, as in Impeller.
        The radial equilibrium reads the hub and tip radii off these contours (_r_hub, _r_tip).

        Blades are not modeled yet. Their leading and trailing edges would be built here, hub to tip, like
        Impeller.LeadingEdge and Impeller.TrailingEdge.
        """
        r_t = self.D_tip.to('m').magnitude / 2
        r_h1 = self.R_hub_1.to('m').magnitude
        r_h2 = self.R_hub_2.to('m').magnitude
        l_ax = self.L_ax.to('m').magnitude

        # hub is a linear taper for now
        self.ShroudContour: Polyline = Polyline([[0, r_t], [l_ax, r_t]])
        self.HubContour: Polyline = Polyline([[0, r_h1], [l_ax, r_h2]])
        self.Inlet: Polyline = Polyline([[0, r_h1], [0, r_t]])
        self.Outlet: Polyline = Polyline([[l_ax, r_h2], [l_ax, r_t]])

        z = np.linspace(0, l_ax, 51)
        self.MeanContour: Polyline = Polyline(np.column_stack([z, self._r_ms(z)]))

    # ---------------------------------------------------------
    # Radial equilibrium (NASA SP-8052 sec. 2.2.1, eqs. 33-41, 74): swirl and axial velocity behind the constant-tip,
    # linearly tapered-hub inducer. Internals are SI floats, z in m from the inducer inlet, b2 a trial beta_2 in rad.
    # ---------------------------------------------------------
    def _r_hub(self, z):
        return np.interp(z, *self.HubContour.points.T)

    def _r_tip(self, z):
        return np.interp(z, *self.ShroudContour.points.T)

    def _r_ms(self, z):
        # rms = root-mean-square radius sqrt((r_tip^2 + r_hub^2)/2), the radius splitting the annulus into equal areas (Glossary)
        return np.sqrt((self._r_tip(z) ** 2 + self._r_hub(z) ** 2) / 2.0)

    def _beta(self, z, b2):
        # Parabolic camber distribution, SP-8052 eq. 74 (specified at the rms radius, sec. 3.1.12)
        b1 = self.beta_1.to('rad').m
        return b1 + (b2 - b1) * (z / self.L_ax.to('m').m) ** 2

    def _dA(self, r, lam, strict=True):
        # Open area per dr: 2 pi r minus Z_blade blades of tangential width thk/sin(beta), with sin(beta) = lam/sqrt(r^2 + lam^2)
        dA = 2.0 * np.pi * r - self.Geometry.Z_blade * self.Geometry.thk.to('m').m * np.sqrt(r**2 + lam**2) / lam
        if strict and np.any(dA <= 0.0):
            raise ValueError("Blade blockage closes the passage; reduce thickness or blade number.")
        return dA

    @cache
    def _lam_flow(self, z):
        # Zero-swirl (psi = 0) lead: fixed point of Q = omega * lam * int dA dr. It does not depend on the blade angles, so
        # it is cached per z: _solve_beta2 asks for it at the outlet for every trial b2
        q, omega = self.DP.Q.to('m^3/s').m, self.DP.N.to('rad/s').m
        r_h, r_t = self._r_hub(z), self._r_tip(z)
        r = np.linspace(r_h, r_t, 400)
        lam = q / (np.pi * omega * (r_t**2 - r_h**2))
        for _ in range(30):
            area = np.trapezoid(self._dA(r, lam, strict=False), r)
            if area <= 0.0:
                raise ValueError("Blade blockage closes the passage; reduce thickness or blade number.")
            lam = q / (omega * area)
        self._dA(r, lam)
        return lam

    def _lead(self, z, b2):
        gamma = self._beta(z, b2) - self.alpha_1.to('rad').m * (1.0 - z / self.L_ax.to('m').m)
        lam_sched = self._r_ms(z) * np.tan(gamma)  # eq. 33: lambda = r tan(beta)
        return max(lam_sched, self._lam_flow(z))  # lam below the zero-swirl lead would give negative swirl

    def _C(self, z, b2):
        # C of eq. 38 at station z (not the design constraints self.C)
        q, omega = self.DP.Q.to('m^3/s').m, self.DP.N.to('rad/s').m
        lam = self._lead(z, b2)
        r = np.linspace(self._r_hub(z), self._r_tip(z), 400)
        dA = self._dA(r, lam)
        # Own derivation: Q = omega*lam*(A0 - C*A1) from eqs. 38-39; SP-8052 eq. 43 is the closed form for a linearly tapered blade
        A0 = np.trapezoid(dA, r)
        A1 = np.trapezoid(dA / (r**2 + lam**2), r)
        return max((A0 - q / (omega * lam)) / A1, 0.0)

    def _psi(self, r, z, b2):
        # Eq. 38: psi_l = C/(r^2 + lambda^2), from radial equilibrium (eq. 34) with perfect guidance (eq. 35) and zero loss (eq. 36)
        return self._C(z, b2) / (np.asarray(r) ** 2 + self._lead(z, b2) ** 2)

    def _outlet_head(self, b2):
        l_ax, omega = self.L_ax.to('m').m, self.DP.N.to('rad/s').m
        r = np.linspace(self._r_hub(l_ax), self._r_tip(l_ax), 400)
        lam = self._lead(l_ax, b2)
        psi = self._psi(r, l_ax, b2)
        dA = self._dA(r, lam)
        c_u = psi * r * omega  # eq. 40
        c_a = (1.0 - psi) * lam * omega  # eq. 39
        # Mass-averaged Euler head g*H = omega * <r c_u> (eq. 36 with c_u = 0 at inlet, eq. 47)
        return omega * np.trapezoid(r * c_u * c_a * dA, r) / np.trapezoid(c_a * dA, r)

    def _solve_beta2(self):
        # SP-8052 eq. 47: g*dH_net = eta_bl * d(u*c_u), eta_bl = 0.85 (sec. 2.2.6)
        target_euler_head = (g * self.Delta_H).to('m^2/s^2').m / self.Geometry.eta_bl
        grid = np.linspace(self.beta_1.to('rad').m, 1.5, 600)
        f = np.array([self._outlet_head(b) - target_euler_head for b in grid])
        idx = np.flatnonzero((f[1:] >= 0) & (f[:-1] < 0))
        if idx.size == 0:
            raise ValueError("Target outlet head is not reachable with this geometry.")
        i = idx[0]
        return brentq(lambda b: self._outlet_head(b) - target_euler_head, grid[i], grid[i + 1])

    def _check(self, r, z):
        if not 0.0 <= z <= self.L_ax.to('m').m:
            raise ValueError("z must lie in [0, L_ax].")
        if np.any(np.asarray(r) < self._r_hub(z) - 1e-9) or np.any(np.asarray(r) > self._r_tip(z) + 1e-9):
            raise ValueError("r must lie in [r_hub(z), r_tip].")

    # ---- Radial equilibrium results: pint in, pint out; z is a scalar, r may be an array ----
    def GetHubRadius(self, z: Q_[float]) -> Q_[float]:
        return Q_(self._r_hub(z.to('m').m), 'm')

    def GetBladeAngle(self, z: Q_[float]) -> Q_[float]:
        """Blade angle at the rms radius, from the tangential direction (SP-8052 eq. 74)"""
        return Q_(np.degrees(self._beta(z.to('m').m, self.beta_2.to('rad').m)), 'deg')

    def GetSwirlRatio(self, r: Q_, z: Q_[float]) -> np.ndarray:
        """psi = c_u / (omega r) of SP-8052 eqs. 38-40, not the head coefficients psi_1 and psi_2"""
        r_m, z_m = r.to('m').m, float(z.to('m').m)
        self._check(r_m, z_m)
        return self._psi(r_m, z_m, self.beta_2.to('rad').m)

    def GetC_u(self, r: Q_, z: Q_[float]) -> Q_:
        """Swirl velocity c_u = psi omega r (SP-8052 eq. 40)"""
        r_m, z_m = r.to('m').m, float(z.to('m').m)
        self._check(r_m, z_m)
        return Q_(self._psi(r_m, z_m, self.beta_2.to('rad').m) * np.asarray(r_m) * self.DP.N.to('rad/s').m, 'm/s')

    def GetC_m(self, r: Q_, z: Q_[float]) -> Q_:
        """
        Axial velocity c_a = (1 - psi) lambda omega (SP-8052 eq. 39), the meridional velocity: the radial equilibrium
        neglects the radial velocity
        """
        r_m, z_m = r.to('m').m, float(z.to('m').m)
        self._check(r_m, z_m)
        b2 = self.beta_2.to('rad').m
        psi = self._psi(r_m, z_m, b2)
        return Q_((1.0 - psi) * self._lead(z_m, b2) * self.DP.N.to('rad/s').m, 'm/s')

    def MeridionalTraces(self, l: str = 'in', color: str = 'black') -> list[go.Scatter]:
        """
        Meridional cross-section as plotly traces in length unit l, z = 0 at the inducer inlet: hub and shroud
        contours, inlet and outlet, rms meanline (see BuildMeridian). Drawn by PlotMeridional, and with the impeller by
        Pump.PlotMeridional. Blades are not drawn yet; their leading and trailing edges would be added here, as in
        Impeller.MeridionalTraces.
        """
        return [
            self.Inlet.Trace(l, name='Inlet', line=dict(dash='dot', color='grey'), showlegend=False),
            self.Outlet.Trace(l, name='Outlet', line=dict(dash='dot', color='grey'), showlegend=False),
            self.ShroudContour.Trace(l, name='Shroud Contour', line=dict(color=color), showlegend=False),
            self.HubContour.Trace(l, name='Hub Contour', line=dict(color=color), showlegend=False),
            self.MeanContour.Trace(l, name='RMS Meanline', line=dict(color='black', dash='dashdot'))
        ]

    def PlotMeridional(self, l: str = 'in', color: str = 'black', name: str = 'Inducer') -> go.Figure:
        """Meridional cross-section (MeridionalTraces) and the axis, in its own figure."""
        fig = go.Figure(self.MeridionalTraces(l, color))

        fig.add_trace(go.Scatter(
            x=[0, self.L_ax.to(l).magnitude], y=[0, 0],
            mode='lines', name='Meridional axis', line=dict(dash='dashdot', color='grey')
        ))

        fig.update_xaxes(title_text=f"Meridional position ({l})")
        fig.update_yaxes(title_text=f"Radial position ({l})", scaleanchor='x', scaleratio=1)

        # legend in plotly's default place, outside the plot: the passage fills the top of the view
        fig.update_layout(
            title=f"{name} Meridional Cross-Section",
            template="plotly_white"
        )

        fig.show()

        return fig

    def PlotSwirl(self, l: str = 'in', name: str = 'Inducer') -> go.Figure:
        """Swirl velocity c_u (GetC_u) from hub to tip at z/L_ax = 0, 0.25, 0.5, 0.75 and 1, against radius in length unit l."""
        fig = go.Figure()

        for frac in (0.0, 0.25, 0.5, 0.75, 1.0):
            z = self.L_ax * frac
            r = Q_(np.linspace(self._r_hub(z.to('m').m), self._r_tip(z.to('m').m), 50), 'm')
            fig.add_trace(go.Scatter(
                x=r.to(l).magnitude, y=self.GetC_u(r, z).to('m/s').magnitude,
                mode='lines', name=f'z/L = {frac:.2f}'
            ))

        fig.update_xaxes(title_text=f"Radial position ({l})")
        fig.update_yaxes(title_text="Swirl velocity c_u (m/s)")

        fig.update_layout(
            title=f"{name} Swirl Velocity",
            template="plotly_white"
        )

        fig.show()

        return fig


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
    print_group('Design geometry', vars(inducer.Geometry).items())
    print_group('Sizing results', (
        (label, getattr(inducer, label)) for label in (
            'Delta_H', 'phi', 'D_tip', 'D_hub', 'N_ss', 'NPSH_3',
            'psi_2', 'R_hub_1', 'R_hub_2', 'L_ax', 'alpha_taper', 'beta_1', 'alpha_1', 'beta_2',
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
