# Radial equilibrium c_theta(r, z) for a tapered-hub inducer (NASA SP-8052 sec. 2.2.1, eqs. 33-41, 74).
import numpy as np
import matplotlib.pyplot as plt
from pint import Quantity as Q_
from scipy.optimize import brentq

from Inducer import Inducer, InducerDesignPoint, InducerInletDesignConstraints

_trapz = np.trapezoid


class InducerRadialEquilibrium:
    """Swirl and axial velocity fields behind a constant-tip, linearly tapered-hub inducer.

    Lead follows beta(z) = beta1 + (beta2 - beta1) (z/L)^2 at the rms radius; beta2 is solved so the
    mass-averaged outlet Euler head equals H_rise / eta_bl. Blade blockage uses N blades of constant
    normal thickness; blade_thickness=None disables it.
    incidence_ratio (alpha/beta, ~0.425) offsets the flow angle from the blade angle at the inlet,
    fading linearly to zero at the outlet; None disables it.

    eta_bl = 0.85  # Blade efficiency, sec. 2.2.6
    """

    def __init__(self, inducer: Inducer, eta_bl: float = 0.85, incidence_ratio: float | None = None,
                 blade_number: int = 3, blade_thickness: Q_ | None = None):
        self.inducer = inducer
        self.eta_bl = eta_bl
        self.incidence_ratio = incidence_ratio
        self.N_blades = blade_number
        self.t = None if blade_thickness is None else blade_thickness.to('m').magnitude

        self.R_t = inducer.D_tip.to('m').magnitude / 2.0
        self.R_h1 = inducer.R_hub_1.to('m').magnitude
        self.R_h2 = inducer.R_hub_2.to('m').magnitude
        self.L = inducer.L_ax.to('m').magnitude
        self.Q = inducer.DP.Q.to('m^3/s').magnitude
        self.omega = inducer.DP.N_shaft.to('rad/s').magnitude
        g_dh = 9.81 * inducer.DP.H_rise.to('m').magnitude
        self.target_euler_head = g_dh / eta_bl  # SP-8052 eq. 47: g*dH_net = eta_bl * d(u*c_u), eta_bl = 0.85 (sec. 2.2.6)

        # Zero-swirl inlet (c_theta = 0), uniform c_a = Q/A1 (eq. 23); flow lead lambda = c_a/omega (eq. 39 with psi=0)
        lam0 = self._lam_flow(0.0)
        gamma1 = np.arctan(lam0 / self._r_ms(0.0))
        if incidence_ratio is None:
            self.beta1, self.alpha1 = gamma1, 0.0
        else:
            self.beta1 = gamma1 / (1.0 - incidence_ratio)  # gamma = beta - alpha, alpha/beta from sec. 2.1.9 / 3.1.9
            self.alpha1 = incidence_ratio * self.beta1
        self.beta2 = self._solve_beta2()

    # ---- SI-float internals -------------------------------------------------
    def _r_hub(self, z):
        return self.R_h1 + (self.R_h2 - self.R_h1) * z / self.L

    def _r_ms(self, z):
        # rms = root-mean-square radius sqrt((r_tip^2 + r_hub^2)/2), the radius splitting the annulus into equal areas (Glossary)
        return np.sqrt((self.R_t**2 + self._r_hub(z) ** 2) / 2.0)

    def _beta(self, z, beta2):
        # Parabolic camber distribution, SP-8052 eq. 74 (specified at the rms radius, sec. 3.1.12)
        return self.beta1 + (beta2 - self.beta1) * (z / self.L) ** 2

    def _dA(self, r, lam, strict=True):
        # Open area per dr: 2 pi r minus N blades of tangential width t/sin(beta), with sin(beta) = lam/sqrt(r^2 + lam^2)
        dA = 2.0 * np.pi * r
        if self.t is not None:
            dA = dA - self.N_blades * self.t * np.sqrt(r**2 + lam**2) / lam
            if strict and np.any(dA <= 0.0):
                raise ValueError("Blade blockage closes the passage; reduce thickness or blade number.")
        return dA

    def _lam_flow(self, z):
        # Zero-swirl (psi = 0) lead: fixed point of Q = omega * lam * int dA dr
        r = np.linspace(self._r_hub(z), self.R_t, 400)
        lam = self.Q / (np.pi * self.omega * (self.R_t**2 - self._r_hub(z) ** 2))
        for _ in range(30):
            area = _trapz(self._dA(r, lam, strict=False), r)
            if area <= 0.0:
                raise ValueError("Blade blockage closes the passage; reduce thickness or blade number.")
            lam = self.Q / (self.omega * area)
        self._dA(r, lam)
        return lam

    def _lead(self, z, beta2):
        gamma = self._beta(z, beta2) - self.alpha1 * (1.0 - z / self.L)
        lam_sched = self._r_ms(z) * np.tan(gamma)  # eq. 33: lambda = r tan(beta)
        return max(lam_sched, self._lam_flow(z))  # lam below the zero-swirl lead would give negative swirl

    def _C(self, z, beta2):
        lam = self._lead(z, beta2)
        r = np.linspace(self._r_hub(z), self.R_t, 400)
        dA = self._dA(r, lam)
        # Own derivation: Q = omega*lam*(A0 - C*A1) from eqs. 38-39; SP-8052 eq. 43 is the closed form for a linearly tapered blade
        A0 = _trapz(dA, r)
        A1 = _trapz(dA / (r**2 + lam**2), r)
        return max((A0 - self.Q / (self.omega * lam)) / A1, 0.0)

    def _psi(self, r, z, beta2):
        # Eq. 38: psi_l = C/(r^2 + lambda^2), from radial equilibrium (eq. 34) with perfect guidance (eq. 35) and zero loss (eq. 36)
        return self._C(z, beta2) / (np.asarray(r) ** 2 + self._lead(z, beta2) ** 2)

    def _avg_angular_momentum(self, beta2):
        # Mass-averaged r c_u at the outlet = angular momentum flux / mass flux; conserved across a torque-free gap
        r = np.linspace(self._r_hub(self.L), self.R_t, 400)
        lam = self._lead(self.L, beta2)
        psi = self._psi(r, self.L, beta2)
        dA = self._dA(r, lam)
        c_u = psi * r * self.omega  # eq. 40
        c_a = (1.0 - psi) * lam * self.omega  # eq. 39
        return _trapz(r * c_u * c_a * dA, r) / _trapz(c_a * dA, r)
    
    def _outlet_head(self, beta2):
        # Mass-averaged Euler head g*H = omega * <r c_u> (eq. 36 with c_u = 0 at inlet, eq. 47)
        return self.omega * self._avg_angular_momentum(beta2)

    def _solve_beta2(self):
        grid = np.linspace(self.beta1, 1.5, 600)
        f = np.array([self._outlet_head(b) - self.target_euler_head for b in grid])
        idx = np.flatnonzero((f[1:] >= 0) & (f[:-1] < 0))
        if idx.size == 0:
            raise ValueError("Target outlet head is not reachable with this geometry.")
        i = idx[0]
        return brentq(lambda b: self._outlet_head(b) - self.target_euler_head, grid[i], grid[i + 1])

    def _check(self, r, z):
        if not 0.0 <= z <= self.L:
            raise ValueError("z must lie in [0, L_ax].")
        rh = self._r_hub(z)
        if np.any(np.asarray(r) < rh - 1e-9) or np.any(np.asarray(r) > self.R_t + 1e-9):
            raise ValueError("r must lie in [r_hub(z), r_tip].")

    # ---- Public API (pint in / pint out; z is scalar, r may be an array) -----
    def hub_radius(self, z: Q_) -> Q_:
        return Q_(self._r_hub(z.to('m').magnitude), 'm')

    def blade_angle(self, z: Q_) -> Q_:
        return Q_(np.degrees(self._beta(z.to('m').magnitude, self.beta2)), 'deg')

    def psi(self, r: Q_, z: Q_) -> np.ndarray:
        r_m, z_m = r.to('m').magnitude, z.to('m').magnitude
        self._check(r_m, z_m)
        return self._psi(r_m, z_m, self.beta2)

    def c_theta(self, r: Q_, z: Q_) -> Q_:
        r_m, z_m = r.to('m').magnitude, z.to('m').magnitude
        self._check(r_m, z_m)
        return Q_(self._psi(r_m, z_m, self.beta2) * np.asarray(r_m) * self.omega, 'm/s')

    def c_axial(self, r: Q_, z: Q_) -> Q_:
        r_m, z_m = r.to('m').magnitude, z.to('m').magnitude
        self._check(r_m, z_m)
        psi = self._psi(r_m, z_m, self.beta2)
        return Q_((1.0 - psi) * self._lead(z_m, self.beta2) * self.omega, 'm/s')

    def GetOutletAngularMomentum(self) -> Q_:
        """Mass-averaged r*c_theta at the outlet; equals g*H_rise / (eta_bl*omega) through the choice of beta2."""
        return Q_(self._avg_angular_momentum(self.beta2), 'm^2/s')


if __name__ == "__main__":
    g = Q_(9.81, 'm/s^2')
    import Main
    lox_pump = Main.CreatePump('LOX')
    ipa_pump = Main.CreatePump('IPA')
    lox_pump.AddInducer(Q_(65, 'm'))
    ipa_pump.AddInducer(Q_(80, 'm'))
    
    fig, axes = plt.subplots(2, 2, figsize=(12, 9), sharex='col', sharey='row')
    for col, inducer in zip(axes.T, (lox_pump.Inducer, ipa_pump.Inducer)):

        eq = InducerRadialEquilibrium(inducer, incidence_ratio=None, blade_number=3, blade_thickness=Q_(0.5, 'mm')) # Assuming 3 blades with 0.5 mm constant thickness
        print(f"beta1 = {np.degrees(eq.beta1):.2f} deg, beta2 = {np.degrees(eq.beta2):.2f} deg")
        for frac in (0.0, 0.25, 0.5, 0.75, 1.0):
            z = inducer.L_ax * frac
            r = np.linspace(eq.hub_radius(z).to('m').magnitude, eq.R_t, 50)
            col[0].plot(r * 1e3, eq.c_theta(Q_(r, 'm'), z).magnitude, label=f"z/L = {frac:.2f}")
            col[1].plot(r * 1e3, eq.c_axial(Q_(r, 'm'), z).magnitude, label=f"z/L = {frac:.2f}")
        col[0].set_title("Uhhh inducer?")
        col[1].set_xlabel("r [mm]")
        for ax in col:
            ax.grid(True, linestyle='--', alpha=0.5)
            ax.legend()
    axes[0, 0].set_ylabel("$c_\\theta$ [m/s]")
    axes[1, 0].set_ylabel("$c_a$ [m/s]")
    plt.tight_layout()
    plt.show()
