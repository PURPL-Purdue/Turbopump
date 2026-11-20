# Radial equilibrium c_theta(r, z) for a tapered-hub inducer (NASA SP-8052 sec. 2.2.1, eqs. 33-41, 74).
import sys
from pathlib import Path

import numpy as np
import matplotlib.pyplot as plt
import plotly.graph_objects as go
from pint import Quantity as Q_
from scipy.optimize import brentq
from scipy import integrate
from warnings import warn

if not __package__:
	# Run directly as a script: put Code/ first and drop this folder so the Inducer package (not Inducer.py) is imported
	_here = Path(__file__).resolve().parent
	sys.path[:] = [p for p in sys.path if Path(p or '.').resolve() != _here]
	sys.path.insert(0, str(_here.parent))

from Inducer.Inducer import Inducer
from Inducer.MeanlineTriangle import INCIDENCE, MeanlineTriangle
from VelocityTriangle import PlotTriangles, VelocityTriangle

"""
@Andrew Fontanetta:
- core functionality/radial equilibrium implementation
Modified by @Charlie Zhu:
- to expose more API and for CAD generation (InducerCAD.py)
- outlet slip and hydraulic efficiency (Gulich) implementation
"""

def _tangent_transition_spiral(r_0: float, r_max: float, wrap_angle: float, n: float = 2):
	"""
	Leading edge theta(r) from theta = 0 at r_0 to theta = wrap_angle at r_max, where it is tangent to the circle r_max (n > 1)
	"""
	# r = r_max - (r_max - r0) * (1 - theta / wrap)^n
	# dr/dtheta = n (r_max - r0) (1 - theta / wrap)^(n - 1) / wrap = 0 at theta = wrap

	# theta = wrap * (1 - ((r_max - r) / (r_max - r0))^(1/n))
	def func(r: float) -> float:
		return wrap_angle * (1 - np.maximum((r_max - r) / (r_max - r_0), 0)**(1 / n))

	return func

class InducerRadialEquilibrium:
	"""Swirl and axial velocity fields behind a constant-tip, linearly tapered-hub inducer.

	Lead follows beta(z) = beta1b + (beta2b - beta1b) (z/L)^2 at the rms radius; beta2b is solved so the
	mass-averaged outlet Euler head equals H_rise / eta_bl, with the flow leaving at beta2b - delta2 (Gülich deviation,
	see _deviation); no deviation would assume SP-8052's perfect guidance. Blade blockage uses N blades of constant
	normal thickness at the tip, tapering linearly to hub_thickness_ratio times that at the hub (sec. 3.1.11,
	eq. 42; SP-8052 recommends thicker at the hub); blade_thickness=None disables it.
	Incidence is an angle quantity, following MeanlineTriangle: gamma1 = beta1b - incidence.
	It defaults to MeanlineTriangle.py's INCIDENCE value and blends linearly into the deviation at the outlet.

	eta_bl = 1 - 0.11 L/t  # Hydraulic efficiency, Gülich eq. 7.6.13, with the tip solidity L/t (solidity), see _eta_hyd
	"""

	def __init__(self, inducer: Inducer, incidence: Q_ | None = None, blade_thickness: Q_=Q_(1, 'mm'),
		hub_thickness_ratio: float = 1.0, meanline: MeanlineTriangle | None = None):
		self.inducer = inducer

		self.sweepback = np.deg2rad(60)

		self.hub_thickness_ratio = hub_thickness_ratio
		if meanline is not None and incidence is not None:
			raise ValueError("Give either incidence or meanline, not both.")
		if meanline is not None:
			self.meanline = meanline
		else:
			if incidence is None:
				incidence = Q_(INCIDENCE, 'deg')
			self.meanline = MeanlineTriangle(inducer, incidence=incidence)
		self.incidence = self.meanline.incidence
		self.Z_blade = self.inducer.Z_blade
		self.blade_thick: float = blade_thickness.to('m').m

		self.R_t = inducer.R_tip.to('m').magnitude
		self.R_h1 = inducer.R_hub_1.to('m').magnitude
		self.R_h2 = inducer.R_hub_2.to('m').magnitude
		self.L = inducer.L_ax.to('m').magnitude
		self.Q = inducer.DP.Q.to('m^3/s').magnitude
		self.omega = inducer.DP.N_shaft.to('rad/s').magnitude
		self._g_dh = 9.81 * inducer.DP.H_rise.to('m').magnitude  # SP-8052 eq. 47: g*dH_net = eta_bl * d(u*c_u)

		# Flow angle from the flow lead at z = 0, including blade blockage (eq. 39 with psi = 0).
		gamma1 = np.arctan(self._lam_flow(0.0) / self.r_ms(0.0))
		self.beta1b = gamma1 + self.incidence
		self.beta2b = self._solve_beta2b()
		self.delta2 = self._deviation(self.beta2b)		# outlet deviation at the rms radius, the flow leaves at beta2b - delta2
		self.solidity = self._solidity(self.beta2b)		# blade length over pitch at the tip
		self.eta_hyd = self._eta_hyd(self.beta2b)

	# ---- SI-float geometry (public, used by InducerCAD) ----------------------
	def r_hub(self, z):
		return self.R_h1 + (self.R_h2 - self.R_h1) * z / self.L

	def r_ms(self, z):
		# rms = root-mean-square radius sqrt((r_tip^2 + r_hub^2)/2), the radius splitting the annulus into equal areas (Glossary)
		return np.sqrt((self.R_t**2 + self.r_hub(z) ** 2) / 2.0)

	def beta(self, z, beta2b):
		# Parabolic camber distribution, SP-8052 eq. 74 (specified at the rms radius, sec. 3.1.12)
		return self.beta1b + (beta2b - self.beta1b) * (z / self.L) ** 2

	# ---- SI-float internals -------------------------------------------------
	def _deviation(self, beta2b):
		# Outlet deviation on the rms line (Gülich, inducer design table), angles in degrees:
		# delta_2 = (2 + (beta_2B - beta_1B) / 3) (t / L)^(1/3), with the tip solidity L/t
		return np.radians((2 + np.degrees(beta2b - self.beta1b) / 3) / self._solidity(beta2b)**(1 / 3))

	def _solidity(self, beta2b):
		# Solidity L/t at the tip (SP-8052 sec. 2.1.15: chord over blade spacing): blade length from the swept leading edge,
		# which reaches the tip at theta = sweepback, to the outlet, over the tip pitch 2 pi R_t / Z
		z = np.linspace(0, self.L, 200)
		lead = self.r_ms(z) * np.tan(self.beta(z, beta2b))		# the same lead at every radius, as in InducerCAD._GetBladeLines
		theta = integrate.cumulative_trapezoid(1 / lead, z, initial=0)
		arclength = integrate.cumulative_trapezoid(np.sqrt(1 + (self.R_t / lead)**2), z, initial=0)		# dl^2 = dz^2 + (R_t dtheta)^2
		length = arclength[-1] - np.interp(self.sweepback, theta, arclength)
		return length / (2 * np.pi * self.R_t / self.Z_blade)

	def _eta_hyd(self, beta2b):
		# Hydraulic efficiency, Gülich eq. 7.6.13: eta = 1 - 0.11 L/t, with the tip solidity L/t
		return 1 - 0.11 * self._solidity(beta2b)

	def _dA(self, r, lam, strict=True):
		# Open area per dr: 2 pi r minus N blades of tangential width t/sin(beta), with sin(beta) = lam/sqrt(r^2 + lam^2)
		dA = 2.0 * np.pi * r
		if self.blade_thick is not None:# Assumed constant thickness for now
			dA = dA - self.Z_blade * self.blade_thick * np.sqrt(r**2 + lam**2) / lam
			if strict and np.any(dA <= 0.0):
				raise ValueError("Blade blockage closes the passage; reduce thickness or blade number.")
		return dA

	def _lam_flow(self, z):
		# Zero-swirl (psi = 0) lead: fixed point of Q = omega * lam * int dA dr
		r = np.linspace(self.r_hub(z), self.R_t, 400)
		lam = self.Q / (np.pi * self.omega * (self.R_t**2 - self.r_hub(z) ** 2))
		for _ in range(30):
			area = np.trapezoid(self._dA(r, lam, strict=False), r)
			if area <= 0.0:
				raise ValueError("Blade blockage closes the passage; reduce thickness or blade number.")
			lam = self.Q / (self.omega * area)
		self._dA(r, lam)
		return lam

	def _lead(self, z, beta2b):
		# Flow angle: the blade angle less the incidence at the inlet and the deviation at the outlet, linearly in between (own choice)
		gamma = self.beta(z, beta2b) - self.incidence * (1.0 - z / self.L) - self._deviation(beta2b) * z / self.L
		lam_sched = self.r_ms(z) * np.tan(gamma)  # eq. 33: lambda = r tan(beta)
		return max(lam_sched, self._lam_flow(z))  # lam below the zero-swirl lead would give negative swirl

	def _C(self, z, beta2b):
		lam = self._lead(z, beta2b)
		r = np.linspace(self.r_hub(z), self.R_t, 400)
		dA = self._dA(r, lam)
		# Own derivation: Q = omega*lam*(A0 - C*A1) from eqs. 38-39; SP-8052 eq. 43 is the closed form for a linearly tapered blade
		A0 = np.trapezoid(dA, r)
		A1 = np.trapezoid(dA / (r**2 + lam**2), r)
		return max((A0 - self.Q / (self.omega * lam)) / A1, 0.0)

	def _psi(self, r, z, beta2b):
		# Eq. 38: psi_l = C/(r^2 + lambda^2), from radial equilibrium (eq. 34) with guidance along the flow lead of _lead (eq. 35) and zero loss (eq. 36)
		return self._C(z, beta2b) / (np.asarray(r) ** 2 + self._lead(z, beta2b) ** 2)

	def _avg_angular_momentum(self, beta2b):
		# Mass-averaged r c_u at the outlet = angular momentum flux / mass flux; conserved across a torque-free gap
		r = np.linspace(self.r_hub(self.L), self.R_t, 400)
		lam = self._lead(self.L, beta2b)
		psi = self._psi(r, self.L, beta2b)
		dA = self._dA(r, lam)
		c_u = psi * r * self.omega  # eq. 40
		c_a = (1.0 - psi) * lam * self.omega  # eq. 39
		return np.trapezoid(r * c_u * c_a * dA, r) / np.trapezoid(c_a * dA, r)
	
	def _outlet_head(self, beta2b):
		# Mass-averaged Euler head g*H = omega * <r c_u> (eq. 36 with c_u = 0 at inlet, eq. 47)
		return self.omega * self._avg_angular_momentum(beta2b)

	def _solve_beta2b(self):
		grid = np.linspace(self.beta1b, 1.5, 600)
		f = np.array([self._outlet_head(b) - self._g_dh / self._eta_hyd(b) for b in grid])
		idx = np.flatnonzero((f[1:] >= 0) & (f[:-1] < 0))
		if idx.size == 0:
			raise ValueError("Target outlet head is not reachable with this geometry.")
		i = idx[0]
		return brentq(lambda b: self._outlet_head(b) - self._g_dh / self._eta_hyd(b), grid[i], grid[i + 1])

	def _check(self, r, z):
		if not 0.0 <= z <= self.L:
			raise ValueError("z must lie in [0, L_ax].")
		rh = self.r_hub(z)
		if np.any(np.asarray(r) < rh - 1e-9) or np.any(np.asarray(r) > self.R_t + 1e-9):
			raise ValueError("r must lie in [r_hub(z), r_tip].")

	# ---- Public API (pint in / pint out, except the SI floats of hub_radius, mean_radius and blade_angle; z is scalar, r may be an array) -----
	def hub_radius(self, z: Q_) -> float:
		"""Hub radius in m."""
		return self.r_hub(z.to('m').magnitude)

	def mean_radius(self, z: Q_) -> float:
		"""Root-mean-square radius in m, which splits the annulus into equal areas."""
		return self.r_ms(z.to('m').magnitude)

	def blade_angle(self, z: Q_) -> float:
		"""Blade angle at the rms radius in rad, from tangential."""
		return self.beta(z.to('m').magnitude, self.beta2b)

	def psi(self, r: Q_, z: Q_) -> np.ndarray:
		r_m, z_m = r.to('m').magnitude, z.to('m').magnitude
		self._check(r_m, z_m)
		return self._psi(r_m, z_m, self.beta2b)

	def c_theta(self, r: Q_, z: Q_) -> Q_:
		r_m, z_m = r.to('m').magnitude, z.to('m').magnitude
		self._check(r_m, z_m)
		return Q_(self._psi(r_m, z_m, self.beta2b) * np.asarray(r_m) * self.omega, 'm/s')

	def c_axial(self, r: Q_, z: Q_, blockage: bool=True) -> Q_:
		r_m, z_m = r.to('m').magnitude, z.to('m').magnitude
		self._check(r_m, z_m)
		psi = self._psi(r_m, z_m, self.beta2b)
		blocked_c_axial = Q_((1.0 - psi) * self._lead(z_m, self.beta2b) * self.omega, 'm/s')

		if blockage:
			return blocked_c_axial
		
		lam = self._lead(z_m, self.beta2b)
		tau = self._dA(r_m, lam) / (2 * np.pi * r_m)    # open-area fraction, 1 = no blockage
		return blocked_c_axial * tau        # circumferentially averaged c_axial

	def GetOutletAngularMomentum(self) -> Q_:
		"""Mass-averaged r*c_theta at the outlet; equals g*H_rise / (eta_bl*omega) through the choice of beta2b."""
		return Q_(self._avg_angular_momentum(self.beta2b), 'm^2/s')

	def MeanlineTriangle(self, z: Q_, notation: str = '') -> VelocityTriangle:
		"""Absolute-frame velocity triangle on the rms-radius meanline at axial position z."""
		r_mean = Q_(self.mean_radius(z), 'm')
		return VelocityTriangle(
			u=Q_(self.omega, 'rad/s') * r_mean,
			c_m=self.c_axial(r_mean, z),
			c_u=self.c_theta(r_mean, z),
			notation=notation,
		)

	def TipTriangle(self, z: Q_, notation: str = '') -> VelocityTriangle:
		"""Absolute-frame velocity triangle at the tip radius at axial position z."""
		R_t = Q_(self.R_t, 'm')
		return VelocityTriangle(
			u=Q_(self.omega, 'rad/s') * R_t,
			c_m=self.c_axial(R_t, z),
			c_u=self.c_theta(R_t, z),
			notation=notation,
		)

	def RelativeVelocity(self, r: Q_, z: Q_) -> tuple[Q_, Q_, Q_]:
		"""Relative velocity W = V - U at (r, z) from the converged radial equilibrium: (W_u, W_m, |W|)."""
		U = Q_(self.omega, 'rad/s') * r.to('m')
		W_u = U - self.c_theta(r, z)
		W_m = self.c_axial(r, z)
		return W_u, W_m, np.sqrt(W_u**2 + W_m**2)

	def FlowAngle(self, r: Q_, z: Q_) -> Q_:
		"""Relative flow angle beta from the tangential direction: tan(beta) = c_m / (U - c_u)."""
		W_u, W_m, _ = self.RelativeVelocity(r, z)
		return Q_(np.degrees(np.arctan2(W_m.magnitude, W_u.to(W_m.units).magnitude)), 'deg')

	def BladeAngle(self, r: Q_, z: Q_) -> Q_:
		"""Blade angle from the parabolic camber (eq. 74, rms radius) with constant lead, r tan(beta) = lambda_blade (eq. 33)."""
		z_m = z.to('m').magnitude
		lam_blade = self.r_ms(z_m) * np.tan(self.beta(z_m, self.beta2b))
		return Q_(np.degrees(np.arctan(lam_blade / r.to('m').magnitude)), 'deg')

	def PlotBetaDistribution(self, z_fracs=(0.0, 0.25, 0.5, 0.75, 1.0), n: int = 50, l: str = 'mm',
							 title: str = 'Inducer relative flow and blade angle, hub to tip') -> go.Figure:
		"""beta against radius at several axial stations; solid = relative flow angle, dashed = blade angle."""
		fig = go.Figure()
		colors = [f'hsl({int(240 * (1 - f))},70%,45%)' for f in z_fracs]
		for frac, color in zip(z_fracs, colors):
			z = Q_(self.L * frac, 'm')
			r = Q_(np.linspace(self.r_hub(z.magnitude), self.R_t, n), 'm')
			x = r.to(l).magnitude
			fig.add_trace(go.Scatter(x=x, y=self.FlowAngle(r, z).magnitude, mode='lines', line=dict(color=color),
									 name=f'z/L = {frac:.2f} flow', legendgroup=f'{frac}'))
			fig.add_trace(go.Scatter(x=x, y=self.BladeAngle(r, z).magnitude, mode='lines', line=dict(color=color, dash='dash'),
									 name=f'z/L = {frac:.2f} blade', legendgroup=f'{frac}'))
		fig.update_layout(title=title, xaxis_title=f'r [{l}]', yaxis_title='β from tangential [deg]', template='plotly_white')
		return fig

	def PlotMeanlineTriangles(self, title: str = 'Inducer meanline velocity triangles') -> go.Figure:
		"""Inlet (z = 0) and outlet (z = L_ax) meanline triangles overlaid."""
		triangles = [
			self.MeanlineTriangle(Q_(0, 'm'), '<sub>1</sub>'),
			self.MeanlineTriangle(Q_(self.L, 'm'), '<sub>2</sub>'),
		]
		return PlotTriangles(triangles, names=['Inlet', 'Outlet'], title=title)


if __name__ == "__main__":
	g = Q_(9.81, 'm/s^2')
	import Main

	lox_pump = Main.CreatePump(Main.Propellant.LOX)
	ipa_pump = Main.CreatePump(Main.Propellant.IPA)
	lox_pump.AddInducer(Q_(65, 'm'))
	ipa_pump.AddInducer(Q_(80, 'm'))
	
	fig, axes = plt.subplots(2, 2, figsize=(12, 9), sharex='col', sharey='row')
	print("--- Meanline Angles after Radial equilibrium ---")
	for col, propellant, inducer in zip(axes.T, ('LOX', 'IPA'), (lox_pump.Inducer, ipa_pump.Inducer)):

		eq = InducerRadialEquilibrium(inducer, blade_thickness=Q_(1, 'mm'))
		print(f"{propellant}: beta1b = {np.degrees(eq.beta1b):.2f} deg, beta2b = {np.degrees(eq.beta2b):.2f} deg, delta2 = {np.degrees(eq.delta2):.2f} deg")
		for frac in (0.0, 0.25, 0.5, 0.75, 1.0):
			z = inducer.L_ax * frac
			r = np.linspace(eq.hub_radius(z), eq.R_t, 50)
			col[0].plot(r * 1e3, eq.c_theta(Q_(r, 'm'), z).magnitude, label=f"z/L = {frac:.2f}")
			col[1].plot(r * 1e3, eq.c_axial(Q_(r, 'm'), z).magnitude, label=f"z/L = {frac:.2f}")
		col[0].set_title(f"{propellant} inducer")
		col[1].set_xlabel("r [mm]")
		for ax in col:
			ax.grid(True, linestyle='--', alpha=0.5)
			ax.legend()

		z_in = Q_(0, 'm')
		r_in = Q_(np.array([eq.R_h1, eq.r_ms(0.0), eq.R_t]), 'm')
		print(f"{propellant} inlet beta (hub, rms, tip) = {np.round(eq.FlowAngle(r_in, z_in).magnitude, 2)} deg")
		eq.PlotBetaDistribution(title=f"{propellant} inducer relative flow and blade angle").show()
		PlotMeanlineTriangles = eq.PlotMeanlineTriangles(
				f"{propellant} inducer meanline velocity triangles"
			).show()
		# InducerCAD(eq).ShowCAD()  # blocks until the viewer is closed

	axes[0, 0].set_ylabel("$c_\\theta$ [m/s]")
	axes[1, 0].set_ylabel("$c_a$ [m/s]")
	plt.tight_layout()
	plt.show()
