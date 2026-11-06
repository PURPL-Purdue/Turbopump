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

# temporary, will move later
import cadquery.func
from cadquery.vis import show

if not __package__:
	# Run directly as a script: put Code/ first and drop this folder so the Inducer package (not Inducer.py) is imported
	_here = Path(__file__).resolve().parent
	sys.path[:] = [p for p in sys.path if Path(p or '.').resolve() != _here]
	sys.path.insert(0, str(_here.parent))

from Inducer.Inducer import Inducer
from Inducer.MeanlineTriangle import MeanlineTriangle
from VelocityTriangle import PlotTriangles, VelocityTriangle

"""
@Andrew Fontanetta:
- core functionality/radial equilibrium implementation
Modified by @Charlie Zhu:
- to expose more API and for CAD generation
- outlet slip and hydraulic efficiency (Gulich) implementation
"""

def _archimedian_spiral(r_0: float, r_max: float, wrap_angle: float) -> float:
	"""
	@Charlie Zhu
	geninuely don't know what i did here
	use Desmos for visualization
	"""
	# r = a*theta - r0

	# r = r0 = a * theta1 - r0
	# r = r_max = a * theta2 - r0
	# wrap = theta2 - theta1

	# theta2 * a = r_max + r0
	# theta1 * a = 2 * r0
	# a * (theta2 - theta1) = a * wrap = r_max - r0

	a = (r_max - r_0) / wrap_angle
	theta_1 = 2 * r_0 / a

	def func(r: float) -> float:
		return (r + r_0) / a - theta_1

	return func

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

def _leading_edge_theta(r_0, r_max, wrap_angle, p=2.0):
	"""theta = wrap_angle * t^p, t = (r - r_0) / (r_max - r_0). Radial at r_0 for p > 1.
	Sweep from radial at r_max: tan(sweep) = p * wrap_angle * r_max / (r_max - r_0)."""
	def func(r):
		t = np.clip((r - r_0) / (r_max - r_0), 0, 1)
		return wrap_angle * t**p
	return func

_leading_edge_theta = _archimedian_spiral

class InducerRadialEquilibrium:
	"""Swirl and axial velocity fields behind a constant-tip, linearly tapered-hub inducer.

	Lead follows beta(z) = beta1 + (beta2b - beta1) (z/L)^2 at the rms radius; beta2b is solved so the
	mass-averaged outlet Euler head equals H_rise / eta_bl, with the flow leaving at beta2b - delta2 (Gülich deviation,
	see _deviation); no deviation would assume SP-8052's perfect guidance. Blade blockage uses N blades of constant
	normal thickness at the tip, tapering linearly to hub_thickness_ratio times that at the hub (sec. 3.1.11,
	eq. 42; SP-8052 recommends thicker at the hub); blade_thickness=None disables it.
	incidence_ratio (alpha/beta, ~0.425) offsets the flow angle from the blade angle at the inlet,
	blending linearly into the deviation at the outlet; None disables it.

	eta_bl = 1 - 0.11 L/t  # Hydraulic efficiency, Gülich eq. 7.6.13, see _eta_bl
	"""

	def __init__(self, inducer: Inducer, incidence_ratio: float | None = None, blade_thickness: Q_ | None = None,
				 hub_thickness_ratio: float = 1.0, meanline: MeanlineTriangle | None = None):
		self.inducer = inducer

		self._sweepback = np.deg2rad(60)

		self.hub_thickness_ratio = hub_thickness_ratio
		# The blade inlet angle comes from the zero-blockage inlet meanline triangle
		self.meanline = meanline if meanline is not None else MeanlineTriangle(inducer, incidence_ratio)
		self.incidence_ratio = self.meanline.incidence_ratio
		self.Z_blade = self.inducer.Z_blade
		self.t: float = None if blade_thickness is None else blade_thickness.to('m').magnitude

		self.R_t = inducer.R_tip.to('m').magnitude
		self.R_h1 = inducer.R_hub_1.to('m').magnitude
		self.R_h2 = inducer.R_hub_2.to('m').magnitude
		self.L = inducer.L_ax.to('m').magnitude
		self.Q = inducer.DP.Q.to('m^3/s').magnitude
		self.omega = inducer.DP.N_shaft.to('rad/s').magnitude
		self._g_dh = 9.81 * inducer.DP.H_rise.to('m').magnitude  # SP-8052 eq. 47: g*dH_net = eta_bl * d(u*c_u)

		# Blade inlet angle from the flow lead at z = 0, including blade blockage (eq. 39 with psi = 0);
		# the meanline triangle is the zero-blockage version and only supplies the incidence ratio
		gamma1 = np.arctan(self._lam_flow(0.0) / self._r_ms(0.0))
		if self.incidence_ratio is None:
			self.beta1, self.alpha1 = gamma1, 0.0
		else:
			self.beta1 = gamma1 / (1.0 - self.incidence_ratio)
			self.alpha1 = self.incidence_ratio * self.beta1
		self.beta2b = self._solve_beta2b()
		self.delta2 = self._deviation(self.beta2b)		# outlet deviation at the rms radius, the flow leaves at beta2b - delta2
		self.eta_hyd = self._eta_hyd(self.beta2b)

	# ---- SI-float internals -------------------------------------------------
	def _r_hub(self, z):
		return self.R_h1 + (self.R_h2 - self.R_h1) * z / self.L

	def _r_ms(self, z):
		# rms = root-mean-square radius sqrt((r_tip^2 + r_hub^2)/2), the radius splitting the annulus into equal areas (Glossary)
		return np.sqrt((self.R_t**2 + self._r_hub(z) ** 2) / 2.0)

	def _beta(self, z, beta2b):
		# Parabolic camber distribution, SP-8052 eq. 74 (specified at the rms radius, sec. 3.1.12)
		return self.beta1 + (beta2b - self.beta1) * (z / self.L) ** 2

	def _deviation(self, beta2b):
		# Outlet deviation on the rms line (Gülich, inducer design table), angles in degrees:
		# delta_2 = (2 + (beta_2B - beta_1B) / 3) (t / L)^(1/3), with the outlet pitch t and the blade length L from the swept leading edge
		z = np.linspace(0, self.L, 200)
		r = self._r_ms(z)
		tan_beta = np.tan(self._beta(z, beta2b))
		theta = integrate.cumulative_trapezoid(1 / (r * tan_beta), z, initial=0)		# wrap, as in GetBladeLines
		get_LE_theta = _leading_edge_theta(self.R_h1, self.R_t, wrap_angle=self._sweepback)
		z_LE = np.interp(0, theta - get_LE_theta(r), z)
		# blade arclength: dl^2 = dz^2 + dr^2 + (r dtheta)^2, with r dtheta = dz / tan(beta)
		arclength = integrate.cumulative_trapezoid(np.sqrt(1 + np.gradient(r, z)**2 + 1 / tan_beta**2), z, initial=0)
		length = arclength[-1] - np.interp(z_LE, z, arclength)
		pitch = 2 * np.pi * r[-1] / self.Z_blade
		return np.radians((2 + np.degrees(beta2b - self.beta1) / 3) * (pitch / length)**(1 / 3))

	def _eta_hyd(self, beta2b):
		# Hydraulic efficiency, Gülich eq. 7.6.13: eta = 1 - 0.11 L/t. Initial estimate: L = r * wrap and t = 2 pi r / Z, so r cancels; no leading edge sweep
		z = np.linspace(0, self.L, 200)
		wrap = np.trapezoid(1 / (self._r_ms(z) * np.tan(self._beta(z, beta2b))), z)		# as in GetBladeLines
		return 1 - 0.11 * self.Z_blade * wrap / (2 * np.pi)

	def _dA(self, r, lam, strict=True):
		# Open area per dr: 2 pi r minus N blades of tangential width t/sin(beta), with sin(beta) = lam/sqrt(r^2 + lam^2)
		dA = 2.0 * np.pi * r
		if self.t is not None:# Assumed constant thickness for now
			dA = dA - self.Z_blade * self.t * np.sqrt(r**2 + lam**2) / lam
			if strict and np.any(dA <= 0.0):
				raise ValueError("Blade blockage closes the passage; reduce thickness or blade number.")
		return dA

	def _lam_flow(self, z):
		# Zero-swirl (psi = 0) lead: fixed point of Q = omega * lam * int dA dr
		r = np.linspace(self._r_hub(z), self.R_t, 400)
		lam = self.Q / (np.pi * self.omega * (self.R_t**2 - self._r_hub(z) ** 2))
		for _ in range(30):
			area = np.trapezoid(self._dA(r, lam, strict=False), r)
			if area <= 0.0:
				raise ValueError("Blade blockage closes the passage; reduce thickness or blade number.")
			lam = self.Q / (self.omega * area)
		self._dA(r, lam)
		return lam

	def _lead(self, z, beta2b):
		# Flow angle: the blade angle less the incidence at the inlet and the deviation at the outlet, linearly in between (own choice)
		gamma = self._beta(z, beta2b) - self.alpha1 * (1.0 - z / self.L) - self._deviation(beta2b) * z / self.L
		lam_sched = self._r_ms(z) * np.tan(gamma)  # eq. 33: lambda = r tan(beta)
		return max(lam_sched, self._lam_flow(z))  # lam below the zero-swirl lead would give negative swirl

	def _C(self, z, beta2b):
		lam = self._lead(z, beta2b)
		r = np.linspace(self._r_hub(z), self.R_t, 400)
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
		r = np.linspace(self._r_hub(self.L), self.R_t, 400)
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
		grid = np.linspace(self.beta1, 1.5, 600)
		f = np.array([self._outlet_head(b) - self._g_dh / self._eta_hyd(b) for b in grid])
		idx = np.flatnonzero((f[1:] >= 0) & (f[:-1] < 0))
		if idx.size == 0:
			raise ValueError("Target outlet head is not reachable with this geometry.")
		i = idx[0]
		return brentq(lambda b: self._outlet_head(b) - self._g_dh / self._eta_hyd(b), grid[i], grid[i + 1])

	def _check(self, r, z):
		if not 0.0 <= z <= self.L:
			raise ValueError("z must lie in [0, L_ax].")
		rh = self._r_hub(z)
		if np.any(np.asarray(r) < rh - 1e-9) or np.any(np.asarray(r) > self.R_t + 1e-9):
			raise ValueError("r must lie in [r_hub(z), r_tip].")

	# ---- Public API (pint in / pint out; z is scalar, r may be an array) -----
	def hub_radius(self, z: Q_) -> Q_:
		return Q_(self._r_hub(z.to('m').magnitude), 'm')

	def mean_radius(self, z: Q_) -> Q_:
		"""Root-mean-square radius, which splits the annulus into equal areas."""
		return Q_(self._r_ms(z.to('m').magnitude), 'm')

	def blade_angle(self, z: Q_) -> Q_:
		return Q_(np.degrees(self._beta(z.to('m').magnitude, self.beta2b)), 'deg')

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
		r_mean = self.mean_radius(z)
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

	def PlotMeanlineTriangles(self, title: str = 'Inducer meanline velocity triangles') -> go.Figure:
		"""Inlet (z = 0) and outlet (z = L_ax) meanline triangles overlaid."""
		triangles = [
			self.MeanlineTriangle(Q_(0, 'm'), '<sub>1</sub>'),
			self.MeanlineTriangle(Q_(self.L, 'm'), '<sub>2</sub>'),
		]
		return PlotTriangles(triangles, names=['Inlet', 'Outlet'], title=title)

	def GetBladeLines(self, num_lines: int=7, resolution: int = 200) -> list[np.ndarray]:
		"""
		@Charlie zhu
		Blade lines from hub to tip, a list of (resolution, 3) arrays of (z, r, theta) in m, m and rad, from the inlet.
		"""
		z_ref = np.linspace(0, self.L, resolution)

		r_hub: np.ndarray[float] = self._r_hub(z_ref)
		r_tip: np.ndarray[float] = np.full(resolution, self.R_t)
		r_diff: np.ndarray[float] = r_tip - r_hub

		points_zrt = []
		
		radial_lead = self._r_ms(z_ref) * np.tan(self._beta(z_ref, self.beta2b))
		theta_ref = integrate.cumulative_trapezoid(y=1 / radial_lead, x=z_ref, initial=0)

		# radii clustered toward the tip
		radii_LE = self.R_h1 + (self.R_t - self.R_h1) * np.sin(np.linspace(0, np.pi / 2, resolution))

		get_LE_theta = _leading_edge_theta(self.R_h1, self.R_t, wrap_angle=self._sweepback)
		theta_LE = get_LE_theta(radii_LE)
		assert theta_LE[-1] <= theta_ref[-1], 'The leading edge wraps further than the blade'
		z_LE_ref = np.interp(theta_LE, theta_ref, z_ref)		# theta_ref increases with z_ref
		
		hub_taper_slope = (self.R_h2 - self.R_h1) / self.L
		z_LE_canted = z_LE_ref - (radii_LE - self._r_hub(z_LE_ref)) * hub_taper_slope

		arclength_LE = np.concatenate(([0], np.cumsum(np.sqrt(np.diff(z_LE_canted)**2 + np.diff(radii_LE)**2 + (radii_LE[1:] * np.diff(theta_LE))**2))))
		s = np.linspace(0, arclength_LE[-1], num_lines)
		# keep z and r values at the arclengths specified s (linearly interpolated)
		z_LE_ref, radii_LE = np.interp(s, arclength_LE, z_LE_ref), np.interp(s, arclength_LE, radii_LE)

		LE_span_fractions: np.ndarray[float] = (radii_LE - self._r_hub(z_LE_ref)) / (self.R_t - self._r_hub(z_LE_ref))

		for alpha in LE_span_fractions:
			
			r_offset = alpha * r_diff
			r = r_hub + r_offset

			tan_beta = radial_lead/r # more necessary than you think lol
			
			theta = integrate.cumulative_trapezoid(
				y=1 / (r * tan_beta), x=z_ref, initial=0
				)
			# adjust z such that meridional sweep is perpendicular to hub
			z = z_ref - r_offset * hub_taper_slope

			points = np.column_stack((z, r, theta)) # (resolution, 3)

			# Angle past the leading edge curve, at each point's own radius: the blade line and the leading edge,
			# projected on the polar (r, theta) plane, cross where it is zero
			past_LE = theta - get_LE_theta(r)
			i = np.argmax(past_LE > 0)		# first point behind the leading edge
			# linear interpolation to the crossing
			t = past_LE[i-1] / (past_LE[i-1] - past_LE[i])
			leading_edge = points[i-1] + t * (points[i] - points[i-1])

			points_zrt.append(np.vstack((leading_edge, points[i:])))

		return points_zrt

	def CreateCAD(self, num_lines: int=7, resolution: int=300, flip_z: bool=True, flip_theta: bool=True) -> cadquery.Assembly:
		"""
		@Charlie Zhu
		Temporarily borrowing from impeller CAD API
		Blade lines from hub to tip as cadquery splines (x = r cos theta, y = r sin theta, z, in mm) and the revolved hub, in the cadquery viewer."""
		from Helpers import CAD
		blade_lines = self.GetBladeLines(num_lines=num_lines, resolution=resolution)
		blade_edges = [cadquery.Edge.makeSpline(CAD.CylindricalToCartesian(p, flip_z, flip_theta)) for p in blade_lines]
		camber = cadquery.func.loft(blade_edges)
		blade_locations = [cadquery.Location((0, 0, 0), (0, 0, 1),  360 * n / self.Z_blade) for n in range(self.Z_blade)]
		blades = camber.moved(*blade_locations)

		# Linearly tapered hub: a line from (z, r) = (0, R_h1) to (L, R_h2) at theta = 0, revolved about the z axis
		hub_edge = cadquery.Edge.makeLine(*CAD.CylindricalToCartesian(np.array([[0, self.R_h1, 0], [self.L, self.R_h2, 0]]), flip_z, flip_theta))
		hub_surface = cadquery.func.revolve(hub_edge, (0, 0, 0), (0, 0, 1))
		assm = (cadquery.Assembly()
			.add(hub_surface, name='hub', color=cadquery.Color(0.6, 0.6, 0.6))
			.add(blades, name='blades', color=cadquery.Color(0.2, 0.3, 0.7))
			#.add(blade_edges, name='blade_edges', color=cadquery.Color(0.9, 0.1, 0.1))
		)
		self._cad = [assm, blade_edges]
		return assm

	def ShowCAD(self):
		if self._cad is None: self.CreateCAD()
		show(*self._cad)


if __name__ == "__main__":
	g = Q_(9.81, 'm/s^2')
	import Main

	lox_pump = Main.CreatePump(Main.Propellant.LOX)
	ipa_pump = Main.CreatePump(Main.Propellant.IPA)
	lox_pump.AddInducer(Q_(65, 'm'))
	ipa_pump.AddInducer(Q_(80, 'm'))
	
	fig, axes = plt.subplots(2, 2, figsize=(12, 9), sharex='col', sharey='row')
	for col, propellant, inducer in zip(axes.T, ('LOX', 'IPA'), (lox_pump.Inducer, ipa_pump.Inducer)):

		eq = InducerRadialEquilibrium(inducer, incidence_ratio=0, blade_thickness=Q_(0.5, 'mm'))
		print(f"{propellant}: beta1 =  {np.degrees(eq.beta1):.2f} deg, beta2b = {np.degrees(eq.beta2b):.2f} deg, delta2 = {np.degrees(eq.delta2):.2f} deg")
		for frac in (0.0, 0.25, 0.5, 0.75, 1.0):
			z = inducer.L_ax * frac
			r = np.linspace(eq.hub_radius(z).to('m').magnitude, eq.R_t, 50)
			col[0].plot(r * 1e3, eq.c_theta(Q_(r, 'm'), z).magnitude, label=f"z/L = {frac:.2f}")
			col[1].plot(r * 1e3, eq.c_axial(Q_(r, 'm'), z).magnitude, label=f"z/L = {frac:.2f}")
		col[0].set_title(f"{propellant} inducer")
		col[1].set_xlabel("r [mm]")
		for ax in col:
			ax.grid(True, linestyle='--', alpha=0.5)
			ax.legend()

		PlotMeanlineTriangles = eq.PlotMeanlineTriangles(
				f"{propellant} inducer meanline velocity triangles"
			).show()

	axes[0, 0].set_ylabel("$c_\\theta$ [m/s]")
	axes[1, 0].set_ylabel("$c_a$ [m/s]")
	plt.tight_layout()
	plt.show()
