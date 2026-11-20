import numpy as np
from scipy import integrate
import cadquery.func
from cadquery.vis import show
from Helpers import CAD
from Inducer.RadialEquilibrium import InducerRadialEquilibrium, leading_edge_theta

"""
@Charlie Zhu
Inducer CAD generation, moved out of RadialEquilibrium.py
"""

class InducerCAD:

	def __init__(self, radial_equilibrium: InducerRadialEquilibrium):
		self.RadialEquilibrium = radial_equilibrium
		self._cad = None

	def GetBladeLines(self, num_lines: int=7, resolution: int = 200) -> list[np.ndarray]:
		"""
		@Charlie zhu
		Blade lines from hub to tip, a list of (resolution, 3) arrays of (z, r, theta) in m, m and rad, from the inlet.
		"""
		eq = self.RadialEquilibrium
		z_ref = np.linspace(0, eq.L, resolution)

		r_hub: np.ndarray[float] = eq.r_hub(z_ref)
		r_tip: np.ndarray[float] = np.full(resolution, eq.R_t)
		r_diff: np.ndarray[float] = r_tip - r_hub

		points_zrt = []

		radial_lead = eq.r_ms(z_ref) * np.tan(eq.beta(z_ref, eq.beta2b))
		theta_ref = integrate.cumulative_trapezoid(y=1 / radial_lead, x=z_ref, initial=0)

		# radii clustered toward the tip
		radii_LE = eq.R_h1 + (eq.R_t - eq.R_h1) * np.sin(np.linspace(0, np.pi / 2, resolution))

		get_LE_theta = leading_edge_theta(eq.R_h1, eq.R_t, wrap_angle=eq.sweepback)
		theta_LE = get_LE_theta(radii_LE)
		assert theta_LE[-1] <= theta_ref[-1], 'The leading edge wraps further than the blade'
		z_LE_ref = np.interp(theta_LE, theta_ref, z_ref)		# theta_ref increases with z_ref

		hub_taper_slope = (eq.R_h2 - eq.R_h1) / eq.L
		z_LE_canted = z_LE_ref - (radii_LE - eq.r_hub(z_LE_ref)) * hub_taper_slope

		arclength_LE = np.concatenate(([0], np.cumsum(np.sqrt(np.diff(z_LE_canted)**2 + np.diff(radii_LE)**2 + (radii_LE[1:] * np.diff(theta_LE))**2))))
		s = np.linspace(0, arclength_LE[-1], num_lines)
		# keep z and r values at the arclengths specified s (linearly interpolated)
		z_LE_ref, radii_LE = np.interp(s, arclength_LE, z_LE_ref), np.interp(s, arclength_LE, radii_LE)

		LE_span_fractions: np.ndarray[float] = (radii_LE - eq.r_hub(z_LE_ref)) / (eq.R_t - eq.r_hub(z_LE_ref))

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
		eq = self.RadialEquilibrium
		blade_lines = self.GetBladeLines(num_lines=num_lines, resolution=resolution)
		blade_edges = [cadquery.Edge.makeSpline(CAD.CylindricalToCartesian(p, flip_z, flip_theta)) for p in blade_lines]
		camber = cadquery.func.loft(blade_edges)
		blade_locations = [cadquery.Location((0, 0, 0), (0, 0, 1),  360 * n / eq.Z_blade) for n in range(eq.Z_blade)]
		blades = camber.moved(*blade_locations)

		# Linearly tapered hub: a line from (z, r) = (0, R_h1) to (L, R_h2) at theta = 0, revolved about the z axis
		hub_edge = cadquery.Edge.makeLine(*CAD.CylindricalToCartesian(np.array([[0, eq.R_h1, 0], [eq.L, eq.R_h2, 0]]), flip_z, flip_theta))
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
