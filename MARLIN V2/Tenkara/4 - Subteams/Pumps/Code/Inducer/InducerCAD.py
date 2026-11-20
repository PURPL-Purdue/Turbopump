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

# TODO: these should be attributes instantiated in constructor
WEDGE_ANGLE = np.deg2rad(6) # radians
END_WEDGE_ANGLE = np.deg2rad(5) # radians

def _Arclength(points_zrt: np.ndarray) -> np.ndarray:
	"""Cumulative arclength (m) along a (z, r, theta) line, dl^2 = dz^2 + dr^2 + (r dtheta)^2."""
	z, r, theta = points_zrt.T
	return np.concatenate(([0], np.cumsum(np.sqrt(np.diff(z)**2 + np.diff(r)**2 + (r[1:] * np.diff(theta))**2))))

def _WedgeLengths(thickness: float) -> tuple[float, float]:
	"""Camber arclength (m) of the leading edge wedge (one-sided, WEDGE_ANGLE) and the trailing edge taper (both sides at END_WEDGE_ANGLE to the camber), for a blade thickness in m."""
	return thickness / np.tan(WEDGE_ANGLE), thickness / (2 * np.tan(END_WEDGE_ANGLE))

def _Offsets(arclength: np.ndarray, thickness: float) -> tuple[np.ndarray, np.ndarray]:
	"""
	Axial offsets (m) of the suction and pressure sides from the camber, along camber arclength, for a blade thickness t in m:
	- leading edge: one-sided wedge, the suction side flat at +t/2 and the pressure side opening from +t/2 at WEDGE_ANGLE, thickness = tan(WEDGE_ANGLE) * s;
	- middle: +-t/2;
	- trailing edge: symmetric taper, both sides closing onto the camber at END_WEDGE_ANGLE, half thickness = tan(END_WEDGE_ANGLE) * distance from the trailing edge.
	"""
	t_LE = np.minimum(arclength * np.tan(WEDGE_ANGLE), thickness)
	half_TE = np.minimum((arclength[-1] - arclength) * np.tan(END_WEDGE_ANGLE), thickness / 2)
	return half_TE, np.maximum(thickness / 2 - t_LE, -half_TE)

def _EdgeBetween(face: cadquery.Shape, a: cadquery.Vector, b: cadquery.Vector, tol: float=1e-4) -> cadquery.Edge:
	"""The edge of face that runs between points a and b (mm), in either direction."""
	near = lambda v, w: (v - w).Length < tol
	return next(e for e in face.Edges()
		if (near(e.startPoint(), a) and near(e.endPoint(), b)) or (near(e.startPoint(), b) and near(e.endPoint(), a)))

class InducerCAD:

	def __init__(self, radial_equilibrium: InducerRadialEquilibrium):
		self.RadialEquilibrium = radial_equilibrium
		self._cad = None

	def _MoveAxially(self, points_zrt: np.ndarray, dz: np.ndarray) -> np.ndarray:
		"""Points moved by dz along z at the same span fraction, so hub points stay on the tapered hub."""
		eq = self.RadialEquilibrium
		z, r, theta = points_zrt.T
		span_fraction = (r - eq.r_hub(z)) / (eq.R_t - eq.r_hub(z))
		z_new = z + dz
		r_new = eq.r_hub(z_new) + span_fraction * (eq.R_t - eq.r_hub(z_new))
		return np.column_stack((z_new, r_new, theta))

	def GetBladeLines(self, num_lines: int=7, resolution: int = 200, num_sections: int=60, num_wedge_sections: int=10) -> tuple[list[np.ndarray], list[np.ndarray]]:
		"""
		@Charlie zhu
		Suction and pressure lines from hub to tip, lists of (n, 3) arrays of (z, r, theta) in m, m and rad. Point j of every line is on section j.
		"""
		eq = self.RadialEquilibrium
		z_ref = np.linspace(0, eq.L, resolution)

		r_hub: np.ndarray[float] = eq.r_hub(z_ref)
		r_tip: np.ndarray[float] = np.full(resolution, eq.R_t)
		r_diff: np.ndarray[float] = r_tip - r_hub

		#points_zrt = []

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
		ss_points_zrt = []
		ps_points_zrt = []

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

			camber_points = np.vstack((leading_edge, points[i:]))
			#points_zrt.append(np.vstack((leading_edge, points[i:])))
			camber_arclen = _Arclength(camber_points)

			# resample at the sections: the leading edge wedge, the middle and the trailing edge taper
			length = camber_arclen[-1]
			LE_length, TE_length = _WedgeLengths(eq.t)
			s_LE_wedge = np.linspace(0, LE_length, num_wedge_sections)
			s_middle = np.linspace(LE_length, length - TE_length, num_sections)
			s_TE_taper = np.linspace(length - TE_length, length, num_wedge_sections)
			s = np.concatenate((s_LE_wedge, s_middle[1:], s_TE_taper[1:]))		# shared ends counted once
			camber_points = np.column_stack([np.interp(s, camber_arclen, c) for c in camber_points.T])

			dz_ss, dz_ps = _Offsets(s, eq.t)
			ss_points_zrt.append(self._MoveAxially(camber_points, dz_ss))
			ps_points_zrt.append(self._MoveAxially(camber_points, dz_ps))

		return ss_points_zrt, ps_points_zrt

	def CreateCAD(self, num_lines: int=7, resolution: int=300, flip_z: bool=True, flip_theta: bool=True, num_sections: int=60,
				  num_wedge_sections: int=10) -> cadquery.Assembly:
		"""
		@Charlie Zhu
		Blade solids and the revolved hub (x = r cos theta, y = r sin theta, z, in mm). Each side is lofted in three stretches between shared spanwise edges,
		the two sides share the leading and trailing edges, and the hub and tip ends are ruled faces."""
		eq = self.RadialEquilibrium
		suction_bladelines, pressure_bladelines = self.GetBladeLines(num_lines, resolution, num_sections, num_wedge_sections)
		ToCQ = lambda points_zrt: CAD.CylindricalToCartesian(points_zrt, flip_z, flip_theta)

		S, P = np.array(suction_bladelines), np.array(pressure_bladelines)		# (num_lines, number of sections, 3)
		assert np.allclose(S[:, [0, -1]], P[:, [0, -1]]), 'Suction and pressure lines should share their leading and trailing edge points'
		# section index of the leading edge, the wedge ends and the trailing edge, as laid out in GetBladeLines
		breaks = [0, num_wedge_sections - 1, num_wedge_sections + num_sections - 2, S.shape[1] - 1]

		# boundary edges of the stretches, the leading and trailing edges shared by both sides
		Spanwise = lambda points, j: cadquery.Edge.makeSpline(ToCQ(points[:, j]))
		leading_edge, trailing_edge = Spanwise(S, 0), Spanwise(S, -1)
		suction_boundaries = [leading_edge, Spanwise(S, breaks[1]), Spanwise(S, breaks[2]), trailing_edge]
		pressure_boundaries = [leading_edge, Spanwise(P, breaks[1]), Spanwise(P, breaks[2]), trailing_edge]

		faces = []
		for k in range(3):
			a, b = breaks[k], breaks[k + 1]
			# C1: the default C2 loft overshoots between sections near the leading edge (one face bulged ~5 mm out, 50 % too much volume)
			Side = lambda points, boundaries: cadquery.func.loft([boundaries[k], *[Spanwise(points, j) for j in range(a + 1, b)], boundaries[k + 1]], continuity='C1')
			suction_face, pressure_face = Side(S, suction_boundaries), Side(P, pressure_boundaries)
			# hub (line 0) and tip (line -1) ends, ruled between the suction and pressure edges; they close to a point at the leading and trailing edges
			# (a filled surface fails on these long, thin wedge-shaped ends)
			end_faces = [cadquery.Face.makeRuledSurface(_EdgeBetween(suction_face, *ToCQ(S[line, [a, b]])), _EdgeBetween(pressure_face, *ToCQ(P[line, [a, b]])))
						 for line in (0, -1)]
			faces += [suction_face, pressure_face, *end_faces]
		blade = cadquery.func.solid(cadquery.func.shell(*faces))

		blade_locations = [cadquery.Location((0, 0, 0), (0, 0, 1),  360 * n / eq.Z_blade) for n in range(eq.Z_blade)]
		blades = blade.moved(*blade_locations)
		# suction and pressure blade lines of the first blade, as splines through the points of GetBladeLines
		blade_line_edges = [cadquery.Edge.makeSpline(ToCQ(line)) for line in (*suction_bladelines, *pressure_bladelines)]
		# wedge ends, hub to tip: the boundary edges between the stretches where a side has a corner, the leading edge wedge on the pressure side and
		# the start of the trailing edge taper on both sides
		wedge_end_edges = [pressure_boundaries[1], suction_boundaries[2], pressure_boundaries[2]]

		# Linearly tapered hub: a line from (z, r) = (0, R_h1) to (L, R_h2) at theta = 0, revolved about the z axis
		hub_edge = cadquery.Edge.makeLine(*CAD.CylindricalToCartesian(np.array([[0, eq.R_h1, 0], [eq.L, eq.R_h2, 0]]), flip_z, flip_theta))
		hub_surface = cadquery.func.revolve(hub_edge, (0, 0, 0), (0, 0, 1))
		assm = (cadquery.Assembly()
			.add(hub_surface, name='hub', color=cadquery.Color(0.6, 0.6, 0.6))
			.add(blades, name='blades', color=cadquery.Color(0.2, 0.3, 0.7, a=0.5))
		)
		self._cad = [assm, *blade_line_edges, *wedge_end_edges]		# the blade line and wedge end edges are shown with debug
		return assm

	def ShowCAD(self, debug: bool=True):
		if self._cad is None: self.CreateCAD()
		show(*(self._cad if debug else self._cad[:1]))
