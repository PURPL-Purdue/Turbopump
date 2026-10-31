import numpy as np
from pint import Quantity as Q_
import plotly.graph_objects as go
from . import Bezier
# (z, r, s_1, s_2) of each intersection of two curves: the point, and its normalized arclength on the curve FindIntersect is
# called on (s_1) and on the other curve (s_2). Ordered by s_1, empty when the curves do not intersect.
IntersectResults = list[tuple[Q_[float], Q_[float], float, float]]

def Cross(a: np.ndarray, b: np.ndarray) -> np.ndarray:
	"""2D cross product a x b = a_z b_r - a_r b_z of (..., 2) arrays of (z, r), broadcast over the leading axes"""
	return a[..., 0] * b[..., 1] - a[..., 1] * b[..., 0]

class MeridionalCurve:
	"""
		Curve in the meridional (z, r) plane evaluated by arclength. Subclasses set length_m and implement Points, all
		lengths in meters:
			Polyline:		straight lines through points: the geometric meanline, the potential flow streamlines, and
							any straight line (two points), such as the trailing edge
			BezierCurve:	cubic Bezier: the bends of the hub and shroud contours, the leading edge
			SegmentChain:	curves joined end to end: the hub and shroud contours (MeridionalContour), the FlowDomain walls
	"""

	length_m: float

	@property
	def Length(self) -> Q_[float]:
		return Q_(self.length_m, 'm')

	def Points(self, s_m: np.ndarray) -> np.ndarray:
		"""(n, 2) array of (z, r) at arclengths s_m (n,) in [0, length_m] from the start of the curve"""
		raise NotImplementedError

	def __call__(self, s: float | np.ndarray) -> tuple[Q_, Q_]:
		"""(z, r) at normalized arclength s in [0, 1]"""
		s_true = np.clip(np.atleast_1d(np.asarray(s, float)), 0, 1) * self.length_m
		pts = self.Points(s_true.ravel())

		z = pts[:, 0].reshape(s_true.shape)
		r = pts[:, 1].reshape(s_true.shape)

		if np.ndim(s) == 0:
			return Q_(z[0], 'm'), Q_(r[0], 'm')
		return Q_(z, 'm'), Q_(r, 'm')

	def Tangent(self, s: float, ds: float = 1e-4) -> np.ndarray:
		"""
			Unit tangent (z, r) at normalized arclength s, toward increasing s, by central differences over ds (one-sided
			at the ends of the curve)
		"""
		a, b = self.Points(np.clip([s - ds, s + ds], 0, 1) * self.length_m)
		return (b - a) / np.linalg.norm(b - a)

	def StaticMoment(self, s_1: float = 0, s_2: float = 1, n: int = 1001) -> Q_[float]:
		"""
			Static moment about the axis between normalized arclengths s_1 and s_2, S = integral(r dx) with x the
			arclength, by the trapezoid rule over n points (Milos eqs. 34-35)
		"""
		x = np.linspace(s_1, s_2, n) * self.length_m
		r = self.Points(x)[:, 1]
		return Q_(np.sum((r[1:] + r[:-1]) / 2 * np.diff(x)), 'm^2')

	def Closest(self, p: Q_[np.ndarray] | np.ndarray, n: int = 1001) -> tuple[Q_[float], float]:
		"""
			Distance from the point p, (z, r) as a Quantity or in meters, to this curve, and the normalized arclength of
			the closest point: the closest of n samples, then 9 times the closest of 21 samples spanning its two
			neighbours, 10 times finer each time. The distance is to the curve itself, not to the chords between samples.
		"""
		p = p.to('m').magnitude if isinstance(p, Q_) else np.asarray(p, float)
		lo, hi, m = 0.0, 1.0, n
		for _ in range(9):
			t = np.linspace(lo, hi, m)
			dist = np.linalg.norm(self.Points(t * self.length_m) - p, axis=1)
			k = np.argmin(dist)
			step = t[1] - t[0]
			lo, hi, m = max(t[k] - step, 0.0), min(t[k] + step, 1.0), 21
		return Q_(dist[k], 'm'), float(t[k])

	def FindIntersect(self, other: 'MeridionalCurve', n: int=1001) -> IntersectResults:
		"""
			Intersections of this curve with other, ordered by normalized arclength along this curve. Both curves are
			sampled at n points, uniform in normalized arclength, and the straight pieces between samples are
			intersected pairwise. Writing a piece of this curve as p + u d and a piece of other as q + v e, they cross at
				u = (q - p) x e / (d x e),	v = (q - p) x d / (d x e)
			with x the 2D cross product, and the crossing lies on both pieces when u and v are in [0, 1]. The point is
			evaluated on this curve, within about ds^2 * curvature / 8 of the other (ds: sample spacing), under 1e-7 m
			for the hub and shroud contours at the default n.

			A curve that starts or ends on the other (within 1 um) is an intersection at that start or end, with s_1 or s_2
			exactly 0 or 1. The pieces alone can miss these: where the curves join tangent-continuously the pieces are
			parallel, and where the other curve is curved its chord can fall just beyond the end. Crossings within one
			sample spacing of each other on both curves are reported once.

			Not handled: curves that overlap along a stretch (their pieces are parallel) beyond the ends of the stretch,
			and a curve with itself.
		"""
		tol = 1e-9					# on u and v, so a crossing at a sample point is not lost to round-off
		tol_end = Q_(1, 'um')		# distance within which a start or end of a curve is on the other curve

		s = np.linspace(0, 1, n)
		A = self.Points(s * self.length_m)
		B = other.Points(s * other.length_m)

		d = np.diff(A, axis=0)[:, None]			# pieces of this curve, (n - 1, 1, 2)
		e = np.diff(B, axis=0)[None]			# pieces of other, (1, n - 1, 2)
		w = B[None, :-1] - A[:-1, None]			# q - p for every pair of pieces, (n - 1, n - 1, 2)
		with np.errstate(divide='ignore', invalid='ignore'):		# parallel pieces: d x e = 0, u and v are inf or nan, never a hit
			den = Cross(d, e)
			u = Cross(w, e) / den
			v = Cross(w, d) / den
			i, j = np.nonzero((u >= -tol) & (u <= 1 + tol) & (v >= -tol) & (v <= 1 + tol))
		s_1 = np.clip((i + u[i, j]) / (n - 1), 0, 1)
		s_2 = np.clip((j + v[i, j]) / (n - 1), 0, 1)

		# starts and ends of one curve on the other, as (s_1, s_2)
		ends = []
		for s_end in (0, 1):
			dist, s_other = other.Closest(self.Points(np.array([s_end * self.length_m]))[0], n)
			if dist < tol_end:
				ends.append((s_end, s_other))
			dist, s_self = self.Closest(other.Points(np.array([s_end * other.length_m]))[0], n)
			if dist < tol_end:
				ends.append((s_self, s_end))

		# a crossing at a sample point is found by up to four neighbouring pairs of pieces, and a crossing at a start or
		# end by its pieces as well, keep one: the start or end first, it is the exact one
		found = []
		for s_a, s_b in ends + sorted(zip(s_1, s_2)):
			if not any(abs(s_a - a) < 1 / (n - 1) and abs(s_b - b) < 1 / (n - 1) for a, b in found):
				found.append((s_a, s_b))

		return [(*self(a), float(a), float(b)) for a, b in sorted(found)]

	def Crop(self, s_1: float, s_2: float) -> 'SubCurve':
		"""The part of this curve between normalized arclengths s_1 < s_2, as a curve of its own"""
		return SubCurve(self, s_1, s_2)

	def Between(self, a: 'MeridionalCurve', b: 'MeridionalCurve') -> 'SubCurve':
		"""
			The part of this curve between where it crosses the curves a and b, once each: e.g. a streamline between
			the leading and trailing edges, streamline.Between(LeadingEdge, TrailingEdge)
		"""
		s = []
		for other in (a, b):
			hits = self.FindIntersect(other)
			assert len(hits) == 1, 'The curve must cross each bound exactly once'
			s.append(hits[0][2])
		return self.Crop(*s)

	def Trace(self, l: str = 'in', n: int = 200, **scatter) -> go.Scatter:
		"""
			The curve as a plotly line in length unit l, through n points evenly spaced in arclength. scatter is passed
			on to go.Scatter (name, line, legendgroup, ...), with mode='lines' unless it says otherwise.
		"""
		z, r = self(np.linspace(0, 1, n))
		return go.Scatter(x=z.to(l).magnitude, y=r.to(l).magnitude, **({'mode': 'lines'} | scatter))

class Polyline(MeridionalCurve):
	"""
		Straight lines through points, an (n, 2) array of (z, r) in meters, evaluated by arclength. Two points make a
		straight line. s_m holds the arclength of each point from the first one.
	"""

	def __init__(self, points: np.ndarray):
		self.points = np.asarray(points, float)
		self.s_m = np.concatenate(([0], np.cumsum(np.linalg.norm(np.diff(self.points, axis=0), axis=1))))
		self.length_m = self.s_m[-1]

	def _Interp(self, values: np.ndarray, s_m: np.ndarray) -> np.ndarray:
		"""values given at each point, (n_points,) or (n_points, k), linearly interpolated at arclengths s_m"""
		if values.ndim == 1:
			return np.interp(s_m, self.s_m, values)
		return np.stack([np.interp(s_m, self.s_m, column) for column in values.T], axis=-1)

	def Points(self, s_m: np.ndarray) -> np.ndarray:
		"""(n, 2) array of (z, r) at arclengths s_m (n,) in [0, length_m] from the first point"""
		return self._Interp(self.points, s_m)

	def Interpolate(self, values: Q_[np.ndarray] | np.ndarray, s: float | np.ndarray) -> Q_[np.ndarray] | np.ndarray:
		"""
			values given at each point, (n_points,) or (n_points, k), at normalized arclengths s in [0, 1], interpolated
			linearly in arclength like the points themselves. A Quantity keeps its units, e.g. the meridional velocity
			along a streamline: streamline.Interpolate(streamline.v_m, s).
		"""
		s_m = np.clip(np.asarray(s, float), 0, 1) * self.length_m
		if isinstance(values, Q_):
			return Q_(self._Interp(values.magnitude, s_m), values.units)
		return self._Interp(np.asarray(values, float), s_m)

class BezierCurve(MeridionalCurve):
	"""
		Cubic Bezier with control points P, a (4, 2) array of (z, r) in meters (Bezier.bez), evaluated by arclength
		through a lookup table of the Bezier parameter t against arclength: n_table points, trapezoidal integration of
		|dB/dt|.
	"""

	def __init__(self, P: np.ndarray, n_table: int = 2001):
		self.P = np.asarray(P, float)

		t = np.linspace(0, 1, n_table)
		speed = np.linalg.norm(Bezier.dbez(self.P, t), axis=1)
		self._t = t
		self._s = np.concatenate(([0], np.cumsum((speed[1:] + speed[:-1]) / 2 * np.diff(t))))
		self.length_m = self._s[-1]

	@classmethod
	def ThroughPoints(cls, P_a: np.ndarray, P_m: np.ndarray, P_b: np.ndarray, n_table: int = 2001) -> 'BezierCurve':
		"""
			Quadratic Bezier from P_a to P_b through P_m, (z, r) in meters, with P_m at the fraction t of the polygon
			length P_a-P_m-P_b. Written as a cubic: the inner control points are 2/3 of the way from each end to the
			quadratic's middle control point P1.
		"""
		P_a, P_m, P_b = (np.asarray(P, float) for P in (P_a, P_m, P_b))
		l_a, l_b = np.linalg.norm(P_m - P_a), np.linalg.norm(P_b - P_m)
		assert l_a > 0 and l_b > 0, 'The middle point must differ from both ends'
		t = l_a / (l_a + l_b)
		P1 = (P_m - (1 - t)**2 * P_a - t**2 * P_b) / (2 * t * (1 - t))
		return cls(np.array([P_a, P_a + 2/3 * (P1 - P_a), P_b + 2/3 * (P1 - P_b), P_b]), n_table)

	@property
	def ControlPoints(self) -> Q_[np.ndarray]:
		return Q_(self.P, 'm')

	def Points(self, s_m: np.ndarray) -> np.ndarray:
		"""(n, 2) array of (z, r) at arclengths s_m (n,) in [0, length_m] from P[0]"""
		return Bezier.bez(self.P, np.interp(s_m, self._s, self._t))

class SegmentChain(MeridionalCurve):
	"""
		Curves joined end to end, evaluated by arclength along the whole chain. The segments can be any
		MeridionalCurves, chains included. Zero-length segments are dropped. s_m holds the arclength of each joint from
		the start: the start of every segment, then the end of the chain. All lengths in meters.
	"""

	def __init__(self, segments: list[MeridionalCurve]):
		self.Segments = [seg for seg in segments if seg.length_m > 0]
		self.s_m = np.concatenate(([0], np.cumsum([seg.length_m for seg in self.Segments])))
		self.length_m = self.s_m[-1]

	def Points(self, s_m: np.ndarray) -> np.ndarray:
		"""(n, 2) array of (z, r) at arclengths s_m (n,) in [0, length_m] from the start of the chain"""
		# segment containing each point, the end of the chain belongs to the last one
		i_seg = np.clip(np.searchsorted(self.s_m, s_m, side='right') - 1, 0, len(self.Segments) - 1)
		pts = np.empty((s_m.size, 2))
		for i, seg in enumerate(self.Segments):
			on_seg = i_seg == i
			pts[on_seg] = seg.Points(s_m[on_seg] - self.s_m[i])
		return pts

class SubCurve(MeridionalCurve):
	"""
		Part of the curve Parent between its normalized arclengths s_1 < s_2, evaluated by arclength from s_1 (see
		MeridionalCurve.Crop and Between). ToParent maps back to the parent, e.g. to interpolate a streamline's v_m.
	"""

	def __init__(self, parent: MeridionalCurve, s_1: float, s_2: float):
		assert 0 <= s_1 < s_2 <= 1, 'Need 0 <= s_1 < s_2 <= 1'
		self.Parent = parent
		self.s_1, self.s_2 = s_1, s_2
		self.length_m = (s_2 - s_1) * parent.length_m

	def Points(self, s_m: np.ndarray) -> np.ndarray:
		"""(n, 2) array of (z, r) at arclengths s_m (n,) in [0, length_m] from s_1 on the parent"""
		return self.Parent.Points(self.s_1 * self.Parent.length_m + np.asarray(s_m, float))

	def ToParent(self, s: float | np.ndarray) -> float | np.ndarray:
		"""Normalized arclength on the parent of normalized arclength s on this curve"""
		return self.s_1 + np.asarray(s, float) * (self.s_2 - self.s_1)

class MeridionalContour(SegmentChain):
	"""
		Meridional contour (z: axial, r: radial) made of tangent-continuous segments:
			1. flat:	straight axial section at r_in, of length flat (Polyline)
			2. bend:	cubic Bezier from axial to the outlet angle, ending at (z_bend, r_bend) (BezierCurve, Bend)
			3. outlet:	straight section at the outlet angle, from (z_bend, r_bend) to (z_out, r_out) (Polyline).
						Absent when r_bend = r_out (default), the bend then ends at the outlet.

		outlet_angle is the angle of the contour from the radial direction at the outlet, positive when
		z increases with r (contour leaning downstream). 0 gives a purely radial outlet.

		The Bezier's end tangents meet at the corner point C, where the outlet tangent line crosses the
		inlet line r = r_in, and both inner control points lie on those tangent lines:
			P1 = P0 + fullness[0] * (C - P0)
			P2 = P3 + fullness[1] * (C - P3)
		fullness -> 1 pulls the curve into the corner (tight bend), fullness -> 0 flattens it
		toward the straight chord P0-P3. ~0.55 approximates a circular arc.

		The contour is evaluated by normalized arclength s in [0, 1] (__call__).
		FlowDomain extends the hub and shroud contours past the inlet and outlet for the FEM solution.
	"""

	def __init__(self, r_in: Q_[float], z_out: Q_[float], r_out: Q_[float], flat: Q_[float],
			outlet_angle: Q_[float] = Q_(0, 'deg'), r_bend: Q_[float] = None,
			fullness: tuple[float, float] = (0.55, 0.55), z_in: Q_[float] = Q_(0, 'm'),
			n_table: int = 2001):

		r_bend = r_bend if r_bend is not None else r_out
		assert r_in < r_bend <= r_out, 'Radii must satisfy r_in < r_bend <= r_out'

		# all lengths stored as magnitudes in meters, angles in radians
		self.z_in = z_in.to('m').magnitude
		self.r_in = r_in.to('m').magnitude
		self.flat = flat.to('m').magnitude
		self.z_out = z_out.to('m').magnitude
		self.r_bend = r_bend.to('m').magnitude
		self.r_out = r_out.to('m').magnitude
		self.outlet_angle = outlet_angle.to('rad').magnitude
		assert abs(self.outlet_angle) < np.pi / 2, 'Outlet angle must be within 90 deg of radial'

		tan_out = np.tan(self.outlet_angle)
		self.z_bend = self.z_out - (self.r_out - self.r_bend) * tan_out

		P0 = np.array([self.z_in + self.flat, self.r_in])
		P3 = np.array([self.z_bend, self.r_bend])
		C = np.array([self.z_bend - (self.r_bend - self.r_in) * tan_out, self.r_in])
		assert C[0] > P0[0], 'Outlet tangent must cross the inlet line after the flat section, reduce flat or outlet_angle'

		self.Bend = BezierCurve(np.array([P0, P0 + fullness[0] * (C - P0), P3 + fullness[1] * (C - P3), P3]), n_table)

		super().__init__([
			Polyline([[self.z_in, self.r_in], P0]),		# flat
			self.Bend,									# bend
			Polyline([P3, [self.z_out, self.r_out]])	# outlet
		])

	@property
	def ControlPoints(self) -> Q_[np.ndarray]:
		"""Control points of the bend"""
		return self.Bend.ControlPoints

def MeanCurve(a: MeridionalCurve, b: MeridionalCurve, n: int = 201) -> Polyline:
	"""Midpoints of a and b at equal normalized arclength, e.g. the geometric meanline between the hub and shroud"""
	s = np.linspace(0, 1, n)
	return Polyline((a.Points(s * a.length_m) + b.Points(s * b.length_m)) / 2)

def QuasiOrthogonalFlow(hub: np.ndarray, shroud: np.ndarray, phi: float | np.ndarray, n_pts: int = 51
		) -> tuple[np.ndarray, np.ndarray]:
	"""
		Continuity through quasi-orthogonals (QOs), the straight lines joining homologous wall points
		hub[i] and shroud[i], given as (n_qo, 2) arrays of (z, r) in meters ordered inlet to outlet.
		Assumes c_m uniform along each QO, no blade blockage. Along each QO, the flow direction t is blended
		linearly from the hub wall tangent to the shroud wall tangent, and continuity through the QO gives
			Q = c_m * integral( (t . n) 2 pi r dl )
		where n is the QO normal. (t . n) corrects for QOs not being perpendicular to the flow.

		Returns
		1. lam: position on each QO (0: hub, 1: shroud) such that the flow fraction phi passes between the
			hub and lam, shape (n_qo,) + shape of phi
		2. area: (n_qo,) effective flow area of each QO, c_m = Q / area
	"""
	n_qo = len(hub)

	# wall flow directions
	t_hub = np.gradient(hub, axis=0)
	t_hub /= np.linalg.norm(t_hub, axis=1, keepdims=True)
	t_shroud = np.gradient(shroud, axis=0)
	t_shroud /= np.linalg.norm(t_shroud, axis=1, keepdims=True)

	chord = shroud - hub
	L = np.linalg.norm(chord, axis=1)
	n = np.column_stack([chord[:, 1], -chord[:, 0]]) / L[:, None]	# QO normal, pointing downstream

	# points along each QO, lam = 0 at hub, 1 at shroud: (n_qo, n_pts, 2)
	lam = np.linspace(0, 1, n_pts)
	pts = hub[:, None, :] + lam[None, :, None] * chord[:, None, :]
	t = (1 - lam[None, :, None]) * t_hub[:, None, :] + lam[None, :, None] * t_shroud[:, None, :]
	t /= np.linalg.norm(t, axis=2, keepdims=True)

	# effective flow area per unit lam, and its cumulative integral from the hub
	dA = np.sum(t * n[:, None, :], axis=2) * 2 * np.pi * pts[:, :, 1] * L[:, None]
	A = np.concatenate((np.zeros((n_qo, 1)), np.cumsum((dA[:, 1:] + dA[:, :-1]) / 2 * np.diff(lam), axis=1)), axis=1)

	# position of the phi flow fraction on each QO
	lam_phi = np.array([np.interp(phi, A_qo / A_qo[-1], lam) for A_qo in A])

	return lam_phi, A[:, -1]

class FlowDomain:
	"""
		Meridional flow domain for the FEM potential flow solution (Milos, Impeller Design Using CAD Techniques
		and Conformal Mapping Method, sections 2-7). The shroud and hub contours are extended by straight lines:
			inlet:	axial, upstream to z_min, inlet_ratio inlet passage heights before the contour inlets.
					The inlet boundary is then a radial line where the flow can be taken as uniform and axial.
					Disturbances from the bend decay upstream roughly like exp(-pi * distance / height).
			outlet:	on at the contour's outlet angle to r_corner = corner_ratio * r_out, then radial out to
					r_max = outlet_ratio * r_out. The outlet boundary is then an axial line where the flow can be
					taken as uniform and radial (Milos uses D_max = 1.25 D_2). A tilted contour has a corner at
					r_corner, where the potential flow velocity is singular. Milos puts it at 1.05 D_2 (point G) so
					this does not distort the impeller outlet velocities. corner_ratio = 1 turns radial at the outlet.

		Along the flow the domain has three zones, ZONES. The extended walls Shroud and Hub are SegmentChains with
		the contour as one segment, and ZoneBounds[wall] holds the arclengths (m) along each wall of: start, impeller
		inlet, impeller outlet, end.
	"""

	ZONES = ('inlet', 'impeller', 'outlet')

	def __init__(self, shroud: MeridionalContour, hub: MeridionalContour,
			inlet_ratio: float = 2.0, outlet_ratio: float = 1.25, corner_ratio: float = 1.05):

		assert shroud.r_in > hub.r_in, 'Shroud must be above the hub at the inlet'
		assert inlet_ratio > 0 and 1 <= corner_ratio < outlet_ratio, 'Extensions must have a positive length, corner before the end'

		self.z_min = min(shroud.z_in, hub.z_in) - inlet_ratio * (shroud.r_in - hub.r_in)
		self.r_corner = corner_ratio * max(shroud.r_out, hub.r_out)
		self.r_max = outlet_ratio * max(shroud.r_out, hub.r_out)

		self.Shroud = self._Extend(shroud)
		self.Hub = self._Extend(hub)

		# wall segments: inlet extension, the contour, outlet extensions. The first two always have a length, so their
		# joints are the impeller inlet and outlet.
		self.ZoneBounds = {
			'shroud': self.Shroud.s_m[[0, 1, 2, -1]],
			'hub': self.Hub.s_m[[0, 1, 2, -1]]
		}

	def _Extend(self, contour: MeridionalContour) -> SegmentChain:
		corner = np.array([contour.z_out + (self.r_corner - contour.r_out) * np.tan(contour.outlet_angle), self.r_corner])
		return SegmentChain([
			Polyline([[self.z_min, contour.r_in], [contour.z_in, contour.r_in]]),	# axial inlet
			contour,
			Polyline([[contour.z_out, contour.r_out], corner]),		# on at the outlet angle
			Polyline([corner, [corner[0], self.r_max]])				# radial outlet
		])

	def Stations(self, n_stream: int) -> tuple[np.ndarray, np.ndarray, dict[str, tuple[int, int]]]:
		"""
			Homologous hub and shroud points, (n_stream + 1, 2) arrays of (z, r), and the first and last station
			of each zone. Within a zone, stations are equally spaced in arclength along each wall, so the zone
			boundaries (impeller inlet and outlet) are QOs. The n_stream intervals are shared between zones in
			proportion to their mean hub and shroud length.
		"""
		lengths = (np.diff(self.ZoneBounds['hub']) + np.diff(self.ZoneBounds['shroud'])) / 2
		n_zone = np.maximum(1, np.round(n_stream * lengths / lengths.sum()).astype(int))
		n_zone[1] += n_stream - n_zone.sum()		# rounding remainder goes to the impeller
		assert np.all(n_zone >= 1), 'n_stream is too small to mesh every zone'

		walls = []
		for wall, bounds in ((self.Hub, self.ZoneBounds['hub']), (self.Shroud, self.ZoneBounds['shroud'])):
			s = [np.linspace(bounds[k], bounds[k + 1], n_zone[k] + 1) for k in range(len(self.ZONES))]
			s = np.concatenate([s[0]] + [zone[1:] for zone in s[1:]])	# neighbouring zones share their boundary station
			walls.append(wall.Points(s))

		first = np.concatenate(([0], np.cumsum(n_zone)))
		zones = {zone: (int(first[k]), int(first[k + 1])) for k, zone in enumerate(self.ZONES)}

		return walls[0], walls[1], zones

	def Mesh(self, n_stream: int = 100, n_span: int = 10, by_flow: bool = True, n_pts: int = 201) -> 'MeridionalMesh':
		"""
			Structured mesh of n_stream x n_span quadrilaterals (Milos section 3, 101 x 11 nodes by default):
			homologous hub and shroud stations (see Stations) are joined by straight QOs, and each QO is divided
			into n_span intervals by the streamlines.
			by_flow:	streamlines at equal flow fractions phi = j / n_span, estimated with QuasiOrthogonalFlow
						(n_pts points per QO). These are exact streamlines on the inlet and outlet boundaries,
						where the flow is uniform, so there the stream function is psi_max * phi.
			otherwise:	QOs divided into equal lengths, as Milos does.
		"""
		hub, shroud, zones = self.Stations(n_stream)
		phi = np.linspace(0, 1, n_span + 1)

		if by_flow:
			lam, _ = QuasiOrthogonalFlow(hub, shroud, phi, n_pts)	# (n_stations, n_lines)
		else:
			lam = np.broadcast_to(phi, (len(hub), n_span + 1))
		grid = hub[:, None, :] + lam[:, :, None] * (shroud - hub)[:, None, :]

		return MeridionalMesh(list(grid.transpose(1, 0, 2)), phi=phi if by_flow else None, zones=zones)

def Polylines(curves: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
	"""(n_curves, n_pts, 2) curves as the x, y of a single plotly trace, with gaps between curves"""
	gaps = np.full((len(curves), 1, 2), np.nan)
	xy = np.concatenate((curves, gaps), axis=1).reshape(-1, 2)
	return xy[:, 0], xy[:, 1]

class MeridionalMesh:
	"""
		Structured mesh of 4-node quadrilateral elements between streamlines, for the FEM solution of the
		meridional stream function and velocity potential (Milos, sections 3-7). All lengths in meters.

		Built from streamlines ordered hub to shroud, each a (n_stations, 2) array of (z, r) sampled at the same
		stations, inlet to outlet. Homologous stations are joined by straight quasi-orthogonals (QOs), so element
		(i, j) lies between QOs i, i + 1 and streamlines j, j + 1.
			grid:		(n_stations, n_lines, 2), grid[i, j] is station i of streamline j
			nodes:		(n_nodes, 2), node i * n_lines + j is grid[i, j]. Numbering along the QOs keeps the
						bandwidth of the global FEM matrix to about n_lines.
			elements:	(n_elements, 4) node indices, counterclockwise in the (z, r) plane:
						(i, j) -> (i + 1, j) -> (i + 1, j + 1) -> (i, j + 1). Element i * (n_lines - 1) + j.
			Boundary:	node indices of the hub (first streamline), shroud (last streamline), inlet (first QO)
						and outlet (last QO), corner nodes belong to two boundaries
			phi:		flow fraction of each streamline (hub 0, shroud 1), None if unknown
			Zones:		first and last QO of each zone, e.g. {'impeller': (i_inlet, i_outlet)}
	"""

	def __init__(self, streamlines: list[np.ndarray], phi: np.ndarray = None, zones: dict[str, tuple[int, int]] = None):

		self.grid = np.stack(streamlines, axis=1)
		self.n_stations, self.n_lines, _ = self.grid.shape
		assert self.n_stations >= 2 and self.n_lines >= 2, 'Need at least 2 streamlines of 2 stations'

		self.nodes = self.grid.reshape(-1, 2)
		node = np.arange(self.n_stations * self.n_lines).reshape(self.n_stations, self.n_lines)
		self.elements = np.column_stack([
			node[:-1, :-1].ravel(), node[1:, :-1].ravel(), node[1:, 1:].ravel(), node[:-1, 1:].ravel()
		])
		self.n_nodes = len(self.nodes)
		self.n_elements = len(self.elements)

		self.Boundary = {'hub': node[:, 0], 'shroud': node[:, -1], 'inlet': node[0], 'outlet': node[-1]}
		self.phi = phi
		self.Zones = zones if zones is not None else {}

		# bilinear elements need a positive Jacobian everywhere: convex, counterclockwise corners
		quad = self.nodes[self.elements]
		to_next = np.roll(quad, -1, axis=1) - quad
		to_prev = np.roll(quad, 1, axis=1) - quad
		corner = Cross(to_next, to_prev)
		assert np.all(corner > 0), 'Elements must be convex with counterclockwise nodes, streamlines or QOs cross'

	@property
	def Areas(self) -> np.ndarray:
		"""(n_elements,) area of each element in the (z, r) plane"""
		quad = self.nodes[self.elements]
		d13 = quad[:, 2] - quad[:, 0]
		d24 = quad[:, 3] - quad[:, 1]
		return Cross(d13, d24) / 2

	@property
	def Centers(self) -> np.ndarray:
		"""(n_elements, 2) element centers (xi = eta = 0), where Milos evaluates the velocities"""
		return self.nodes[self.elements].mean(axis=1)

	@property
	def AspectRatios(self) -> np.ndarray:
		"""(n_elements,) longest over shortest side, Milos keeps the sides of similar length for accuracy"""
		quad = self.nodes[self.elements]
		sides = np.linalg.norm(np.roll(quad, -1, axis=1) - quad, axis=2)
		return sides.max(axis=1) / sides.min(axis=1)

	def ZoneElements(self, zone: str) -> np.ndarray:
		"""Indices of the elements between the first and last QO of a zone"""
		first, last = self.Zones[zone]
		return np.arange(first * (self.n_lines - 1), last * (self.n_lines - 1))

	def Traces(self, l: str = 'in', color: str = 'black', qo_color: str = 'grey', zone_bounds: bool = True
			) -> list[go.Scatter]:
		"""
			Mesh QOs and streamlines (one legend entry), and the zone boundaries (QOs between zones) if zone_bounds,
			as plotly traces in length unit l
		"""
		scale = Q_(1, 'm').to(l).magnitude

		qo_x, qo_y = Polylines(self.grid * scale)
		line_x, line_y = Polylines(self.grid.transpose(1, 0, 2) * scale)
		traces = [
			go.Scatter(x=qo_x, y=qo_y, mode='lines', name='FEM mesh', legendgroup='mesh', showlegend=False,
				line=dict(color=qo_color, width=0.5)),
			go.Scatter(x=line_x, y=line_y, mode='lines', name='FEM mesh', legendgroup='mesh',
				line=dict(color=color, width=1))
		]

		inner_bounds = sorted({i for bounds in self.Zones.values() for i in bounds} - {0, self.n_stations - 1})
		if zone_bounds and inner_bounds:
			x, y = Polylines(self.grid[inner_bounds] * scale)
			traces.append(go.Scatter(x=x, y=y, mode='lines', name='Zone boundaries', line=dict(color='red', width=2)))

		return traces

	def Plot(self, l: str = 'in', color: str = 'black') -> go.Figure:
		"""Streamlines and QOs of the mesh, with the zone boundaries highlighted"""
		fig = go.Figure(self.Traces(l, color))

		fig.update_xaxes(title_text=f"Axial position ({l})")
		fig.update_yaxes(title_text=f"Radial position ({l})", scaleanchor='x', scaleratio=1)
		fig.update_layout(
			title=f"Meridional Mesh: {self.n_elements} elements, {self.n_nodes} nodes",
			template="plotly_white"
		)

		fig.show()

		return fig
