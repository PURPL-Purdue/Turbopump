import Bezier
import numpy as np
from pint import Quantity as Q_

class _LineSegment:
	"""Straight segment from P0 to P1, evaluated by arclength"""

	def __init__(self, P0: np.ndarray, P1: np.ndarray):
		self.P0 = P0
		self.P1 = P1
		self.length = np.linalg.norm(P1 - P0)

	def __call__(self, s: np.ndarray) -> np.ndarray:
		"""(n, 2) points at arclengths s from P0"""
		return self.P0 + (s / self.length)[:, None] * (self.P1 - self.P0)

class _BezierSegment:
	"""Cubic Bezier with control points P, evaluated by arclength"""

	def __init__(self, P: np.ndarray, n_table: int):
		self.P = P

		# lookup table of Bezier parameter t vs arclength, trapezoidal integration of |dB/dt|
		t = np.linspace(0, 1, n_table)
		speed = np.linalg.norm(Bezier.dbez(P, t), axis=1)
		self._t = t
		self._s = np.concatenate(([0], np.cumsum((speed[1:] + speed[:-1]) / 2 * np.diff(t))))
		self.length = self._s[-1]

	def __call__(self, s: np.ndarray) -> np.ndarray:
		"""(n, 2) points at arclengths s from P[0]"""
		return Bezier.bez(self.P, np.interp(s, self._s, self._t))

class MeridionalContour:
	"""
		Meridional contour (z: axial, r: radial) made of tangent-continuous segments:
			1. flat:	straight axial section at r_in, of length flat
			2. bend:	cubic Bezier from axial to the outlet angle, ending at (z_bend, r_bend)
			3. outlet:	straight section at the outlet angle, from (z_bend, r_bend) to (z_out, r_out).
						Absent when r_bend = r_out (default), the bend then ends at the outlet.

		outlet_angle is the angle of the contour from the radial direction at the outlet, positive when
		z increases with r (contour leaning downstream). 0 gives a purely radial outlet.

		The Bezier's end tangents meet at the corner point C, where the outlet tangent line crosses the
		inlet line r = r_in, and both inner control points lie on those tangent lines:
			P1 = P0 + fullness[0] * (C - P0)
			P2 = P3 + fullness[1] * (C - P3)
		fullness -> 1 pulls the curve into the corner (tight bend), fullness -> 0 flattens it
		toward the straight chord P0-P3. ~0.55 approximates a circular arc.

		The contour can be evaluated by normalized arclength s in [0, 1] (__call__), or as r(z)
		(GetRadius) when z is monotonic along the contour, i.e. outlet_angle >= 0.
		Segments are evaluated generically, so the contour can be extended upstream of the inlet and
		past the outlet (e.g. for an FEM domain) by adding segments.
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

		self.P = np.array([P0, P0 + fullness[0] * (C - P0), P3 + fullness[1] * (C - P3), P3])

		segments = [
			_LineSegment(np.array([self.z_in, self.r_in]), P0),		# flat
			_BezierSegment(self.P, n_table),						# bend
			_LineSegment(P3, np.array([self.z_out, self.r_out]))	# outlet
		]
		self._segments = [seg for seg in segments if seg.length > 0]
		self._s_start = np.concatenate(([0], np.cumsum([seg.length for seg in self._segments])))

		self.length_m = self._s_start[-1]
		self.Length: Q_[float] = Q_(self.length_m, 'm')

	@property
	def ControlPoints(self) -> Q_[np.ndarray]:
		return Q_(self.P, 'm')

	def __call__(self, s: float | np.ndarray) -> tuple[Q_, Q_]:
		"""(z, r) at normalized arclength s in [0, 1]"""
		s_true = np.clip(np.atleast_1d(np.asarray(s, float)), 0, 1) * self.length_m
		s_flat = s_true.ravel()

		# segment containing each point, s = 1 belongs to the last one
		i_seg = np.minimum(np.searchsorted(self._s_start, s_flat, side='right') - 1, len(self._segments) - 1)
		pts = np.empty((s_flat.size, 2))
		for i, seg in enumerate(self._segments):
			on_seg = i_seg == i
			pts[on_seg] = seg(s_flat[on_seg] - self._s_start[i])

		z = pts[:, 0].reshape(s_true.shape)
		r = pts[:, 1].reshape(s_true.shape)

		if np.ndim(s) == 0:
			return Q_(z[0], 'm'), Q_(r[0], 'm')
		return Q_(z, 'm'), Q_(r, 'm')

	def GetRadius(self, z: Q_) -> Q_:
		"""
			r at axial position z, NaN outside [z_in, z_out]. Needs outlet_angle >= 0 so z is monotonic.
			z is constant on a radial outlet section, so there it returns r_bend at z_out.
		"""
		assert self.outlet_angle >= 0, 'r(z) is multivalued when the contour leans upstream at the outlet'
		z_m = np.atleast_1d(z.to('m').magnitude).astype(float)

		# sample each segment end to end so the bend's end point is in the table, then keep only points
		# that advance in z, which drops the constant z of a radial outlet section
		zr = np.vstack([seg(np.linspace(0, seg.length, 2001)) for seg in self._segments])
		rising = np.concatenate(([True], zr[1:, 0] > np.maximum.accumulate(zr[:-1, 0])))
		r = np.interp(z_m, zr[rising, 0], zr[rising, 1])
		r[(z_m < self.z_in) | (z_m > self.z_out)] = np.nan

		return Q_(r[0] if np.ndim(z.magnitude) == 0 else r, 'm')
