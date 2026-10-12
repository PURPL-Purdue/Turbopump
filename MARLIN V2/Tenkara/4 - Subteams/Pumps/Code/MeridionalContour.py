import Bezier
import numpy as np
from pint import Quantity as Q_

class MeridionalContour:
	"""
		Meridional contour (z: axial, r: radial) made of three tangent-continuous segments:
			1. flat:	straight axial section at r_in, of length flat
			2. bend:	cubic Bezier from axial to purely radial, ending at (z_out, r_bend)
			3. radial:	straight radial section at z_out, from r_bend to r_out

		The Bezier's end tangents meet at the corner point C = (z_out, r_in), and both inner
		control points lie on those tangent lines:
			P1 = P0 + fullness[0] * (C - P0)
			P2 = P3 + fullness[1] * (C - P3)
		fullness -> 1 pulls the curve into the corner (tight bend), fullness -> 0 flattens it
		toward the straight chord P0-P3. ~0.55 approximates a circular arc.

		The contour can be evaluated by normalized arclength s in [0, 1] (__call__), or as r(z)
		(GetRadius) since z is monotonic along the contour. z is constant on the radial section,
		so r(z) only describes the contour up to r_bend.
	"""

	def __init__(self, r_in: Q_[float], z_out: Q_[float], r_bend: Q_[float], r_out: Q_[float],
			flat: Q_[float], fullness: tuple[float, float] = (0.55, 0.55), z_in: Q_[float] = Q_(0, 'm'),
			n_table: int = 2001):

		assert r_in < r_bend <= r_out, 'Radii must satisfy r_in < r_bend <= r_out'
		assert z_in + flat < z_out, 'Flat section must end before the outlet'

		# all lengths stored as magnitudes in meters
		self.z_in = z_in.to('m').magnitude
		self.r_in = r_in.to('m').magnitude
		self.flat = flat.to('m').magnitude
		self.z_out = z_out.to('m').magnitude
		self.r_bend = r_bend.to('m').magnitude
		self.r_out = r_out.to('m').magnitude

		P0 = np.array([self.z_in + self.flat, self.r_in])
		P3 = np.array([self.z_out, self.r_bend])
		C = np.array([self.z_out, self.r_in])

		self.P = np.array([P0, P0 + fullness[0] * (C - P0), P3 + fullness[1] * (C - P3), P3])

		# lookup table of Bezier parameter t vs arclength, trapezoidal integration of |dB/dt|
		t = np.linspace(0, 1, n_table)
		speed = np.linalg.norm(Bezier.dbez(self.P, t), axis=1)
		self._t = t
		self._s = np.concatenate(([0], np.cumsum((speed[1:] + speed[:-1]) / 2 * np.diff(t))))
		self._pts = Bezier.bez(self.P, t)

		self.s_bend_m = self.flat + self._s[-1]		# arclength at start of radial section
		self.length_m = self.s_bend_m + self.r_out - self.r_bend
		self.Length: Q_[float] = Q_(self.length_m, 'm')

	@property
	def ControlPoints(self) -> Q_[np.ndarray]:
		return Q_(self.P, 'm')

	def __call__(self, s: float | np.ndarray) -> tuple[Q_, Q_]:
		"""(z, r) at normalized arclength s in [0, 1]"""
		s_true = np.clip(np.atleast_1d(np.asarray(s, float)), 0, 1) * self.length_m
		in_bez = (s_true > self.flat) & (s_true < self.s_bend_m)
		in_rad = s_true >= self.s_bend_m

		z = self.z_in + s_true
		r = np.full_like(s_true, self.r_in)

		t = np.interp(s_true[in_bez] - self.flat, self._s, self._t)
		pts = Bezier.bez(self.P, t)
		z[in_bez] = pts[:, 0]
		r[in_bez] = pts[:, 1]

		z[in_rad] = self.z_out
		r[in_rad] = self.r_bend + s_true[in_rad] - self.s_bend_m

		if np.ndim(s) == 0:
			return Q_(z[0], 'm'), Q_(r[0], 'm')
		return Q_(z, 'm'), Q_(r, 'm')

	def GetRadius(self, z: Q_) -> Q_:
		"""r at axial position z, NaN outside [z_in, z_out]. Returns r_bend at z_out (start of radial section)"""
		z_m = np.atleast_1d(z.to('m').magnitude).astype(float)

		r = np.interp(z_m, self._pts[:, 0], self._pts[:, 1])
		r[z_m <= self.P[0, 0]] = self.r_in
		r[(z_m < self.z_in) | (z_m > self.z_out)] = np.nan

		return Q_(r[0] if np.ndim(z.magnitude) == 0 else r, 'm')
