import numpy as np
from collections.abc import Callable
import warnings
import scipy.integrate
from dataclasses import dataclass
from plotly import graph_objects as go

class BetaCMM:
	Get: Callable[[float], float]
	Theta: np.ndarray[float]
	Domain: tuple[float, float]

	def __init__(self, get, wrap, domain):
		self.Get = get
		self.Domain = domain
		self.Theta = wrap
		self.TotalWrap = wrap[-1]

class ConformalMapping:

	def __init__(self, beta_1: float, beta_2: float, arclen: float, target_wrap: float, radii: np.ndarray) -> Callable[[float | np.ndarray], float | np.ndarray]:
		self.beta_1 = beta_1
		self.beta_2 = beta_2
		self._arclen = arclen
		self._target_wrap = target_wrap
		self._radii = radii
		self._resolution = len(radii)

	def _SolveCoefficients(self, A_f: float) -> np.ndarray[float]:
		"""
		A(x) = a x^3 + b x^2 + c x + d, coefficients highest power first (np.polyval order)
		"""
		x_f = self._arclen
		lhs = np.array([
			[0       ,      0,      0, 1],	# A(0)
			[x_f**3  , x_f**2, x_f**1, 1],	# A(x_f)
			[0       ,      0,      1, 0],	# A'(0)
			[3*x_f**2,  2*x_f,		1, 0]	# A'(x_f)
		])
		# A(0) = 0, A(x_f) = A_f, A' = 1/tan(beta) at both ends
		rhs = np.array([0, A_f, 1/np.tan(self.beta_1), 1/np.tan(self.beta_2)])
		return np.linalg.solve(lhs, rhs)


	def _GetDistribution(self, A_f: float) -> tuple[BetaCMM, np.ndarray]:

		coeffs = self._SolveCoefficients(A_f)
		A_array = np.polyval(coeffs, np.linspace(0, self._arclen, len(self._radii)))
		theta = scipy.integrate.cumulative_trapezoid(1/self._radii, A_array, dx=self._arclen / self._resolution, initial=0)

		def GetBeta(x: float) -> float:
			return np.atan(1 / np.polyval(np.polyder(coeffs), x))

		return BetaCMM(GetBeta, theta, (0, self._arclen)), coeffs

	def SolveDistribution(self, target_wrap: float=None) -> BetaCMM:

		if target_wrap is None:
			target_wrap = self._target_wrap
		c_1, c_2 = 1/np.tan(self.beta_1), 1/np.tan(self.beta_2)
		# beta rises monotonically (A'' <= 0 on the whole blade) only for A_f in the middle third
		A_f_max = self._arclen * (2*c_2 + c_1) / 3	# A''(x_f) = 0: beta flat at the outlet, least wrap
		A_f_min = self._arclen * (c_2 + 2*c_1) / 3	# A''(0) = 0: beta flat at the inlet, most wrap

		beta_high, _ = self._GetDistribution(A_f_max)
		beta_low, _ = self._GetDistribution(A_f_min)

		A_f = A_f_max + (target_wrap - beta_high.TotalWrap) * (A_f_min - A_f_max) / (beta_low.TotalWrap - beta_high.TotalWrap)
		A_f_clamped = np.clip(A_f, min(A_f_max, A_f_min), max(A_f_max, A_f_min))

		beta_dis, coeffs = self._GetDistribution(A_f_clamped)
		if A_f_clamped != A_f:
			warnings.warn(f'Cannot achieve target wrap angle {np.rad2deg(target_wrap):.1f} deg with a monotonic beta; resorting to {np.rad2deg(beta_dis.TotalWrap):.1f} deg')
		
		fig = go.Figure()
		y=np.linspace(0, self._arclen, 100)
		fig.add_trace(go.Scatter(
			x = np.polyval(coeffs, y),
			y = y
		))
		#fig.show()

		return beta_dis

		