import numpy as np
from collections.abc import Callable
import warnings
import scipy.integrate
from dataclasses import dataclass

class BetaCMM:
	Get: Callable[[float], float]
	Theta: np.ndarray[float]
	Domain: tuple[float, float]

	def __init__(self, get, wrap, domain):
		self.Get = get
		self.Domain = domain
		self.Theta = wrap
		self.TotalWrap = wrap[-1]

def AofX(coeffs: np.ndarray[float], x: float | np.ndarray) -> float | np.ndarray:
	return coeffs[0] * x**3 + coeffs[1] * x**2 + coeffs[2] * x + coeffs[3]

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
		Returns a, b, c, and d wrapped in a 1D numpy array
		A(x) = a * x^3 + b * x^2 + c * x + d
		A'(x) = 3*a * x^2 + 2*b * x + c
		"""
		x_f = self._arclen
		lhs = np.array([
			[0       ,      0,   0, 1],	# A(x_0 = 0)	= d
			[x_f**3  , x_f**2, x_f, 1], # A(x_f)		= a * x_f^3 + b * x_f^2 + c * x_f + d
			[0       ,      0,   1, 0],	# A'(x0 = 0) 	= c
			[3*x_f**2,  2*x_f,   1, 0]	# A'(x_f) 		= 3*a * x_f^2 + 2 * b * x_f + c
		])

		# A(x_0 = 0) = 0
		# A(x_f) = A_f
		# A'(x_0) = tan(beta_1)
		# A'(x_f) = tan(beta_2)
		rhs = np.array([0, A_f, np.tan(self.beta_1), np.tan(self.beta_2)])
		return np.linalg.solve(lhs, rhs)

	def _GetDistribution(self, A_f: float) -> BetaCMM:

		coeffs = self._SolveCoefficients(A_f)
		A_array = AofX(coeffs, np.linspace(0, self._arclen, len(self._radii)))
		theta = scipy.integrate.cumulative_trapezoid(1/self._radii, A_array, dx=self._arclen / self._resolution, initial=0)

		def GetBeta(x: float) -> float:
			return np.atan(1 / (coeffs[0] * x**2 + coeffs[1] * x + coeffs[2]))

		return BetaCMM(GetBeta, theta, (0, self._arclen))

	def SolveDistribution(self, target_wrap: float=None) -> BetaCMM:

		if target_wrap is None:
			target_wrap = self._target_wrap
		# arcx / A_max = tan(Beta2)
		A_f_max = self._arclen / np.tan(self.beta_2)
		A_f_min = self._arclen / np.tan(self.beta_1)

		beta_high = self._GetDistribution(A_f_max)
		beta_low = self._GetDistribution(A_f_min)

		residual_high = beta_high.TotalWrap - target_wrap
		residual_low = beta_low.TotalWrap - target_wrap

		if residual_high * residual_low > 0:
			
			beta_dis: BetaCMM = beta_low if np.abs(residual_high) > np.abs(residual_low) else beta_high
			warnings.warn(f'Cannot achieve target wrap angle {np.rad2deg(target_wrap):.1f} deg; resorting to {np.rad2deg(beta_dis.TotalWrap):.1f} deg')
			return beta_dis
		
		residual_new = float('inf')
		guess_high, guess_low = A_f_max, A_f_min

		while abs(residual_new) > 1e-4:
			new_guess: float = (guess_high + guess_low) / 2
			beta_dis: BetaCMM = self._GetDistribution(new_guess)
			residual_new = beta_dis.TotalWrap - target_wrap

			if residual_new * residual_low > 0:
				guess_low, residual_low = new_guess, residual_new
			else:
				guess_high, residual_high = new_guess, residual_new

		return beta_dis

		