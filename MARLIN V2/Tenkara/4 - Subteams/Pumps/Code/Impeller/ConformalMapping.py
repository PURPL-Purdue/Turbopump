import numpy as np
from collections.abc import Callable
from warnings import warn
from scipy.integrate import cumulative_trapezoid
from dataclasses import dataclass
from plotly import graph_objects as go
from typing import Literal

"""
@Charlie Zhu
Minimal AI
Based on Milo's paper
Conformal mapping for impeller blade design
"""

class ParabolicBeta:
	def __init__(self, coeffs_1: tuple, coeffs_2: tuple, x_j: float):
		self.x_j = x_j
		self.coeffs_1 = coeffs_1
		self.coeffs_2 = coeffs_2

	def GetBeta(self, x: np.ndarray[float]) -> np.ndarray[float]:
		return np.where(
			x < self.x_j,
			np.polyval(self.coeffs_1, x),
			np.polyval(self.coeffs_2, x)
			)

	def GetTheta(self, radii: np.ndarray, arc_x: np.ndarray) -> np.ndarray:
		return cumulative_trapezoid(1/ (radii * np.tan(self.GetBeta(arc_x))), x=arc_x, dx=arc_x[-1] / len(arc_x), initial=0)

def _solveParabolas(y0: float, yf: float, xj: float, xf: float) -> tuple[tuple[float], tuple[float]]:
	"""
		Creates two parabolas that when pieced together are continuous and tangent. Returns the piecewise curve as a function of ``x``.
		1. ``y0`` is the start value of the first parabola
		2. ``yf`` is the end value of the second parabola
		3. ``xj`` specifies the ``x`` value at which the parabolas meet.
		4. ``xf`` specifies the ``x`` value at which the second parabola ends.
	"""

	# f1(x) = a1 * x^2 + b1 * x + c1
	# f2(x) = a2 * x^2 + b2 * x + c2

	# f1'(x) = 2 * a1 * x + b1
	# f2'(x) = 2 * a2 * x + b2

	A = np.array([
		[2*xj,	-2*xj,	1,	-1,  0,  0], # f1'(xj) = f2'(xj) => 	(2 * xj) * a1 - (2 * xj) * a2 + b1 - b2 = 0
		[xj**2,	-xj**2, xj, -xj, 1,	-1], # f1(xj) = f2(xj) =>		xj^2 * a1 - xj^2 * a2 + xj * b1 - xj * b2 + c1 - c2 = 0
		[0,		0,		1,	0, 	 0,	 0], # f1'(0) = 0 =>			b1 = 0
		[0,		0,		0,	0,   1,	 0], # f1(0) = y0 =>			c1 = y0
		[0,		2*xf,	0,	1,   0,  0], # f2'(xf) = 0 =>			(2 * xf) * a2 + b2 = 0
		[0,		xf**2,	0,	xf,  0,  1]	 # f2(xf) = yf =>			xf^2 * a2 + xf * b2 + c2 = yf
	], dtype=float)

	b = np.array([0, 0, 0, y0, 0, yf])	# see right hand side of equations above

	# Solve Ax = b directly using NumPy's built-in solver
	sol = np.linalg.solve(A, b)
	a_1, a_2, b_1, b_2, c_1, c_2 = sol

	return ((a_1, b_1, c_1), (a_2, b_2, c_2))

def GetBetaDistribution(radii: np.ndarray, x_f: float, beta_1: float, beta_2: float, s_j: float, resolution: int|None=None) -> tuple[ParabolicBeta, float]:
	x_j = s_j * x_f
	coeffs_1, coeffs_2 = _solveParabolas(y0=beta_1, yf=beta_2, xj=x_j, xf=x_f)
	distr = ParabolicBeta(coeffs_1, coeffs_2, x_j)
	total_wrap = ParabolicBeta.GetTheta(radii=radii, arc_x=np.linspace(0, x_f, len(radii)))[-1]
	return distr, total_wrap


class CamberLine:

	def __init__(self, radii: np.ndarray, x_f: float, beta_1: float, beta_2: float):
		self.beta_1 = beta_1
		self.beta_2 = beta_2
		self.x_f = x_f
		self.radii = radii

	def GetBetaDistributionFromSj(self, s_j: float) -> tuple[ParabolicBeta, float]:
		x_j = s_j * self.x_f
		coeffs_1, coeffs_2 = _solveParabolas(y0=self.beta_1, yf=self.beta_2, xj=x_j, xf=self.x_f)
		distr = ParabolicBeta(coeffs_1, coeffs_2, x_j)
		total_wrap = distr.GetTheta(radii=self.radii, arc_x=np.linspace(0, self.x_f, len(self.radii)))[-1]
		return distr, total_wrap

	def SolveBetaDistributionFromWrap(self, target_wrap) -> tuple[ParabolicBeta, float]:
		# sj = 0 or 1 makes JoinParabola's system singular, so stay just inside
		guess_high = 0.99
		guess_low = 0.01

		beta_high, wrap_high = self.GetBetaDistributionFromSj(s_j=guess_high)
		beta_low, wrap_low = self.GetBetaDistributionFromSj(s_j=guess_low)

		residual_new = float('inf')
		residual_high, residual_low = wrap_high - target_wrap, wrap_low - target_wrap
		# same sign at both ends: the target is not between them, take the closer end
		if residual_high * residual_low > 0:
			if np.abs(residual_high) > np.abs(residual_low):
				beta, wrap = beta_low, wrap_low
			else:
				beta, wrap = beta_high, wrap_high
			warn(f'The target wrap {np.rad2deg(target_wrap):.1f} deg is out of reach by moving sj alone, resorting to {np.rad2deg(wrap)} deg')
			return beta
		# bisection: keep the half of the bracket where the residual changes sign
		while abs(residual_new) > 1e-4: # approximately 0.006 degrees
			new_guess: float = (guess_high + guess_low) / 2
			beta_dis, total_wrap = self.GetBetaDistributionFromSj(s_j=new_guess)
			residual_new = total_wrap - target_wrap

			if residual_new * residual_low > 0:
				guess_low, residual_low = new_guess, residual_new
			else:
				guess_high, residual_high = new_guess, residual_new
		return beta_dis, total_wrap

class TransposedBladeLine():

	def __init__(self, camber: CamberLine, thick_m: float):
		self.camber = camber
		self.s = thick_m

	def GetA_Sides(self, parabolic_beta: ParabolicBeta, resolution: int=None):

		if resolution is None: resolution = len(self.camber.radii)
		L_arr = np.linspace(0, self.camber.x_f, resolution)

		Beta_i = parabolic_beta.GetBeta(L_arr)

		sin_B = np.sin(Beta_i)
		cos_B = np.cos(Beta_i)

		A_i: np.ndarray = cumulative_trapezoid(1/np.tan(Beta_i), L_arr, initial=0)
		
		A_s = A_i + self.s/2 * sin_B
		L_s = L_arr - self.s/2 * cos_B
		A_p = A_i - self.s/2 * cos_B
		L_p = L_arr + self.s/2 * sin_B

		return A_s, L_s, A_p, L_p