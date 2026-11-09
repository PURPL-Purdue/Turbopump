import numpy as np
from warnings import warn
from scipy.integrate import cumulative_trapezoid
from dataclasses import dataclass
from typing import Literal
from plotly import graph_objects as go
from .PotentialFlow import Streamline
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

	def GetTheta(self, radii: np.ndarray, x: np.ndarray) -> np.ndarray:
		return cumulative_trapezoid(1/ (radii * np.tan(self.GetBeta(x))), x=x, initial=0)

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

class CamberLine:

	def __init__(self, radii: np.ndarray, beta_1: float, beta_2: float, x_f: float):
		self.beta_1 = beta_1
		self.beta_2 = beta_2
		self.x_f = x_f
		self.radii = radii

	def GetBetaDistributionFromSj(self, s_j: float) -> tuple[ParabolicBeta, float]:
		x_j = s_j * self.x_f
		coeffs_1, coeffs_2 = _solveParabolas(y0=self.beta_1, yf=self.beta_2, xj=x_j, xf=self.x_f)
		distr = ParabolicBeta(coeffs_1, coeffs_2, x_j)
		
		uh = np.linspace(0, self.x_f, len(self.radii))
		total_wrap = distr.GetTheta(radii=self.radii, x=uh)[-1]
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

def _CamberThetaAt(A: np.ndarray, A_c: np.ndarray, theta_c: np.ndarray, r_c: np.ndarray) -> np.ndarray:
	"""
		Wrap angle of the camber line at abscissa A of the conformal plane (Milos Fig. 30). A_c increases strictly
		(beta in (0, pi/2)). Past either end of the camber, dtheta = dA / r with the radius of that end.
	"""
	theta = np.interp(A, A_c, theta_c)
	theta = np.where(A < A_c[0], theta_c[0] + (A - A_c[0]) / r_c[0], theta)
	return np.where(A > A_c[-1], theta_c[-1] + (A - A_c[-1]) / r_c[-1], theta)

class TransposedBladeLine():

	def __init__(self, camber_length: float, thick_m: float, le_ratio: float=2.0):
		"""
			le_ratio is k_e = a / b of the leading-edge ellipse (Milos eq. 61, Fig. 26): semi-axis a = k_e s/2 along
			the camber line, b = s/2 across it
		"""
		assert le_ratio > 0, 'The leading-edge ellipse needs a positive a/b'
		self.L_max = camber_length
		self.s = thick_m
		self.k_e = le_ratio

	def _GetA_Sides(self, parabolic_beta: ParabolicBeta, resolution: int, le_resolution: int, length_ratio: float|None=None):
		"""
			Suction side, pressure side and camber line in the conformal plane, (A, L) each, in m. The sides are offset
			by the semi-thickness normal to the camber line (Milos eqs. 59-60). The leading edge is a quarter ellipse
			(eq. 61, Fig. 26) on le_resolution points, then resolution points of constant thickness s up to the
			trailing edge, which is cut off normal to the camber line.
		"""

		L_start = (
			0 if (length_ratio is None or length_ratio == 0)
		    else self.L_max * (1 - length_ratio)
		)

		# x_cl of eq. 61 is the arclength along the camber line in the conformal plane, dsigma = dL / sin(beta).
		# sigma >= L - L_start, so the first a of L always reaches sigma = a
		a = self.k_e * self.s / 2
		L_fine = np.linspace(L_start, L_start + a, 100)
		sigma_fine = cumulative_trapezoid(1 / np.sin(parabolic_beta.GetBeta(L_fine)), L_fine, initial=0)

		# even steps of the ellipse angle t (x_M = a cos t) cluster the points at the nose
		x_cl = a * (1 - np.cos(np.linspace(0, np.pi / 2, le_resolution)))
		x_M = self.k_e * self.s / 2 - x_cl
		y_M = np.sqrt((self.k_e * self.s / 2)**2 - x_M**2) / self.k_e	# semi-thickness, eq. 61

		L_le = np.interp(x_cl, sigma_fine, L_fine)
		assert L_le[-1] < self.L_max, 'The leading-edge ellipse is longer than the blade'

		L_i = np.concatenate((L_le, np.linspace(L_le[-1], self.L_max, resolution)[1:]))
		half_thick = np.concatenate((y_M, np.full(resolution - 1, self.s / 2)))

		Beta_i = parabolic_beta.GetBeta(L_i)
		
		sin_B = np.sin(Beta_i)
		cos_B = np.cos(Beta_i) # np.sqrt(1 - sin_B*2) blud is not the trig master

		A_i: np.ndarray = cumulative_trapezoid(1/np.tan(Beta_i), L_i, initial=0)

		if length_ratio:

			pre_res = int(resolution * (1 - length_ratio) + 0.5)
			L_pre = np.linspace(0, L_start, pre_res)
			Beta_pre = parabolic_beta.GetBeta(L_pre)

			A_pre = np.trapezoid(1/np.tan(Beta_pre), L_pre)

			A_i += A_pre
		
		A_s = A_i + half_thick * sin_B
		L_s = L_i - half_thick * cos_B
		A_p = A_i - half_thick * sin_B
		L_p = L_i + half_thick * cos_B

		return A_s, L_s, A_p, L_p, A_i, L_i

	def GetSides(self, parabolic_beta: ParabolicBeta, stream: Streamline, x_0_mean: float, resolution: int, le_resolution: int,
			length_ratio: float|None=None, side_theta: Literal['offset', 'milos']='offset'):
		"""
			(z, r) and wrap angle theta of the suction side, pressure side and camber line, then the conformal-plane
			lines of _GetA_Sides. Each line has le_resolution + resolution - 1 points and starts at theta = 0 on the
			camber leading edge, where the sides meet it. The camber theta is integral of dA / r (Milos eqs. 57-58).
			side_theta sets the theta of the sides:
				'offset':	theta of its own camber point + (A_side - A_camber) / r_camber. The side point lies the
							semi-thickness off the camber point, normal to it, and dA = r dtheta locally
				'milos':	theta of the camber line at the same A (Milos Fig. 30)
		"""
		conformals = self._GetA_Sides(parabolic_beta, resolution, le_resolution, length_ratio)
		A_s, L_s, A_p, L_p, A_i, L_i = conformals

		zr_suction = stream.Points(L_s + x_0_mean)
		zr_pressure = stream.Points(L_p + x_0_mean)
		zr_camber = stream.Points(L_i + x_0_mean)

		radii_camber = zr_camber[:, 1]
		theta_camber = cumulative_trapezoid(1/radii_camber, A_i, initial=0)

		if side_theta == 'offset':
			theta_suction = theta_camber + (A_s - A_i) / radii_camber
			theta_pressure = theta_camber + (A_p - A_i) / radii_camber
		elif side_theta == 'milos':
			theta_suction = _CamberThetaAt(A_s, A_i, theta_camber, radii_camber)
			theta_pressure = _CamberThetaAt(A_p, A_i, theta_camber, radii_camber)
		else:
			raise ValueError(f"side_theta must be 'offset' or 'milos', not {side_theta!r}")

		return (zr_suction, theta_suction, zr_pressure, theta_pressure, zr_camber, theta_camber, conformals)