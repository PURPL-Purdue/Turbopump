import numpy as np
from collections.abc import Callable

def GetBetaDistribution(beta_1_rad: float, beta_2_rad: float, x_j: float, x_f: float) -> Callable[[float | np.ndarray], float | np.ndarray]:
	"""
		Creates two parabolas that when pieced together are continuous and tangent. Returns the piecewise curve as a function of ``x``.
		1. ``beta_1_rad`` is the start value of the first parabola
		2. ``beta_2_rad`` is the end value of the second parabola
		3. ``x_j`` specifies the ``x`` value at which the parabolas meet.
		4. ``x_f`` specifies the ``x`` value at which the second parabola ends.

		f1(x) = a1 * x^2 + b1 * x + c1
		f2(x) = a2 * x^2 + b2 * x + c2

		f1'(x) = 2 * a1 * x + b1
		f2'(x) = 2 * a2 * x + b2

		f1'(x_j) = f2'(x_j)
		=> 	(2 * x_j) * a1 - (2 * x_j) * a2 + b1 - b2 = 0

		f1(x_j) = f2(x_j)
		=>	x_j^2 * a1 - x_j^2 * a2 + x_j * b1 - x_j * b2 + c1 - c2 = 0

		f1'(0) = 0
		=> 	b1 = 0

		f1(0) = beta_1
		=>	c1 = beta_1

		f2'(x_f) = 0
		=>	(2 * x_f) * a2 + b2 = 0

		f2(x_f) = beta_2
		=> 	x_f^2 * a2 + x_f * b2 + c2 = beta_2
	"""

	A = np.array([
		[2*x_j, -2*x_j, 1, -1, 0, 0],
		[x_j**2, -x_j**2, x_j, -x_j, 1, -1],
		[0, 0, 1, 0, 0, 0],
		[0, 0, 0, 0, 1, 0],
		[0, 2*x_f, 0, 1, 0, 0],
		[0, x_f**2, 0, x_f, 0, 1]
	], dtype=float)

	b = np.array([0, 0, 0, beta_1_rad, 0, beta_2_rad])

	# Solve Ax = b directly using NumPy's built-in solver
	sol = np.linalg.solve(A, b)
	a_1, a_2, b_1, b_2, c_1, c_2 = sol

	def func(x: float | np.ndarray[float]) -> float | np.ndarray[float]:
		return np.where(
			x < x_j,
		    a_1*x**2 + b_1*x + c_1,
			a_2*x**2 + b_2*x + c_2
			)

	return func

if __name__ == '__main__':
	from plotly import graph_objects as go

	x_f = 10
	Distribution = GetBetaDistribution(beta_1_rad=1, beta_2_rad=10, x_j=3, x_f=x_f)

	x = np.linspace(0, x_f, 100)
	y = Distribution(x)

	fig = go.Figure(go.Scatter(
		x=x, y=y, mode='lines'
	))

	fig.show()