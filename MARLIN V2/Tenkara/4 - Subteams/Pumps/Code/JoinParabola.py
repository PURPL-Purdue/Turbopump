import numpy as np
from collections.abc import Callable

def JoinParabola(y0: float, yf: float, xj: float, xf: float) -> Callable[[float | np.ndarray], float | np.ndarray]:
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

	def func(x: float | np.ndarray[float]) -> float | np.ndarray[float]:
		return np.where(
			x < xj,
		    a_1*x**2 + b_1*x + c_1,
			a_2*x**2 + b_2*x + c_2
			)

	return func

if __name__ == '__main__':
	from plotly import graph_objects as go

	xf = 10
	xj = 3
	func = JoinParabola(y0=1, yf=10, xj=xj, xf=xf)

	x1, x2 = np.linspace(0, xj, 100), np.linspace(xj, xf, 100)
	y1, y2 = func(x1), func(x2)

	fig = go.Figure([go.Scatter(
		x=x1, y=y1, mode='lines', line=dict(color='firebrick', width=2)
		), go.Scatter(
			x=x2, y=y2, mode='lines', line=dict(color='green', width=2)
		)])

	fig.update_layout(template='plotly_white', xaxis_title=r'$x$', yaxis_title=r'$y$')

	fig.show()