"""
thing
"""
import numpy as np
import plotly.graph_objects as go
from pint import Quantity as Q_

class VelocityTriangle:
	u: Q_[float]
	c: Q_[float]
	c_m: Q_[float]
	c_u: Q_[float]
	w: Q_[float]

	def __init__(self, u: Q_[float], c_m: Q_[float], c_u: Q_[float], w: Q_[float]):
		self.u = u
		self.c_m = c_m
		self.c_u = c_u
		self.w = w
		self.c = np.sqrt(c_m**2 + c_u**2)

	def Plot(self, title='Velocity Triangle', unit='m/s', station=0) -> go.Figure:

		fig = go.Figure()
		c_u, c_m, c = self.c_u.to(unit).magnitude, self.c_m.to(unit).magnitude, self.c.to(unit).magnitude
		u, w = self.u.to(unit).magnitude, self.w.to(unit).magnitude

		# C: absolute velocity, origin -> (C_u, C_m)
		fig.add_trace(go.Scatter(x=[0, c_u], y=[0, c_m], mode='lines', line=dict(width=2),
			name=f'C<sub>{station}</sub> = {c:.1f} {unit}'))

		# U: blade speed, origin -> (U, 0)
		fig.add_trace(go.Scatter(x=[0, u], y=[0, 0], mode='lines', line=dict(width=2),
			name=f'U<sub>{station}</sub> = {u:.1f} {unit}'))

		# W: relative velocity, (C_u, C_m) -> (U, 0)
		fig.add_trace(go.Scatter(x=[c_u, u], y=[c_m, 0], mode='lines', line=dict(width=2),
			name=f'W<sub>{station}</sub> = {w:.1f} {unit}'))

		# C_m: meridional leg, (C_u, 0) -> (C_u, C_m)
		fig.add_trace(go.Scatter(x=[c_u, c_u], y=[0, c_m], mode='lines', line=dict(width=2),
			name=f'C<sub>m{station}</sub> = {c_m:.1f} {unit}'))

		# Labels
		fig.add_annotation(x=c_u / 2, y=c_m / 2, text=f'C<sub>{station}</sub>',
			showarrow=False, font=dict(size=14), xanchor='left')
		fig.add_annotation(x=u / 2, y=-0.04 * c_m, text=f'U<sub>{station}</sub>',
			showarrow=False, font=dict(size=14), xanchor='center')
		fig.add_annotation(x=(u + c_u) / 2, y=c_m / 2, text=f'W<sub>{station}</sub>',
			showarrow=False, font=dict(size=14), xanchor='left')
		fig.add_annotation(x=c_u, y=c_m / 2, text=f'C<sub>m{station}</sub>',
			showarrow=False, font=dict(size=14), xanchor='left', yanchor='middle')

		# Axes through the origin
		fig.add_hline(y=0, line_width=0.8)
		fig.add_vline(x=0, line_width=0.8)

		fig.update_layout(
			title=title,
			xaxis=dict(title=f'Tangential velocity [{unit}]', showgrid=True,
				gridcolor='rgba(128,128,128,0.2)', zeroline=False),
			yaxis=dict(title=f'Meridional velocity [{unit}]', showgrid=True,
				gridcolor='rgba(128,128,128,0.2)', zeroline=False,
				scaleanchor='x', scaleratio=1),  # equal axis scaling
			template='plotly_white',
		)

		return fig