"""
thing
"""
import numpy as np
import plotly.graph_objects as go
from pint import Quantity as Q_

# one colour per vector: its line, arrowhead and label share it
COLORS = dict(
	c="#233ccc",
	u="#14A606",
	w="#E62F0E",
	c_m="#E4940a",
	c_u="#a10da6"
	)

class VelocityTriangle:
	u: Q_[float]
	c: Q_[float]
	c_m: Q_[float]
	c_u: Q_[float]
	w: Q_[float]
	beta_b: Q_[float] | None	# blade (metal) angle from the tangent, as Beta2B; None if not applicable

	def __init__(self, u: Q_[float], c_m: Q_[float], c_u: Q_[float], w: Q_[float], beta_b: Q_[float] | None = None):
		self.u = u
		self.c_m = c_m
		self.c_u = c_u
		self.w = w
		self.beta_b = beta_b
		self.c = np.sqrt(c_m**2 + c_u**2)

	@property
	def alpha(self) -> Q_[float]:
		"""Absolute flow angle, from the tangential axis (+u) towards +c_m (centrifugal convention)."""
		return np.arctan2(self.c_m, self.c_u)

	@property
	def beta(self) -> Q_[float]:
		"""Relative flow angle, from the tangent (-u) towards +c_m (centrifugal convention)."""
		return np.arctan2(self.c_m, self.u - self.c_u)

	@property
	def Ideal(self) -> 'VelocityTriangle | None':
		"""
			Triangle without slip: same u and c_m, but the relative flow follows the blade angle exactly,
			w_u = c_m / tan(beta_b) and w = c_m / sin(beta_b). None if beta_b is not set.
		"""
		if self.beta_b is None:
			return None
		return VelocityTriangle(
			u=self.u, c_m=self.c_m,
			c_u=self.u - self.c_m / np.tan(self.beta_b.to('rad').m),
			w=self.c_m / np.sin(self.beta_b.to('rad').m)
		)

	def Plot(self, title='Velocity Triangle', unit='m/s', station: int|str = '', convention: str = 'centrifugal') -> go.Figure:
		"""
			Plots the velocity triangle as the vector sum c = u + w: U and C start at the origin, W runs from
			the tip of U to the tip of C, and C_m is the dashed meridional component.

			If beta_b is set, the ideal triangle (relative flow along the blade angle, no slip) is overlaid with
			dashed C and W. U and C_m are the same in both, so they are drawn once.

			convention sets how the angles are measured, and is printed under the title:
			- 'centrifugal': from the tangential axis towards the meridional direction, positive for normal pump flow.
			  alpha from +u, beta and beta_b from -u (the blade tangent).
			- 'axial': from the meridional axis (+c_m), positive towards +u, the direction of blade motion.
			  alpha = atan(c_u/c_m), beta = atan((c_u-u)/c_m), which is negative when the flow lags the blade.
			  beta_b is converted from its tangent reference.
		"""
		if convention not in ('centrifugal', 'axial'):
			raise ValueError(f"convention must be 'centrifugal' or 'axial', not {convention!r}")
		axial = convention == 'axial'

		fig = go.Figure()
		c_u, c_m, c = self.c_u.to(unit).magnitude, self.c_m.to(unit).magnitude, self.c.to(unit).magnitude
		u, w = self.u.to(unit).magnitude, self.w.to(unit).magnitude
		ideal = self.Ideal
		if ideal is not None:
			c_u_id, c_id, w_id = ideal.c_u.to(unit).magnitude, ideal.c.to(unit).magnitude, ideal.w.to(unit).magnitude
		st, st_id = f'{station}', (f'{station},id' if station != '' else 'id')

		def Arrow(p0, p1, name, color, dash='solid'):
			# Shaft and head are separate traces. The head is an 'arrow' marker with its tip on p1, turned along the
			# previous segment on screen (size 0 hides the tail). The shaft is drawn from the tip backwards, so its
			# first dash lies under the head, and a dash always reaches the head base (the head is ~11 px long, the
			# pattern is 6 px on / 4 px off). marker.standoff would pull the head off p1 and leave a gap.
			fig.add_trace(go.Scatter(x=[p1[0], p0[0]], y=[p1[1], p0[1]], mode='lines', name=name,
				line=dict(width=(2 if dash=='solid' else 1), color=color, dash='6px,4px' if dash == 'dash' else dash)))
			fig.add_trace(go.Scatter(x=[p0[0], p1[0]], y=[p0[1], p1[1]], mode='markers', showlegend=False, hoverinfo='skip',
				marker=dict(symbol='arrow', size=[0, 14], angleref='previous', color=color, opacity=1, line=dict(width=0))))

		def Label(x, y, nx, ny, text, pad=6):
			# text sits at (x, y) pushed by pad px along the unit direction (nx, ny)
			fig.add_annotation(x=x, y=y, text=text, showarrow=False, font=dict(size=14),
				xanchor='left' if nx > 0.3 else 'right' if nx < -0.3 else 'center',
				yanchor='bottom' if ny > 0.3 else 'top' if ny < -0.3 else 'middle',
				xshift=pad * nx, yshift=pad * ny)

		def SideLabel(p0, p1, centroid, outward, text):
			# label on the outer or inner side of the segment p0-p1, relative to the centroid of its triangle
			d = np.subtract(p1, p0)
			n = np.array([-d[1], d[0]]) / np.hypot(*d)
			mid = np.add(p0, p1) / 2
			if (np.dot(n, mid - centroid) < 0) == outward:
				n = -n
			Label(mid[0], mid[1], n[0], n[1], text)

		def Arc(vertex, th_ref, th_vec, R, dash='solid'):
			# arc from the reference ray to the vector, returns its midpoint for the label
			th = np.linspace(th_ref, th_vec, 40)
			fig.add_trace(go.Scatter(x=vertex[0] + R * np.cos(th), y=vertex[1] + R * np.sin(th), mode='lines',
				line=dict(width=1.2, color='black', dash=dash), showlegend=False, hoverinfo='skip'))
			th_mid = (th_ref + th_vec) / 2
			return vertex[0] + R * np.cos(th_mid), vertex[1] + R * np.sin(th_mid)

		# --- real triangle ---
		O, U_tip, C_tip = (0, 0), (u, 0), (c_u, c_m)
		centroid = np.array([u + c_u, c_m]) / 3
		swirl = abs(c_u) > 1e-6 * c		# without swirl C and C_m are the same vector

		# on the side where both triangles have a vector, the outer one is labelled outside and the inner one inside
		ideal_right = ideal is not None and c_u_id > c_u
		C_out = ideal is None or ideal_right
		W_out = ideal is None or not ideal_right

		c_name = f'C<sub>{st}</sub> = ' if swirl else f'C<sub>{st}</sub> = C<sub>m{st}</sub> = '
		Arrow(O, U_tip, f'U<sub>{st}</sub> = {u:.1f} {unit}', COLORS['u'])
		Arrow(O, C_tip, f'{c_name}{c:.1f} {unit}', COLORS['c'])
		Arrow(U_tip, C_tip, f'W<sub>{st}</sub> = {w:.1f} {unit}', COLORS['w'])
		Label(u / 2, 0, 0, -1, f'U<sub>{st}</sub> = {u:.1f} {unit}')
		SideLabel(O, C_tip, centroid, C_out, f'{c_name}{c:.1f} {unit}')
		SideLabel(U_tip, C_tip, centroid, W_out, f'W<sub>{st}</sub> = {w:.1f} {unit}')
		if swirl:
			Arrow((c_u, 0), C_tip, f'C<sub>m{st}</sub> = {c_m:.1f} {unit}', COLORS['c_u'])
			Label(c_u, c_m / 2, 1, 0, f'C<sub>m{st}</sub> = {c_m:.1f} {unit}')

		# --- ideal triangle (relative flow along the blade angle), dashed on top of the real one ---
		if ideal is not None:
			C_tip_id = (c_u_id, c_m)
			centroid_id = np.array([u + c_u_id, c_m]) / 3
			Arrow(O, C_tip_id, f'C<sub>{st_id}</sub> = {c_id:.1f} {unit}', COLORS['c'], 'dash')
			Arrow(U_tip, C_tip_id, f'W<sub>{st_id}</sub> = {w_id:.1f} {unit}', COLORS['w'], 'dash')
			Arrow((c_u_id, 0), C_tip_id, f'C<sub>m{st},id</sub> = {c_m:.1f} {unit}', COLORS['c_u'], 'dash')
			SideLabel(O, C_tip_id, centroid_id, not C_out, f'C<sub>{st_id}</sub> = {c_id:.1f} {unit}')
			SideLabel(U_tip, C_tip_id, centroid_id, not W_out, f'W<sub>{st_id}</sub> = {w_id:.1f} {unit}')

		# --- angles ---
		# vertex 0 is the origin (alpha), vertex 1 is the tip of U (beta, beta_b). `ref` is the angle of the
		# reference ray in the plot, `sense` is the direction in which the angle is positive
		ref = (np.pi / 2, np.pi / 2) if axial else (0, np.pi)
		sense = (-1, -1) if axial else (1, -1)
		# label offsets in px from the arc midpoint (x right, y up), outside the triangle
		offsets = dict(alpha=(-45, 30), beta=(45, 30), beta_b=(45, 55)) if axial \
			else dict(alpha=(-45, -30), beta=(45, -30), beta_b=(45, 35))
		R = 0.25 * (u - c_u)

		def Angle(k, th_vec, R, name, offset, dash='solid'):
			vertex = O if k == 0 else U_tip
			x, y = Arc(vertex, ref[k], th_vec, R, dash)
			value = np.degrees(sense[k] * (th_vec - ref[k]))
			# ax/ay only apply with showarrow=True, so the offset is a shift of the text itself
			fig.add_annotation(x=x, y=y, xshift=offset[0], yshift=offset[1],
				text=f'{name} = {value:.1f}°', showarrow=False,
				font=dict(size=14))

		Angle(0, np.arctan2(c_m, c_u), 0.25 * c_u, f'α<sub>{st}</sub>', offsets['alpha'])
		Angle(1, np.arctan2(c_m, c_u - u), R, f'β<sub>{st}</sub>', offsets['beta'])
		
		if ideal is not None:
			# the ideal W is along the blade, so its angle is beta_b
			R_ideal = 0.25 * (u - c_u_id)
			Angle(1, np.arctan2(c_m, c_u_id - u), R_ideal, f'β<sub>{st}B</sub>', offsets['beta_b'], 'dash')
		if axial:
			for x0 in (0, u):
				fig.add_trace(go.Scatter(x=[x0, x0], y=[0, 2 * R], mode='lines', showlegend=False, hoverinfo='skip',
					line=dict(width=1, color='gray', dash='dot')))

		# Axes through the origin
		fig.add_hline(y=0, line_width=0.8, layer='below')
		fig.add_vline(x=0, line_width=0.8, layer='below')

		if axial:
			convention_note = 'Angles from the meridional axis, positive towards +U (blade motion)'
		else:
			convention_note = 'Angles from the tangential axis towards the meridional direction (α from +U, β and blade angle from −U)'

		# room on both sides for the angle labels
		xs = [0, u, c_u] + ([c_u_id] if ideal is not None else [])
		pad = 0.2 * (max(xs) - min(xs))

		fig.update_layout(
			title=f'{title}<br><sup>{convention_note}</sup>',
			xaxis=dict(title=f'Tangential velocity [{unit}]', showgrid=True,
				gridcolor='rgba(128,128,128,0.2)', zeroline=False,
				range=[min(xs) - pad, max(xs) + pad]),
			yaxis=dict(title=f'Meridional velocity [{unit}]', showgrid=True,
				gridcolor='rgba(128,128,128,0.2)', zeroline=False,
				scaleanchor='x', scaleratio=1),  # equal axis scaling
			template='plotly_white',
		)
		fig.show()
		return fig
