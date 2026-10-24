# Centrifugal pump impeller and inducer turbomachinery velocity triangle plotting
import numpy as np
import plotly.graph_objects as go
from pint import Quantity as Q_
from typing import Optional

# one colour per vector: its line, arrowhead and label share it
COLORS = dict(
	c="#233ccc",
	u="#14A606",
	w="#E62F0E",
	c_m="#E4940a",
	c_u="#a10da6"
	)

# line dash of each superimposed triangle, in order. '6px,4px' is a short dash, so one always reaches the arrowhead
DASHES = ['solid', '6px,4px', 'dot', 'dashdot']

# the blade direction of a triangle with an incidence: thin, black and long dash-dot, unlike any vector or triangle
BLADE_LINE = dict(width=1, color='black', dash='longdashdot')

# how the angles are measured, printed under the title (see PlotTriangles)
CONVENTION_NOTES = dict(
	centrifugal='Angles from the tangential axis towards the meridional direction (α from +U, β from −U)',
	axial='Angles from the meridional velocity C<sub>m</sub> to C and W, positive towards +U (blade motion)'
	)

class VelocityTriangle:
	"""
		Velocity triangle c = u + w at one station. Angles from the tangent towards +c_m (centrifugal convention):
		alpha from +u to c, beta from -u to w, so tan(beta) = c_m / (u - c_u)
		The optional incidence i is the blade angle minus the flow angle, beta_B = beta + i (Gülich Fig. 3.1), positive
		when the blade is steeper than w. PlotTriangles then draws the blade direction
	"""
	u: Q_[float]
	c: Q_[float]
	c_m: Q_[float]
	c_u: Q_[float]
	w: Q_[float]
	beta: Q_[float]
	alpha: Q_[float]
	notation: str		# HTML after each symbol in the plot labels, e.g. '<sub>2B</sub>' or "<sub>1</sub><sup>'</sup>"
	incidence: Optional[Q_[float]]		# blade angle minus beta, None to draw no blade

	def __init__(self, u: Q_[float], c_m: Q_[float], c_u: Q_[float], notation: str='', incidence: Optional[Q_[float]]=None):
		self.notation = notation
		self.u = u
		self.c_m = c_m
		self.c_u = c_u

		self.beta = np.arctan2(self.c_m, u - c_u)
		self.alpha = np.arctan2(self.c_m, self.c_u)
		self.w = np.sqrt(c_m**2 + (u - self.c_u)**2)
		self.c = np.sqrt(c_m**2 + self.c_u**2)
		self.incidence = incidence

	@classmethod
	def FromBeta(cls, u: Q_[float], c_m: Q_[float], beta: Q_[float], notation: str='', incidence: Optional[Q_[float]]=None) -> 'VelocityTriangle':
		"""Triangle whose relative flow leaves at beta from the tangent: c_u = u - c_m / tan(beta)"""
		return cls(u, c_m, u - c_m / np.tan(beta.to('rad')), notation, incidence)

	def _Vectors(self, unit: str) -> dict[str, tuple[np.ndarray, np.ndarray]]:
		"""
			Start and end of U, C, W and C_m in the plot plane (tangential, meridional), in unit. C_m is the meridional
			leg under the tip of C, left out without swirl, where it is C itself
		"""
		u, c_u, c_m = (v.to(unit).m for v in (self.u, self.c_u, self.c_m))
		O, U_tip, C_tip = np.array([0.0, 0.0]), np.array([u, 0.0]), np.array([c_u, c_m])
		vectors = dict(U=(O, U_tip), C=(O, C_tip), W=(U_tip, C_tip))
		if abs(c_u) > 1e-6 * np.hypot(c_u, c_m):
			vectors['C_m'] = (np.array([c_u, 0.0]), C_tip)
		return vectors

def _Arrow(fig: go.Figure, p0, p1, color: str, dash: str, text: str, group: str, row: int|None, col: int|None):
	# Shaft and head are separate traces. The head is an 'arrow' marker with its tip on p1, turned along the
	# previous segment on screen (size 0 hides the tail). The shaft is drawn from the tip backwards, so its
	# first dash lies under the head, and a dash always reaches the head base (the head is ~11 px long, longer than
	# any gap in DASHES). marker.standoff would pull the head off p1 and leave a gap.
	fig.add_trace(go.Scatter(x=[p1[0], p0[0]], y=[p1[1], p0[1]], mode='lines', legendgroup=group, showlegend=False,
		hovertemplate=f'{text}<extra></extra>', line=dict(width=2 if dash == 'solid' else 1, color=color, dash=dash)),
		row=row, col=col)
	fig.add_trace(go.Scatter(x=[p0[0], p1[0]], y=[p0[1], p1[1]], mode='markers', legendgroup=group, showlegend=False,
		hoverinfo='skip', marker=dict(symbol='arrow', size=[0, 14], angleref='previous', color=color, opacity=1,
		line=dict(width=0))), row=row, col=col)

def _Label(fig: go.Figure, x, y, nx, ny, text: str, color: str, row: int|None, col: int|None, pad: float=6):
	# text sits at (x, y) pushed by pad px along the unit direction (nx, ny). It is centred only beside a nearly level
	# line (within ~6°): half a wide label beside a steeper line would reach back across it
	fig.add_annotation(x=x, y=y, text=text, showarrow=False, font=dict(size=14, color=color),
		xanchor='left' if nx > 0.1 else 'right' if nx < -0.1 else 'center',
		yanchor='bottom' if ny > 0.3 else 'top' if ny < -0.3 else 'middle',
		xshift=pad * nx, yshift=pad * ny, row=row, col=col)

def _Angle(fig: go.Figure, vertex, th_ref: float, th_vec: float, sense: int, R: float, name: str, dash: str, group: str,
		level: int, side: int, row: int|None, col: int|None):
	# Arc at vertex from the reference ray th_ref to the vector th_vec (plot angles, rad). The value is positive in the
	# direction sense (+1 counterclockwise), and printed as line `level` of a stack above the vertex, on the left
	# (side -1) or right (+1)
	th = np.linspace(th_ref, th_vec, 40)
	fig.add_trace(go.Scatter(x=vertex[0] + R * np.cos(th), y=vertex[1] + R * np.sin(th), mode='lines', legendgroup=group,
		line=dict(width=1.2, color='black', dash=dash), showlegend=False, hoverinfo='skip'), row=row, col=col)
	value = np.degrees(sense * (th_vec - th_ref)) + 0.0		# + 0.0 turns -0.0 into 0.0
	fig.add_annotation(x=vertex[0], y=vertex[1], text=f'{name} = {value:.1f}°', showarrow=False, font=dict(size=14),
		xanchor='left' if side > 0 else 'right', yanchor='bottom', xshift=6 * side, yshift=4 + 20 * level,
		row=row, col=col)

def _LabelPlaces(group: list[tuple[np.ndarray, np.ndarray]], defaults: list[np.ndarray], tol: float) -> list[tuple[float, np.ndarray]]:
	"""
		Where to label each vector (p0, p1) in group: the fraction t along it, and the unit normal towards the label side.
		The side faces away from the mean midpoint of the group, or is the default side when the vector has no clear
		side (alone, or crossing the middle). A vector in line with an earlier one takes the opposite side to it.
		Labels sharing a side are spread along their vectors, the outermost nearest the start: vectors from a common
		point are furthest apart near their ends, where the inner labels need the room.
	"""
	mids = [(p0 + p1) / 2 for p0, p1 in group]
	centre = np.mean(mids, axis=0)
	lefts = [np.array([p0[1] - p1[1], p1[0] - p0[0]]) / np.hypot(*(p1 - p0)) for p0, p1 in group]
	offsets = [n @ (mid - centre) for n, mid in zip(lefts, mids)]
	normals = [n * np.sign(offset) if abs(offset) > tol else default for n, offset, default in zip(lefts, offsets, defaults)]

	for k, (q0, q1) in enumerate(group):
		for j, (p0, _) in enumerate(group[:k]):
			if abs(lefts[j] @ (q0 - p0)) < tol and abs(lefts[j] @ (q1 - p0)) < tol:
				normals[k] = -normals[j]

	t = [0.5] * len(group)
	for side in (1, -1):
		on_side = sorted((k for k in range(len(group)) if np.sign(normals[k] @ lefts[k]) == side), key=lambda k: -side * offsets[k])
		for j, k in enumerate(on_side):
			t[k] = (j + 1) / (len(on_side) + 1)
	return list(zip(t, normals))

def PlotTriangles(triangles: list[VelocityTriangle], fig: go.Figure|None=None, row: int|None=None, col: int|None=None,
		names: list[str]|None=None, dashes: list[str]|None=None, unit: str='m/s', convention: str='centrifugal',
		title: str='Velocity Triangles') -> go.Figure:
	"""
		Plots velocity triangles on top of each other, each as the vector sum c = u + w: U and C start at the origin,
		W runs from the tip of U to the tip of C, and C_m is the meridional leg under the tip of C. Colour marks the
		vector (COLORS) and dash the triangle (DASHES), so the triangles may differ in anything: c_m, beta, even u.
		A vector shared by several triangles, such as a common U, is drawn once.

		Each vector's magnitude is printed beside it, labelled with the triangle's notation (e.g. W<sub>1</sub>), on the
		side facing away from the same vector of the other triangles. The angles are stacked above their vertex, the
		first triangle on top. The legend shows the names (or else the notations) once per figure, and clicking an
		entry toggles that triangle in every subplot.

		A triangle with an incidence also gets its blade direction, as in Gülich Fig. 3.1: a line (BLADE_LINE) from the
		vertex of beta at beta + i, which is W turned by i towards the meridional for i > 0, labelled with beta + i at
		its end. Its arc i from W is stacked under the triangle's beta, and a W label on the blade's side goes beyond the
		blade line.

		convention sets how the angles are measured (CONVENTION_NOTES):
		- 'centrifugal': from the tangential axis towards the meridional direction, positive for normal pump flow.
		  alpha from +u at the origin, beta from -u (the blade tangent) at the tip of U.
		- 'axial': from the meridional velocity c_m to c (alpha) and to w (beta), drawn at the tip of C,
		  positive towards +u, the direction of blade motion.
		  alpha = atan(c_u/c_m), beta = atan((c_u-u)/c_m), which is negative when the flow lags the blade.

		Without fig a new figure is made, titled with the convention note under title. With fig, and row and col for a
		subplot of make_subplots, the caller sets the layout, e.g. the note once in the figure title.
	"""
	if convention not in CONVENTION_NOTES:
		raise ValueError(f"convention must be one of {list(CONVENTION_NOTES)}, not {convention!r}")
	axial = convention == 'axial'

	new_fig = fig is None
	if new_fig:
		fig = go.Figure()

	n = len(triangles)
	if dashes is None:
		dashes = [DASHES[i % len(DASHES)] for i in range(n)]
	groups = [name or tri.notation for tri, name in zip(triangles, names or [''] * n)]
	vectors = [tri._Vectors(unit) for tri in triangles]
	scale = max(np.hypot(*(p1 - p0)) for v in vectors for p0, p1 in v.values())
	tol = 1e-6 * scale

	# beta at its vertex, clockwise from the reference ray th_ref to W in both conventions (plot angles, rad), and the
	# blade direction at beta + i, which is clockwise from W by i
	betas = []
	for tri, v in zip(triangles, vectors):
		U_tip, C_tip = v['W']
		if axial:
			vertex, th_ref, (dx, dy) = C_tip, -np.pi / 2, U_tip - C_tip		# W reversed, from the tip of C
		else:
			vertex, th_ref, (dx, dy) = U_tip, np.pi, C_tip - U_tip
		th_w = np.arctan2(dy, dx)
		betas.append((vertex, th_ref, th_w, None if tri.incidence is None else th_w - tri.incidence.to('rad').m))

	# --- legend: one entry per triangle, once per figure ---
	if n > 1 or names is not None:
		shown = {trace.legendgroup for trace in fig.data if trace.showlegend}
		for group, dash in zip(groups, dashes):
			if group not in shown:
				fig.add_trace(go.Scatter(x=[None], y=[None], mode='lines', name=group, legendgroup=group, showlegend=True,
					line=dict(color='black', dash=dash)), row=row, col=col)

	# --- vectors, each kind as a group whose labels keep away from each other ---
	for kind in ('U', 'C', 'W', 'C_m'):
		group, owners = [], []		# each vector once, with the triangles that have it
		for i, v in enumerate(vectors):
			if kind not in v:
				continue
			p0, p1 = v[kind]
			for (q0, q1), sharing in zip(group, owners):
				if np.allclose(p0, q0, atol=tol) and np.allclose(p1, q1, atol=tol):
					sharing.append(i)
					break
			else:
				group.append((p0, p1))
				owners.append([i])

		# without anything to avoid: U below the axis, C_m to the right of its leg, C and W outside their triangle
		defaults = []
		for (p0, p1), (i, *_) in zip(group, owners):
			if kind == 'U':
				defaults.append(np.array([0.0, -1.0]))
			elif kind == 'C_m':
				defaults.append(np.array([1.0, 0.0]))
			else:
				n_out = np.array([p0[1] - p1[1], p1[0] - p0[0]]) / np.hypot(*(p1 - p0))
				centroid = (vectors[i]['U'][1] + vectors[i]['C'][1]) / 3
				defaults.append(n_out if n_out @ ((p0 + p1) / 2 - centroid) > 0 else -n_out)

		for (p0, p1), sharing, (t, normal) in zip(group, owners, _LabelPlaces(group, defaults, 1e-3 * scale)):
			i = sharing[0]
			notation = triangles[i].notation
			if kind == 'C_m':
				sym = f'C<sub>m</sub>{notation}'
			elif kind == 'C' and 'C_m' not in vectors[i]:
				sym = f'C{notation} = C<sub>m</sub>{notation}'
			else:
				sym = kind + notation
			text = f'{sym} = {np.hypot(*(p1 - p0)):.1f} {unit}'
			color = COLORS['c_u' if kind == 'C_m' else kind.lower()]
			# a shared vector stays when one of its triangles is toggled off
			_Arrow(fig, p0, p1, color, dashes[i], text, groups[i] if len(sharing) == 1 else '', row, col)
			at = p0 + t * (p1 - p0)
			# a W label on the side of its blade line goes beyond it, at the same distance from the vertex
			if kind == 'W':
				for vertex, _, _, th_blade in (betas[k] for k in sharing):
					if th_blade is not None:
						on_blade = vertex + np.hypot(*(at - vertex)) * np.array([np.cos(th_blade), np.sin(th_blade)])
						if normal @ (on_blade - at) > 0:
							at = on_blade
			_Label(fig, at[0], at[1], normal[0], normal[1], text, color, row, col)

	# --- angles: nested arcs, values stacked above the vertex, the first triangle on top ---
	if axial:
		R_0 = 0.25 * min(min(np.hypot(*(p1 - p0)) for k, (p0, p1) in v.items() if k != 'U') for v in vectors)
	else:
		R_alpha = 0.2 * min(min(np.hypot(*v['U'][1]), np.hypot(*v['C'][1])) for v in vectors)
		R_beta = 0.2 * min(min(np.hypot(*v['U'][1]), np.hypot(*(v['W'][1] - v['W'][0]))) for v in vectors)
	n_beta = n + sum(tri.incidence is not None for tri in triangles)		# lines in the beta stack: each beta, then its i
	line = 0
	blade_ends = []
	for i, (tri, v, (vertex, th_ref, th_w, th_blade)) in enumerate(zip(triangles, vectors, betas)):
		O = v['U'][0]
		c_u, c_m = v['C'][1]
		grow = 1 + 0.3 * i
		alpha, beta = 'α' + tri.notation, 'β' + tri.notation
		args = (dashes[i], groups[i])
		if axial:
			# At the tip of C, from the C_m leg (pointing down) to C and W reversed, which is the same angle as from
			# +C_m to C and W
			_Angle(fig, vertex, -np.pi / 2, np.arctan2(-c_m, -c_u), -1, grow * R_0, alpha, *args, n - 1 - i, -1, row, col)
			R = grow * R_0
		else:
			_Angle(fig, O, 0, np.arctan2(c_m, c_u), 1, grow * R_alpha, alpha, *args, n - 1 - i, -1, row, col)
			R = grow * R_beta
		_Angle(fig, vertex, th_ref, th_w, -1, R, beta, *args, n_beta - 1 - line, 1, row, col)
		line += 1

		if th_blade is not None:
			# The line overshoots W. Its arc is solid, as it is too short for a dash, and near the far end, where the
			# small angle opens widest
			W = np.hypot(*(v['W'][1] - v['W'][0]))
			end = vertex + 1.15 * W * np.array([np.cos(th_blade), np.sin(th_blade)])
			text = f'blade, {beta} + i{tri.notation} = {np.degrees(th_ref - th_blade):.1f}°'
			fig.add_trace(go.Scatter(x=[vertex[0], end[0]], y=[vertex[1], end[1]], mode='lines', legendgroup=groups[i],
				showlegend=False, hovertemplate=f'{text}<extra></extra>', line=BLADE_LINE), row=row, col=col)
			_Label(fig, end[0], end[1], np.cos(th_blade), np.sin(th_blade), text, 'black', row, col)
			_Angle(fig, vertex, th_w, th_blade, -1, 0.75 * W, 'i' + tri.notation, 'solid', groups[i], n_beta - 1 - line, 1,
				row, col)
			line += 1
			blade_ends.append(end)

	# --- axes through the origin, equal scaling ---
	fig.add_hline(y=0, line_width=0.8, layer='below', row=row, col=col)
	fig.add_vline(x=0, line_width=0.8, layer='below', row=row, col=col)
	xs = [0] + [p[0] for v in vectors for p in (v['U'][1], v['C'][1])] + [end[0] for end in blade_ends]
	pad = 0.2 * (max(xs) - min(xs))		# room on both sides for the angle labels
	fig.update_xaxes(title_text=f'Tangential velocity [{unit}]', range=[min(xs) - pad, max(xs) + pad], row=row, col=col)
	fig.update_yaxes(title_text=f'Meridional velocity [{unit}]', row=row, col=col,
		scaleanchor='x' if row is None else fig.get_subplot(row, col).yaxis.anchor, scaleratio=1)

	if new_fig:
		fig.update_layout(title=f'{title}<br><sup>{CONVENTION_NOTES[convention]}</sup>', template='plotly_white')
		fig.update_xaxes(showgrid=True, gridcolor='rgba(128,128,128,0.2)', zeroline=False)
		fig.update_yaxes(showgrid=True, gridcolor='rgba(128,128,128,0.2)', zeroline=False)

	return fig
