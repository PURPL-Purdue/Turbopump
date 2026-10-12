from . import Bezier
import numpy as np
from pint import Quantity as Q_
import plotly.graph_objects as go

class _LineSegment:
	"""Straight segment from P0 to P1, evaluated by arclength"""

	def __init__(self, P0: np.ndarray, P1: np.ndarray):
		self.P0 = P0
		self.P1 = P1
		self.length = np.linalg.norm(P1 - P0)

	def __call__(self, s: np.ndarray) -> np.ndarray:
		"""(n, 2) points at arclengths s from P0"""
		return self.P0 + (s / self.length)[:, None] * (self.P1 - self.P0)

class _BezierSegment:
	"""Cubic Bezier with control points P, evaluated by arclength"""

	def __init__(self, P: np.ndarray, n_table: int):
		self.P = P

		# lookup table of Bezier parameter t vs arclength, trapezoidal integration of |dB/dt|
		t = np.linspace(0, 1, n_table)
		speed = np.linalg.norm(Bezier.dbez(P, t), axis=1)
		self._t = t
		self._s = np.concatenate(([0], np.cumsum((speed[1:] + speed[:-1]) / 2 * np.diff(t))))
		self.length = self._s[-1]

	def __call__(self, s: np.ndarray) -> np.ndarray:
		"""(n, 2) points at arclengths s from P[0]"""
		return Bezier.bez(self.P, np.interp(s, self._s, self._t))

class SegmentChain:
	"""
		Segments joined end to end, evaluated by arclength along the whole chain. Zero-length segments
		are dropped. All lengths in meters.
	"""

	def __init__(self, segments: list):
		self._segments = [seg for seg in segments if seg.length > 0]
		self._s_start = np.concatenate(([0], np.cumsum([seg.length for seg in self._segments])))

		self.length_m = self._s_start[-1]
		self.Length: Q_[float] = Q_(self.length_m, 'm')

	def Points(self, s_m: np.ndarray) -> np.ndarray:
		"""(n, 2) array of (z, r) at arclengths s_m (n,) in [0, length_m] from the start of the chain"""
		# segment containing each point, the end of the chain belongs to the last one
		i_seg = np.clip(np.searchsorted(self._s_start, s_m, side='right') - 1, 0, len(self._segments) - 1)
		pts = np.empty((s_m.size, 2))
		for i, seg in enumerate(self._segments):
			on_seg = i_seg == i
			pts[on_seg] = seg(s_m[on_seg] - self._s_start[i])
		return pts

	def __call__(self, s: float | np.ndarray) -> tuple[Q_, Q_]:
		"""(z, r) at normalized arclength s in [0, 1]"""
		s_true = np.clip(np.atleast_1d(np.asarray(s, float)), 0, 1) * self.length_m
		pts = self.Points(s_true.ravel())

		z = pts[:, 0].reshape(s_true.shape)
		r = pts[:, 1].reshape(s_true.shape)

		if np.ndim(s) == 0:
			return Q_(z[0], 'm'), Q_(r[0], 'm')
		return Q_(z, 'm'), Q_(r, 'm')

class MeridionalContour(SegmentChain):
	"""
		Meridional contour (z: axial, r: radial) made of tangent-continuous segments:
			1. flat:	straight axial section at r_in, of length flat
			2. bend:	cubic Bezier from axial to the outlet angle, ending at (z_bend, r_bend)
			3. outlet:	straight section at the outlet angle, from (z_bend, r_bend) to (z_out, r_out).
						Absent when r_bend = r_out (default), the bend then ends at the outlet.

		outlet_angle is the angle of the contour from the radial direction at the outlet, positive when
		z increases with r (contour leaning downstream). 0 gives a purely radial outlet.

		The Bezier's end tangents meet at the corner point C, where the outlet tangent line crosses the
		inlet line r = r_in, and both inner control points lie on those tangent lines:
			P1 = P0 + fullness[0] * (C - P0)
			P2 = P3 + fullness[1] * (C - P3)
		fullness -> 1 pulls the curve into the corner (tight bend), fullness -> 0 flattens it
		toward the straight chord P0-P3. ~0.55 approximates a circular arc.

		The contour is evaluated by normalized arclength s in [0, 1] (__call__).
		FlowDomain extends the hub and shroud contours past the inlet and outlet for the FEM solution.
	"""

	def __init__(self, r_in: Q_[float], z_out: Q_[float], r_out: Q_[float], flat: Q_[float],
			outlet_angle: Q_[float] = Q_(0, 'deg'), r_bend: Q_[float] = None,
			fullness: tuple[float, float] = (0.55, 0.55), z_in: Q_[float] = Q_(0, 'm'),
			n_table: int = 2001):

		r_bend = r_bend if r_bend is not None else r_out
		assert r_in < r_bend <= r_out, 'Radii must satisfy r_in < r_bend <= r_out'

		# all lengths stored as magnitudes in meters, angles in radians
		self.z_in = z_in.to('m').magnitude
		self.r_in = r_in.to('m').magnitude
		self.flat = flat.to('m').magnitude
		self.z_out = z_out.to('m').magnitude
		self.r_bend = r_bend.to('m').magnitude
		self.r_out = r_out.to('m').magnitude
		self.outlet_angle = outlet_angle.to('rad').magnitude
		assert abs(self.outlet_angle) < np.pi / 2, 'Outlet angle must be within 90 deg of radial'

		tan_out = np.tan(self.outlet_angle)
		self.z_bend = self.z_out - (self.r_out - self.r_bend) * tan_out

		P0 = np.array([self.z_in + self.flat, self.r_in])
		P3 = np.array([self.z_bend, self.r_bend])
		C = np.array([self.z_bend - (self.r_bend - self.r_in) * tan_out, self.r_in])
		assert C[0] > P0[0], 'Outlet tangent must cross the inlet line after the flat section, reduce flat or outlet_angle'

		self.P = np.array([P0, P0 + fullness[0] * (C - P0), P3 + fullness[1] * (C - P3), P3])

		super().__init__([
			_LineSegment(np.array([self.z_in, self.r_in]), P0),		# flat
			_BezierSegment(self.P, n_table),						# bend
			_LineSegment(P3, np.array([self.z_out, self.r_out]))	# outlet
		])

	@property
	def ControlPoints(self) -> Q_[np.ndarray]:
		return Q_(self.P, 'm')

def QuasiOrthogonalFlow(hub: np.ndarray, shroud: np.ndarray, phi: float | np.ndarray, n_pts: int = 51
		) -> tuple[np.ndarray, np.ndarray]:
	"""
		Continuity through quasi-orthogonals (QOs), the straight lines joining homologous wall points
		hub[i] and shroud[i], given as (n_qo, 2) arrays of (z, r) in meters ordered inlet to outlet.
		Assumes c_m uniform along each QO, no blade blockage. Along each QO, the flow direction t is blended
		linearly from the hub wall tangent to the shroud wall tangent, and continuity through the QO gives
			Q = c_m * integral( (t . n) 2 pi r dl )
		where n is the QO normal. (t . n) corrects for QOs not being perpendicular to the flow.

		Returns
		1. lam: position on each QO (0: hub, 1: shroud) such that the flow fraction phi passes between the
			hub and lam, shape (n_qo,) + shape of phi
		2. area: (n_qo,) effective flow area of each QO, c_m = Q / area
	"""
	n_qo = len(hub)

	# wall flow directions
	t_hub = np.gradient(hub, axis=0)
	t_hub /= np.linalg.norm(t_hub, axis=1, keepdims=True)
	t_shroud = np.gradient(shroud, axis=0)
	t_shroud /= np.linalg.norm(t_shroud, axis=1, keepdims=True)

	chord = shroud - hub
	L = np.linalg.norm(chord, axis=1)
	n = np.column_stack([chord[:, 1], -chord[:, 0]]) / L[:, None]	# QO normal, pointing downstream

	# points along each QO, lam = 0 at hub, 1 at shroud: (n_qo, n_pts, 2)
	lam = np.linspace(0, 1, n_pts)
	pts = hub[:, None, :] + lam[None, :, None] * chord[:, None, :]
	t = (1 - lam[None, :, None]) * t_hub[:, None, :] + lam[None, :, None] * t_shroud[:, None, :]
	t /= np.linalg.norm(t, axis=2, keepdims=True)

	# effective flow area per unit lam, and its cumulative integral from the hub
	dA = np.sum(t * n[:, None, :], axis=2) * 2 * np.pi * pts[:, :, 1] * L[:, None]
	A = np.concatenate((np.zeros((n_qo, 1)), np.cumsum((dA[:, 1:] + dA[:, :-1]) / 2 * np.diff(lam), axis=1)), axis=1)

	# position of the phi flow fraction on each QO
	lam_phi = np.array([np.interp(phi, A_qo / A_qo[-1], lam) for A_qo in A])

	return lam_phi, A[:, -1]

class FlowDomain:
	"""
		Meridional flow domain for the FEM potential flow solution (Milos, Impeller Design Using CAD Techniques
		and Conformal Mapping Method, sections 2-7). The shroud and hub contours are extended by straight lines:
			inlet:	axial, upstream to z_min, inlet_ratio inlet passage heights before the contour inlets.
					The inlet boundary is then a radial line where the flow can be taken as uniform and axial.
					Disturbances from the bend decay upstream roughly like exp(-pi * distance / height).
			outlet:	on at the contour's outlet angle to r_corner = corner_ratio * r_out, then radial out to
					r_max = outlet_ratio * r_out. The outlet boundary is then an axial line where the flow can be
					taken as uniform and radial (Milos uses D_max = 1.25 D_2). A tilted contour has a corner at
					r_corner, where the potential flow velocity is singular. Milos puts it at 1.05 D_2 (point G) so
					this does not distort the impeller outlet velocities. corner_ratio = 1 turns radial at the outlet.

		Along the flow the domain has three zones, ZONES. The extended walls Shroud and Hub are SegmentChains,
		and ZoneBounds[wall] holds the arclengths (m) along each wall of: start, impeller inlet, impeller
		outlet, end.
	"""

	ZONES = ('inlet', 'impeller', 'outlet')

	def __init__(self, shroud: MeridionalContour, hub: MeridionalContour,
			inlet_ratio: float = 2.0, outlet_ratio: float = 1.25, corner_ratio: float = 1.05):

		assert shroud.r_in > hub.r_in, 'Shroud must be above the hub at the inlet'
		assert inlet_ratio > 0 and 1 <= corner_ratio < outlet_ratio, 'Extensions must have a positive length, corner before the end'

		self.z_min = min(shroud.z_in, hub.z_in) - inlet_ratio * (shroud.r_in - hub.r_in)
		self.r_corner = corner_ratio * max(shroud.r_out, hub.r_out)
		self.r_max = outlet_ratio * max(shroud.r_out, hub.r_out)

		self.Shroud = self._Extend(shroud)
		self.Hub = self._Extend(hub)

		# wall segments: inlet extension, the contour's segments, outlet extension
		self.ZoneBounds = {
			'shroud': self.Shroud._s_start[[0, 1, 1 + len(shroud._segments), -1]],
			'hub': self.Hub._s_start[[0, 1, 1 + len(hub._segments), -1]]
		}

	def _Extend(self, contour: MeridionalContour) -> SegmentChain:
		corner = np.array([contour.z_out + (self.r_corner - contour.r_out) * np.tan(contour.outlet_angle), self.r_corner])
		return SegmentChain([
			_LineSegment(np.array([self.z_min, contour.r_in]), np.array([contour.z_in, contour.r_in])),	# axial inlet
			*contour._segments,
			_LineSegment(np.array([contour.z_out, contour.r_out]), corner),		# on at the outlet angle
			_LineSegment(corner, np.array([corner[0], self.r_max]))				# radial outlet
		])

	def Stations(self, n_stream: int) -> tuple[np.ndarray, np.ndarray, dict[str, tuple[int, int]]]:
		"""
			Homologous hub and shroud points, (n_stream + 1, 2) arrays of (z, r), and the first and last station
			of each zone. Within a zone, stations are equally spaced in arclength along each wall, so the zone
			boundaries (impeller inlet and outlet) are QOs. The n_stream intervals are shared between zones in
			proportion to their mean hub and shroud length.
		"""
		lengths = (np.diff(self.ZoneBounds['hub']) + np.diff(self.ZoneBounds['shroud'])) / 2
		n_zone = np.maximum(1, np.round(n_stream * lengths / lengths.sum()).astype(int))
		n_zone[1] += n_stream - n_zone.sum()		# rounding remainder goes to the impeller
		assert np.all(n_zone >= 1), 'n_stream is too small to mesh every zone'

		walls = []
		for wall, bounds in ((self.Hub, self.ZoneBounds['hub']), (self.Shroud, self.ZoneBounds['shroud'])):
			s = [np.linspace(bounds[k], bounds[k + 1], n_zone[k] + 1) for k in range(len(self.ZONES))]
			s = np.concatenate([s[0]] + [zone[1:] for zone in s[1:]])	# neighbouring zones share their boundary station
			walls.append(wall.Points(s))

		first = np.concatenate(([0], np.cumsum(n_zone)))
		zones = {zone: (int(first[k]), int(first[k + 1])) for k, zone in enumerate(self.ZONES)}

		return walls[0], walls[1], zones

	def Mesh(self, n_stream: int = 100, n_span: int = 10, by_flow: bool = True, n_pts: int = 201) -> 'MeridionalMesh':
		"""
			Structured mesh of n_stream x n_span quadrilaterals (Milos section 3, 101 x 11 nodes by default):
			homologous hub and shroud stations (see Stations) are joined by straight QOs, and each QO is divided
			into n_span intervals by the streamlines.
			by_flow:	streamlines at equal flow fractions phi = j / n_span, estimated with QuasiOrthogonalFlow
						(n_pts points per QO). These are exact streamlines on the inlet and outlet boundaries,
						where the flow is uniform, so there the stream function is psi_max * phi.
			otherwise:	QOs divided into equal lengths, as Milos does.
		"""
		hub, shroud, zones = self.Stations(n_stream)
		phi = np.linspace(0, 1, n_span + 1)

		if by_flow:
			lam, _ = QuasiOrthogonalFlow(hub, shroud, phi, n_pts)	# (n_stations, n_lines)
		else:
			lam = np.broadcast_to(phi, (len(hub), n_span + 1))
		grid = hub[:, None, :] + lam[:, :, None] * (shroud - hub)[:, None, :]

		return MeridionalMesh(list(grid.transpose(1, 0, 2)), phi=phi if by_flow else None, zones=zones)

def Polylines(curves: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
	"""(n_curves, n_pts, 2) curves as the x, y of a single plotly trace, with gaps between curves"""
	gaps = np.full((len(curves), 1, 2), np.nan)
	xy = np.concatenate((curves, gaps), axis=1).reshape(-1, 2)
	return xy[:, 0], xy[:, 1]

class MeridionalMesh:
	"""
		Structured mesh of 4-node quadrilateral elements between streamlines, for the FEM solution of the
		meridional stream function and velocity potential (Milos, sections 3-7). All lengths in meters.

		Built from streamlines ordered hub to shroud, each a (n_stations, 2) array of (z, r) sampled at the same
		stations, inlet to outlet. Homologous stations are joined by straight quasi-orthogonals (QOs), so element
		(i, j) lies between QOs i, i + 1 and streamlines j, j + 1.
			grid:		(n_stations, n_lines, 2), grid[i, j] is station i of streamline j
			nodes:		(n_nodes, 2), node i * n_lines + j is grid[i, j]. Numbering along the QOs keeps the
						bandwidth of the global FEM matrix to about n_lines.
			elements:	(n_elements, 4) node indices, counterclockwise in the (z, r) plane:
						(i, j) -> (i + 1, j) -> (i + 1, j + 1) -> (i, j + 1). Element i * (n_lines - 1) + j.
			Boundary:	node indices of the hub (first streamline), shroud (last streamline), inlet (first QO)
						and outlet (last QO), corner nodes belong to two boundaries
			phi:		flow fraction of each streamline (hub 0, shroud 1), None if unknown
			Zones:		first and last QO of each zone, e.g. {'impeller': (i_inlet, i_outlet)}
	"""

	def __init__(self, streamlines: list[np.ndarray], phi: np.ndarray = None, zones: dict[str, tuple[int, int]] = None):

		self.grid = np.stack(streamlines, axis=1)
		self.n_stations, self.n_lines, _ = self.grid.shape
		assert self.n_stations >= 2 and self.n_lines >= 2, 'Need at least 2 streamlines of 2 stations'

		self.nodes = self.grid.reshape(-1, 2)
		node = np.arange(self.n_stations * self.n_lines).reshape(self.n_stations, self.n_lines)
		self.elements = np.column_stack([
			node[:-1, :-1].ravel(), node[1:, :-1].ravel(), node[1:, 1:].ravel(), node[:-1, 1:].ravel()
		])
		self.n_nodes = len(self.nodes)
		self.n_elements = len(self.elements)

		self.Boundary = {'hub': node[:, 0], 'shroud': node[:, -1], 'inlet': node[0], 'outlet': node[-1]}
		self.phi = phi
		self.Zones = zones if zones is not None else {}

		# bilinear elements need a positive Jacobian everywhere: convex, counterclockwise corners
		quad = self.nodes[self.elements]
		to_next = np.roll(quad, -1, axis=1) - quad
		to_prev = np.roll(quad, 1, axis=1) - quad
		corner = to_next[:, :, 0] * to_prev[:, :, 1] - to_next[:, :, 1] * to_prev[:, :, 0]
		assert np.all(corner > 0), 'Elements must be convex with counterclockwise nodes, streamlines or QOs cross'

	@property
	def Areas(self) -> np.ndarray:
		"""(n_elements,) area of each element in the (z, r) plane"""
		quad = self.nodes[self.elements]
		d13 = quad[:, 2] - quad[:, 0]
		d24 = quad[:, 3] - quad[:, 1]
		return (d13[:, 0] * d24[:, 1] - d13[:, 1] * d24[:, 0]) / 2

	@property
	def Centers(self) -> np.ndarray:
		"""(n_elements, 2) element centers (xi = eta = 0), where Milos evaluates the velocities"""
		return self.nodes[self.elements].mean(axis=1)

	@property
	def AspectRatios(self) -> np.ndarray:
		"""(n_elements,) longest over shortest side, Milos keeps the sides of similar length for accuracy"""
		quad = self.nodes[self.elements]
		sides = np.linalg.norm(np.roll(quad, -1, axis=1) - quad, axis=2)
		return sides.max(axis=1) / sides.min(axis=1)

	def ZoneElements(self, zone: str) -> np.ndarray:
		"""Indices of the elements between the first and last QO of a zone"""
		first, last = self.Zones[zone]
		return np.arange(first * (self.n_lines - 1), last * (self.n_lines - 1))

	def Traces(self, l: str = 'in', color: str = 'black', qo_color: str = 'grey', zone_bounds: bool = True
			) -> list[go.Scatter]:
		"""
			Mesh QOs and streamlines (one legend entry), and the zone boundaries (QOs between zones) if zone_bounds,
			as plotly traces in length unit l
		"""
		scale = Q_(1, 'm').to(l).magnitude

		qo_x, qo_y = Polylines(self.grid * scale)
		line_x, line_y = Polylines(self.grid.transpose(1, 0, 2) * scale)
		traces = [
			go.Scatter(x=qo_x, y=qo_y, mode='lines', name='FEM mesh', legendgroup='mesh', showlegend=False,
				line=dict(color=qo_color, width=0.5)),
			go.Scatter(x=line_x, y=line_y, mode='lines', name='FEM mesh', legendgroup='mesh',
				line=dict(color=color, width=1))
		]

		inner_bounds = sorted({i for bounds in self.Zones.values() for i in bounds} - {0, self.n_stations - 1})
		if zone_bounds and inner_bounds:
			x, y = Polylines(self.grid[inner_bounds] * scale)
			traces.append(go.Scatter(x=x, y=y, mode='lines', name='Zone boundaries', line=dict(color='red', width=2)))

		return traces

	def Plot(self, l: str = 'in', color: str = 'black') -> go.Figure:
		"""Streamlines and QOs of the mesh, with the zone boundaries highlighted"""
		fig = go.Figure(self.Traces(l, color))

		fig.update_xaxes(title_text=f"Axial position ({l})")
		fig.update_yaxes(title_text=f"Radial position ({l})", scaleanchor='x', scaleratio=1)
		fig.update_layout(
			title=f"Meridional Mesh: {self.n_elements} elements, {self.n_nodes} nodes",
			template="plotly_white"
		)

		fig.show()

		return fig
