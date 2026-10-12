import numpy as np
from pint import Quantity as Q_
import plotly.graph_objects as go
import plotly.subplots as subplots
from plotly.colors import sample_colorscale
import scipy.sparse as sparse
from scipy.sparse.linalg import spsolve
from scipy.interpolate import PchipInterpolator
from MeridionalContour import MeridionalMesh, QuasiOrthogonalFlow, Polylines

# corners of the bilinear quadrilateral in (xi, eta), in element node order, and its 2 x 2 Gauss points (weights 1)
_CORNERS = np.array([[-1, -1], [1, -1], [1, 1], [-1, 1]])
_GAUSS = _CORNERS / np.sqrt(3)

def UniformFlowBoundary(mesh: MeridionalMesh) -> tuple[np.ndarray, np.ndarray]:
	"""
		Boundary nodes and their stream function, normalized by its shroud value so it is the flow fraction
		(Milos eqs. 7-13): 0 on the hub, 1 on the shroud, and on the inlet and outlet the flow is taken as uniform
		and normal to the boundary, so psi grows with integral(r dl) from the hub. That is (r^2 - r_hub^2) on a
		radial inlet line and linear on an axial outlet line.
	"""
	psi = np.full(mesh.n_nodes, np.nan)

	for name in ('inlet', 'outlet'):
		nodes = mesh.Boundary[name]				# hub to shroud
		pts = mesh.nodes[nodes]
		r_dl = np.linalg.norm(np.diff(pts, axis=0), axis=1) * (pts[1:, 1] + pts[:-1, 1]) / 2
		psi[nodes] = np.concatenate(([0], np.cumsum(r_dl))) / r_dl.sum()

	psi[mesh.Boundary['hub']] = 0
	psi[mesh.Boundary['shroud']] = 1

	fixed = np.flatnonzero(~np.isnan(psi))
	return fixed, psi[fixed]

def SolveStreamFunction(mesh: MeridionalMesh, fixed: np.ndarray, values: np.ndarray) -> np.ndarray:
	"""
		Stokes stream function psi (n_nodes,) of axisymmetric irrotational flow, from its values on the fixed nodes
		(Dirichlet problem, Milos eqs. 4 and 14-16). Zero azimuthal vorticity gives
			d2psi/dz2 + d2psi/dr2 - (1/r) dpsi/dr = r div((1/r) grad psi) = 0
		solved in the divergence form, which makes the global matrix symmetric, with bilinear elements and
		2 x 2 Gauss quadrature.
	"""
	X = mesh.nodes[mesh.elements]										# (n_elements, 4, 2) corner (z, r)

	# shape functions N_k = (1 + xi xi_k)(1 + eta eta_k) / 4 and their (xi, eta) derivatives at the Gauss points
	N = (1 + _GAUSS[:, None, 0] * _CORNERS[:, 0]) * (1 + _GAUSS[:, None, 1] * _CORNERS[:, 1]) / 4	# (gauss, corner)
	dN = np.stack([
		_CORNERS[:, 0] * (1 + _GAUSS[:, None, 1] * _CORNERS[:, 1]),
		_CORNERS[:, 1] * (1 + _GAUSS[:, None, 0] * _CORNERS[:, 0])
	], axis=1) / 4														# (gauss, xi or eta, corner)

	J = np.einsum('gak,ekb->egab', dN, X)								# d(z, r) / d(xi, eta)
	grad = np.linalg.solve(J, np.broadcast_to(dN, J.shape[:2] + dN.shape[1:]))	# dN / d(z, r)
	weight = np.linalg.det(J) / np.einsum('gk,ek->eg', N, X[:, :, 1])	# det(J) / r
	K_e = np.einsum('egak,egal,eg->ekl', grad, grad, weight)

	rows = np.broadcast_to(mesh.elements[:, :, None], K_e.shape)
	cols = np.broadcast_to(mesh.elements[:, None, :], K_e.shape)
	K = sparse.coo_matrix((K_e.ravel(), (rows.ravel(), cols.ravel())), shape=(mesh.n_nodes, mesh.n_nodes)).tocsr()

	psi = np.zeros(mesh.n_nodes)
	psi[fixed] = values
	free = np.setdiff1d(np.arange(mesh.n_nodes), fixed)
	psi[free] = spsolve(K[free][:, free].tocsc(), -K[free][:, fixed] @ values)

	return psi

def StreamFunctionGradient(mesh: MeridionalMesh, psi: np.ndarray) -> np.ndarray:
	"""
		(n_nodes, 2) gradient (dpsi/dz, dpsi/dr) at the nodes. psi and the node coordinates are differentiated
		along the mesh lines with the same second order differences (central inside, one-sided on the boundaries)
		and the chain rule maps them to (z, r). On the walls this is more accurate than the element center values
		Milos extrapolates from (eq. 17). The differences are taken against the mean distance between stations
		and between streamlines, and along the streamlines within each zone, where the mesh is smooth, averaging
		the two one-sided values on the stations neighbouring zones share. So they stay second order at the
		impeller inlet and outlet, where the station spacing changes.
	"""
	f = psi.reshape(mesh.n_stations, mesh.n_lines)
	z, r = mesh.grid[..., 0], mesh.grid[..., 1]
	i = np.concatenate(([0], np.cumsum(np.linalg.norm(np.diff(mesh.grid, axis=0), axis=2).mean(axis=1))))
	j = np.concatenate(([0], np.cumsum(np.linalg.norm(np.diff(mesh.grid, axis=1), axis=2).mean(axis=0))))

	def AlongStreamlines(x: np.ndarray) -> np.ndarray:
		dx = np.zeros_like(x)
		count = np.zeros(len(x))
		for first, last in sorted(mesh.Zones.values()) or [(0, mesh.n_stations - 1)]:
			zone = slice(first, last + 1)
			dx[zone] += np.gradient(x[zone], i[zone], axis=0, edge_order=2 if last - first > 1 else 1)
			count[zone] += 1
		return dx / count[:, None]

	f_i, z_i, r_i = (AlongStreamlines(x) for x in (f, z, r))
	f_j, z_j, r_j = (np.gradient(x, j, axis=1, edge_order=2) for x in (f, z, r))

	det = z_i * r_j - z_j * r_i
	return np.stack([(f_i * r_j - f_j * r_i) / det, (z_i * f_j - z_j * f_i) / det], axis=-1).reshape(-1, 2)

class PotentialFlow:
	"""
		Inviscid, incompressible, irrotational axisymmetric flow through the meridional flow domain of a
		MeridionalMesh, by FEM (Milos, Impeller Design Using CAD Techniques and Conformal Mapping Method, sections
		4-7). There are no blades: this is the through-flow infinitely many, infinitely thin blades would see.
			psi:		(n_nodes,) stream function normalized by its shroud value, i.e. the fraction of the flow
						passing between the hub and each node. Streamline phi is the contour psi = phi.
			v_z, v_r:	(n_nodes,) axial and radial velocity in m/s (Milos eq. 6)
							v_z = (Q / 2 pi) (1/r) dpsi/dr,  v_r = -(Q / 2 pi) (1/r) dpsi/dz
			v_m:		(n_nodes,) meridional velocity magnitude in m/s

		Swirl of the incoming flow does not change this solution as long as it is a free vortex (r c_u the same on
		every streamline) with uniform total head. Otherwise the meridional flow has vorticity and psi follows the
		Bragg-Hawthorne equation instead, with r^2 dH/dpsi - C dC/dpsi (C = r c_u) on the right hand side.
	"""

	def __init__(self, mesh: MeridionalMesh, flow_Q: Q_[float]):
		assert mesh.n_stations >= 3 and mesh.n_lines >= 3, 'Need at least 3 nodes along every mesh line'

		self.Mesh = mesh
		self.Q = flow_Q

		self.psi = SolveStreamFunction(mesh, *UniformFlowBoundary(mesh))
		assert np.all(np.diff(self.psi.reshape(mesh.n_stations, mesh.n_lines), axis=1) > 0), \
			'psi must increase from hub to shroud along every QO'

		grad = StreamFunctionGradient(mesh, self.psi) * flow_Q.to('m^3/s').magnitude / (2 * np.pi) / mesh.nodes[:, 1:]
		self.v_z = grad[:, 1]
		self.v_r = -grad[:, 0]
		self.v_m = np.hypot(self.v_z, self.v_r)

	def GetStreamline(self, phi: float) -> tuple[Q_, Q_, Q_, Q_]:
		"""
			Streamline psi = phi (0: hub, 1: shroud), one point per QO, located by monotone cubic interpolation of psi
			along each QO (Milos uses cubic splines), with the velocity components interpolated the same way.
			Returns z, r, arclength from the impeller inlet (negative upstream), and the meridional velocity v_m
		"""
		assert 0 <= phi <= 1, 'phi is a flow fraction, 0 at the hub and 1 at the shroud'
		m = self.Mesh
		psi = self.psi.reshape(m.n_stations, m.n_lines)
		vel = np.column_stack([self.v_z, self.v_r]).reshape(m.n_stations, m.n_lines, 2)

		pts = np.empty((m.n_stations, 2))
		v = np.empty((m.n_stations, 2))
		for i, qo in enumerate(m.grid):
			t = np.concatenate(([0], np.cumsum(np.linalg.norm(np.diff(qo, axis=0), axis=1))))	# distance from the hub
			t_phi = PchipInterpolator(psi[i], t)(phi)
			pts[i] = [np.interp(t_phi, t, qo[:, 0]), np.interp(t_phi, t, qo[:, 1])]
			v[i] = PchipInterpolator(t, vel[i])(t_phi)

		arclen = np.concatenate(([0], np.cumsum(np.linalg.norm(np.diff(pts, axis=0), axis=1))))
		arclen -= arclen[m.Zones.get('impeller', (0, 0))[0]]

		return Q_(pts[:, 0], 'm'), Q_(pts[:, 1], 'm'), Q_(arclen, 'm'), Q_(np.hypot(v[:, 0], v[:, 1]), 'm/s')

	def Plot(self, phis: np.ndarray = np.linspace(0, 1, 11), l: str = 'in') -> go.Figure:
		"""
			Left: meridional velocity over the flow domain, with the FEM mesh it is solved on, the streamlines phis
			and the impeller inlet and outlet.
			Right: meridional velocity along each streamline against arclength from the impeller inlet, colored from
			hub to shroud and bold over the impeller (Milos Fig. 8), and Q / A of each QO (c_m uniform across the
			passage) for comparison.
		"""
		m = self.Mesh
		scale = Q_(1, 'm').to(l).magnitude
		first, last = m.Zones.get('impeller', (0, m.n_stations - 1))
		lines = [self.GetStreamline(phi) for phi in phis]

		fig = subplots.make_subplots(
			rows=1, cols=2,
			column_widths=[0.5, 0.5],
			horizontal_spacing=0.15,
			subplot_titles=("Meridional Velocity Field and FEM Mesh", "Meridional Velocity Along Streamlines")
		)

		# ----- velocity field, subplot (1, 1) -----
		# carpet over the structured mesh, flattened with the station index a varying fastest
		zr = m.grid.transpose(1, 0, 2).reshape(-1, 2) * scale
		a = np.tile(np.arange(m.n_stations), m.n_lines)
		b = np.repeat(np.arange(m.n_lines), m.n_stations)
		carpet_axis = dict(smoothing=0, showgrid=False, minorgridcount=0, showticklabels='none')
		fig.add_trace(go.Carpet(
			a=a, b=b, x=zr[:, 0], y=zr[:, 1], carpet='domain', aaxis=carpet_axis, baxis=carpet_axis
		), row=1, col=1)
		fig.add_trace(go.Contourcarpet(
			a=a, b=b, z=self.v_m.reshape(m.n_stations, m.n_lines).T.ravel(), carpet='domain',
			ncontours=40, contours=dict(showlines=False), colorscale='Jet',
			colorbar=dict(title='v_m (m/s)', x=fig.layout.xaxis.domain[1] + 0.01, len=0.9, thickness=15)
		), row=1, col=1)

		# mesh over the velocity field, translucent so the colors show through
		mesh_color = 'rgba(80, 80, 80, 0.45)'
		for trace in m.Traces(l, color=mesh_color, qo_color=mesh_color, zone_bounds=False):
			fig.add_trace(trace, row=1, col=1)

		x, y = Polylines(np.stack([np.column_stack([z.to(l).m, r.to(l).m]) for z, r, _, _ in lines]))
		fig.add_trace(go.Scatter(
			x=x, y=y,
			mode='lines', name='Streamlines', line=dict(color='black', width=1),
			visible="legendonly"
		), row=1, col=1)
		x, y = Polylines(m.grid[[first, last]] * scale)
		fig.add_trace(go.Scatter(
			x=x, y=y, mode='lines', name='Impeller inlet and outlet', line=dict(color='red', dash='dash', width=2)
		), row=1, col=1)

		fig.update_xaxes(title_text=f"Meridional position ({l})", row=1, col=1)
		fig.update_yaxes(title_text=f"Radial position ({l})", scaleanchor='x', scaleratio=1, row=1, col=1)

		# ----- velocity along the streamlines, subplot (1, 2) -----
		colors = sample_colorscale('Viridis', np.linspace(0, 0.85, len(phis)))
		for phi, (_, _, s, v_m), color in zip(phis, lines, colors):
			s, v_m = s.to(l).m, v_m.to('m/s').m
			name = f"ψ = {phi:.2f}" + (" (hub)" if phi == 0 else " (shroud)" if phi == 1 else "")
			hover = f"ψ = {phi:g}<br>s = %{{x:.3f}} {l}<br>v_m = %{{y:.2f}} m/s<extra></extra>"

			fig.add_trace(go.Scatter(
				x=s, y=v_m, mode='lines', name=name, legendgroup=name, showlegend=False, hovertemplate=hover,
				line=dict(color=color, width=1), opacity=0.5
			), row=1, col=2)
			fig.add_trace(go.Scatter(
				x=s[first:last + 1], y=v_m[first:last + 1], mode='lines', name=name, legendgroup=name,
				hovertemplate=hover, line=dict(color=color, width=2)
			), row=1, col=2)

		_, area = QuasiOrthogonalFlow(m.grid[:, 0], m.grid[:, -1], 0.5)
		_, _, s_mean, _ = self.GetStreamline(0.5)
		fig.add_trace(go.Scatter(
			x=s_mean.to(l).m, y=self.Q.to('m^3/s').m / area, mode='lines',
			name='Q / QO area (uniform c_m)',
			line=dict(color='black', dash='dashdot', width=2)
		), row=1, col=2)
		fig.add_vline(x=0, line=dict(color='grey', dash='dot'), row=1, col=2)

		fig.update_xaxes(title_text=f"Streamline arclength from impeller inlet ({l})", showgrid=True,
			gridcolor='rgba(128,128,128,0.3)', row=1, col=2)
		fig.update_yaxes(title_text="Meridional velocity v_m (m/s)", rangemode='tozero', showgrid=True,
			gridcolor='rgba(128,128,128,0.3)', row=1, col=2)

		# a separate legend for each subplot, see layout below
		fig.update_traces(legend='legend2', row=1, col=2)

		fig.update_layout(
			title=f"Meridional Potential Flow at Q = {self.Q.to('L/s').m:.2f} L/s ({m.n_elements} elements)",
			# field legend in its empty upper-left corner, streamline legend right of its plot
			legend=dict(x=0.01, y=0.99, xanchor='left', yanchor='top', bgcolor='rgba(255,255,255,0.8)'),
			legend2=dict(x=1.02, y=1, xanchor='left', yanchor='top'),
			template="plotly_white"
		)

		fig.show()

		return fig
