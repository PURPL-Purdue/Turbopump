# 3D view of the inducer built from the radial-equilibrium blade-angle distribution.
# Name starts with a digit, so run it as a script (like 1DPerformance.py).
import sys
from pathlib import Path
from typing import TYPE_CHECKING

import numpy as np
import plotly.graph_objects as go
from pint import Quantity as Q_

if TYPE_CHECKING:
	from Inducer.RadialEquilibrium import InducerRadialEquilibrium


def GetBladeLines(eq: 'InducerRadialEquilibrium', n_lines: int = 11, resolution: int = 100) -> list[np.ndarray]:
	"""
	Blade camber lines from hub (first) to tip, each a (resolution, 3) array of (z, r, theta) in m, m and rad.
	Along a line at constant span fraction, the wrap follows d(theta) = dz / (r tan(beta_blade)) (NASA TN D-1562 eq. 9, with beta
	from the tangential direction and z as the axial coordinate); the leading edge is straight at theta = 0.
	"""
	z = np.linspace(0.0, eq.L, resolution)
	r_hub = eq.hub_radius(Q_(z, 'm')).to('m').magnitude
	lines = []
	for frac in np.linspace(0.0, 1.0, n_lines):
		r = r_hub + frac * (eq.R_t - r_hub)
		beta = np.radians(eq.BladeAngle(Q_(r, 'm'), Q_(z, 'm')).magnitude)
		integrand = 1.0 / (r * np.tan(beta))
		theta = np.concatenate(([0.0], np.cumsum(0.5 * (integrand[1:] + integrand[:-1]) * np.diff(z))))
		lines.append(np.column_stack((z, r, theta)))
	return lines


def _GridMesh(grid: np.ndarray, offset: int) -> tuple[np.ndarray, np.ndarray]:
	"""Vertices (m*n, 3) and triangle indices (2(m-1)(n-1), 3) of a quad grid of shape (m, n, 3)."""
	m, n = grid.shape[:2]
	idx = np.arange(m * n).reshape(m, n) + offset
	a, b, c, d = idx[:-1, :-1].ravel(), idx[:-1, 1:].ravel(), idx[1:, 1:].ravel(), idx[1:, :-1].ravel()
	return grid.reshape(-1, 3), np.concatenate([np.column_stack((a, b, c)), np.column_stack((a, c, d))])


def _Cartesian(z: np.ndarray, r: np.ndarray, theta: np.ndarray, scale: float) -> np.ndarray:
	return np.stack((r * np.cos(theta), r * np.sin(theta), z), axis=-1) * scale


def PlotInducer3D(eq: 'InducerRadialEquilibrium', l: str = 'mm', n_lines: int = 11, resolution: int = 100,
				  title: str = 'Inducer 3D geometry') -> go.Figure:
	"""All Z_blade blades (thickened by eq.t if set), the tapered hub, and the tip circles at the inlet and outlet."""
	scale = Q_(1, 'm').to(l).magnitude
	lines = GetBladeLines(eq, n_lines, resolution)
	z, r, theta = (np.stack([line[:, i] for line in lines]) for i in range(3))  # each (n_lines, resolution)

	# Tangential half-width of the blade, t / (2 r sin(beta)), with beta the blade angle at each point
	if eq.t is None:
		half = np.zeros_like(theta)
	else:
		beta = np.radians(np.stack([eq.BladeAngle(Q_(r[i], 'm'), Q_(z[i], 'm')).magnitude for i in range(n_lines)]))
		half = eq.t / (2.0 * r * np.sin(beta))

	fig = go.Figure()
	for k in range(eq.Z_blade):
		shift = 2.0 * np.pi * k / eq.Z_blade
		sides = [_Cartesian(z, r, theta + s * half + shift, scale) for s in (-1.0, 1.0)]  # (n_lines, resolution, 3)
		# Two faces plus the hub, tip, leading and trailing edge strips close the blade
		grids = [sides[0], sides[1],
				 np.stack((sides[0][0], sides[1][0])), np.stack((sides[0][-1], sides[1][-1])),
				 np.stack((sides[0][:, 0], sides[1][:, 0])), np.stack((sides[0][:, -1], sides[1][:, -1]))]
		verts, tris, offset = [], [], 0
		for grid in grids:
			v, t = _GridMesh(grid, offset)
			verts.append(v)
			tris.append(t)
			offset += len(v)
		verts, tris = np.concatenate(verts), np.concatenate(tris)
		fig.add_trace(go.Mesh3d(x=verts[:, 0], y=verts[:, 1], z=verts[:, 2], i=tris[:, 0], j=tris[:, 1], k=tris[:, 2],
								color='royalblue', flatshading=True, name='Blades', legendgroup='blades', showlegend=k == 0))

	# Hub surface of revolution
	z_hub = np.linspace(0.0, eq.L, resolution)
	r_hub = eq.hub_radius(Q_(z_hub, 'm')).to('m').magnitude
	phi = np.linspace(0.0, 2.0 * np.pi, 73)
	hub = _Cartesian(z_hub[:, None] * np.ones_like(phi), r_hub[:, None] * np.ones_like(phi), phi[None, :] * np.ones((resolution, 1)), scale)
	fig.add_trace(go.Surface(x=hub[..., 0], y=hub[..., 1], z=hub[..., 2], colorscale=[[0, 'lightgray'], [1, 'lightgray']],
							 showscale=False, name='Hub', showlegend=True, opacity=1.0))

	# Tip circles at the inlet and outlet, as housing reference
	for z_tip in (0.0, eq.L):
		ring = _Cartesian(np.full_like(phi, z_tip), np.full_like(phi, eq.R_t), phi, scale)
		fig.add_trace(go.Scatter3d(x=ring[:, 0], y=ring[:, 1], z=ring[:, 2], mode='lines', line=dict(color='black', width=3),
								   name='Tip', legendgroup='tip', showlegend=z_tip == 0.0))

	fig.update_layout(title=title, template='plotly_white',
					  scene=dict(aspectmode='data', xaxis_title=f'x [{l}]', yaxis_title=f'y [{l}]', zaxis_title=f'z (axial) [{l}]'))
	return fig


if __name__ == '__main__':
	# Put Code/ first and drop this folder so the Inducer package (not Inducer.py) is imported
	_here = Path(__file__).resolve().parent
	sys.path[:] = [p for p in sys.path if Path(p or '.').resolve() != _here]
	sys.path.insert(0, str(_here.parent))

	import Main
	from Inducer.RadialEquilibrium import InducerRadialEquilibrium

	lox_pump = Main.CreatePump(Main.Propellant.LOX)
	ipa_pump = Main.CreatePump(Main.Propellant.IPA)
	lox_pump.AddInducer(Q_(65, 'm'))
	ipa_pump.AddInducer(Q_(80, 'm'))

	for propellant, pump in (('LOX', lox_pump), ('IPA', ipa_pump)):
		eq = InducerRadialEquilibrium(pump.Inducer, blade_thickness=Q_(0.5, 'mm'))
		lines = GetBladeLines(eq)
		print(f"{propellant}: wrap hub = {np.degrees(lines[0][-1, 2]):.1f} deg, tip = {np.degrees(lines[-1][-1, 2]):.1f} deg")
		PlotInducer3D(eq, title=f'{propellant} inducer 3D geometry').show()
