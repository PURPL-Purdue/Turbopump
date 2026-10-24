import numpy as np
from pint import Quantity as Q_
from plotly import graph_objects as go, subplots
from plotly.colors import sample_colorscale
from InducerRadialEquilibrium import InducerRadialEquilibrium
from MeridionalDesign import MeridionalDesign
from MeridionalContour import MeridionalCurve, Polylines
from PotentialFlow import Streamline
from typing import Optional
from VelocityTriangle import VelocityTriangle, PlotTriangles, CONVENTION_NOTES
from dataclasses import dataclass

"""Notes for myself
	- Large wrap angle θ = ∫ dm / (r tan β) largely caused by low Beta1 and Beta2
	- Can also decrease total meridional length by decreasing axial length
	- Use blade solidity as a check for blade length
"""

g = Q_(9.81, 'm/s^2')

@dataclass
class _TriInfo:
	c_m: Q_
	Beta: Q_

@dataclass
class BladeLine:
	"""Mean blade line (camber line) on streamline phi, stations from the leading edge to the trailing edge"""
	phi: float
	m: Q_			# meridional length from the leading edge
	z: Q_
	r: Q_
	Beta: Q_		# blade angle from the tangent
	Theta: Q_		# angle ahead of the trailing edge, in the direction of rotation
	v_m: Q_			# meridional velocity of the potential flow (no blades)

	@property
	def Wrap(self) -> Q_:
		"""Wrap angle: how far around the axis the blade line turns from the leading to the trailing edge"""
		return self.Theta[0]

@dataclass
class BladeLoading:
	"""Relative velocity and static pressure along a BladeLine: mean blade surface, suction side, pressure side"""
	w: Q_
	w_ss: Q_
	w_ps: Q_
	p: Q_
	p_ss: Q_
	p_ps: Q_

class ImpellerBladeDesign:

	def __init__(self, meridional_design: MeridionalDesign, inducer_radial_equilibrium: Optional[InducerRadialEquilibrium]=None):

		self.Pump = meridional_design.Pump
		self._mer = meridional_design
		self._rad_eq = inducer_radial_equilibrium
		self.K_inlet = inducer_radial_equilibrium.GetOutletAngularMomentum() if inducer_radial_equilibrium is not None else Q_(0, 'm^2/s')
		
		# it's better to have slight positive incidence than slight negative incidence
		self._INCIDENCE_PADDING = np.pi/180 * 2	# degrees of extra angle

	@property
	def Meridional(self) -> MeridionalDesign:
		"""Meridional design the blades are built on: contours, potential flow, leading and trailing edges"""
		return self._mer

	def GetInletC_u(self, r: Q_):
		return self.K_inlet / r

	def _SolveInletFlow(self, phi: float) -> tuple[_TriInfo, _TriInfo, Q_, float]:
		"""
			Radius r_1 (m), flow angle beta_1' just inside the blades (rad, = blade angle at zero incidence) and degree of
			blockage rho_1 where the leading edge crosses streamline phi. Tangential blockage (CFturbo; Hassan et al.,
			Sci. Rep. 2024, eq. 17) with the meridional leading edge angle lambda_1, which slightly overestimates it:
				rho_1 = 1 - a / sin(beta_1'),		a = z s / (pi d_1 sin(lambda_1)),	s = blade_thick
				tan(beta_1') = tan(beta_1) / rho_1		(Milos eq. 33)
			Solved explicitly (own derivation): sin(beta_1' - beta_1) = a cos(beta_1)
		"""
		streamline = self._mer.PotentialFlow.GetStreamline(phi)

		# leading edge angle to the exact contours at the hub and shroud, to the polyline streamline in between
		if phi == 0:
			lambda_1 = self._GetLeadingEdgeAngle(self._mer.HubContour)
		elif phi == 1:
			lambda_1 = self._GetLeadingEdgeAngle(self._mer.ShroudContour)
		else:
			lambda_1 = self._GetLeadingEdgeAngle(streamline)

		r_1, c_m1, w_u1 = self._GetLeadingEdgeFlow(streamline)
		impeller = self._mer.Pump.Impeller

		# blockage coefficient: blade thickness over blade pitch, widened by the leading edge sweep
		pitch = 2 * np.pi * r_1 / impeller.Z_blade
		a: Q_[float] = impeller.blade_thick / (pitch * np.sin(lambda_1.to('rad').m))
		assert a < 1, 'The blades fill the inlet circumference, reduce the blade number or thickness'

		beta_1: Q_ = np.arctan2(c_m1, w_u1)							# without blockage
		beta_1_prime = beta_1 + np.arcsin(a * np.cos(beta_1))		# with blockage
		blockage_factor: float = (1 - a / np.sin(beta_1_prime)).to('dimensionless').m
		c_m1_prime: Q_[float] = c_m1 / blockage_factor

		return _TriInfo(c_m1_prime, beta_1_prime), _TriInfo(c_m1, beta_1), r_1, blockage_factor

	def _GetLeadingEdgeFlow(self, streamline: Streamline) -> tuple[Q_, Q_, Q_]:
		"""
			Radius r_1 (m), meridional velocity c_m1 and relative swirl w_u1 = u_1 - c_u1 (m/s) where the leading edge
			crosses streamline, without blade blockage
		"""
		s = self._GetEdgeCrossing(streamline, self._mer.LeadingEdge)
		_, r_1 = streamline(s)

		omega = self._mer.Pump.Impeller.DP.N_shaft.to('rad/s')
		c_m1: Q_[float] = streamline.Interpolate(streamline.v_m, s).to('m/s')
		w_u1 = omega * r_1 - self.GetInletC_u(r_1).to('m/s')
		assert w_u1 > 0, 'The preswirl exceeds the blade speed at the leading edge'

		return r_1, c_m1, w_u1

	def _GetEdgeCrossing(self, streamline: Streamline, edge: MeridionalCurve) -> float:
		"""Normalized arclength along streamline where it crosses edge, the leading or trailing edge"""
		hits = streamline.FindIntersect(edge)
		assert len(hits) == 1, 'The edge must cross the streamline exactly once'
		return hits[0][2]

	def _GetLeadingEdgeAngle(self, curve: MeridionalCurve) -> Q_[float]:
		"""
			Acute angle in the meridional plane between the leading edge and curve (the hub or shroud contour, or a
			streamline) where they meet, from the unit tangents of both: arccos|t_LE . t_curve|
		"""
		hits = curve.FindIntersect(self._mer.LeadingEdge)
		assert len(hits) == 1, 'The leading edge must meet the curve exactly once'
		_, _, s_curve, s_LE = hits[0]
		cos_angle = abs(curve.Tangent(s_curve) @ self._mer.LeadingEdge.Tangent(s_LE))
		return Q_(np.arccos(min(cos_angle, 1.0)), 'rad').to('deg')

	def PlotInletTriangles(self, phis: list[float]=[0, 0.5, 1]) -> go.Figure:
		"""
			Inlet triangles where the leading edge crosses each streamline phi: upstream of the blades (1), just inside
			them with the blade blockage (1'), and along the blade angle (1B), whose difference from 1' is the incidence
		"""
		# at most 3 rows, further triangles fill new columns
		n_rows = min(len(phis), 3)
		n_cols = -(-len(phis) // n_rows)

		# plotly takes the titles row by row, the triangles fill column by column
		titles = [''] * (n_rows * n_cols)
		for i, phi in enumerate(phis):
			titles[i % n_rows * n_cols + i // n_rows] = f'Flow Fraction ϕ = {phi}'

		fig = subplots.make_subplots(rows=n_rows, cols=n_cols, subplot_titles=titles)
		for i, phi in enumerate(phis):
			prime, no_block, r, _ = self._SolveInletFlow(phi)
			u = self.Pump.DP.N_shaft * r
			triangles = [
				VelocityTriangle.FromBeta(u, no_block.c_m, no_block.Beta, '<sub>1</sub>'),
				VelocityTriangle.FromBeta(u, prime.c_m, prime.Beta, "<sub>1</sub>'"),
				VelocityTriangle.FromBeta(u, prime.c_m, prime.Beta + Q_(self._INCIDENCE_PADDING, 'rad'), '<sub>1B</sub>'),
				]
			PlotTriangles(triangles, fig, row=i % n_rows + 1, col=i // n_rows + 1,
				names=['upstream', 'inside blades', 'blade angle'])

		fig.update_layout(template='simple_white',
			title=f'Inlet Velocity Triangles<br><sup>{CONVENTION_NOTES["centrifugal"]}</sup>')

		return fig

	def GetBladeLine(self, phi: float, n: int=201) -> BladeLine:
		"""
			Mean blade line on streamline phi, at n stations equally spaced in meridional length m from the leading to
			the trailing edge. The blade angle goes linearly in m from beta_1B at the leading edge (beta_1' plus the
			incidence padding, see _SolveInletFlow) to Beta2B at the trailing edge.
			Along a streamline the blade turns r dtheta = dm / tan(beta_B) around the axis (NASA TN D-1562, eq. 9, whose
			beta is from the meridional direction), so its angle ahead of the trailing edge, in the direction of rotation, is
				theta(m) = integral from m to the trailing edge of dm / (r tan(beta_B))
			by the trapezoid rule. theta = 0 at the trailing edge on every streamline, so the trailing edge lies in an
			axial plane, and the leading edge is ahead of it since the blades sweep back.
		"""
		streamline = self._mer.PotentialFlow.GetStreamline(phi)
		s_LE, s_TE = (self._GetEdgeCrossing(streamline, edge) for edge in (self._mer.LeadingEdge, self._mer.TrailingEdge))
		s = np.linspace(s_LE, s_TE, n)
		zr = streamline.Points(s * streamline.length_m)
		m = (s - s_LE) * streamline.length_m

		beta_1B = self._SolveInletFlow(phi)[0].Beta.to('rad').m + self._INCIDENCE_PADDING
		beta = np.interp(m, [0, m[-1]], [beta_1B, self.Pump.Impeller.Beta2B.to('rad').m])

		dtheta_dm = 1 / (zr[:, 1] * np.tan(beta))
		theta = np.concatenate(([0], np.cumsum((dtheta_dm[1:] + dtheta_dm[:-1]) / 2 * np.diff(m))))	# from the LE

		return BladeLine(phi, Q_(m, 'm'), Q_(zr[:, 0], 'm'), Q_(zr[:, 1], 'm'), Q_(beta, 'rad'),
			Q_(theta[-1] - theta, 'rad'), streamline.Interpolate(streamline.v_m, s))

	@property
	def InletTotalPressure(self) -> Q_[float]:
		"""Total pressure ahead of the impeller, the same from hub to shroud: feed pressure plus the inducer head"""
		H_inducer = self._rad_eq.inducer.DP.H_rise if self._rad_eq is not None else Q_(0, 'm')
		return (self.Pump.DP.feed_pressure + self.Pump.DP.rho * g * H_inducer).to('psi')

	def GetBladeLoading(self, line: BladeLine) -> BladeLoading:
		"""
			Velocities and static pressures on both blade faces along line, with the static pressure taken to vary
			linearly from blade to blade (NASA TN D-1562, App. C). The flow follows the blades exactly (no incidence, no
			slip) and is inviscid, so it is least accurate near both edges and has no losses.
				t_u = t / sin(beta_B)								circumferential blade thickness (eq. 10, no lean)
				b_u = 2 pi r / Z - t_u								open pitch between two blades
				w_m = v_m (2 pi r / Z) / b_u							potential flow v_m with the blade blockage
				w_u = w_m / tan(beta_B),	c_u = u - w_u
				dp / rho = b_u w_m (1/r) d(r c_u)/dm					pressure side minus suction side (eq. C3, rewritten)
				w_ps = sqrt(w^2 - dp / rho),	w_ss = sqrt(w^2 + dp / rho)		(eqs. C6, C7)
			The blades must add angular momentum r c_u to the flow, and the pressure difference across them is the force
			that does it. A negative w^2 - dp / rho means an eddy on the pressure side: w_ps is then shown negative, as in
			TN D-1562 (p. 34).
			Static pressure on the mean surface from the rothalpy p / rho + (w^2 - u^2) / 2, the same everywhere (eq. B22),
			with the inlet total pressure and swirl r c_u = K_inlet. p_ps, p_ss = p +- dp / 2 (eqs. C4, C5).
		"""
		impeller = self.Pump.Impeller
		Z = impeller.Z_blade
		omega = impeller.DP.N_shaft.to('rad/s').m
		rho = self.Pump.DP.rho.to('kg/m^3').m
		m, r, beta, v_m = line.m.to('m').m, line.r.to('m').m, line.Beta.to('rad').m, line.v_m.to('m/s').m

		b_u = 2 * np.pi * r / Z - impeller.blade_thick.to('m').m / np.sin(beta)
		w_m = v_m * (2 * np.pi * r / Z) / b_u
		w_u = w_m / np.tan(beta)
		w = np.hypot(w_m, w_u)
		c_u = omega * r - w_u
		dp_rho = b_u * w_m * np.gradient(r * c_u, m) / r

		def SignedSqrt(x: np.ndarray) -> Q_:
			return Q_(np.sign(x) * np.sqrt(np.abs(x)), 'm/s')

		rothalpy = self.InletTotalPressure.to('Pa').m / rho - omega * self.K_inlet.to('m^2/s').m
		p = rho * (rothalpy + ((omega * r)**2 - w**2) / 2)

		return BladeLoading(
			w=Q_(w, 'm/s'), w_ss=SignedSqrt(w**2 + dp_rho), w_ps=SignedSqrt(w**2 - dp_rho),
			p=Q_(p, 'Pa').to('psi'), p_ss=Q_(p - rho * dp_rho / 2, 'Pa').to('psi'), p_ps=Q_(p + rho * dp_rho / 2, 'Pa').to('psi')
		)

	def PlotPlanView(self, phis: list[float]=[0, 0.5, 1], l: str='in') -> go.Figure:
		"""
			Plan view, looking along the axis: the mean blade lines (GetBladeLine) on streamlines phis for all the blades,
			with the hub, inlet and outlet circles. Rotation is counterclockwise.
		"""
		impeller = self.Pump.Impeller
		Z = impeller.Z_blade
		fig = go.Figure()

		circle = np.linspace(0, 2 * np.pi, 361)
		for d, name in ((impeller.d_hub, 'Hub d<sub>hub</sub>'), (impeller.d_1, 'Inlet d<sub>1</sub>'), (impeller.d_2, 'Outlet d<sub>2</sub>')):
			R = d.to(l).m / 2
			fig.add_trace(go.Scatter(
				x=R * np.cos(circle), y=R * np.sin(circle), mode='lines', name=name, hoverinfo='skip',
				line=dict(color='grey', dash='dot', width=1)
			))

		colors = sample_colorscale('Viridis', np.linspace(0, 0.85, len(phis)))
		for phi, color in zip(phis, colors):
			line = self.GetBladeLine(phi)
			r, theta = line.r.to(l).m, line.Theta.to('rad').m
			blades = np.stack([
				np.column_stack([r * np.cos(theta + 2 * np.pi * k / Z), r * np.sin(theta + 2 * np.pi * k / Z)]) for k in range(Z)
				])
			x, y = Polylines(blades)
			name = f'ϕ = {phi:g}' + (' (hub)' if phi == 0 else ' (shroud)' if phi == 1 else '') + \
				f', wrap {line.Wrap.to("deg").m:.0f}°'
			fig.add_trace(go.Scatter(x=x, y=y, mode='lines', name=name, line=dict(color=color, width=2)))

		fig.update_xaxes(title_text=f'x ({l})')
		fig.update_yaxes(title_text=f'y ({l})', scaleanchor='x', scaleratio=1)
		fig.update_layout(
			title=f'Impeller Plan View, {Z} Blades<br><sup>β<sub>B</sub> linear in meridional length from the leading '
				f'to the trailing edge. The blades overlap when the wrap exceeds 360°/Z = {360 / Z:.0f}°. '
				'Rotation counterclockwise.</sup>',
			template='plotly_white'
		)

		return fig

	def PlotBladeLoading(self, phis: list[float]=[0, 0.5, 1]) -> go.Figure:
		"""
			Blade loading (GetBladeLoading) along streamlines phis, one column each, against the meridional length from
			the leading edge as a fraction of the blade length. Top: relative velocities. Bottom: static pressures, with
			the vapor pressure.
		"""
		fig = subplots.make_subplots(
			rows=2, cols=len(phis), shared_xaxes=True, shared_yaxes='rows', vertical_spacing=0.08,
			subplot_titles=[f'Flow Fraction ϕ = {phi:g}' for phi in phis]
		)
		styles = {
			'Mean blade surface': dict(color='black', dash='dashdot', width=1),
			'Suction side': dict(color='#d62728', width=2),
			'Pressure side': dict(color='#1f77b4', dash='dash', width=2),
		}
		p_vapor = self.Pump.DP.inlet_vapor_pressure.to('psi').m

		for col, phi in enumerate(phis, start=1):
			line = self.GetBladeLine(phi)
			load = self.GetBladeLoading(line)
			m_bar = (line.m / line.m[-1]).m

			for (name, style), w, p in zip(styles.items(), (load.w, load.w_ss, load.w_ps), (load.p, load.p_ss, load.p_ps)):
				fig.add_trace(go.Scatter(
					x=m_bar, y=w.to('m/s').m, mode='lines', name=name, legendgroup=name, showlegend=col == 1, line=style,
					hovertemplate=f'{name}<br>m/m<sub>TE</sub> = %{{x:.2f}}<br>w = %{{y:.1f}} m/s<extra></extra>'
				), row=1, col=col)
				fig.add_trace(go.Scatter(
					x=m_bar, y=p.to('psi').m, mode='lines', name=name, legendgroup=name, showlegend=False, line=style,
					hovertemplate=f'{name}<br>m/m<sub>TE</sub> = %{{x:.2f}}<br>p = %{{y:.0f}} psi<extra></extra>'
				), row=2, col=col)

			fig.add_hline(y=0, line=dict(color='grey', width=1), row=1, col=col)
			fig.add_hline(y=p_vapor, line=dict(color='grey', dash='dot'), row=2, col=col,
				annotation_text='Vapor pressure' if col == 1 else None, annotation_position='bottom right')
			fig.update_xaxes(title_text='Meridional length from LE, m / m<sub>TE</sub>', row=2, col=col)

		fig.update_yaxes(title_text='Relative velocity w (m/s)', row=1, col=1)
		fig.update_yaxes(title_text='Static pressure (psi)', row=2, col=1)
		fig.update_layout(
			title=f'Blade Loading (NASA TN D-1562, App. C)<br><sup>Inviscid, the flow follows the blades (no incidence, '
				f'no slip): unreliable near both edges. Inlet total pressure {self.InletTotalPressure.m:.0f} psi.</sup>',
			template='plotly_white'
		)

		return fig
