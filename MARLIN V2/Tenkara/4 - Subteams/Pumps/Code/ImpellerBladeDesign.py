import numpy as np
from pint import Quantity as Q_
from plotly import graph_objects as go, subplots, express
from Inducer.RadialEquilibrium import InducerRadialEquilibrium
from MeridionalDesign import MeridionalDesign
from MeridionalContour import MeridionalCurve
from PotentialFlow import Streamline
from typing import Optional
from VelocityTriangle import VelocityTriangle, PlotTriangles, CONVENTION_NOTES
from dataclasses import dataclass
from JoinParabola import JoinParabola
from collections.abc import Callable
from plotly import graph_objects as go
from scipy import integrate
from warnings import warn
from ConformalMapping import ConformalMapping, BetaCMM

@dataclass
class _TriInfo:
	c_m: Q_
	Beta: Q_

@dataclass
class BetaDistribution:
	Get: Callable[[float], float]
	Domain: tuple[float, float]
	TotalWrap: float

def _TraceCircle(r: float, t: list[float] | np.ndarray, name: str=None, color: str='grey', dash: str='dot') -> go.Scatter:
	return go.Scatter(
		x=r * np.cos(t), y=r * np.sin(t),
		mode='lines', line=dict(dash=dash, color=color), name=name,
	)

@dataclass
class BladeParams:
	sj: float = 0.05
	lead_sweep: Q_[float] = Q_(10, 'deg'),
	resolution: int = 200

class ImpellerBladeDesign:

	def __init__(self, meridional_design: MeridionalDesign, inducer_radial_equilibrium: Optional[InducerRadialEquilibrium]=None, blade_params: BladeParams = BladeParams()):

		self._params = blade_params
		self.Pump = meridional_design.Pump
		self.MeridionalDesign = meridional_design
		self._rad_eq = inducer_radial_equilibrium
		self.K_inlet = inducer_radial_equilibrium.GetOutletAngularMomentum() if inducer_radial_equilibrium is not None else Q_(0, 'm^2/s')
		
		# it's better to have slight positive incidence than slight negative incidence
		self._INCIDENCE_PADDING = np.pi/180 * 2	# degrees of extra angle

		self._CreateBladeLines2(resolution=self._params.resolution)

	def _GetInletC_u(self, r: Q_, free_vortex: bool=True) -> Q_:

		if free_vortex:
			# free vortex:
			return self.K_inlet/r

		# superimpose inducer outlet velocity distribution
		if r < self._rad_eq.inducer.R_hub_2:
			r = self._rad_eq.inducer.R_hub_2
		elif r > self._rad_eq.inducer.D_tip/2:
			r = self._rad_eq.inducer.D_tip/2
		
		return self._rad_eq.c_theta(r, self._rad_eq.inducer.L_ax)

	def _GetInletC_m(self, r: Q_) -> Q_:
		if r < self._rad_eq.inducer.R_hub_2:
			r = self._rad_eq.inducer.R_hub_2
		elif r > self._rad_eq.inducer.D_tip/2:
			r = self._rad_eq.inducer.D_tip/2

		return self._rad_eq.c_axial(r, self._rad_eq.inducer.L_ax)

	def _SolveInletFlow(self, streamline: MeridionalCurve) -> tuple[_TriInfo, _TriInfo, Q_, float]:
		"""
			Radius r_1 (m), flow angle beta_1' just inside the blades (rad, = blade angle at zero incidence) and degree of
			blockage rho_1 where the leading edge crosses streamline phi. Tangential blockage (CFturbo; Hassan et al.,
			Sci. Rep. 2024, eq. 17) with the meridional leading edge angle lambda_1, which slightly overestimates it:
				rho_1 = 1 - a / sin(beta_1'),		a = z s / (pi d_1 sin(lambda_1)),	s = blade_thick
				tan(beta_1') = tan(beta_1) / rho_1		(Milos eq. 33)
			Solved explicitly (own derivation): sin(beta_1' - beta_1) = a cos(beta_1)
		"""
		lambda_1 = self._GetLeadingEdgeAngle(streamline)

		r_1, c_m1, w_u1 = self._GetLeadingEdgeFlow(streamline)
		impeller = self.MeridionalDesign.Pump.Impeller

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
		hits = streamline.FindIntersect(self.MeridionalDesign.LeadingEdge)
		assert len(hits) == 1, 'The leading edge must cross the streamline exactly once'
		_, r_1, s, _ = hits[0]

		omega = self.MeridionalDesign.Pump.Impeller.DP.N_shaft.to('rad/s')

		#c_m1: Q_[float] = self._GetInletC_m(r_1)
		c_m1: Q_[float] = streamline.Interpolate(streamline.v_m, s).to('m/s')
		
		w_u1 = omega * r_1 - self._GetInletC_u(r_1)
		assert w_u1 > 0, 'The preswirl exceeds the blade speed at the leading edge'

		return r_1, c_m1, w_u1

	def _GetLeadingEdgeAngle(self, curve: MeridionalCurve) -> Q_[float]:
		"""
			Acute angle in the meridional plane between the leading edge and curve (the hub or shroud contour, or a
			streamline) where they meet, from the unit tangents of both: arccos|t_LE . t_curve|
		"""
		hits = curve.FindIntersect(self.MeridionalDesign.LeadingEdge)
		assert len(hits) == 1, 'The leading edge must meet the curve exactly once'
		_, _, s_curve, s_LE = hits[0]
		cos_angle = abs(curve.Tangent(s_curve) @ self.MeridionalDesign.LeadingEdge.Tangent(s_LE))
		return Q_(np.arccos(min(cos_angle, 1.0)), 'rad').to('deg')

	def PlotInletTriangles(self, phis: list[float]=[0, 0.5, 1]) -> go.Figure:
		"""
			Inlet triangles where the leading edge crosses each streamline phi: upstream of the blades (1) and just inside
			them with the blade blockage (1'), with the blade direction at the incidence i = beta_1B - beta_1' from W_1'
			(Gülich Fig. 3.1)
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
			prime, no_block, r, _ = self._SolveInletFlow(self.MeridionalDesign.PotentialFlow.GetStreamline(phi))
			u = self.Pump.DP.N_shaft * r
			triangles = [
				VelocityTriangle.FromBeta(u, no_block.c_m, no_block.Beta, '<sub>1</sub>'),
				VelocityTriangle.FromBeta(u, prime.c_m, prime.Beta, "<sub>1</sub>'", Q_(self._INCIDENCE_PADDING, 'rad')),
				]
			PlotTriangles(triangles, fig, row=i % n_rows + 1, col=i // n_rows + 1, names=['upstream', 'inside blades'])

		fig.update_layout(template='simple_white',
			title=f'Inlet Velocity Triangles<br><sup>{CONVENTION_NOTES["centrifugal"]}</sup>')

		return fig

	def _GetBetaDistribution(self, r: np.ndarray[float], beta_1: float, x_f: float, sj: float=None) -> BetaDistribution:
		"""
			Returns the wrap angle (in radians) of a blade streamline at the specified flow_fraction, from the leading
			to the trailing edge. xj (m from the leading edge) is where the two beta parabolas join, Milos's x_3;
			defaults to mid-blade
		"""
		beta_2: float = self.Pump.Impeller.Beta2B.to('rad').m

		x_j = x_f * sj if sj is not None else x_f * 0.5
		GetBeta = JoinParabola(y0=beta_1, yf=beta_2, xj=x_j, xf=x_f)

		# approximate wrap angle using trapezoidal Riemann sum
		# Milo's pg. 51, eq 51 and 52
		x = np.linspace(0, x_f, len(r))
		Beta = GetBeta(x)

		total_wrap = np.trapezoid(y=1/(r * np.tan(Beta)), x=x)

		return BetaDistribution(
			Get=GetBeta, Domain=(0, x_f), TotalWrap=total_wrap
		)

	def _SolveBetaDistribution(self, beta_1: float, r: np.ndarray[float], x_f: float, target_wrap: float) -> BetaDistribution:
		# sj = 0 or 1 makes JoinParabola's system singular, so stay just inside
		guess_high = 0.99
		guess_low = 0.01

		beta_high = self._GetBetaDistribution(r=r, beta_1=beta_1, x_f=x_f, sj=guess_high)
		beta_low = self._GetBetaDistribution(r=r, beta_1=beta_1, x_f=x_f, sj=guess_low)

		residual_new = float('inf')
		residual_high, residual_low = beta_high.TotalWrap - target_wrap, beta_low.TotalWrap - target_wrap
		# same sign at both ends: the target is not between them, take the closer end
		if residual_high * residual_low > 0:
			beta=beta_low if np.abs(residual_high) > np.abs(residual_low) else beta_high
			warn(f'The target wrap {np.rad2deg(target_wrap):.1f} deg is out of reach by moving sj alone, resorting to {np.rad2deg(beta.TotalWrap)} deg')
			return beta
		# bisection: keep the half of the bracket where the residual changes sign
		while abs(residual_new) > 1e-4:
			new_guess: float = (guess_high + guess_low) / 2
			beta_dis = self._GetBetaDistribution(r=r, beta_1=beta_1, x_f=x_f, sj=new_guess)
			residual_new = beta_dis.TotalWrap - target_wrap

			if residual_new * residual_low > 0:
				guess_low, residual_low = new_guess, residual_new
			else:
				guess_high, residual_high = new_guess, residual_new

		return beta_dis

	def _CreateBladeLines2(self, num_streams: int=10, resolution: int=200) -> tuple[list[np.ndarray[float]], list[BetaDistribution]]:
		betas: list[BetaDistribution] = []
		points_3d: list[np.ndarray[float]] = []
		lead_sweep = self._params.lead_sweep.to('rad').m

		for i in range(num_streams):

			phi = i/num_streams

			streamline = self.MeridionalDesign.PotentialFlow.GetStreamline(phi)
			beta_1 = self._SolveInletFlow(streamline)[0].Beta.to('rad').m
			blade_stream = streamline.Between(self.MeridionalDesign.LeadingEdge, self.MeridionalDesign.TrailingEdge)

			stream_zr = blade_stream.Points(np.linspace(0, blade_stream.length_m, resolution))
			r = stream_zr[:, 1]

			if i == 0: r_0 = r[0]
			arc_lead = np.sin(lead_sweep) * (r[0] - r_0)
			theta_offset = arc_lead / r[0]

			cm = ConformalMapping(
				beta_1=beta_1, beta_2=self.Pump.Impeller.Beta2B.to('rad').m,
				arclen=blade_stream.length_m, target_wrap=np.deg2rad(95) + theta_offset, radii=r
				)
			
			beta_dis: BetaCMM = cm.SolveDistribution()

			betas.append(beta_dis)
			points_3d.append(np.column_stack((stream_zr, beta_dis.Theta)))

		self._bladelines = (points_3d, betas)
		return self._bladelines

	def _CreateBladeLines(self, num_streams: int=10, resolution: int=200) -> tuple[list[np.ndarray[float]], list[BetaDistribution]]:
		"""
			Blade lines (z, r, theta) along num_streams streamlines, hub (phi = 0) to shroud (phi = 1), each sampled at
			resolution points from the leading edge (x = 0) to the trailing edge (Milos p. 51, eqs. 51-52):
				beta(x)  = JoinParabola from beta_1' to Beta2B, joined at x_j = sj x_f
				theta(x) = theta_1 + integral_0^x dx / (r tan(beta))
			The shroud is the reference (sj = BladeParams.sj, theta_1 = 0). On the other streamlines sj is bisected so
			every trailing edge ends at the shroud's wrap (radial in plan view) and the leading edge sweeps:
				theta_1  = sin(lead_sweep) (r_1,shroud - r_1) / r_1
			Returns two lists, one entry per streamline, hub first:
				points: (resolution, 3) arrays of (z, r, theta), in m, m and rad
				betas:  the BetaDistribution, whose Get(x) is beta (rad) at x (m) on its Domain
		"""

		lead_sweep: float = self._params.lead_sweep.to('rad').m
		ref_phi = 1 # shroud streamline
		ref_stream = self.MeridionalDesign.PotentialFlow.GetStreamline(ref_phi)

		ref_beta_1 = self._SolveInletFlow(ref_stream)[0].Beta.to('rad').m
		ref_bladestream = ref_stream.Between(self.MeridionalDesign.LeadingEdge, self.MeridionalDesign.TrailingEdge)
		ref_stream_zr = ref_bladestream.Points(np.linspace(0, ref_bladestream.length_m, resolution))
		r_ref = ref_stream_zr[:, 1]

		ref_beta_dis = self._GetBetaDistribution(
			r=r_ref, beta_1=ref_beta_1, x_f=ref_bladestream.length_m, sj=self._params.sj
			)

		betas: list[BetaDistribution] = []
		points_3d: list[np.ndarray[float]] = []

		for phi in np.linspace(0, 1, num_streams):

			if phi != ref_phi:

				streamline = self.MeridionalDesign.PotentialFlow.GetStreamline(phi)
				beta_1 = self._SolveInletFlow(streamline)[0].Beta.to('rad').m
				blade_stream = streamline.Between(self.MeridionalDesign.LeadingEdge, self.MeridionalDesign.TrailingEdge)

				stream_zr = blade_stream.Points(np.linspace(0, blade_stream.length_m, len(r_ref)))
				r = stream_zr[:, 1]
				
				arc_lead = np.sin(lead_sweep) * (r[0] - r_ref[0])
				theta_offset = arc_lead / r[0]

				beta_dis = self._SolveBetaDistribution(
					beta_1=beta_1, r=r, x_f=blade_stream.length_m,
					target_wrap=ref_beta_dis.TotalWrap + theta_offset
				)
			else:
				beta_dis, stream_zr, r = ref_beta_dis, ref_stream_zr, r_ref
			
			arc_m=np.linspace(*beta_dis.Domain, len(r))
			beta = beta_dis.Get(arc_m)

			theta: np.ndarray = integrate.cumulative_trapezoid(1/(r * np.tan(beta)), arc_m, initial=0)
			theta += ref_beta_dis.TotalWrap - beta_dis.TotalWrap

			betas.append(beta_dis)
			points_3d.append(np.column_stack((stream_zr, theta)))

		self._bladelines = (points_3d, betas)
		return self._bladelines

	def GetBladeLines(self):
		return self._bladelines

	def PlotBladeLines(self, unit:str='in.', propellant: str='') -> go.Figure:

		fig = subplots.make_subplots(rows=1, cols=2, column_widths=[0.6, 0.4],
				subplot_titles=(f'{propellant} Impeller Plan View (Bottom Up)', r'Variation of <i>β</i>')
			)

		points, betas = self.GetBladeLines()
		num_streams = len(points)
		colors = express.colors.qualitative.Plotly

		for i in range(num_streams):

			phi: float = i / num_streams
			r, theta, beta_dis = points[i][:, 1], points[i][:, 2], betas[i]

			arc_m=np.linspace(*beta_dis.Domain, len(r))
			beta = beta_dis.Get(arc_m)

			trace_name = f'Phi = {phi:.2f}, Wrap = {(beta_dis.TotalWrap * 180/np.pi):.1f} °'

			theta_blade: float = 2*np.pi / self.Pump.Impeller.Z_blade
			Z = self.Pump.Impeller.Z_blade

			# all Z blades in one trace, so each streamline keeps one color and legend entry; NaN breaks the line between blades
			theta_n = np.concatenate([np.append(theta + n*theta_blade, np.nan) for n in range(Z)])
			r_n = np.tile(np.append(r, np.nan), Z)

			color = colors[i % len(colors)]
			fig.add_trace(go.Scatter(
				x=Q_(r_n * np.cos(theta_n), 'm').to(unit).m,
				y=Q_(r_n * np.sin(theta_n), 'm').to(unit).m,
				name=trace_name, legendgroup=f'{i}', line=dict(color=color)
			), row=1, col=1)

			fig.add_trace(go.Scatter(
				x=Q_(arc_m, 'm').to(unit).m, y=np.rad2deg(beta),
				name=trace_name, legendgroup=f'{i}', line=dict(color=color), showlegend=False
			), row=1, col=2)

		fig.update_xaxes(title_text=rf'Meridional arclength <i>m</i> ({unit})', row=1, col=2)
		fig.update_yaxes(title_text=r'Blade angle <i>β</i> (°)', row=1, col=2)

		R_MAX = self.Pump.Impeller.d_2.to(unit).m / 2 * 1.2

		t = np.linspace(0, 2*np.pi, 100)
		fig.add_traces(
			(
				_TraceCircle(
					r=self.Pump.Impeller.d_1.to(unit).m / 2,
					t=t, name=r'd<sub>1</sub>', dash='longdash'
				),
				_TraceCircle(
					r=self.Pump.Impeller.d_hub.to(unit).m / 2,
					t=t, name=r'd<sub>hub</sub>', dash='dash'
				),
				_TraceCircle(
					r=self.Pump.Impeller.d_2.to(unit).m / 2,
					t=t, name=r'd<sub>2</sub>', dash='longdash'
				)
			), rows=1, cols=1
		)

		fig.update_xaxes(
			title_text=rf'<i>x</i> ({unit})', row=1, col=1,
			range=[-R_MAX, R_MAX], zeroline=True, constrain='domain'
			)
		fig.update_yaxes(
			title_text=rf'<i>y</i> ({unit})', row=1, col=1,
			range=[-R_MAX, R_MAX], zeroline=True, constrain='domain', scaleanchor='x', scaleratio=1
			)

		fig.update_layout(dict(
			template='plotly_white', title_text=f'{propellant} Impeller Blade Layout'
		))
		return fig