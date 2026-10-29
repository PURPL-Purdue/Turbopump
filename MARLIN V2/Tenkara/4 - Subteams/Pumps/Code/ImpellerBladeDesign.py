import numpy as np
from pint import Quantity as Q_
from plotly import graph_objects as go, subplots
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

@dataclass
class _TriInfo:
	c_m: Q_
	Beta: Q_

@dataclass
class BetaDistribution:
	Get: Callable[[float], float]
	Domain: tuple[float, float]

def _TraceCircle(r: float, t: list[float] | np.ndarray, name: str=None, color: str='grey', dash: str='dot') -> go.Scatter:
	return go.Scatter(
		x=r * np.cos(t), y=r * np.sin(t),
		mode='lines', line=dict(dash=dash, color=color), name=name
	)

@dataclass
class BladeParams:
	sj: float = 0.05
	lead_sweep: Q_[float] = Q_(10, 'deg')

class ImpellerBladeDesign:

	def __init__(self, meridional_design: MeridionalDesign, inducer_radial_equilibrium: Optional[InducerRadialEquilibrium]=None, blade_params: BladeParams = BladeParams()):

		self._params = blade_params
		self.Pump = meridional_design.Pump
		self._mer = meridional_design
		self._rad_eq = inducer_radial_equilibrium
		self.K_inlet = inducer_radial_equilibrium.GetOutletAngularMomentum() if inducer_radial_equilibrium is not None else Q_(0, 'm^2/s')
		
		# it's better to have slight positive incidence than slight negative incidence
		self._INCIDENCE_PADDING = np.pi/180 * 2	# degrees of extra angle

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
		hits = streamline.FindIntersect(self._mer.LeadingEdge)
		assert len(hits) == 1, 'The leading edge must cross the streamline exactly once'
		_, r_1, s, _ = hits[0]

		omega = self._mer.Pump.Impeller.DP.N_shaft.to('rad/s')

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
		hits = curve.FindIntersect(self._mer.LeadingEdge)
		assert len(hits) == 1, 'The leading edge must meet the curve exactly once'
		_, _, s_curve, s_LE = hits[0]
		cos_angle = abs(curve.Tangent(s_curve) @ self._mer.LeadingEdge.Tangent(s_LE))
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
			prime, no_block, r, _ = self._SolveInletFlow(self._mer.PotentialFlow.GetStreamline(phi))
			u = self.Pump.DP.N_shaft * r
			triangles = [
				VelocityTriangle.FromBeta(u, no_block.c_m, no_block.Beta, '<sub>1</sub>'),
				VelocityTriangle.FromBeta(u, prime.c_m, prime.Beta, "<sub>1</sub>'", Q_(self._INCIDENCE_PADDING, 'rad')),
				]
			PlotTriangles(triangles, fig, row=i % n_rows + 1, col=i // n_rows + 1, names=['upstream', 'inside blades'])

		fig.update_layout(template='simple_white',
			title=f'Inlet Velocity Triangles<br><sup>{CONVENTION_NOTES["centrifugal"]}</sup>')

		return fig

	def _GetWrap(self, r: np.ndarray[float], beta_1: float, x_f: float, sj: float=None) -> tuple[float, BetaDistribution]:
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

		return np.trapezoid(y=1/(r * np.tan(Beta)), x=x), BetaDistribution(
			Get=GetBeta, Domain=(0, x_f)
		)

	def _SetWrap(self, beta_1: float, r: np.ndarray[float], x_f: float, target_wrap: float) -> tuple[BetaDistribution, float]:
		# sj = 0 or 1 makes JoinParabola's system singular, so stay just inside
		guess_high = 0.99
		guess_low = 0.01

		wrap_high, beta_high = self._GetWrap(r=r, beta_1=beta_1, x_f=x_f, sj=guess_high)
		wrap_low, beta_low = self._GetWrap(r=r, beta_1=beta_1, x_f=x_f, sj=guess_low)

		residual_new = float('inf')
		residual_high, residual_low = wrap_high - target_wrap, wrap_low - target_wrap
		# same sign at both ends: the target is not between them, take the closer end
		if residual_high * residual_low > 0:
			print(f'The target wrap {target_wrap*180/np.pi:.1f} is out of reach by moving sj alone')
			if np.abs(residual_high) > np.abs(residual_low):
				return beta_low, wrap_low
			return beta_high, wrap_high

		# bisection: keep the half of the bracket where the residual changes sign
		while abs(residual_new) > 1e-4:
			new_guess: float = (guess_high + guess_low) / 2
			wrap_new, beta_dis = self._GetWrap(r=r, beta_1=beta_1, x_f=x_f, sj=new_guess)
			residual_new = wrap_new - target_wrap

			if residual_new * residual_low > 0:
				guess_low, residual_low = new_guess, residual_new
			else:
				guess_high, residual_high = new_guess, residual_new

		return beta_dis, wrap_new

	def Uh(self):
		lead_sweep: float = self._params.lead_sweep.to('rad').m
		ref_phi = 1
		ref_stream = self._mer.PotentialFlow.GetStreamline(ref_phi)

		ref_beta_1 = self._SolveInletFlow(ref_stream)[0].Beta.to('rad').m
		ref_bladestream = ref_stream.Between(self._mer.LeadingEdge, self._mer.TrailingEdge)
		r_ref = ref_bladestream.Points(np.linspace(0, ref_bladestream.length_m, 200))[:, 1]

		ref_wrap, mean_beta_dis = self._GetWrap(
			r=r_ref, beta_1=ref_beta_1, x_f=ref_bladestream.length_m, sj=self._params.sj
			)

		fig = go.Figure()
		pfig = go.Figure()

		for phi in np.linspace(0, 1, 10):

			if phi != ref_phi:

				streamline = self._mer.PotentialFlow.GetStreamline(phi)
				beta_1 = self._SolveInletFlow(streamline)[0].Beta.to('rad').m
				blade_stream = streamline.Between(self._mer.LeadingEdge, self._mer.TrailingEdge)
				
				r = blade_stream.Points(np.linspace(0, blade_stream.length_m, len(r_ref)))[:, 1]
				
				arc_lead = np.sin(lead_sweep) * (r[0] - r_ref[0])
				theta_offset = arc_lead / r[0]

				beta_dis, wrap = self._SetWrap(
					beta_1=beta_1, r=r, x_f=blade_stream.length_m,
					target_wrap=ref_wrap + theta_offset
				)
			else:
				beta_dis, wrap, r = mean_beta_dis, ref_wrap, r_ref
			
			arc_x=np.linspace(*beta_dis.Domain, len(r))
			beta = beta_dis.Get(arc_x)

			trace_name = f'Phi = {phi:.2f}, Wrap = {(wrap * 180/np.pi):.1f} °'
			fig.add_trace(go.Scatter(
				x=Q_(arc_x, 'm').to('in').m,
				y=np.rad2deg(beta),
				name=trace_name
			))

			theta: np.ndarray = integrate.cumulative_trapezoid(1/(r * np.tan(beta)), arc_x, initial=0)
			theta += ref_wrap - wrap

			theta_blade: float = 2*np.pi / self.Pump.Impeller.Z_blade
			Z = self.Pump.Impeller.Z_blade

			# all Z blades in one trace, so each streamline keeps one color and legend entry; NaN breaks the line between blades
			theta_n = np.concatenate([np.append(theta + n*theta_blade, np.nan) for n in range(Z)])
			r_n = np.tile(np.append(r, np.nan), Z)

			pfig.add_trace(go.Scatter(
				x=Q_(r_n * np.cos(theta_n), 'm').to('in').m,
				y=Q_(r_n * np.sin(theta_n), 'm').to('in').m,
				name=trace_name
			))

		fig.update_layout(dict(
			template='plotly_white',
			xaxis_title=r'Arclength <i>x</i> (in.)',
			yaxis_title=r'Blade angle <i>β</i> (°)'
		))

		R_MAX = self.Pump.Impeller.d_2.to('in').m / 2 * 1.2

		t = np.linspace(0, 2*np.pi, 100)

		pfig.add_trace(_TraceCircle(
			r=self.Pump.Impeller.d_1.to('in').m / 2,
			t=t, name=r'd<sub>1</sub>', dash='longdash'
		))

		pfig.add_trace(_TraceCircle(
			r=self.Pump.Impeller.d_hub.to('in').m / 2,
			t=t, name=r'd<sub>hub</sub>', dash='dash'
		))

		pfig.add_trace(_TraceCircle(
			r=self.Pump.Impeller.d_2.to('in').m / 2,
			t=t, name=r'd<sub>2</sub>', dash='longdash'
		))


		pfig.update_layout(dict(
			template='plotly_white',
			xaxis_title=r'<i>x</i> (in.)',
			yaxis_title=r'<i>y</i> (in.)',
			yaxis=dict(range=[-R_MAX, R_MAX], zeroline=True, constrain='domain', scaleanchor='x', scaleratio=1),
			xaxis=dict(range=[-R_MAX, R_MAX], zeroline=True, constrain='domain')
		))
		
		pfig.show()
		fig.show()