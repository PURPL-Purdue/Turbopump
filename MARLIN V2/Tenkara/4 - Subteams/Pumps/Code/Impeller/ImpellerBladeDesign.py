import numpy as np
from pint import Quantity as Q_
from plotly import graph_objects as go, subplots, express
from Inducer.RadialEquilibrium import InducerRadialEquilibrium
from typing import Optional, Literal
from VelocityTriangle import VelocityTriangle, PlotTriangles, CONVENTION_NOTES
from dataclasses import dataclass
from collections.abc import Callable
from plotly import graph_objects as go
from scipy import integrate
from warnings import warn
from .MeridionalDesign import MeridionalDesign
from .MeridionalContour import MeridionalCurve
from .PotentialFlow import Streamline
from . import ConformalMapping

"""
@Charlie Zhu
Minimal AI
"""

cum_trapz = integrate.cumulative_trapezoid

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
	lead_sweep: Q_[float] = Q_(10, 'deg')
	resolution: int = 200
	target_wrap: Q_[float] = Q_(120, 'deg')
	use_cubic: bool = True
	splitter_length_ratio: float = 0.6
	le_ratio: float = 2.0			# a / b of the leading-edge ellipse (Milos eq. 61, a/b = 2 in Fig. 27)
	le_resolution: int = 25			# points on the leading-edge ellipse, nose to full thickness
	side_theta: Literal['offset', 'milos'] = 'offset'	# wrap angle of the sides, see ConformalMapping.TransposedBladeLine.GetSides

class ImpellerBladeDesign:

	def __init__(self, meridional_design: MeridionalDesign, inducer_radial_equilibrium: Optional[InducerRadialEquilibrium]=None, blade_params: BladeParams = BladeParams()):

		self._params = blade_params
		self.Pump = meridional_design.Pump
		self.MeridionalDesign = meridional_design
		self._rad_eq = inducer_radial_equilibrium
		self.K_inlet = inducer_radial_equilibrium.GetOutletAngularMomentum() if inducer_radial_equilibrium is not None else Q_(0, 'm^2/s')

		self.HasSplitters: bool = self.Pump.Impeller.Z_blade_2 != self.Pump.Impeller.Z_blade_1
		# it's better to have slight positive incidence than slight negative incidence
		self._INCIDENCE_PADDING = np.pi/180 * 0	# degrees of extra angle

		self.__conformal_fig = None

		self._bladelines = self._CreateBladeLines(resolution=self._params.resolution)
		if self.HasSplitters:
			self._splitter_bladelines = self._CreateBladeLines(resolution=self._params.resolution, length_ratio=blade_params.splitter_length_ratio)

	def _GetInducerOutletC_uC_m(self, radius: Q_) -> tuple[Q_, Q_]:
		if radius < self._rad_eq.inducer.R_hub_2:
			radius = self._rad_eq.inducer.R_hub_2
			warn(f'Querying radius value less than inducer hub radius')
		elif radius > self._rad_eq.inducer.R_tip:
			radius = self._rad_eq.inducer.R_tip
			warn(f'Querying radius value greater than inducer tip radius')

		return (self._rad_eq.c_theta(radius, self._rad_eq.inducer.L_ax),
		        self._rad_eq.c_axial(radius, self._rad_eq.inducer.L_ax, blockage=False))

	def _SolveInletFlow(self, streamline: MeridionalCurve, free_vortex: bool=False) -> tuple[_TriInfo, _TriInfo, Q_, float]:
		"""
			Radius r_1 (m), flow angle beta_1' just inside the blades (rad, = blade angle at zero incidence) and degree of
			blockage rho_1 = 1 / tau_1 where the leading edge crosses streamline phi. Tangential blockage (CFturbo; Hassan
			et al., Sci. Rep. 2024, eq. 17):
				rho_1 = 1 - a / sin(beta_1'),		a = z s / (pi d_1),		s = blade_thick
				tan(beta_1') = tan(beta_1) / rho_1		(Milos eq. 33)
			Solved explicitly (own derivation): sin(beta_1' - beta_1) = a cos(beta_1)
		"""
		
		if self.Pump.Inducer is None: free_vortex = True
		
		r_1, c_m1, w_u1 = self._GetLeadingEdgeFlow(streamline, free_vortex=free_vortex)
		impeller = self.MeridionalDesign.Pump.Impeller

		# blockage coefficient: blade thickness over blade pitch
		pitch = 2 * np.pi * r_1 / impeller.Z_blade_1
		a: Q_[float] = impeller.blade_thick / pitch
		assert a < 1, 'The blades fill the inlet circumference, reduce the blade number or thickness'

		beta_1: Q_ = np.arctan2(c_m1, w_u1)							# without blockage
		beta_1_prime = beta_1 + np.arcsin(a * np.cos(beta_1))		# with blockage
		blockage_factor: float = (1 - a / np.sin(beta_1_prime)).to('dimensionless').m
		c_m1_prime: Q_[float] = c_m1 / blockage_factor

		return _TriInfo(c_m1_prime, beta_1_prime), _TriInfo(c_m1, beta_1), r_1, blockage_factor

	def _GetLeadingEdgeFlow(self, streamline: Streamline, free_vortex: bool=True) -> tuple[Q_, Q_, Q_]:
		"""
			Radius (m), meridional velocity c_m1 and relative swirl w_u1 = u_1 - c_u1 (m/s) where the leading edge
			crosses streamline, without blade blockage
		"""
		hits = streamline.FindIntersect(self.MeridionalDesign.LeadingEdge)
		assert len(hits) == 1, 'The leading edge must cross the streamline exactly once'
		axial_pos_z, radius, s, _ = hits[0]

		omega = self.MeridionalDesign.Pump.Impeller.DP.N_shaft.to('rad/s')
		u_1 = omega * radius # leading edge tangential speed

		if free_vortex:
			# free vortex defines: c_u * r = k	=> constant specific angular momentum
			c_u1 = self.K_inlet / radius
			c_m1: Q_[float] = streamline.Interpolate(streamline.v_m, s).to('m/s')
		else:
			c_u1, c_m1 = self._GetInducerOutletC_uC_m(radius)
		
		w_u1 = u_1 - c_u1
		assert w_u1 > 0, 'The preswirl exceeds the blade speed at the leading edge'

		return radius, c_m1, w_u1

	def PlotInletTriangles(self, phis: list[float]=[0, 0.5, 1], free_vortex: bool=True) -> go.Figure:
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
			prime, no_block, r, _ = self._SolveInletFlow(self.MeridionalDesign.PotentialFlow.GetStreamline(phi), free_vortex)
			u = self.Pump.DP.N_shaft * r
			triangles = [
				VelocityTriangle.FromBeta(u, no_block.c_m, no_block.Beta, '<sub>1</sub>'),
				VelocityTriangle.FromBeta(u, prime.c_m, prime.Beta, "<sub>1</sub>'", Q_(self._INCIDENCE_PADDING, 'rad')),
				]
			PlotTriangles(triangles, fig, row=i % n_rows + 1, col=i // n_rows + 1, names=['upstream', 'inside blades'])

		inlet_model = f'Inlet velocity model: {'free vortex' if free_vortex else 'superimposed inducer outlet'}'
		fig.update_layout(template='simple_white',
			title=f'{self.Pump.Propellant.name} Impeller Inlet Velocity Triangles<br><sup>{inlet_model}<br>{CONVENTION_NOTES["centrifugal"]}</sup>')

		return fig
	
	def __make_conformal_mapping_fig(self) -> go.Figure:

		if self.HasSplitters:
			fig = subplots.make_subplots(rows=1, cols=2, shared_yaxes=True, subplot_titles=('Full blades', 'Splitter blades'))
		else:
			fig = go.Figure()
		fig.update_layout(template='plotly_white', title_text='Blade Conformal Mapping')
		fig.update_xaxes(title_text='A(x) (meters)')
		fig.update_yaxes(title_text='Meridional arclength (meters)')
		return fig

	_SurfaceSide = Literal['suction', 'pressure', 'camber']

	def GetBladeSurfaceLines(self, side: _SurfaceSide | tuple[_SurfaceSide, ...]='camber', splitter: bool=False) -> list[np.ndarray] | tuple[list[np.ndarray], ...]:
		"""
			Blade lines of one surface (side a string), or a tuple of them in the same order (side a tuple of strings)
		"""
		assert not splitter or self.HasSplitters, 'No splitter blades'

		bladeline_surfaces, _ = self._splitter_bladelines if splitter else self._bladelines
		index = {'suction': 0, 'pressure': 1, 'camber': 2}
		if isinstance(side, str):
			return bladeline_surfaces[index[side]]
		return tuple(bladeline_surfaces[index[s]] for s in side)

	def _CreateBladeLines(self, num_streams: int=11, resolution: int=200, length_ratio: float=1.0) -> tuple[tuple[list[np.ndarray], list[np.ndarray], list[np.ndarray]], list[BetaDistribution]]:
		assert length_ratio > 0 and length_ratio <= 1.0, 'Bro'
		splitter: bool = length_ratio != 1.0
		
		lead_sweep: float = self._params.lead_sweep.to('rad').m
		ref_phi = 1 # shroud streamline
		ref_stream = self.MeridionalDesign.PotentialFlow.GetStreamline(ref_phi)

		def GetFullBladeArclengthArray(stream: Streamline, le: MeridionalCurve, te: MeridionalCurve):
			le_intersect = stream.FindIntersect(le, n=resolution)
			te_intersect = stream.FindIntersect(te, n=resolution)
			x_0 = le_intersect[0][2] * stream.length_m
			x_f = te_intersect[0][2] * stream.length_m

			return np.linspace(x_0, x_f, resolution)

		ref_beta_1 = self._SolveInletFlow(ref_stream)[0].Beta.to('rad').m + self._INCIDENCE_PADDING
		ref_arcx = GetFullBladeArclengthArray(ref_stream, self.MeridionalDesign.LeadingEdge, self.MeridionalDesign.TrailingEdge)

		radii_ref = ref_stream.Points(ref_arcx)[:, 1]
		beta_2 = self.Pump.Impeller.Beta2B.to('rad').m

		_, REF_TOTAL_WRAP = (
			ConformalMapping.CamberLine(radii=radii_ref, x_f=ref_arcx[-1] - ref_arcx[0], beta_1=ref_beta_1, beta_2=beta_2)
			.GetBetaDistributionFromSj(self._params.sj)
		)
		
		betas: list[BetaDistribution] = []
		surface_points: tuple[list[np.ndarray], ...] = ([], [], [])	# suction, pressure, camber

		if not self.__conformal_fig:
			self.__conformal_fig = self.__make_conformal_mapping_fig()

		for phi in np.linspace(0, 1, num_streams):

			streamline = self.MeridionalDesign.PotentialFlow.GetStreamline(phi)
			beta_1 = self._SolveInletFlow(streamline)[0].Beta.to('rad').m + self._INCIDENCE_PADDING
			arc_x = GetFullBladeArclengthArray(streamline, self.MeridionalDesign.LeadingEdge, self.MeridionalDesign.TrailingEdge)
			
			stream_zr = streamline.Points(arc_x)
			radii = stream_zr[:, 1]
			
			arc_lead = np.sin(lead_sweep) * (radii[0] - radii_ref[0])
			sweep_offset = arc_lead / radii[0]
			camber = ConformalMapping.CamberLine(radii=radii, x_f=arc_x[-1] - arc_x[0], beta_1=beta_1, beta_2=beta_2)
			parab_beta, mean_wrap = (camber
				.SolveBetaDistributionFromWrap(REF_TOTAL_WRAP + sweep_offset)
			)

			transposed_bl = ConformalMapping.TransposedBladeLine(
				camber.x_f, self.Pump.Impeller.blade_thick.to('m').m, le_ratio=self._params.le_ratio
			)

			if splitter:
				splitter_offset: float = 2*np.pi / self.Pump.Impeller.Z_blade_2
				length_ratio = self._params.splitter_length_ratio
			else:
				splitter_offset = 0
				length_ratio = None

			*sides, conformals = transposed_bl.GetSides(
				parabolic_beta=parab_beta, stream=streamline,
				x_0_mean=arc_x[0], resolution=resolution, le_resolution=self._params.le_resolution,
				length_ratio=length_ratio, side_theta=self._params.side_theta
				)

			# GetSides returns (zr, theta) for suction, pressure, camber in turn. All three turn with the camber, whose
			# trailing edge goes to the reference wrap, so the sides keep their offsets from it
			theta_shift = REF_TOTAL_WRAP - sides[-1][-1] + splitter_offset
			for points, zr, theta in zip(surface_points, sides[0::2], sides[1::2]):
				points.append(np.column_stack((zr, theta + theta_shift)))

			betas.append(
				BetaDistribution(Get=parab_beta.GetBeta, Domain=(0, arc_x[-1] - arc_x[0]), TotalWrap=mean_wrap)
			)

			A_s, L_s, A_p, L_p, _, _ = conformals
			self.__conformal_fig.add_trace(go.Scatter(
				x = A_s,
				y = L_s, name=f'phi = {phi}', showlegend=False
			), row=1, col=2 if splitter else 1)

			self.__conformal_fig.add_trace(go.Scatter(
				x = A_p,
				y = L_p, name=f'phi = {phi}', showlegend=False
			), row=1, col=2 if splitter else 1)

		return (surface_points, betas)

	def PlotBladeLines(self, unit:str='in') -> go.Figure:

		fig = subplots.make_subplots(rows=1, cols=2, column_widths=[0.6, 0.4],
				subplot_titles=(f'{self.Pump.Propellant.name} Impeller Plan View (Bottom Up)', r'Variation of <i>β</i>')
			)
		
		colors = express.colors.qualitative.Plotly

		for is_splitter in (False, True) if self.HasSplitters else (False,) :

			bladelines = self._splitter_bladelines if is_splitter else self._bladelines

			blade_surface_points, betas = bladelines
			num_streams = len(betas)

			for i in range(num_streams):

				phi: float = i / (num_streams - 1)
				beta_dis = betas[i]

				points = blade_surface_points[2] # camber points

				r, theta = points[i][:, 1], points[i][:, 2]

				arc_m=np.linspace(*beta_dis.Domain, len(r))
				beta = beta_dis.Get(arc_m)

				theta_blade: float = 2*np.pi / self.Pump.Impeller.Z_blade_1

				# all Z blades in one trace, so each streamline keeps one color and legend entry; NaN breaks the line between blades
				theta_n = np.concatenate([np.append(theta + n*theta_blade, np.nan) for n in range(self.Pump.Impeller.Z_blade_1)])
				r_n = np.tile(np.append(r, np.nan), self.Pump.Impeller.Z_blade_1)

				color = colors[i % len(colors)]
				trace_name = f'Phi = {phi:.2f}, Wrap = {(beta_dis.TotalWrap * 180/np.pi):.1f} °'

				# plot plan view
				fig.add_trace(go.Scatter(
					x=Q_(r_n * np.cos(theta_n), 'm').to(unit).m,
					y=Q_(r_n * np.sin(theta_n), 'm').to(unit).m,
					name=trace_name, legendgroup=f'{i}', line=dict(color=color), showlegend=False
				), row=1, col=1)

				if is_splitter: continue

				# plot Beta distribution, if blade lines are regular blades (not splitter)
				fig.add_trace(go.Scatter(
					x=Q_(arc_m, 'm').to(unit).m, y=np.rad2deg(beta),
					name=trace_name, legendgroup=f'{i}', line=dict(color=color), showlegend=True
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
			template='plotly_white', title_text=f'{self.Pump.Propellant.name} Impeller Blade Layout'
		))
		return fig

	def ShowPlotConformalMappingFunction(self):
		self.__conformal_fig.show()