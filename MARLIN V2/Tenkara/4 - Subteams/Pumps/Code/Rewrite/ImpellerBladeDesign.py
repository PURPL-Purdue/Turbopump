import numpy as np
from pint import Quantity as Q_
from plotly import graph_objects as go, subplots
from InducerRadialEquilibrium import InducerRadialEquilibrium
from MeridionalDesign import MeridionalDesign
from MeridionalContour import MeridionalCurve
from PotentialFlow import Streamline
from typing import Optional
from VelocityTriangle import VelocityTriangle
from dataclasses import dataclass

@dataclass
class _TriInfo:
	c_m: Q_
	Beta: Q_

class ImpellerBladeDesign:

	def __init__(self, meridional_design: MeridionalDesign, inducer_radial_equilibrium: Optional[InducerRadialEquilibrium]=None,
			n_streamlines: int = 11):

		self.Pump = meridional_design.Pump
		self._mer = meridional_design
		self._rad_eq = inducer_radial_equilibrium
		self.K_inlet = inducer_radial_equilibrium.GetOutletAngularMomentum() if inducer_radial_equilibrium is not None else Q_(0, 'm^2/s')
		
		# it's better to have slight positive incidence than slight negative incidence
		self._INCIDENCE_PADDING = np.pi/180 * 2	# degrees of extra angle

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
		hits = streamline.FindIntersect(self._mer.LeadingEdge)
		assert len(hits) == 1, 'The leading edge must cross the streamline exactly once'
		_, r_1, s, _ = hits[0]

		omega = self._mer.Pump.Impeller.DP.N_shaft.to('rad/s')
		c_m1: Q_[float] = streamline.Interpolate(streamline.v_m, s).to('m/s')
		w_u1 = omega * r_1 - self.GetInletC_u(r_1).to('m/s')
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
		# at most 3 rows, further triangles fill new columns
		n_rows = min(len(phis), 3)
		n_cols = -(-len(phis) // n_rows)

		# plotly takes the titles row by row, the triangles fill column by column
		titles = [''] * (n_rows * n_cols)
		for i, phi in enumerate(phis):
			titles[i % n_rows * n_cols + i // n_rows] = f'Inlet Velocity Triangles at Flow Fraction ϕ = {phi}'

		fig = subplots.make_subplots(rows=n_rows, cols=n_cols, subplot_titles=titles)
		for i, phi in enumerate(phis):
			prime, no_block, r, _ = self._SolveInletFlow(phi)
			u = self.Pump.DP.N_shaft * r

			no_block_w_u = no_block.c_m / np.tan(no_block.Beta)
			
			triangle = VelocityTriangle(
				u=u, c_m=no_block.c_m,
				c_u = u - no_block_w_u,	w=np.sqrt(no_block_w_u**2 + no_block.c_m**2),
				beta_b=prime.Beta, station=1
				)

			triangle.Plot(
				fig=fig,
				row=i % n_rows + 1,
				col=i // n_rows + 1
				)

		fig.update_layout(template='simple_white')
		fig.update_xaxes(title_text='Tangential velocity (m/s)')
		fig.update_yaxes(title_text='Meridional velocity (m/s)')

		return fig
