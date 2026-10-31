from plotly import graph_objects as go
from Pump import Pump
from pint import Quantity as Q_
from .MeridionalContour import MeridionalContour, MeridionalCurve, FlowDomain, MeanCurve, BezierCurve, Polyline
from .PotentialFlow import PotentialFlow
import numpy as np
from dataclasses import dataclass

@dataclass
class MeridionalGeometry:
	shroud_fullness: tuple[float, float] = (1, 0.3)
	hub_fullness: tuple[float, float] = (0.5, 0.7)
	# contour angles from radial at the outlet, positive leaning downstream (+z), see BuildMeridian
	shroud_outlet_angle: Q_[float] = Q_(5, 'deg')
	hub_outlet_angle: Q_[float] = Q_(0, 'deg')
	hub_taper_ratio: float = 0.3

	Z_e_factor: float = 1
	shroud_flat: float = 0.1

class MeridionalDesign:

	def __init__(self, pump: Pump, geometry: MeridionalGeometry=MeridionalGeometry()):
		"""
			Shroud: inlet radius d_1/2 -> outlet at (Z_e, d_2/2)
			Hub: inlet radius d_hub/2 -> outlet at (Z_e + b_2, d_2/2)
			Each contour reaches the outlet at its outlet angle from radial, positive leaning downstream (+z).
			Near the outlet the passage width is b_2 + (r_2 - r) * (tan(shroud angle) - tan(hub angle)),
			so the passage area 2 pi r b does not increase toward the outlet when
				tan(shroud_outlet_angle) - tan(hub_outlet_angle) >= b_2 / r_2
			Shrouds are typically tilted 2-8 deg from radial (Milos, Impeller Design Using CAD Techniques
			and Conformal Mapping Method).

			The hub reaches its outlet angle early and runs straight over the last taper_ratio * (d_2 - d_1)/2
			while the shroud is still bending, so the passage narrows down to b_2.
			Flat lengths are given as a fraction of the inlet blade height (d_1 - d_hub)/2.
			See MeridionalContour for the meaning of fullness.
		"""

		self.Pump = pump
		impeller = pump.Impeller

		# pseudo b_1 not true b_1
		b_1 = (impeller.d_1 - impeller.d_hub) / 2	# inlet blade height
		r_2 = impeller.d_2 / 2
		dr = (impeller.d_2 - impeller.d_1) / 2		# radial extent of the shroud

		self.Z_e: Q_[float] = (		# axial extension, Gulich, pg. 348
			geometry.Z_e_factor * ((impeller.d_2 - impeller.d_1)*((impeller.n_q/74)**(1.07)))
			)

		self.ShroudContour = MeridionalContour(
			r_in=impeller.d_1/2, z_out=self.Z_e, r_out=r_2, outlet_angle=geometry.shroud_outlet_angle,
			flat=geometry.shroud_flat * b_1, fullness=geometry.shroud_fullness
		)
		self.HubContour = MeridionalContour(
			r_in=impeller.d_hub/2, z_out=self.Z_e + impeller.b_2, r_out=r_2, outlet_angle=geometry.hub_outlet_angle,
			r_bend=r_2 - geometry.hub_taper_ratio * dr, flat=0 * b_1, fullness=geometry.hub_fullness
		)
		# geometric meanline, midway between the hub and shroud at equal normalized arclength, an approximation of the
		# mean (phi = 0.5) streamline
		self.MeanContour = MeanCurve(self.HubContour, self.ShroudContour)

		self.PotentialFlow = PotentialFlow(
			mesh=FlowDomain(self.ShroudContour, self.HubContour).Mesh(by_flow=False),
			flow_Q=self.Pump.DesignFlowRate
		)

		# leading edge in the meridional view, on the hub, the meanline and the shroud: all three have the same static
		# moment about the axis, S = integral(r dx) from the leading edge to the impeller outlet. S sets Pfleiderer's
		# slip coefficient p = psi' r_2^2 / (z S) (Milos eqs. 34-35), so the slip is the same from hub to shroud.
		s_LE_shroud = 0									# reference, normalized arclength along the shroud: the impeller inlet
		static_moment_ref = self.ShroudContour.StaticMoment(s_LE_shroud)

		def GetLeadingEdge(contour: MeridionalCurve) -> float:
			"""
				Normalized arclength along contour where its static moment to the impeller outlet matches the reference,
				by bisection. The static moment shrinks as the leading edge moves toward the impeller outlet, so the
				solution always lies between s_a (static moment too large) and s_b (too small), and halving that
				interval converges to it.
			"""
			s_a, s_b = 0, 1
			assert contour.StaticMoment(s_a) >= static_moment_ref, \
				'Contour is too short to match the static moment of the shroud'
			while s_b - s_a > 1e-4:
				s = (s_a + s_b) / 2
				if contour.StaticMoment(s) > static_moment_ref:
					s_a = s
				else:
					s_b = s
			return (s_a + s_b) / 2

		# leading edge points, hub to shroud
		contours = (self.HubContour, self.MeanContour, self.ShroudContour)
		s_LE = (GetLeadingEdge(self.HubContour), GetLeadingEdge(self.MeanContour), s_LE_shroud)
		zr_LE = np.array([[q.m for q in contour(s)] for contour, s in zip(contours, s_LE)])	# (3, 2) (z, r) in m
		
		# quadratic Bezier through the three points, hub to shroud (BezierCurve.ThroughPoints)
		self.LeadingEdge: BezierCurve = BezierCurve.ThroughPoints(*zr_LE)
		
		# trailing edge in the meridional view: a straight line across the impeller outlet, hub to shroud
		zr_TE = np.array([[q.to('m').m for q in contour(1)] for contour in (self.HubContour, self.ShroudContour)])	# (2, 2) (z, r) in m
		self.TrailingEdge: Polyline = Polyline(zr_TE)

	def MeridionalTraces(self, l: str = 'in', color: str = 'black') -> list[go.Scatter]:
		"""
			Meridional cross-section as plotly traces in length unit l, z = 0 at the impeller inlet: hub and shroud
			contours, leading edge, geometric meanline, inlet and outlet. Drawn by PlotMeridional, and with the inducer
			by Pump.PlotMeridional.
		"""

		d_1	= self.Pump.Impeller.d_1.to(l)
		Z_e = self.Z_e.to(l)
		b_2 = self.Pump.Impeller.b_2.to(l)
		d_2 = self.Pump.Impeller.d_2.to(l)
		d_hub = self.Pump.Impeller.d_hub.to(l)

		return [
			go.Scatter(
				x=[0, 0], y=[d_hub.m/2, d_1.m/2],
				mode='lines', name='Inlet', line=dict(dash='dot', color='grey'),
				showlegend=False
			),
			go.Scatter(
				x=[Z_e.m, Z_e.m + b_2.m],
				y=[d_2.m/2, d_2.m/2],
				mode='lines', name='Outlet', line=dict(dash='dot', color='grey'),
				showlegend=False
			),
			self.ShroudContour.Trace(l, name='Shroud Contour', line=dict(color=color), showlegend=False),
			self.HubContour.Trace(l, name='Hub Contour', line=dict(color=color), showlegend=False),
			# leading edge, with markers on its hub, meanline and shroud points (see BuildMeridian)
			self.LeadingEdge.Trace(l, name='Leading Edge', legendgroup='LE', line=dict(color=color)),
			self.MeanContour.Trace(l, name='Geometric Meanline', line=dict(color='black', dash='dashdot'))
		]

	def PlotMeridional(self, l: str = 'in', color: str = 'black') -> go.Figure:
		"""
			Meridional cross-section (MeridionalTraces) and the axis, in its own figure.
		"""

		fig = go.Figure(self.MeridionalTraces(l, color))

		fig.add_trace(go.Scatter(
			x=[0, self.Z_e.to(l).m + self.Pump.Impeller.b_2.to(l).m], y=[0, 0],
			mode='lines', name='Meridional axis', line=dict(dash='dashdot', color='grey')
		))

		fig.update_xaxes(title_text=f"Meridional position ({l})")
		fig.update_yaxes(title_text=f"Radial position ({l})", scaleanchor='x', scaleratio=1)

		fig.update_layout(
			title="Impeller Meridional Cross-Section",
			# legend in the empty upper-left corner
			legend=dict(x=0.01, y=0.99, xanchor='left', yanchor='top', bgcolor='rgba(255,255,255,0.8)'),
			template="plotly_white"
		)

		return fig

