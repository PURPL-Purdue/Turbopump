import cadquery as cq
from cadquery.vis import show
from cadquery.func import spline, revolve
import numpy as np
from pint import Quantity as Q_
from dataclasses import dataclass
import plotly.graph_objects as go
import plotly.subplots
import EmpiricalRelations as emp
from VelocityTriangle import VelocityTriangle
from MeridionalContour import MeridionalContour, MeridionalCurve, FlowDomain, MeanCurve, BezierCurve, Polyline
from PotentialFlow import PotentialFlow

g = Q_(9.81, 'm/s^2')

@dataclass
class InputGeometry:
	Z_blade: int			# number of impeller blades
	Beta2B: Q_[float]		# blade backsweep angle at outlet (from tangent)
	d_1: Q_[float]			# impeller inlet diameter
	d_2: Q_[float]			# impeller outlet diameter
	b_2: Q_[float]			# impeller outlet height
	blade_thick: Q_[float]			# blade thickness at exit
	d_hub: Q_[float]		# diameter of hub, strictly speaking inlet only

	shroud_fullness: tuple[float, float] = (1, 0.3)
	hub_fullness: tuple[float, float] = (0.5, 0.7)
	# contour angles from radial at the outlet, positive leaning downstream (+z), see BuildMeridian
	shroud_outlet_angle: Q_[float] = Q_(5, 'deg')
	hub_outlet_angle: Q_[float] = Q_(0, 'deg')
	hub_taper_ratio: float = 0.3

	Z_e_factor: float = 1

@dataclass
class DesignPoint:
	Q: Q_[float]			# volumetric flowrate
	N_shaft: Q_[float]		# shaft rotating speed
	H: Q_[float]			# headrise, developed head
	n_hyd_BEP: float		# hydraulic efficiency at BEP
	rho: Q_[float]			# average density of the fluid

# Wiesner slip factor: https://manual.cfturbo.com/en/bl_te_wiesner.html
def GetWiesnerSlipRatio(beta2b: Q_[float], Z: int) -> Q_[float]:
	return 1 - np.sqrt(np.sin(beta2b)) / Z**0.7

def GetBladeBlockage(beta2b: Q_[float], Z: int, thick: Q_[float], b2: Q_[float]) -> Q_[float]:
	return thick * b2 * Z / np.sin(beta2b)

class Impeller:

	def __init__(self, geometry: InputGeometry, design_point: DesignPoint):

		assert geometry.d_1 > geometry.d_hub, 'Inlet diameter must be larger than hub diameter'
		
		self.DP: DesignPoint = design_point

		self.SpecificSpeed: float = (	# imperial specific speed
					self.DP.N_shaft.to('rpm') * np.sqrt(self.DP.Q.to('gallon/min')) /
					self.DP.H.to('ft')**0.75
					).magnitude
		
		self.n_q: float = (		# metric specific speed
			self.DP.N_shaft.to('rpm') * np.sqrt(self.DP.Q.to('m^3/s')) /
			self.DP.H.to('m')**0.75
			).magnitude
		
		exit_area: Q_[float] = (np.pi * geometry.d_2 * geometry.b_2 - GetBladeBlockage(
			geometry.Beta2B, geometry.Z_blade, geometry.blade_thick, geometry.b_2
		))

		self.Z_blade = geometry.Z_blade
		self.Beta2B = geometry.Beta2B
		self.d_1 = geometry.d_1
		self.d_2 = geometry.d_2
		self.b_2 = geometry.b_2
		self.blade_thick = geometry.blade_thick
		self.d_hub = geometry.d_hub
		self.WiesnerSlip = GetWiesnerSlipRatio(geometry.Beta2B, geometry.Z_blade)
		self.Area2 = exit_area
		self.Z_e: Q_[float] = (		# axial extension, Gulich, pg. 348
			geometry.Z_e_factor*((self.d_2 - self.d_1)*((self.n_q/74)**(1.07)))
			)

		self.ShroudContour: MeridionalContour
		self.HubContour: MeridionalContour
		self.MeanContour: MeridionalCurve
		self.BuildMeridian(
			hub_fullness=geometry.hub_fullness,
			shroud_fullness=geometry.shroud_fullness,
			shroud_outlet_angle=geometry.shroud_outlet_angle,
			hub_outlet_angle=geometry.hub_outlet_angle,
			taper_ratio=geometry.hub_taper_ratio
		)

		# exit tip velocity at design point
		self.U_2_design = (self.DP.N_shaft * self.d_2/2).to('m/s')
		# meridional exit flow velocity at design point
		self.C_m2_design = (self.DP.Q / self.Area2).to('m/s')

		_, enthalpy, _ = self.GetSpecificWork()
		# AKA Stage loading: Δh / u^2 = g*H / u^2
		self.HeadCoeff: float = (enthalpy/self.U_2_design**2).to('dimensionless').magnitude
		# AKA Flow factor: c / u
		self.FlowCoeff: float = (self.C_m2_design / self.U_2_design).to('dimensionless').magnitude

	def BuildMeridian(self,
			shroud_flat_ratio: float = 0.2, hub_flat_ratio: float = 0.0, taper_ratio: float = 0.3,
			shroud_fullness: tuple[float, float] = (1, 0.3),
			hub_fullness: tuple[float, float] = (0.5, 0.7),
			shroud_outlet_angle: Q_[float] = Q_(5, 'deg'),
			hub_outlet_angle: Q_[float] = Q_(0, 'deg')):
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
		b_1 = (self.d_1 - self.d_hub) / 2	# inlet blade height
		r_2 = self.d_2 / 2
		dr = (self.d_2 - self.d_1) / 2		# radial extent of the shroud

		self.ShroudContour = MeridionalContour(
			r_in=self.d_1/2, z_out=self.Z_e, r_out=r_2, outlet_angle=shroud_outlet_angle,
			flat=shroud_flat_ratio * b_1, fullness=shroud_fullness
		)
		self.HubContour = MeridionalContour(
			r_in=self.d_hub/2, z_out=self.Z_e + self.b_2, r_out=r_2, outlet_angle=hub_outlet_angle,
			r_bend=r_2 - taper_ratio * dr, flat=hub_flat_ratio * b_1, fullness=hub_fullness
		)
		# geometric meanline, midway between the hub and shroud at equal normalized arclength, an approximation of the
		# mean (phi = 0.5) streamline
		self.MeanContour = MeanCurve(self.HubContour, self.ShroudContour)

		self.PotentialFlow = PotentialFlow(
			mesh=FlowDomain(self.ShroudContour, self.HubContour).Mesh(by_flow=False),
			flow_Q=self.DP.Q
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
		
		hub_LE, shroud_LE = zr_LE[0], zr_LE[len(zr_LE) - 1]
		LE_line = hub_LE - shroud_LE
		self.lambdaLE = Q_(np.atan(-LE_line[1] / LE_line[0]), 'rad').to('deg')
		
		# quadratic Bezier through the three points, hub to shroud (BezierCurve.ThroughPoints)
		self.LeadingEdge: BezierCurve = BezierCurve.ThroughPoints(*zr_LE)
		
		# trailing edge in the meridional view: a straight line across the impeller outlet, hub to shroud
		zr_TE = np.array([[q.to('m').m for q in contour(1)] for contour in (self.HubContour, self.ShroudContour)])	# (2, 2) (z, r) in m
		self.TrailingEdge: Polyline = Polyline(zr_TE)

	@property
	def Area1(self) -> Q_[float]:
		return np.pi / 4 * (self.d_1**2 - self.d_hub**2) # TODO: blade blockage?
	@property
	def NPSH_i(self) -> Q_[float]:
		"""
			This method of calculating NPSH_i accounts for inlet velocity and losses.
			Therefore, when comparing with available NSPH, total pressure should be used for NPSH_a calculation
			NPSH_i = (lambda_c * 1/2 * c_m^2 + lambda_w * 1/2 * w^2) / g

			The physical intuition for the NPSH coefficients is that if the coefficients were equal to 1,
			then the drop in static pressure due to suction is merely the dynamic pressure due to
			acceleration of the flow at the inlet due to inlet area reduction (c_m, merdional velocity)
			and the acceleration of the flow over and around the blade leading edge (w, relative velocity)

		"""

		# Gulich, Centrifugal Pumps, chap. 6
		lambda_c = 1.2	# NPSH coefficient for main flow acceleration and losses at inlet
		lambda_w = 1.0	# NPSH coefficient for excess velocity due to flow around leading edge
		
		c_1m = self.DP.Q / self.Area1	# meridional velocity at inlet
		w_1 = np.sqrt(					# relative velocity of flow at leading edge
			(self.DP.N_shaft * self.d_1/2)**2 + c_1m**2
			)
		return (1/2 / g * (lambda_c * c_1m**2 + lambda_w * w_1**2)).to('m')

	def GetMassFlow(self, volumetric_flowrate: Q_[float]) -> Q_[float]:
		volumetric_flowrate = volumetric_flowrate if volumetric_flowrate is not None else self.DP.Q
		return volumetric_flowrate * self.DP.rho
	
	def H_euler(self, Q):
		"""
		https://ntrs.nasa.gov/api/citations/19950013379/downloads/19950013379.pdf
		page 4 and 5
		
		Theoretical Euler head as a function of Q
		H_theoretical(Q) = U_2^2/g - [U_2 / (g*Area2*tan(beta2b))] * Q
		Or
		H_theoretical(Q) = U_2^2/g - [N_shaft*cot(beta2b)/(2*pi*b2*g)] * Q
		https://youtu.be/H7XfYO_-cEg?t=2705
		"""
		sigma = self.WiesnerSlip
		U_2_design = self.U_2_design
		N = self.DP.N_shaft
		beta2b = self.Beta2B
		b2 = self.b_2

		return sigma * U_2_design**2/g - N/np.tan(beta2b)/(2*np.pi*b2*g) * Q

	def SweepInletDiam(self) -> go.Figure:

		npsh = []
		diam = []
		min_NPSH = Q_(np.inf, 'm')
		opt_diam = Q_(-1, 'in')

		old_d_1 = self.d_1

		for d_1 in np.linspace(self.d_hub, self.d_2, 100):
			self.d_1 = d_1
			NPSH_i = self.NPSH_i
			
			if NPSH_i < min_NPSH:
				min_NPSH = NPSH_i.to('m')
				opt_diam = d_1.to('in')
			
			npsh.append(NPSH_i.magnitude)
			diam.append(d_1.to('in').magnitude)
			
		self.d_1 = old_d_1

		fig = go.Figure()
		fig.add_trace(go.Scatter(
			x=diam, y=npsh,
			mode='lines', name='NPSH_i'
		))
		fig.add_trace(go.Scatter(
			x=[opt_diam.magnitude], y=[min_NPSH.magnitude],
			mode='markers',
			name=f'Optimal: {opt_diam:.2f}, {min_NPSH:.0f}',
			marker=dict(size=9, symbol='circle')
		))
		fig.update_layout(
			title='NPSH inception vs inlet diameter',
			xaxis_title='Inlet diameter (in.)',
			yaxis_title='NPSH_i (m)',
			template='plotly_white'
		)
		fig.update_xaxes(showgrid=True)
		fig.update_yaxes(showgrid=True)
		fig.show()

		return fig

	def GetC_u2(self, flow_Q: Q_[float]=None, speed_N: Q_[float]=None) -> Q_[float]:
		flow_Q = flow_Q if not flow_Q is None else self.DP.Q
		speed_N = speed_N if not speed_N is None else self.DP.N_shaft
		
		u_2 = speed_N * self.d_2/2
		c_m2 = flow_Q / self.Area2
		slip = u_2 * (1 - self.WiesnerSlip)
		w_u2 = c_m2 / np.tan(self.Beta2B) + slip
		w_2 = np.sqrt(c_m2**2 + w_u2**2)
		c_u2 = (u_2 - w_u2)
		return c_u2

	def GetC_u1(self, flow_Q: Q_[float]=None, speed_N: Q_[float]=None) -> Q_[float]:
		return Q_(0, 'm/s')

	def GetSpecificWork(self, speed_N: Q_[float]=None, flow_Q: Q_[float]=None) -> tuple[Q_[float], Q_[float], float]:
		"""
			Returns
			1. Brake work, work consumed by the impeller per unit mass of fluid
			2. Fluid work, work done on the fluid per unit mass of fluid
			3. flowspeed factor, non-dimensional number, low confidence of analysis if greater than 2
		"""
		if speed_N is None:
			speed_N = self.DP.N_shaft
		if flow_Q is None:
			flow_Q = self.DP.Q

		#vels: VelocityTriangle
		#f: float
		c_u1, u_1 = self.GetC_u1(flow_Q, speed_N), (self.d_1 + self.d_hub)/4 * speed_N
		c_u2, u_2 = self.GetC_u2(flow_Q, speed_N), self.d_2/2 * speed_N
		# Euler turbomachinery equation
		brake_work = (c_u2 * u_2 - c_u1 * u_1).to('kJ/kg')

		f=emp.FlowSpeedRatio(flow_Q, speed_N, self.DP.Q, self.DP.N_shaft)
		# empirical prediction, valid for 0 < f < 2
		n_hydraulic = emp.HydraulicEfficiency(f) * self.DP.n_hyd_BEP

		fluid_work = brake_work * n_hydraulic

		return (brake_work, fluid_work, f)
		
	def PlotPerformanceHQ(self, rpm_sweep: Q_[list[float]]=None) -> go.Figure:
		N_sweep = rpm_sweep if rpm_sweep is not None else Q_([self.DP.N_shaft])
		Q_sweep = Q_(np.linspace(0, 2.0 * self.DP.Q.to('L/s').magnitude, 50), 'L/s')

		H_predict = []
		Q_predict = []

		for speed_N in N_sweep:

			head_trace = []
			flow_trace = []

			for flow_Q in Q_sweep:
				
				_, fluid_work, f = self.GetSpecificWork(speed_N=speed_N, flow_Q=flow_Q)
				if f > 2: continue

				headrise = fluid_work / g
		
				head_trace.append(headrise.to('m').magnitude)
				flow_trace.append(flow_Q.to('L/s').magnitude)
		
			H_predict.append(head_trace)
			Q_predict.append(flow_trace)

		fig = go.Figure()
		# Plot theoretical Euler (linear)
		H_theoretical = self.H_euler(Q_sweep)
		fig.add_trace(go.Scatter(
			x=Q_sweep.to('L/s').magnitude,
			y=H_theoretical.to('m').magnitude,
			mode='lines', name='Theoretical Euler',
			line=dict(color='gray', dash='dash')
		))

		# Plot RPM swept traces
		for i in range(len(H_predict)):
			fig.add_trace(go.Scatter(
				x=Q_predict[i], y=H_predict[i],
				mode='lines',
				name='Empirical prediction' if i == 0 else f'Empirical prediction ({N_sweep[i].magnitude:.0f} RPM)',
				line=dict(color='red'),
				showlegend=(i == 0)
			))
			if Q_predict[i] and H_predict[i]:
				fig.add_annotation(
					x=Q_predict[i][-1], y=H_predict[i][-1],
					text=f"{N_sweep[i].magnitude:.0f} RPM",
					showarrow=False, xshift=25, yshift=0,
					font=dict(size=10)
				)

		# Plot design point head and flowrate
		fig.add_trace(go.Scatter(
			x=[self.DP.Q.to('L/s').magnitude],
			y=[self.DP.H.to('m').magnitude],
			mode='markers', name='Design point',
			marker=dict(color='black', size=7)
		))

		fig.update_layout(
			title='H-Q Curve',
			xaxis_title='Volumetric flow rate [L/s]',
			yaxis_title='Head [m]',
			yaxis=dict(rangemode='tozero'),
			legend=dict(x=1, y=1, xanchor='right', yanchor='top'),
			template='plotly_white'
		)
		fig.update_xaxes(showgrid=True, gridcolor='rgba(128,128,128,0.3)')
		fig.update_yaxes(showgrid=True, gridcolor='rgba(128,128,128,0.3)')
		fig.show()

		return fig

	def MeridionalTraces(self, l: str = 'in', color: str = 'black') -> list[go.Scatter]:
		"""
			Meridional cross-section as plotly traces in length unit l, z = 0 at the impeller inlet: hub and shroud
			contours, leading edge, geometric meanline, inlet and outlet. Drawn by PlotMeridional, and with the inducer
			by Pump.PlotMeridional.
		"""

		d_1	= self.d_1.to(l)
		Z_e = self.Z_e.to(l)
		b_2 = self.b_2.to(l)
		d_2 = self.d_2.to(l)
		d_hub = self.d_hub.to(l)

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
			The FEM mesh is shown over the potential flow solution by self.PotentialFlow.Plot.
		"""

		fig = go.Figure(self.MeridionalTraces(l, color))

		fig.add_trace(go.Scatter(
			x=[0, self.Z_e.to(l).m + self.b_2.to(l).m], y=[0, 0],
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

		fig.show()

		return fig

	def GetVelocities(self, station: int, phi: float=0.5) -> VelocityTriangle:
		assert (station == 1 or station == 2), 'Station must be 1 or 2'
		streamline = self.PotentialFlow.GetStreamline(phi)
		edge_curve = (self.LeadingEdge if station == 1 else self.TrailingEdge)
		intersects = edge_curve.FindIntersect(streamline)

		assert len(intersects) == 1, 'It should have exactly one intersect'

		_, r, _, s_2 = intersects[0]
		c_m: Q_[float] = streamline.Interpolate(streamline.v_m, s_2)
		u = (r * self.DP.N_shaft).to('m/s')

		if station == 1:
			# assume no inlet swirl (tangential velocity is 0)
			c_u = self.GetC_u1()
		else:
			c_u = self.GetC_u2()
		#BetaB = Q_(np.atan(c_m / (u - c_u)), 'rad').to('deg')
		w = np.sqrt(c_m**2 + (u - c_u)**2)
		
		vel_tri = VelocityTriangle(
			u=u, c_m=c_m, c_u=c_u, w = w, station=station, beta_b=(self.Beta2B if station==2 else None)
		)
		#vel_tri.Plot()
		
		return vel_tri

	def SweepVelocities(self, phis: list[float] | np.ndarray=[0, 0.5, 1]) -> go.Figure:
		grid = plotly.subplots.make_subplots(
			rows=len(phis), cols=2,
			subplot_titles=[f'{label} Velocities at Flow Fraction ϕ = {phi}' for phi in phis for label in ['Inlet', 'Outlet']]
			)
		for i in range(len(phis)):
			phi = phis[i]
			for station in [1, 2]:
				velocity_triangle = self.GetVelocities(station, phi)
				velocity_triangle.Plot(
					convention='centrifugal',
					station=station,
					fig=grid, row=i + 1, col=station
				)

		grid.update_layout(
			title='Inlet & Outlet Velocities at Various Flow Fractions',
			template='plotly_white'
			)
		return grid

	def ShroudHubRevolve(self, resolution: int=200):
		
		s = np.linspace(0, 1, resolution)
		hub_z, hub_r = self.HubContour(s)
		shroud_z, shroud_r = self.ShroudContour(s)
		hub = list(zip(hub_z.to('mm').m, hub_r.to('mm').m))              # [(axial, radial), ...]
		shroud = list(zip(shroud_z.to('mm').m, shroud_r.to('mm').m))

		part = (cq.Workplane("YZ")
			.spline(hub)										# hub: inlet → outlet
			.lineTo(*shroud[-1])								# outlet line across b_2
			.spline(shroud[::-1][1:], includeCurrent=True)		# shroud: outlet → inlet
			.close()											# inlet line back to the hub
			.revolve(360, (1, 0, 0), (0, 0, 0))					# revolve about the sketch x-axis
		)
		show(part, alpha = 0.8)

		# hub_wire = spline(*[cq.Vector(0, z, r) for z, r in hub])   # points in the global YZ plane
		# hub_face = revolve(hub_wire, (0, 0, 0), (0, 1, 0), 360)    # axis = global Y

		# shroud_wire = spline(*[cq.Vector(0, z, r) for z, r in shroud])   # points in the global YZ plane
		# shroud_face = revolve(shroud_wire, (0, 0, 0), (0, 1, 0), 360)    # axis = global Y
		# #show(hub_face, shroud_face, alpha = 0.8)