from Pump import Pump
from plotly import graph_objects as go
import numpy as np
from pint import Quantity as Q_
from typing import Optional
import EmpiricalRelations

g = Q_(9.81, 'm/s^2')

def GetSpecificWork(pump: Pump, Q_flow: Optional[Q_]=None, N_speed: Optional[Q_]=None) -> Q_:
	# assuming no preswirl and perfect matching between inducer and impeller
	# in which case the total headrise across the pump becomes:
	# H = [u_2 * c_u2] / g		(station 2 denotes impeller outlet)
	if Q_flow is None: Q_flow = pump.DesignFlowRate
	if N_speed is None: N_speed = pump.DP.N_shaft
	u_2, _, c_u2, _ = pump.Impeller.GetOutletVelocities(Q_flow, N_speed)
	return (u_2 * c_u2).to('kJ/kg')

def GetPowerConsumption(pump: Pump, Q_flow: Optional[Q_]=None, N_speed: Optional[Q_]=None) -> Q_:
	return (pump.DP.mass_flow * GetSpecificWork(pump, Q_flow, N_speed)).to('kW')

def GetTorque(pump: Pump, Q_flow: Optional[Q_]=None, N_speed: Optional[Q_]=None) -> Q_:
	return (GetPowerConsumption(pump, Q_flow, N_speed) / N_speed).to('kN*m')

def GetEulerHeadRise(pump: Pump, Q_flow: Optional[Q_]=None, N_speed: Optional[Q_]=None) -> Q_:
	return (GetSpecificWork(pump, Q_flow, N_speed) / g).to('m')

def GetActualHeadRise(pump: Pump, Q_flow: Optional[Q_]=None, N_speed: Optional[Q_]=None) -> Q_:

	if Q_flow is None: Q_flow = pump.DesignFlowRate
	if N_speed is None: N_speed = pump.DP.N_shaft

	f = (EmpiricalRelations.FlowSpeedRatio(Q_flow, N_speed, pump.DesignFlowRate, pump.DP.N_shaft)).to('dimensionless').m
	euler_head = GetEulerHeadRise(pump, Q_flow, N_speed)
	eta_hydraulic = pump.Impeller.eta_hyd * EmpiricalRelations.HydraulicEfficiency(f)

	return eta_hydraulic * euler_head

def PlotHQ(pump: Pump, rpm_step_down: float=5000, rpm_sweeps: int=4) -> go.Figure:

	f_max = 2	# for flow speed ratio > 2, the emperical correlation has low accuracy

	fig = go.Figure()

	for i in range(rpm_sweeps):
		rpm = pump.DP.N_shaft - Q_(i * rpm_step_down, 'rpm')
		# DNE flowrate, such that empirical flow speed ratio < 2
		Q_sweep_max = (rpm / pump.DP.N_shaft * f_max * pump.DesignFlowRate).to('L/s')

		flow_sweep = Q_(np.linspace(0, Q_sweep_max.m), 'L/s')
		head_sweep = GetActualHeadRise(pump, flow_sweep, rpm)

		fig.add_trace(go.Scatter(
			x=flow_sweep.m, y=head_sweep.m,
			mode='lines', name='Headrise prediction (empirical correlation)', showlegend=False
		))

		fig.add_annotation(
			x=Q_sweep_max.m, y=head_sweep[len(head_sweep)-1].m, text=f"{rpm.m:d} RPM",
			showarrow=False, xshift=40, yshift=0, font=dict(size=12)
		)

	Q_euler_max = pump.DesignFlowRate.to('L/s') * f_max
	H_euler = GetEulerHeadRise(pump, Q_([0, Q_euler_max.m], 'L/s'))

	fig.add_trace(go.Scatter(
		x=[0, Q_euler_max.m],
		y=H_euler.m,
		mode='lines', name='Theoretical Euler',
		line=dict(color='gray', dash='dash')
	))

	# Plot design point head and flowrate
	fig.add_trace(go.Scatter(
		x=[pump.DesignFlowRate.to('L/s').m],
		y=[pump.DesignHeadRise.to('m').m],
		mode='markers', name='Design point',
		marker=dict(color='black', size=7)
	))

	fig.update_layout(
		title='Head Rise vs Flowrate (H-Q) Curve',
		xaxis_title='Volumetric flow rate [L/s]',
		yaxis_title='Head [m]',
		yaxis=dict(rangemode='tozero'),
		legend=dict(x=1, y=1, xanchor='right', yanchor='top'),
		template='plotly_white'
	)

	return fig

def GuessInducerlessNPSH_r(pump: Pump, lambda_c: float=1.1, lambda_w: float=1.0, blade_blockage_factor: float=0.7) -> Q_:

	"""
		This method of calculating NPSH_i accounts for inlet velocity and losses.
		Therefore, when comparing with available NSPH, total pressure should be used for NPSH_a calculation
		NPSH_i = (lambda_c * 1/2 * c_m^2 + lambda_w * 1/2 * w^2) / g

		The physical intuition for the NPSH coefficients is that if the coefficients were equal to 1,
		then the drop in static pressure due to suction is merely the dynamic pressure due to
		acceleration of the flow at the inlet due to inlet area reduction (c_m, merdional velocity)
		and the acceleration of the flow over and around the blade leading edge (w, relative velocity)

		https://manual.cfturbo.com/en/md_parameters_pumpvent.html
	"""
	impeller = pump.Impeller
	area_1 = (impeller.d_1**2 - impeller.d_hub**2)/4 * np.pi * blade_blockage_factor
	# assume no preswirl, c_m1 = c_1
	c_m1 = impeller.DP.Q / area_1
	u_1 = impeller.DP.N_shaft * impeller.d_1/2
	w_1: Q_[float] = np.sqrt(c_m1**2 + u_1**2)

	return (1/(2 * g) * (lambda_c * c_m1**2 + lambda_w * w_1**2)).to('m')

def GetInducerHeadRequirement(pump: Pump) -> Q_:
	no_inducer_NPSH_r: Q_[float] = GuessInducerlessNPSH_r(pump)
	difference = (no_inducer_NPSH_r - pump.NPSH_a).to('m')
	# round to the nearest 5 m
	# this is to prevent the inducer requirements from changing too much with small changes in design point
	return Q_(5 * round(difference.m / 5), 'm')

if __name__ == '__main__':
	import Main
	lox_pump = Main.LOXPump()
	#lox_pump.Impeller.OutletVelocityTriangleAtDesign.Plot().show()
	#PlotHQ(lox_pump).show()
	GuessInducerlessNPSH_r(lox_pump)