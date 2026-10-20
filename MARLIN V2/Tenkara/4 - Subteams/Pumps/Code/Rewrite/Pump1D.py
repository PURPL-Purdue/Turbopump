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
	print(Q_flow, N_speed)
	u_2, _, c_u2, _ = pump.Impeller.GetOutletVelocities(Q_flow, N_speed)
	return (u_2 * c_u2).to('kJ/kg')

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
		Q_sweep_max = (rpm / pump.DP.N_shaft * f_max * pump.DesignFlowRate).to('L/s')
		flow_sweep = Q_(np.linspace(0, Q_sweep_max.m), 'L/s')
		head_sweep = GetActualHeadRise(pump, flow_sweep, rpm)
		fig.add_trace(go.Scatter(
			x=flow_sweep.m, y=head_sweep.m,
			mode='lines', name='Headrise prediction (empirical correlation)', showlegend=False
		))

		fig.add_annotation(
			x=Q_sweep_max.m, y=head_sweep[len(head_sweep)-1].m,
			text=f"{rpm.m:d} RPM",
			showarrow=False, xshift=40, yshift=0,
			font=dict(size=10)
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

if __name__ == '__main__':
	import ConstructHardware
	lox_pump = ConstructHardware.LOXPump()
	#lox_pump.Impeller.OutletVelocityTriangleAtDesign.Plot().show()
	PlotHQ(lox_pump).show()