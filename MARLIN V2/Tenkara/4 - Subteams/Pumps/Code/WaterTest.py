import Pump1D
import Main
from Pump import Pump
import numpy as np
from pint import Quantity as Q_
from CoolProp import CoolProp as CP
from plotly import graph_objects as go

"""
	Water flow test: flow rate scales linearly with shaft speed (Q/N fixed at the design value), so the velocity
	triangles stay similar (no viscous effects). The inlet is at atmospheric pressure.
"""

T_WATER = Q_(1, 'degC')
P_INLET = Q_(1, 'atm')

def PlotWaterTest(pump: Pump, rho_water: Q_, p_vapor_water: Q_, res: int=100) -> go.Figure:
	speed_sweep: Q_ = np.linspace(0, 1, res) * pump.DP.N_shaft
	flow_sweep: Q_ = np.linspace(0, 1, res) * pump.DesignFlowRate

	power_sweep = Pump1D.GetPowerConsumption(pump, flow_sweep, speed_sweep, rho=rho_water).to('kW')
	NPSH_r_sweep = Pump1D.GetInducerlessNPSH_r(pump, Q_flow=flow_sweep, N_speed=speed_sweep).to('m')
	NPSH_a = ((P_INLET - p_vapor_water) / (rho_water * Pump1D.g)).to('m')

	rpm = speed_sweep.to('rpm').m
	# NPSH_r rises monotonically with speed, so interpolate for the speed where it reaches NPSH_a
	rpm_limit = np.interp(NPSH_a.m, NPSH_r_sweep.m, rpm)
	power_limit = np.interp(rpm_limit, rpm, power_sweep.m)
	print(f'{pump.Propellant.name} water test: NPSH_a = {NPSH_a:.2f~P}, NPSH_r = NPSH_a at {rpm_limit:.0f} rpm, '
		f'{np.interp(rpm_limit, rpm, flow_sweep.to("L/s").m):.2f} L/s, {power_limit:.2f} kW')

	fig = go.Figure()
	fig.add_trace(go.Scatter(x=rpm, y=power_sweep.m, mode='lines', name='Power (Euler)', yaxis='y'))
	fig.add_trace(go.Scatter(x=rpm, y=NPSH_r_sweep.m, mode='lines', name='NPSH_r (inducerless)', yaxis='y2'))
	fig.add_trace(go.Scatter(
		x=[rpm[0], rpm[-1]], y=[NPSH_a.m]*2, mode='lines', name='NPSH_a',
		line=dict(color='black', dash='dash'), yaxis='y2'
	))
	fig.add_trace(go.Scatter(
		x=[rpm_limit], y=[power_limit], mode='markers', name='Power at NPSH_r = NPSH_a',
		marker=dict(color='black', size=7), yaxis='y'
	))

	fig.update_layout(
		title=f'{pump.Propellant.name} pump water test (Q ∝ N, inlet at {P_INLET:~P}, water at {T_WATER:~P})',
		xaxis=dict(title='Shaft speed [rpm]'),
		yaxis=dict(title='Power [kW]', rangemode='tozero'),
		yaxis2=dict(title='NPSH [m]', overlaying='y', side='right', rangemode='tozero', showgrid=False),
		legend=dict(x=0.02, y=0.98, xanchor='left', yanchor='top'),
		template='plotly_white'
	)
	return fig

if __name__ == '__main__':
	T = T_WATER.to('K').m
	rho_water = Q_(CP.PropsSI('D', 'T', T, 'P', P_INLET.to('Pa').m, 'Water'), 'kg/m^3')
	p_vapor_water = Q_(CP.PropsSI('P', 'T', T, 'Q', 0, 'Water'), 'Pa')

	pump = Main.CreatePump(Main.Propellant.IPA)
	PlotWaterTest(pump, rho_water, p_vapor_water).show()