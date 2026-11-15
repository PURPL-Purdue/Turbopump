import numpy as np
import plotly.graph_objects as go
from pint import Quantity as Q_
# Inputs
d_1 = 0.02413
d_h = 0.01016	# [m] tip diameter
lambda_w = 0.06
lambda_c = 1.1
Q = 0.001771	# [m^3/s]
n = 35000		# [rpm]
g = 9.81
H = 35

u1 = np.pi * n * d_1 / 60
cm_1 = (4 * Q) / (np.pi * (d_1**2 - d_h**2))
w_1 = np.sqrt(cm_1**2 + u1**2)
Npsh = lambda_c * (cm_1**2 / (2 * g)) + lambda_w * (w_1**2 / (2 * g))
flow_coeff = cm_1 / u1
nss = n*np.sqrt(Q)/(Npsh**(3/4))
tau = 2*g*Npsh/(u1**2)
K = (tau - flow_coeff**2)/(1 + flow_coeff**2)
psi = H / (u1**2 / g)
print(Npsh, flow_coeff, nss, psi)
#print(K)

def PlotContour(d_tip: Q_, d_hub_1: Q_, d_hub_2: Q_, L_ax: Q_ = None, l: str = 'mm', color: str = 'blue', name: str = 'Inducer') -> go.Figure:
	"""Meridional section from given diameters, z = 0 at the inlet. L_ax defaults to 0.6 * d_tip."""
	if L_ax is None:
		L_ax = 0.6 * d_tip
	r_t, r_1, r_2, length = (v.to(l).magnitude for v in (d_tip / 2, d_hub_1 / 2, d_hub_2 / 2, L_ax))
	line = lambda x, y, **kw: go.Scatter(x=x, y=y, mode='lines', showlegend=False, **kw)
	fig = go.Figure([
		line([0, length], [r_t, r_t], name='Shroud', line=dict(color=color)),
		line([0, length], [r_1, r_2], name='Hub', line=dict(color=color)),
		line([0, 0], [r_1, r_t], name='Inlet', line=dict(dash='dot', color='grey')),
		line([length, length], [r_2, r_t], name='Outlet', line=dict(dash='dot', color='grey')),
		line([0, length], [0, 0], name='Axis', line=dict(dash='dashdot', color='grey')),
	])
	fig.update_xaxes(title_text=f"Meridional position ({l})")
	fig.update_yaxes(title_text=f"Radial position ({l})", scaleanchor='x', scaleratio=1)
	fig.update_layout(title=f"{name} Inducer Meridional Section", template="plotly_white", showlegend=False)
	return fig

#if __name__ == "__main__":
#	d_tip = Q_(1, 'in')
#	d_hub_1 = Q_(0.55, 'in')
#	d_hub_2 = Q_(0.6, 'in')
#	PlotContour(d_tip, d_hub_1, d_hub_2, name='LOX').show()