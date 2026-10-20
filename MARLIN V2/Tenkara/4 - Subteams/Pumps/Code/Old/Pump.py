import numpy as np
import plotly.graph_objects as go
from pint import Quantity as Q_
from Inducer import Inducer
from Impeller import Impeller
import ConstructHardware

class Pump:
	"""
		Inducer followed by the impeller on the same shaft.
		Meridional coordinates are the impeller's, z = 0 at the impeller inlet: the inducer lies upstream at negative z,
		its outlet axial_gap before the impeller inlet.
	"""

	def __init__(self, inducer: Inducer, impeller: Impeller, axial_gap: Q_[float]=Q_(0, 'in')):
		self.Inducer = inducer
		self.Impeller = impeller
		self.axial_gap = axial_gap		# from the inducer outlet to the impeller inlet

	def PlotMeridional(self, l: str = 'in', color: str = 'black') -> go.Figure:
		"""
			Meridional cross-section of the inducer and the impeller in one figure: Inducer.MeridionalTraces moved
			upstream, Impeller.MeridionalTraces as they are, and the axis under both.
		"""

		fig = go.Figure()

		z_inducer = -(self.Inducer.L_ax + self.axial_gap).to(l).m		# inducer inlet
		for trace in self.Inducer.MeridionalTraces(l, color):
			trace.x = np.asarray(trace.x) + z_inducer
			fig.add_trace(trace)
		fig.add_traces(self.Impeller.MeridionalTraces(l, color))

		fig.add_trace(go.Scatter(
			x=[z_inducer, (self.Impeller.Z_e + self.Impeller.b_2).to(l).m], y=[0, 0],
			mode='lines', name='Meridional axis', line=dict(dash='dashdot', color='grey')
		))

		fig.update_xaxes(title_text=f"Meridional position ({l})")
		fig.update_yaxes(title_text=f"Radial position ({l})", scaleanchor='x', scaleratio=1)

		fig.update_layout(
			title="Pump Meridional Cross-Section",
			# legend in the empty upper-left corner, above the inducer
			legend=dict(x=0.01, y=0.99, xanchor='left', yanchor='top', bgcolor='rgba(255,255,255,0.8)'),
			template="plotly_white"
		)

		fig.show()

		return fig

def LOX() -> Pump:
	return Pump(
		inducer=ConstructHardware.InducerLOX(),
		impeller=ConstructHardware.ImpellerLOX()
		)

def IPA() -> Pump:
	return Pump(
		inducer=ConstructHardware.InducerIPA(),
		impeller=ConstructHardware.ImpellerIPA()
		)
