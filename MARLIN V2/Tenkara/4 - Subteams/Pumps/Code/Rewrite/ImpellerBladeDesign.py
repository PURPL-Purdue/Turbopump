import numpy as np
from pint import Quantity as Q_
from Impeller import Impeller
from Inducer import Inducer
from InducerRadialEquilibrium import InducerRadialEquilibrium
from MeridionalDesign import MeridionalDesign
from MeridionalContour import MeridionalCurve
from typing import Optional

class ImpellerBladeDesign:

	def __init__(self, meridional_design: MeridionalDesign, inducer_radial_equilibrium: Optional[InducerRadialEquilibrium]=None):

		self._mer = meridional_design
		self._rad_eq = inducer_radial_equilibrium
		self.K_inlet = inducer_radial_equilibrium.GetOutletAngularMomentum() if inducer_radial_equilibrium is not None else Q_(0, 'm^2/s')

		# angle between the leading edge and the hub contour in the meridional plane. Milos's epsilon_1 (eq. 30), the tilt
		# of the blade surface to the stream surface, follows from it once the blade angle beta_1' is known (own
		# derivation, for a leading edge lying in a meridional plane): tan(epsilon_1) = tan(lambda) / cos(beta_1')
		self._lambda_1_hub = self._GetLeadingEdgeAngle(self._mer.HubContour)

	def GetInletC_u(self, r: Q_):
		return self.K_inlet / r

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

	