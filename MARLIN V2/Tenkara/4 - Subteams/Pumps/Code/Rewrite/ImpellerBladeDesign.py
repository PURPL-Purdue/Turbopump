import numpy as np
from pint import Quantity as Q_
from Impeller import Impeller
from Inducer import Inducer
from InducerRadialEquilibrium import InducerRadialEquilibrium
from MeridionalDesign import MeridionalDesign
from MeridionalContour import MeridionalCurve
from PotentialFlow import Streamline
from typing import Optional

class ImpellerBladeDesign:

	def __init__(self, meridional_design: MeridionalDesign, inducer_radial_equilibrium: Optional[InducerRadialEquilibrium]=None):

		self._mer = meridional_design
		self._rad_eq = inducer_radial_equilibrium
		self.K_inlet = inducer_radial_equilibrium.GetOutletAngularMomentum() if inducer_radial_equilibrium is not None else Q_(0, 'm^2/s')

		# blade angle, tilt of the blade surface to the stream surface and degree of blockage at the hub end of the leading
		# edge, with the leading edge angle taken to the exact hub contour rather than to the polyline hub streamline
		self._beta_1_hub, self._epsilon_1_hub, self._rho_1_hub = self._SolveInletBlockage(
			self._mer.PotentialFlow.GetStreamline(0), self._GetLeadingEdgeAngle(self._mer.HubContour)
		)

	def GetInletC_u(self, r: Q_):
		return self.K_inlet / r

	def _SolveInletBlockage(self, streamline: Streamline, lambda_1: Q_[float]) -> tuple[Q_[float], Q_[float], float]:
		"""
			Blade angle beta_1', tilt epsilon_1 of the blade surface to the stream surface and degree of blockage rho_1
			where the leading edge crosses streamline, with lambda_1 the angle between the leading edge and the streamline
			in the meridional plane (Milos eqs. 31-33):
				tan(beta_1') = v_m1 / (rho_1 (u_1 - c_u1))			eq. 33, with the preswirl c_u1 of eq. 42
				tan(epsilon_1) = tan(lambda_1) / cos(beta_1')		own derivation, leading edge in a meridional plane
				rho_1 = 1 - z s_1r / (2 pi r_1 sin(beta_1') sin(epsilon_1))		eq. 31
			s_1r is the thickness normal to the blade surface, taken as blade_thick: constant from inlet to outlet (Milos p. 44).
			Milos iterates on rho_1 from 0.8, which oscillates without converging when the blockage is heavy. Instead
			beta_1' is found by bisection: as beta_1' grows, the rho_1 the velocity triangle needs falls and the rho_1 the
			blades give (eq. 31) rises, so the two are equal at exactly one beta_1' between 0 and 90 deg.
		"""
		hits = streamline.FindIntersect(self._mer.LeadingEdge)
		assert len(hits) == 1, 'The leading edge must cross the streamline exactly once'
		_, r_1, s, _ = hits[0]

		impeller = self._mer.Pump.Impeller
		v_m1 = streamline.Interpolate(streamline.v_m, s).to('m/s').m
		w_u1 = (impeller.DP.N_shaft.to('rad/s').m * r_1.to('m').m) - self.GetInletC_u(r_1).to('m/s').m
		blade_pitch_ratio = impeller.Z_blade * impeller.blade_thick.to('m').m / (2 * np.pi * r_1.to('m').m)	# z s_1r / (2 pi r_1)
		tan_lambda = np.tan(lambda_1.to('rad').m)
		assert w_u1 > 0, 'The preswirl exceeds the blade speed at the leading edge'
		assert blade_pitch_ratio < 1, 'The blades fill the inlet circumference'

		def Epsilon(beta: float) -> float:
			return np.arctan2(tan_lambda, np.cos(beta))

		def Mismatch(beta: float) -> float:
			"""rho_1 the velocity triangle needs minus rho_1 the blades give, decreasing in beta"""
			return v_m1 / (w_u1 * np.tan(beta)) - (1 - blade_pitch_ratio / (np.sin(beta) * np.sin(Epsilon(beta))))

		beta_a, beta_b = 0, np.pi / 2
		while beta_b - beta_a > 1e-12:
			beta = (beta_a + beta_b) / 2
			if Mismatch(beta) > 0:
				beta_a = beta
			else:
				beta_b = beta
		beta_1 = (beta_a + beta_b) / 2
		rho_1 = v_m1 / (w_u1 * np.tan(beta_1))

		return Q_(beta_1, 'rad').to('deg'), Q_(Epsilon(beta_1), 'rad').to('deg'), rho_1

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

	