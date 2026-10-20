from dataclasses import dataclass
from typing import Optional
from pint import Quantity as Q_
from VelocityTriangle import VelocityTriangle
import numpy as np

g = Q_(9.81, 'm/s^2')

@dataclass
class InputGeometry:

	eta_hyd: float			# impeller hydraulic efficiency at BEP
	Z_blade: int			# number of impeller blades
	Beta2B: Q_[float]		# blade backsweep angle at outlet (from tangent)
	d_1: Q_[float]			# impeller inlet diameter
	d_2: Q_[float]			# impeller outlet diameter
	b_2: Q_[float]			# impeller outlet height
	thk2: Q_[float]			# blade thickness at exit
	d_hub: Q_[float]		# diameter of hub, strictly speaking inlet only

	# shroud_fullness: tuple[float, float] = (1, 0.3)
	# hub_fullness: tuple[float, float] = (0.5, 0.7)
	# # contour angles from radial at the outlet, positive leaning downstream (+z), see BuildMeridian
	# shroud_outlet_angle: Q_[float] = Q_(5, 'deg')
	# hub_outlet_angle: Q_[float] = Q_(0, 'deg')
	# hub_taper_ratio: float = 0.3

	Z_e_factor: float = 1

@dataclass
class DesignPoint:
	Q: Q_[float]			# volumetric flowrate
	N_shaft: Q_[float]		# shaft rotating speed
	H_rise: Q_[float]		# impeller developed head
	#rho: Q_[float]			# average density of the fluid

# Wiesner slip factor: https://manual.cfturbo.com/en/bl_te_wiesner.html
def GetWiesnerSlipRatio(beta2b: Q_[float], Z: int) -> float:
	return 1 - np.sqrt(np.sin(beta2b.to('rad').m)) / Z**0.7

def GetExitBladeBlockage(beta2b: Q_[float], Z: int, thick: Q_[float], b2: Q_[float]) -> Q_[float]:
	return thick * b2 * Z / np.sin(beta2b)

class Impeller:

	def __init__(self, geometry: InputGeometry, design_point: DesignPoint):

		assert geometry.d_1 > geometry.d_hub, 'Inlet diameter must be larger than hub diameter'
		
		self.DP: DesignPoint = design_point

		exit_blockage_area = GetExitBladeBlockage(
			beta2b=geometry.Beta2B,
			Z=geometry.Z_blade,
			thick=geometry.thk2,
			b2=geometry.b_2
		)
		self.Area2: Q_[float] = (np.pi * geometry.d_2 * geometry.b_2) - exit_blockage_area

		self.eta_hyd = geometry.eta_hyd
		self.Z_blade = geometry.Z_blade
		self.Beta2B = geometry.Beta2B
		self.d_1 = geometry.d_1
		self.d_2 = geometry.d_2
		self.b_2 = geometry.b_2
		self.thk2 = geometry.thk2
		self.d_hub = geometry.d_hub
		self.WiesnerSlip: float = GetWiesnerSlipRatio(geometry.Beta2B, geometry.Z_blade)
		self.Z_e: Q_[float] = (		# axial extension, Gulich, pg. 348
			geometry.Z_e_factor*((self.d_2 - self.d_1)*((self.n_q/74)**(1.07)))
			)

		u_2, c_m2, c_u2, w_2 = self.GetOutletVelocities()

		# AKA Flow factor: c / u
		self.FlowCoefficient: float = (c_m2 / u_2).to('dimensionless').m
		
		self.OutletVelocityTriangleAtDesign = VelocityTriangle(
			u=u_2, c_m=c_m2, c_u=c_u2, w=w_2, beta_b=self.Beta2B, station=2
			)

	@property
	def HeadCoefficient(self):
		"""head coefficient: g*H / u^2 | AKA Stage loading: Δh / u^2"""
		u_2 = self.d_2/2 * self.DP.N_shaft
		return (g * self.DP.H_rise / u_2**2).to('dimensionless').m

	@property
	def n_s(self) -> float:
		"""Imperial specific speed"""
		return (
			self.DP.N_shaft.to('rpm') * np.sqrt(self.DP.Q.to('gallon/min')) /
			self.DP.H_rise.to('ft')**0.75
			).m

	@property
	def n_q(self) -> float:
		"""Metric specific speed"""
		return (
			self.DP.N_shaft.to('rpm') * np.sqrt(self.DP.Q.to('m^3/s')) /
			self.DP.H_rise.to('m')**0.75
			).magnitude

	def GetOutletVelocities(self, Q_flow: Optional[Q_]=None, N_speed: Optional[Q_]=None) -> Q_:
		"""
			Returns u, c_m, c_u, and w
			at the outlet of the impeller, and at the specified flowrate and spin-speed
		"""

		if Q_flow is None:
			Q_flow = self.DP.Q
		if N_speed is None:
			N_speed = self.DP.N_shaft
		# Use velocity triangle to visualize (draw it by hand)
		u_2 = (self.d_2/2 * N_speed).to('m/s')
		c_m2 = (Q_flow / self.Area2).to('m/s')

		c_u2_ideal: Q_[float] = u_2 - c_m2 / np.tan(self.Beta2B.to('rad').m)
		# Wiesner slip = 1 - (c_u2_ideal - c_u2)/u_2
		# https://manual.cfturbo.com/en/bl_te_aungier.html
		c_u2 = c_u2_ideal - u_2*(1 - self.WiesnerSlip)
		w_2: Q_[float] = np.sqrt(c_m2**2 + (u_2 - c_u2)**2)

		return u_2, c_m2, c_u2, w_2

	def PrintDimensionlessParameters(self):

		print("\n--- Impeller dimensionless parameters ---")
		print(f"Specific speed (imperial)   = {self.n_s:.0f}")
		print(f"Specific speed (metric)     = {self.n_q:.1f}")

		print(f"Slip factor (σ)             = {self.WiesnerSlip:.4f}")
		print(f"Head coefficient (ψ)        = {self.HeadCoefficient:.4f}")
		print(f"Flow coefficient (ϕ)        = {self.FlowCoefficient:.4f}")
