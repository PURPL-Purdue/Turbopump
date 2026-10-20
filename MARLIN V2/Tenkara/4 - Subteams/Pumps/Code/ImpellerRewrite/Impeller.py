from dataclasses import dataclass
from pint import Quantity as Q_
import numpy as np

g = Q_(9.81, 'm/s^2')

@dataclass
class InputGeometry:
	Z_blade: int			# number of impeller blades
	Beta2B: Q_[float]		# blade backsweep angle at outlet (from tangent)
	d_1: Q_[float]			# impeller inlet diameter
	d_2: Q_[float]			# impeller outlet diameter
	b_2: Q_[float]			# impeller outlet height
	thk2: Q_[float]			# blade thickness at exit
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
			geometry.Beta2B, geometry.Z_blade, geometry.thk2, geometry.b_2
		))

		self.Z_blade = geometry.Z_blade
		self.Beta2B = geometry.Beta2B
		self.d_1 = geometry.d_1
		self.d_2 = geometry.d_2
		self.b_2 = geometry.b_2
		self.thk2 = geometry.thk2
		self.d_hub = geometry.d_hub
		self.WiesnerSlip = GetWiesnerSlipRatio(geometry.Beta2B, geometry.Z_blade)
		self.Area2 = exit_area
		self.Z_e: Q_[float] = (		# axial extension, Gulich, pg. 348
			geometry.Z_e_factor*((self.d_2 - self.d_1)*((self.n_q/74)**(1.07)))
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