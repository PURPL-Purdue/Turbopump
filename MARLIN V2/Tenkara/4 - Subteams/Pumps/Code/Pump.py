import Impeller
import Inducer
from pint import Quantity as Q_
from dataclasses import dataclass

g = Q_(9.81, 'm/s^2')

@dataclass
class PumpDesignPoint:
	feed_pressure: Q_[float]			# feed pressure (tank pressure)
	discharge_pressure: Q_[float]		# desired outlet pressure
	inlet_vapor_pressure: Q_[float]		# vapor pressure of the working fluid at the inlet
	N_shaft: Q_[float]					# shaft spin speed
	mass_flow: Q_[float]				# mass flow rate (strictly speaking inlet, should be almost equal to outlet)
	rho: Q_[float]						# average fluid density

class Pump:

	def __init__(self, design_point: PumpDesignPoint, impeller_geometry: Impeller.InputGeometry):
		self.DP = design_point

		self.DesignHeadRise = (self.DP.discharge_pressure - self.DP.feed_pressure) / self.DP.rho / g
		self.DesignFlowRate = self.DP.mass_flow / self.DP.rho

		self.Impeller = Impeller.Impeller(
			geometry=impeller_geometry,
			design_point=Impeller.DesignPoint(
				Q=design_point.mass_flow/design_point.rho,
				N_shaft=design_point.N_shaft,
				H_rise=self.DesignHeadRise # initial
			)
		)

		self.NPSH_a = ((design_point.feed_pressure - design_point.inlet_vapor_pressure) / design_point.rho / g).to('m')

	def AddInducer(self, inducer_head: Q_[float]) -> Inducer.Inducer:
		self.Inducer = Inducer.Inducer(
			design_point=Inducer.InducerDesignPoint(
				Q=self.DesignFlowRate, N_shaft=self.DP.N_shaft, H_rise=inducer_head
				)
		)
		