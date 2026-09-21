import numpy as np
from pint import Quantity as Q_
from dataclasses import dataclass
from matplotlib import pyplot as plt
import EmpiricalRelations as emp
from VelocityTriangle import VelocityTriangle

g = Q_(9.81, 'm/s^2')

@dataclass
class InputGeometry:
	Z_blade: int			# number of impeller blades
	Beta2B: Q_[float]		# blade backsweep angle at outlet (from tangent)
	d_2: Q_[float]			# impeller outlet diameter
	b_2: Q_[float]			# impeller outlet height
	thk2: Q_[float]			# blade thickness at exit
	d_hub: Q_[float]		# diameter of hub, strictly speaking inlet only

@dataclass
class DesignPoint:
	Q: Q_[float]			# volumetric flowrate
	N_shaft: Q_[float]		# shaft rotating speed
	H: Q_[float]			# headrise, developed head
	n_hyd_BEP: float		# hydraulic efficiency at BEP

# Wiesner slip factor: https://manual.cfturbo.com/en/bl_te_wiesner.html
def GetWiesnerSlipRatio(beta2b: float | Q_[float], Z: int) -> float:
	return 1 - np.sqrt(np.sin(beta2b)) / Z**0.7

def GetBladeBlockage(beta2b: float, Z: int, thick: float, b2: float) -> float:
	return thick * b2 * Z / np.sin(beta2b)

class Impeller:

	def __init__(self, geometry: InputGeometry, design_point: DesignPoint):

		exit_area: Q_[float] = (np.pi * geometry.d_2 * geometry.b_2 - GetBladeBlockage(
			geometry.Beta2B, geometry.Z_blade, geometry.thk2, geometry.b_2
		))

		self.Z_blade = geometry.Z_blade
		self.Beta2B = geometry.Beta2B
		self.d_2 = geometry.d_2
		self.b_2 = geometry.b_2
		self.thk2 = geometry.thk2
		self.d_hub = geometry.d_hub
		self.WiesnerSlip = GetWiesnerSlipRatio(geometry.Beta2B, geometry.Z_blade)
		self.Area2 = exit_area

		self.DP: DesignPoint = design_point

		self.SpecificSpeed: float = (
			self.DP.N_shaft.to('rpm') * np.sqrt(self.DP.Q.to('gallon/min')) /
			self.DP.H.to('ft')**0.75
			).magnitude

		self.n_q: float = (
			self.DP.N_shaft.to('rpm') * np.sqrt(self.DP.Q.to('m^3/s')) /
			self.DP.H.to('m')**0.75
			).magnitude

		# exit tip velocity at design point
		self.U_2_design = (self.DP.N_shaft * self.d_2/2).to('m/s')
		# meridional exit flow velocity at design point
		self.C_m2_design = (self.DP.Q / self.Area2).to('m/s')

		_, enthalpy, _ = self.GetSpecificWork()
		# AKA Stage loading: Δh / u^2 = g*H / u^2
		self.HeadCoeff: float = (enthalpy/self.U_2_design**2).to('dimensionless').magnitude
		# AKA Flow factor: c / u
		self.FlowCoeff: float = (self.C_m2_design / self.U_2_design).to('dimensionless').magnitude
		
		# TODO: Optimize for inlet diameter d_1
		self.d_1 = Q_(1, 'in')

	@property
	def Area1(self) -> Q_[float]:
		return np.pi / 4 * (self.d_1**2 - self.d_hub**2) # TODO: blade blockage?
	@property
	def NPSH_i(self) -> Q_[float]:
		"""
			This method of calculating NPSH_i accounts for inlet velocity and losses.
			Therefore, when comparing with available NSPH, total pressure should be used for NPSH_a calculation
			NPSH_i = (gamma_c * 1/2 * c_m^2 + gamma_w * 1/2 * w^2) / g

			The physical intuition for the NPSH coefficients is that if the coefficients were equal to 1,
			then the drop in static pressure due to suction is merely the dynamic pressure due to
			acceleration of the flow at the inlet due to inlet area reduction (c_m, merdional velocity)
			and the acceleration of the flow over and around the blade leading edge (w, relative velocity)

		"""

		# Gulich, Centrifugal Pumps, chap. 6
		gamma_c = 1.2	# NPSH coefficient for main flow acceleration and losses at inlet
		gamma_w = 1.0	# NPSH coefficient for excess velocity due to flow around leading edge
		
		c_1m = self.DP.Q / self.Area1	# meridional velocity at inlet
		w_1 = np.sqrt(					# relative velocity of flow at leading edge
			(self.DP.N_shaft * self.d_1/2)**2 + c_1m**2
			)
		return (1/2 / g * (gamma_c * c_1m**2 + gamma_w * w_1**2)).to('m')
	
	def H_euler(self, Q):
		"""
		https://ntrs.nasa.gov/api/citations/19950013379/downloads/19950013379.pdf
		page 4 and 5
		
		Theoretical Euler head as a function of Q
		H_theoretical(Q) = U_2^2/g - [U_2 / (g*Area2*tan(beta2b))] * Q
		Or
		H_theoretical(Q) = U_2^2/g - [N_shaft*cot(beta2b)/(2*pi*b2*g)] * Q
		https://youtu.be/H7XfYO_-cEg?t=2705
		"""
		sigma = self.WiesnerSlip
		U_2_design = self.U_2_design
		N = self.DP.N_shaft
		beta2b = self.Beta2B
		b2 = self.b_2

		return sigma * U_2_design**2/g - N/np.tan(beta2b)/(2*np.pi*b2*g) * Q

	def GetInletVelocities(self, flow_Q: Q_[float]=None, speed_N: Q_[float]=None) -> VelocityTriangle:
		# no preswirl
		flow_Q = flow_Q if not flow_Q is None else self.DP.Q
		speed_N = speed_N if not speed_N is None else self.DP.N_shaft

		c_m1 = flow_Q / self.Area1
		c_u1 = Q_(0, 'm/s')
		u_1 = speed_N * self.d_1/2
		w_u1 = u_1 - c_u1
		c = np.sqrt(c_m1**2 + c_u1**2)
		w = np.sqrt(w_u1**2 + c_m1**2)

		return VelocityTriangle(
			u=u_1, c_m=c_m1, c_u=c_u1, w=w
		)

	def GetOutletVelocities(self, flow_Q: Q_[float]=None, speed_N: Q_[float]=None) -> tuple[VelocityTriangle, float]:
		flow_Q = flow_Q if not flow_Q is None else self.DP.Q
		speed_N = speed_N if not speed_N is None else self.DP.N_shaft
		
		f = emp.FlowSpeedRatio(
				flow_Q, speed_N,
				self.DP.Q, self.DP.N_shaft
				).to('dimensionless').magnitude
		
		U_2 = speed_N * self.d_2/2
		C_m2 = flow_Q / self.Area2
		slip = U_2 * (1 - self.WiesnerSlip)
		W_u2 = C_m2 / np.tan(self.Beta2B) + slip
		W_2 = np.sqrt(C_m2**2 + W_u2**2)
		C_u2 = (U_2 - W_u2)
		
		return VelocityTriangle(
			u=U_2, c_m=C_m2, c_u=C_u2, w=W_2
		), f

	def GetSpecificWork(self, speed_N: Q_[float]=None, flow_Q: Q_[float]=None) -> tuple[Q_[float], Q_[float], float]:
		if speed_N is None:
			speed_N = self.DP.N_shaft
		if flow_Q is None:
			flow_Q = self.DP.Q

		#vels: VelocityTriangle
		#f: float
		vels, f = self.GetOutletVelocities(flow_Q, speed_N)
		# Euler turbomachinery equation
		work_consumed = (vels.c_u * vels.u).to('kJ/kg')

		# empirical prediction, valid for 0 < f < 2
		n_hydraulic = emp.HydraulicEfficiency(f) * self.DP.n_hyd_BEP

		fluid_work = work_consumed * n_hydraulic

		return work_consumed, fluid_work, f
		
	def PlotPerformanceHQ(self, rpm_sweep: Q_[list[float]]=None) -> None:
		N_sweep = rpm_sweep if rpm_sweep is not None else Q_([self.DP.N_shaft])
		Q_sweep = Q_(np.linspace(0, 2.0 * self.DP.Q.to('L/s').magnitude, 50), 'L/s')

		H_predict = []
		Q_predict = []

		for speed_N in N_sweep:

			head_trace = []
			flow_trace = []

			for flow_Q in Q_sweep:
				
				_, fluid_work, f = self.GetSpecificWork(speed_N=speed_N, flow_Q=flow_Q)
				if f > 2: continue

				headrise = fluid_work / g
		
				head_trace.append(headrise.to('m').magnitude)
				flow_trace.append(flow_Q.to('L/s').magnitude)
		
			H_predict.append(head_trace)
			Q_predict.append(flow_trace)

		plt.figure()
		# Plot theoretical Euler (linear)
		H_theoretical = self.H_euler(Q_sweep)
		plt.plot(Q_sweep.to('L/s').magnitude, H_theoretical.to('m').magnitude,
				'--', color='gray', label='Theoretical Euler')

		# Plot RPM swept traces
		for i in range(len(H_predict)):

			plt.plot(Q_predict[i], H_predict[i],
					'-', color='red', label='Empirical prediction' if i == 0 else None)
			
			plt.annotate(f"{N_sweep[i].magnitude:.0f} RPM",
						xy=(Q_predict[i][-1], H_predict[i][-1]),
						xytext=(4, 0), textcoords='offset points',
						fontsize=8, va='center')

		# Plot design point head and flowrate
		plt.plot(self.DP.Q.to('L/s').magnitude, self.DP.H.to('m').magnitude,
				'o', color='black', markersize=4, label='Design point',
				zorder=5)

		plt.xlabel('Volumetric flow rate [L/s]')
		plt.ylabel('Head [m]')
		plt.title('H-Q Curve')
		plt.ylim(ymin=0)
		plt.legend(loc='upper right')
		plt.grid(True, alpha=0.3)