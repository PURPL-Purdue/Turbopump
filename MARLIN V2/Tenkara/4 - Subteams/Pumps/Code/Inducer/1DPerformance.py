import numpy as np
from dataclasses import dataclass
from pint import Quantity as Q_

g = Q_(9.81, 'm/s^2')

# ---------------------------------------------------------------------------------------------------------------------
# 1D performance (no blade geometry) of a sized inducer. Takes the geometry from InducerInletSizing.InletSizing (anything with the same
# attributes works): design_point (Q, N_shaft, H_rise), D_tip, D_hub_1, D_hub_2, lambda_w, constants.lambda_c
# ---------------------------------------------------------------------------------------------------------------------

@dataclass
class Inducer1D:
	"""1D performance parameters of the sized inducer (no velocity triangles)"""
	# cavitation
	phi: float					# flow coefficient from lambda_w, sqrt(lambda_w / (2 (lambda_c + lambda_w)))
	N_ss: float					# suction specific speed (rpm, m^3/s, m)
	NPSH_3: Q_[float]			# 3% head drop NPSH, (N sqrt(Q) / N_ss)^(4/3)
	NPSH_a: Q_[float] | None	# available, if given
	NPSH_margin: Q_[float] | None	# NPSH_a - NPSH_3
	# head
	H_rise: Q_[float]
	U_tip: Q_[float]
	psi_2: float				# g H / U_tip^2
	# areas
	A_inlet: Q_[float]
	A_outlet: Q_[float]
	area_ratio: float			# A_outlet / A_inlet

def Calculate1DPerformance(sizing, NPSH_a: Q_[float] | None = None) -> Inducer1D:
	"""1D performance from the sizing result: cavitation (phi, N_ss, NPSH_3), head coefficient and annulus areas"""
	dp = sizing.design_point
	q, n_rpm = dp.Q.to('m^3/s').m, dp.N_shaft.to('rpm').m
	lambda_w, lambda_c = sizing.lambda_w, sizing.constants.lambda_c
	rho = sizing.hub_tip_ratio

	phi = float(np.sqrt(lambda_w / (2 * (lambda_c + lambda_w))))
	N_ss = float((98 / (lambda_c + lambda_w) ** 0.25) * np.sqrt((1 - rho**2) / lambda_w))
	NPSH_3 = Q_((n_rpm * np.sqrt(q) / N_ss) ** (4 / 3), 'm')

	U_tip = (dp.N_shaft.to('rad/s') * sizing.D_tip / 2).to('m/s')
	A_in = (np.pi * (sizing.D_tip**2 - sizing.D_hub_1**2) / 4).to('m^2')
	A_out = (np.pi * (sizing.D_tip**2 - sizing.D_hub_2**2) / 4).to('m^2')

	return Inducer1D(
		phi=phi, N_ss=N_ss, NPSH_3=NPSH_3, NPSH_a=NPSH_a,
		NPSH_margin=(NPSH_a.to('m') - NPSH_3) if NPSH_a is not None else None,
		H_rise=dp.H_rise.to('m'), U_tip=U_tip, psi_2=float((g * dp.H_rise / U_tip**2).to('').m),
		A_inlet=A_in, A_outlet=A_out, area_ratio=float((A_out / A_in).to('').m)
	)

def Print1DPerformance(p: Inducer1D, name: str = 'Inducer') -> None:
	rows = lambda title, items: (print(f"\n{title}"), [print(f"{k:<36} = {v}") for k, v in items])
	print(f"\n===== {name} inducer 1D performance =====")
	rows('Cavitation', [
		('phi (from lambda_w)', f"{p.phi:.4f}"),
		('N_ss', f"{p.N_ss:.1f}"),
		('NPSH_3', f"{p.NPSH_3.to('m'):.2f}"),
		*([('NPSH_a', f"{p.NPSH_a.to('m'):.2f}"), ('NPSH margin (NPSH_a - NPSH_3)', f"{p.NPSH_margin.to('m'):.2f}")] if p.NPSH_a is not None else []),
	])
	rows('Head', [('H_rise', f"{p.H_rise:.2f}"), ('U_tip', f"{p.U_tip:.2f}"), ('psi_2 (g H / U_tip^2)', f"{p.psi_2:.4f}")])
	rows('Areas', [('A_inlet', f"{p.A_inlet.to('mm^2'):.2f}"), ('A_outlet', f"{p.A_outlet.to('mm^2'):.2f}"), ('A_outlet / A_inlet', f"{p.area_ratio:.4f}")])

if __name__ == '__main__':
	import sys
	from pathlib import Path
	sys.path.append(str(Path(__file__).resolve().parent.parent))		# Code folder, where Main.py sits
	import Main
	# Runs the sizing for the cases in Sizing.py and prints the 1D performance of each result
	from Sizing import (constants, ImpellerInlet, InducerDesignPoint, InletTolerance, HeadRequirementModel,
		OptimizeImpellerInlet)

	cases = {
		'LOx': dict(Q=Q_(1.771, 'L/s'), NPSH_a=Q_(52.3, 'm'), baseline=ImpellerInlet(Q_(1, 'in'), Q_(0.6, 'in'))),
		'IPA': dict(Q=Q_(2.735, 'L/s'), NPSH_a=Q_(88.8, 'm'), baseline=ImpellerInlet(Q_(1.2, 'in'), Q_(0.7, 'in'))),
	}
	for name, c in cases.items():
		dp = InducerDesignPoint(Q=c['Q'], N_shaft=Q_(35000, 'rpm'), H_rise=Q_(0, 'm'))
		result = OptimizeImpellerInlet(dp, constants, c['baseline'], HeadRequirementModel(NPSH_a=c['NPSH_a']), InletTolerance(0.20, 0.20))
		Print1DPerformance(Calculate1DPerformance(result.Sizing, NPSH_a=c['NPSH_a']), name)