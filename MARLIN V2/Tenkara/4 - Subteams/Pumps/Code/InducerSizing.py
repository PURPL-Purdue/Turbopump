import numpy as np
from dataclasses import dataclass, field, replace
from pint import Quantity as Q_
from Inducer import Inducer, InducerDesignPoint, InducerInletDesignConstraints, constraints, g

@dataclass
class ImpellerInlet:
	"""The inducer exit has to match these: the impeller inlet tip (shroud) and hub diameters"""
	d_1: Q_[float]			# impeller inlet diameter, the inducer tip diameter
	d_hub: Q_[float]		# impeller hub diameter at the inlet, the inducer exit hub diameter

@dataclass
class InletSizingLimits:
	sigma: tuple[float, float] = (0.03, 0.06)			# allowable design cavitation parameter (sheet range)
	hub_tip_ratio: tuple[float, float] = (0.2, 0.6)		# allowable inlet hub to tip ratio (range swept in InducerSizing)

@dataclass
class InletSizing:
	constraints: InducerInletDesignConstraints		# the solved constraints: Inducer built with them has the matching exit
	Inducer: Inducer
	psi_2: float									# head coefficient the hub taper has to reach
	taper_ratio: float								# sqrt(psi_1 / psi_2)
	Warnings: list[str] = field(default_factory=list)

def SizeInletForImpeller(design_point: InducerDesignPoint, base: InducerInletDesignConstraints, impeller: ImpellerInlet,
		limits: InletSizingLimits = InletSizingLimits()) -> InletSizing:
	"""
		Inducer inlet sizing from the impeller inlet. The exit of the inducer is the impeller inlet, so
			tip radius R_t = d_1 / 2 and exit hub radius R_h2 = d_hub / 2
		are given, and the inlet hub radius R_h1 (hub to tip ratio rho = R_h1 / R_t) and the cavitation parameter sigma
		follow from the two relations Inducer uses, each solved for its unknown:

		1. Hub taper (head coefficient). With U_t from the given R_t, psi_2 = g dH / U_t^2 and
				R_h2^2 = R_t^2 - (R_t^2 - R_h1^2) sqrt(psi_1 / psi_2)
			solved for the inlet hub radius
				R_h1^2 = R_t^2 - (R_t^2 - R_h2^2) / sqrt(psi_1 / psi_2)

		2. Optimal tip diameter (Gulich, sheet methodology) at that rho,
				D_t = 3.25 (Q / (N (1 - rho^2)))^(1/3) ((K + sigma) / sigma)^(1/6)	(m, m^3/s, rpm)
			solved for sigma, with x = (D_t / (3.25 (Q / (N (1 - rho^2)))^(1/3)))^6
				sigma = K / (x - 1)

		Nothing iterates: U_t is fixed by the impeller. The result is checked by building the Inducer with the solved
		rho and sigma, which must reproduce D_t and R_h2. A sigma or rho outside the limits means that impeller inlet
		cannot be matched at this design point, the other levers are psi_1, the head requirement and the impeller inlet
		itself. Those are listed in Warnings.
	"""
	r_t = impeller.d_1.to('m').m / 2
	r_h2 = impeller.d_hub.to('m').m / 2
	assert 0 < r_h2 < r_t, 'The impeller hub diameter must be smaller than its inlet diameter'

	q = design_point.Q.to('m^3/s').m
	n_rpm = design_point.N_shaft.to('rpm').m
	u_t = design_point.N_shaft.to('rad/s').m * r_t
	delta_H = design_point.H_rise.to('m').m

	psi_2 = g.to('m/s^2').m * delta_H / u_t**2		# H_rise: the head the inducer has to supply
	taper_ratio = float(np.sqrt(base.psi_1 / psi_2))
	if taper_ratio >= 1:
		raise ValueError(f"Baseline psi_1 = {base.psi_1} is not below the required psi_2 = {psi_2:.4f}: no hub taper is needed")

	# 1. hub taper
	r_h1_sq = r_t**2 - (r_t**2 - r_h2**2) / taper_ratio
	if r_h1_sq <= 0:
		raise ValueError(
			f"The required head rise needs a hub radius below zero at the inlet (R_h1^2 = {r_h1_sq:.3e} m^2): "
			f"R_h2 / R_t = {r_h2 / r_t:.3f} is too large for sqrt(psi_1 / psi_2) = {taper_ratio:.3f}"
		)
	rho = float(np.sqrt(r_h1_sq) / r_t)

	# 2. Gulich tip diameter
	x = (2 * r_t / (3.25 * (q / (n_rpm * (1 - rho**2))) ** (1 / 3))) ** 6
	if x <= 1:
		raise ValueError(f"The impeller inlet diameter {impeller.d_1.to('mm'):.2f} is smaller than the Gulich tip diameter "
			f"for any sigma at rho = {rho:.3f}")
	sigma = float(base.K / (x - 1))

	solved = replace(base, hub_tip_ratio=rho, sigma=sigma)
	inducer = Inducer(design_point, solved)

	# the inducer built from the solution has to match the impeller inlet
	assert np.isclose(inducer.D_tip.to('m').m, 2 * r_t, rtol=1e-9), 'Inducer tip diameter does not match the impeller inlet'
	assert np.isclose(inducer.R_hub_2.to('m').m, r_h2, rtol=1e-9), 'Inducer exit hub radius does not match the impeller inlet'

	warnings = []
	if not limits.sigma[0] <= sigma <= limits.sigma[1]:
		warnings.append(f"sigma = {sigma:.4f} is outside the limits {limits.sigma[0]} to {limits.sigma[1]}")
	if not limits.hub_tip_ratio[0] <= rho <= limits.hub_tip_ratio[1]:
		warnings.append(f"hub to tip ratio = {rho:.3f} is outside the limits {limits.hub_tip_ratio[0]} to {limits.hub_tip_ratio[1]}")

	return InletSizing(constraints=solved, Inducer=inducer, psi_2=psi_2, taper_ratio=taper_ratio, Warnings=warnings)

def PrintInletSizing(sizing: InletSizing, name: str = 'Inducer') -> None:
	inducer = sizing.Inducer
	print(f"\n===== {name} inducer inlet sized from the impeller inlet =====")
	print(f"Impeller inlet (given)     d_1 = {inducer.D_tip.to('mm'):.3f}, d_hub = {2 * inducer.R_hub_2.to('mm'):.3f}")
	print(f"Required head coefficient  psi_2 = {sizing.psi_2:.4f} (psi_1 = {sizing.constraints.psi_1}, sqrt ratio = {sizing.taper_ratio:.3f})")
	print(f"Inlet hub radius           R_h1 = {inducer.R_hub_1.to('mm'):.3f}")
	print(f"Hub to tip ratio           rho = {sizing.constraints.hub_tip_ratio:.4f}")
	print(f"Cavitation parameter       sigma = {sizing.constraints.sigma:.4f}")
	print(f"Hub taper angle            alpha = {inducer.alpha_taper.to('deg'):.2f}")
	print(f"Axial length               L_a = {inducer.L_ax.to('mm'):.3f}")
	print(f"Suction specific speed     N_ss = {inducer.N_ss:.1f}, NPSH_3 = {inducer.NPSH_3.to('m'):.1f}")
	for warning in sizing.Warnings:
		print(f"WARNING: {warning}")

if __name__ == '__main__':
	cases = {
		'LOx': dict(
			design_point=InducerDesignPoint(Q=Q_(1.69, 'L/s'), N_shaft=Q_(35000, 'rpm'), H_rise=Q_(40, 'm')),
			impeller_inlet=ImpellerInlet(d_1=Q_(1, 'in'), d_hub=Q_(0.6, 'in'))
		),
		'IPA': dict(
			design_point=InducerDesignPoint(Q=Q_(2.739, 'L/s'), N_shaft=Q_(35000, 'rpm'), H_rise=Q_(44.5, 'm')),
			impeller_inlet=ImpellerInlet(d_1=Q_(28.01, 'mm'), d_hub=Q_(19.39, 'mm'))
		),
	}
	# ------------------------------------------------------------------------------------------

	for name, case in cases.items():
		PrintInletSizing(SizeInletForImpeller(case['design_point'], constraints, case['impeller_inlet']), name)