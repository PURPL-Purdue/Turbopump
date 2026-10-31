import numpy as np
import plotly.graph_objects as go
from dataclasses import dataclass, field, replace
from pint import Quantity as Q_

g = Q_(9.81, 'm/s^2')

# ---------------------------------------------------------------------------------------------------------------------
# Inducer inlet sizing from the impeller inlet.
# ---------------------------------------------------------------------------------------------------------------------

@dataclass
class InducerDesignPoint:
	Q: Q_[float]			# volumetric flow rate
	N_shaft: Q_[float]		# shaft speed
	H_rise: Q_[float]		# head the inducer has to supply

@dataclass
class SizingConstants:
	lambda_w: float = 0.03					# nominal design cavitation parameter (sheet uses 0.03 to 0.06)
	lambda_c: float = 1.1						# empirical constant in phi = sqrt(sigma / (2 (K + sigma)))
	psi_1: float = 0.08					# baseline head coefficient of the untapered inducer (NASA SP-8052)
	axial_length_ratio: float = 0.6		# L_a / D_t

constants = SizingConstants()

@dataclass
class ImpellerInlet:
	"""The inducer exit has to match these: the impeller inlet tip (shroud) and hub diameters"""
	d_1: Q_[float]			# impeller inlet diameter, the inducer tip diameter
	d_hub: Q_[float]		# impeller hub diameter at the inlet, the inducer exit hub diameter

@dataclass
class InletSizingLimits:
	lambda_w: tuple[float, float] = (0.03, 0.06)			# allowable design cavitation parameter (sheet range)
	hub_tip_ratio: tuple[float, float] = (0.2, 0.6)		# allowable inlet hub to tip ratio (range swept in InducerSizing)

def GulichTipDiameter(Q_flow: Q_[float], N_shaft: Q_[float], rho: float, lambda_w: float, lambda_c: float) -> Q_[float]:
	"""Optimal tip diameter of the sheet / Gulich relation, D_t = 3.25 (Q / (N (1 - rho^2)))^(1/3) ((K + sigma) / sigma)^(1/6)"""
	q, n_rpm = Q_flow.to('m^3/s').m, N_shaft.to('rpm').m		# the 3.25 is dimensional: m^3/s, rpm, m
	return Q_(3.25 * (q / (n_rpm * (1 - rho**2))) ** (1 / 3) * ((lambda_c + lambda_w) / lambda_w) ** (1 / 6), 'm')

@dataclass
class InletSizing:
	design_point: InducerDesignPoint
	constants: SizingConstants						# the solved constants: rho is stored below, sigma and axial_length_ratio are the solved ones
	D_tip: Q_[float]								# inducer tip diameter, equal to the impeller inlet d_1
	D_hub_1: Q_[float]								# inducer inlet hub diameter
	D_hub_2: Q_[float]								# inducer exit hub diameter, equal to the impeller inlet d_hub
	hub_tip_ratio: float							# rho = D_hub_1 / D_tip
	lambda_w: float									# cavitation parameter that makes D_tip the Gulich optimum at rho
	L_ax: Q_[float]									# axial length
	alpha_taper: Q_[float]							# hub taper angle
	psi_2: float									# head coefficient the hub taper has to reach
	taper_ratio: float								# sqrt(psi_1 / psi_2)
	D_opt_nominal: Q_[float]						# Gulich optimal tip diameter at the solved rho and the nominal sigma (constants.sigma)
	D_opt_band: tuple[Q_[float], Q_[float]]			# the same over the sigma limits, (at the largest sigma, at the smallest)
	tip_error: float								# D_tip / D_opt_nominal - 1: how far the impeller inlet is from the optimum
	Warnings: list[str] = field(default_factory=list)

	@property
	def R_tip(self) -> Q_[float]: return self.D_tip / 2
	@property
	def R_hub_1(self) -> Q_[float]: return self.D_hub_1 / 2
	@property
	def R_hub_2(self) -> Q_[float]: return self.D_hub_2 / 2

def SizeInletForImpeller(design_point: InducerDesignPoint, base: SizingConstants, impeller: ImpellerInlet,
		limits: InletSizingLimits = InletSizingLimits()) -> InletSizing:
	
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
	r_h1 = float(np.sqrt(r_h1_sq))
	rho = r_h1 / r_t

	# 2. Gulich tip diameter
	x = (2 * r_t / (3.25 * (q / (n_rpm * (1 - rho**2))) ** (1 / 3))) ** 6
	if x <= 1:
		raise ValueError(f"The impeller inlet diameter {impeller.d_1.to('mm'):.2f} is smaller than the Gulich tip diameter "
			f"for any lambda_w at rho = {rho:.3f}")
	lambda_w = float(base.lambda_c / (x - 1))

	# the Gulich relation at the solved rho and sigma has to give back the impeller inlet
	assert np.isclose(GulichTipDiameter(design_point.Q, design_point.N_shaft, rho, lambda_w, base.lambda_c).m, 2 * r_t, rtol=1e-9), \
		'Gulich tip diameter does not match the impeller inlet'

	# Tip diameter check. The optimum depends on rho and sigma, so a rho inside its limits says nothing about it. The
	# optimal diameter at the solved rho is compared with the impeller inlet at the nominal sigma and over the sigma limits
	# (a smaller sigma gives a larger optimal tip diameter). The impeller inlet has to fall inside that band.
	D_nominal = GulichTipDiameter(design_point.Q, design_point.N_shaft, rho, base.lambda_w, base.lambda_c)
	D_band = (
		GulichTipDiameter(design_point.Q, design_point.N_shaft, rho, limits.lambda_w[1], base.lambda_c),
		GulichTipDiameter(design_point.Q, design_point.N_shaft, rho, limits.lambda_w[0], base.lambda_c)
	)
	tip_error = float(2 * r_t / D_nominal.m - 1)

	l_ax = base.axial_length_ratio * 2 * r_t
	alpha = float(np.degrees(np.arctan((r_h2 - r_h1) / l_ax)))

	warnings = []
	if not D_band[0].m * (1 - 1e-9) <= 2 * r_t <= D_band[1].m * (1 + 1e-9):
		warnings.append(
			f"impeller inlet d_1 = {impeller.d_1.to('mm'):.2f} is outside the Gulich optimal tip diameter band "
			f"{D_band[0].to('mm'):.2f} to {D_band[1].to('mm'):.2f} at rho = {rho:.3f} "
			f"(lambda_w limits {limits.lambda_w[0]} to {limits.lambda_w[1]}, it needs lambda_w = {lambda_w:.4f})"
		)
	if not limits.hub_tip_ratio[0] <= rho <= limits.hub_tip_ratio[1]:
		warnings.append(f"hub to tip ratio = {rho:.3f} is outside the limits {limits.hub_tip_ratio[0]} to {limits.hub_tip_ratio[1]}")

	return InletSizing(
		design_point=design_point, constants=replace(base, lambda_w=lambda_w),
		D_tip=Q_(2 * r_t, 'm'), D_hub_1=Q_(2 * r_h1, 'm'), D_hub_2=Q_(2 * r_h2, 'm'), hub_tip_ratio=rho, lambda_w=lambda_w,
		L_ax=Q_(l_ax, 'm'), alpha_taper=Q_(alpha, 'deg'), psi_2=psi_2, taper_ratio=taper_ratio,
		D_opt_nominal=D_nominal, D_opt_band=D_band, tip_error=tip_error, Warnings=warnings
	)

def PrintInletSizing(sizing: InletSizing, name: str = 'Inducer') -> None:
	print(f"\n===== {name} inducer inlet sized from the impeller inlet =====")
	print(f"Impeller inlet (given)     d_1 = {sizing.D_tip.to('mm'):.3f}, d_hub = {sizing.D_hub_2.to('mm'):.3f}")
	print(f"Required head coefficient  psi_2 = {sizing.psi_2:.4f} (psi_1 = {sizing.constants.psi_1}, sqrt ratio = {sizing.taper_ratio:.3f})")
	print(f"Inlet hub radius           R_h1 = {sizing.R_hub_1.to('mm'):.3f}")
	print(f"Hub to tip ratio           rho = {sizing.hub_tip_ratio:.4f}")
	print(f"Cavitation parameter       lambda_w = {sizing.lambda_w:.4f} needed to make d_1 the optimum")
	print(f"Optimal tip diameter       {sizing.D_opt_nominal.to('mm'):.3f} at nominal lambda_w, d_1 is {100 * sizing.tip_error:+.1f}% from it")
	print(f"Optimal tip band           {sizing.D_opt_band[0].to('mm'):.3f} to {sizing.D_opt_band[1].to('mm'):.3f} over the lambda_w limits")
	print(f"Hub taper angle            alpha = {sizing.alpha_taper:.2f}")
	print(f"Axial length               L_a = {sizing.L_ax.to('mm'):.3f}, L_a / D_t = {sizing.constants.axial_length_ratio:.3f}")
	for warning in sizing.Warnings:
		print(f"WARNING: {warning}")

def MeridionalTraces(sizing: InletSizing, l: str = 'mm', color: str = 'blue') -> list[go.Scatter]:
	"""
		Meridional section of the sized inducer as plotly traces in length unit l, z = 0 at the inlet: the straight shroud at
		the tip radius, the hub tapering from R_h1 to R_h2 over L_a, the inlet and outlet stations and the axis.
	"""
	scale = Q_(1, 'm').to(l).magnitude
	r_t, r_1, r_2, length = (v.to('m').m * scale for v in (sizing.R_tip, sizing.R_hub_1, sizing.R_hub_2, sizing.L_ax))
	line = lambda x, y, **kw: go.Scatter(x=x, y=y, mode='lines', showlegend=False, **kw)
	return [
		line([0, length], [r_t, r_t], name='Shroud', line=dict(color=color)),
		line([0, length], [r_1, r_2], name='Hub', line=dict(color=color)),
		line([0, 0], [r_1, r_t], name='Inlet', line=dict(dash='dot', color='grey')),
		line([length, length], [r_2, r_t], name='Outlet', line=dict(dash='dot', color='grey')),
		line([0, length], [0, 0], name='Axis', line=dict(dash='dashdot', color='grey')),
	]

def PlotContour(sizing: InletSizing, l: str = 'mm', color: str = 'blue', name: str = 'Inducer') -> go.Figure:
	"""Plain meridional section of the sized inducer (MeridionalTraces)"""
	fig = go.Figure(MeridionalTraces(sizing, l, color))
	fig.update_xaxes(title_text=f"Meridional position ({l})")
	fig.update_yaxes(title_text=f"Radial position ({l})", scaleanchor='x', scaleratio=1)
	fig.update_layout(title=f"{name} Inducer Meridional Section", template="plotly_white", showlegend=False)
	return fig

# ---------------------------------------------------------------------------------------------------------------------
# Iterative search: the largest impeller inlet that satisfies every condition
# ---------------------------------------------------------------------------------------------------------------------

@dataclass
class HeadRequirementModel:
	"""
		Head the inducer has to supply. It depends on the impeller inlet, so it moves when the search changes d_1 and d_hub.
		Same as Pump1D.GetInducerHeadRequirement: NPSH_r of the impeller without an inducer (GuessInducerlessNPSH_r), minus
		NPSH_a, rounded up to a multiple of step.
	"""
	NPSH_a: Q_[float]
	lambda_c: float = 1.1
	lambda_w: float = 1.0
	blockage_factor: float = 0.7
	step: Q_[float] = Q_(5, 'm')

	def Head(self, Q_flow: Q_[float], N_shaft: Q_[float], d_1, d_hub) -> np.ndarray:
		"""Head requirement in m for impeller inlet diameters d_1 and d_hub (m, floats or arrays)"""
		q, omega = Q_flow.to('m^3/s').m, N_shaft.to('rad/s').m
		with np.errstate(divide='ignore'):		# d_hub = d_1 has no inlet area, the point is then never valid
			area_1 = (d_1**2 - d_hub**2) / 4 * np.pi * self.blockage_factor
			c_m1 = q / area_1
		u_1 = omega * d_1 / 2
		w_1_sq = c_m1**2 + u_1**2
		npsh_r = (self.lambda_c * c_m1**2 + self.lambda_w * w_1_sq) / (2 * g.to('m/s^2').m)
		step = self.step.to('m').m
		return step * np.round((npsh_r - self.NPSH_a.to('m').m) / step + 0.5)

@dataclass
class InletTolerance:
	"""How far the impeller inlet may move from its baseline, as a fraction of each diameter, either way"""
	d_1: float = 0.10
	d_hub: float = 0.10

@dataclass
class OptimizedInlet:
	Impeller: ImpellerInlet					# the largest impeller inlet that satisfies every condition (or the closest if not Feasible)
	H_rise: Q_[float]						# inducer head requirement at that impeller inlet
	Sizing: InletSizing
	Feasible: bool
	Binding: list[str]						# the limits the result sits on
	Violations: list[str] = field(default_factory=list)		# only if not Feasible: what is still violated

def _Evaluate(d_1, d_hub, H_rise, design_point: InducerDesignPoint, base: SizingConstants,
		limits: InletSizingLimits) -> dict[str, np.ndarray]:
	"""
		The conditions of SizeInletForImpeller on arrays of impeller inlet diameters d_1, d_hub (m) and the head H_rise (m):
		rho, sigma, and whether each condition holds. Invalid points (no taper, hub radius below zero, no Gulich
		solution) get nan.
	"""
	q, n_rpm = design_point.Q.to('m^3/s').m, design_point.N_shaft.to('rpm').m
	r_t, r_h2 = d_1 / 2, d_hub / 2
	u_t = design_point.N_shaft.to('rad/s').m * r_t
	psi_2 = g.to('m/s^2').m * H_rise / u_t**2

	with np.errstate(all='ignore'):
		s = np.sqrt(base.psi_1 / psi_2)
		r_h1_sq = r_t**2 - (r_t**2 - r_h2**2) / s
		valid = (s < 1) & (r_h1_sq > 0) & (r_h2 < r_t)
		rho = np.where(valid, np.sqrt(r_h1_sq) / r_t, np.nan)
		x = (d_1 / (3.25 * (q / (n_rpm * (1 - rho**2))) ** (1 / 3))) ** 6
		valid &= x > 1
		lambda_w = np.where(valid, base.lambda_c / (x - 1), np.nan)

	return dict(
		rho=rho, lambda_w=lambda_w, valid=valid,
		lambda_w_ok=(lambda_w >= limits.lambda_w[0]) & (lambda_w <= limits.lambda_w[1]),
		rho_ok=(rho >= limits.hub_tip_ratio[0]) & (rho <= limits.hub_tip_ratio[1])
	)

def OptimizeImpellerInlet(design_point: InducerDesignPoint, base: SizingConstants, baseline: ImpellerInlet,
		head: HeadRequirementModel | Q_[float], tolerance: InletTolerance = InletTolerance(),
		limits: InletSizingLimits = InletSizingLimits(), objective: str = 'd_1',
		n_d_1: int = 801, n_d_hub: int = 201) -> OptimizedInlet:
	
	assert objective in ('d_1', 'area'), "objective is 'd_1' or 'area'"
	d1_0, dh_0 = baseline.d_1.to('m').m, baseline.d_hub.to('m').m
	d_1, d_hub = np.meshgrid(
		np.linspace(d1_0 * (1 - tolerance.d_1), d1_0 * (1 + tolerance.d_1), n_d_1),
		np.linspace(dh_0 * (1 - tolerance.d_hub), dh_0 * (1 + tolerance.d_hub), n_d_hub),
		indexing='ij'
	)

	if isinstance(head, HeadRequirementModel):
		H_rise = head.Head(design_point.Q, design_point.N_shaft, d_1, d_hub)
	else:
		H_rise = np.full(d_1.shape, head.to('m').m)

	e = _Evaluate(d_1, d_hub, H_rise, design_point, base, limits)
	ok = e['valid'] & e['lambda_w_ok'] & e['rho_ok']

	if ok.any():
		if objective == 'd_1':
			d1_best = d_1[ok].max()
			ok_best = ok & (d_1 >= d1_best * (1 - 1e-12))
			k = np.argmin(np.where(ok_best, np.abs(d_hub / dh_0 - 1), np.inf))		# fewest changes to the hub
		else:
			k = np.argmax(np.where(ok, d_1**2 - d_hub**2, -np.inf))
	else:
		# closest: the valid points with the smallest total violation of the sigma and rho limits, relative to the limits
		lo_s, hi_s = limits.lambda_w
		lo_r, hi_r = limits.hub_tip_ratio
		viol = (np.maximum(lo_s - e['sigma'], 0) / lo_s + np.maximum(e['sigma'] - hi_s, 0) / hi_s
			+ np.maximum(lo_r - e['rho'], 0) / lo_r + np.maximum(e['rho'] - hi_r, 0) / hi_r)
		viol = np.where(e['valid'], viol, np.inf)
		if not np.isfinite(viol).any():
			raise ValueError("No impeller inlet within the tolerance has a valid hub taper: the head requirement is too "
				"large (or too small, psi_2 below psi_1) for every d_1 and d_hub tried. Widen the tolerance or change psi_1.")
		k = np.argmin(viol)

	k = np.unravel_index(k, d_1.shape)
	impeller = ImpellerInlet(d_1=Q_(d_1[k], 'm').to(baseline.d_1.units), d_hub=Q_(d_hub[k], 'm').to(baseline.d_hub.units))
	H = Q_(H_rise[k], 'm')
	feasible = bool(ok.any())

	sized = replace(design_point, H_rise=H)
	sizing = SizeInletForImpeller(sized, base, impeller, limits)

	# the limits the result sits on, within one grid step
	binding, violations = [], list(sizing.Warnings) if not feasible else []
	d1_step = 2 * d1_0 * tolerance.d_1 / (n_d_1 - 1)
	dh_step = 2 * dh_0 * tolerance.d_hub / (n_d_hub - 1)
	lambda_w_r, rho_r = sizing.lambda_w, sizing.hub_tip_ratio
	if abs(lambda_w_r - limits.lambda_w[0]) < 0.02 * limits.lambda_w[0]: binding.append(f"lambda_w at its lower limit {limits.lambda_w[0]} (largest optimal tip diameter)")
	if abs(lambda_w_r - limits.lambda_w[1]) < 0.02 * limits.lambda_w[1]: binding.append(f"lambda_w at its upper limit {limits.lambda_w[1]}")
	if abs(rho_r - limits.hub_tip_ratio[0]) < 0.02 * limits.hub_tip_ratio[0]: binding.append(f"hub to tip ratio at its lower limit {limits.hub_tip_ratio[0]}")
	if abs(rho_r - limits.hub_tip_ratio[1]) < 0.02 * limits.hub_tip_ratio[1]: binding.append(f"hub to tip ratio at its upper limit {limits.hub_tip_ratio[1]}")
	if abs(d_1[k] - d1_0 * (1 + tolerance.d_1)) < d1_step: binding.append("d_1 at the upper edge of its tolerance")
	if abs(d_1[k] - d1_0 * (1 - tolerance.d_1)) < d1_step: binding.append("d_1 at the lower edge of its tolerance")
	if abs(d_hub[k] - dh_0 * (1 + tolerance.d_hub)) < dh_step: binding.append("d_hub at the upper edge of its tolerance")
	if abs(d_hub[k] - dh_0 * (1 - tolerance.d_hub)) < dh_step: binding.append("d_hub at the lower edge of its tolerance")


	return OptimizedInlet(Impeller=impeller, H_rise=H, Sizing=sizing, Feasible=feasible, Binding=binding,
		Violations=violations)

def PrintOptimizedInlet(result: OptimizedInlet, baseline: ImpellerInlet, name: str = 'Inducer') -> None:
	imp = result.Impeller
	change = lambda new, old: 100 * (new.to('m').m / old.to('m').m - 1)
	print(f"\n#### {name}: {'largest feasible' if result.Feasible else 'NO FEASIBLE'} impeller inlet ####")
	print(f"Impeller d_1               {imp.d_1.to('mm'):.3f} ({change(imp.d_1, baseline.d_1):+.1f}% from baseline {baseline.d_1.to('mm'):.3f})")
	print(f"Impeller d_hub             {imp.d_hub.to('mm'):.3f} ({change(imp.d_hub, baseline.d_hub):+.1f}% from baseline {baseline.d_hub.to('mm'):.3f})")
	print(f"Inducer head requirement   H_rise = {result.H_rise.to('m'):.1f}")
	for item in result.Binding:
		print(f"Binding: {item}")
	for item in result.Violations:
		print(f"VIOLATED: {item}")
	PrintInletSizing(result.Sizing, name)

# ---------------------------------------------------------------------------------------------------------------------
# Pulling the impeller inlet and design point from the impeller code
# ---------------------------------------------------------------------------------------------------------------------

def ImpellerInletFromImpeller(impeller) -> ImpellerInlet:
	"""Impeller.Impeller (pump.Impeller): its inlet diameter d_1 and hub diameter d_hub"""
	return ImpellerInlet(d_1=impeller.d_1, d_hub=impeller.d_hub)

def SizeInducerFromPump(pump, base: SizingConstants = constants, limits: InletSizingLimits = InletSizingLimits(),
		tolerance: InletTolerance | None = None, head_model: HeadRequirementModel | None = None) -> InletSizing:
	"""
		Size the inducer from a Pump.Pump. The impeller inlet is pump.Impeller.d_1 and d_hub, flow and speed come from the pump.
		tolerance None: the impeller is kept as it is and the head requirement is Pump1D.GetInducerHeadRequirement.
		tolerance given: the impeller inlet may move within it (OptimizeImpellerInlet, largest feasible inlet), the new
		impeller inlet is in the result's D_tip and D_hub_2, the impeller itself is not changed.
	"""
	baseline = ImpellerInletFromImpeller(pump.Impeller)
	model = head_model or HeadRequirementModel(NPSH_a=pump.NPSH_a)
	design_point = InducerDesignPoint(Q=pump.DesignFlowRate, N_shaft=pump.DP.N_shaft, H_rise=Q_(0, 'm'))
	if tolerance is not None:
		return OptimizeImpellerInlet(design_point, base, baseline, model, tolerance, limits).Sizing
	H = model.Head(design_point.Q, design_point.N_shaft, baseline.d_1.to('m').m, baseline.d_hub.to('m').m)
	return SizeInletForImpeller(replace(design_point, H_rise=Q_(float(H), 'm')), base, baseline, limits)


if __name__ == '__main__':
	# ------------------------------------------------------------------------------------------
	# EDIT THE NUMBERS BELOW, one entry per inducer, then run:  python InducerInletSizing.py
	#	Q, N_shaft:		flow rate and shaft speed
	#	NPSH_a:			net positive suction head available at the pump inlet
	#	d_1, d_hub:		baseline impeller inlet (tip, shroud) and hub diameters, they may move within the tolerance
	# Any pint unit works, e.g. Q_(0.94, 'in'), Q_(25, 'gal/min'), Q_(3000, 'rad/s').
	# ------------------------------------------------------------------------------------------
	cases = {
		'LOx': dict(
			design_point=InducerDesignPoint(Q=Q_(1.771, 'L/s'), N_shaft=Q_(35000, 'rpm'), H_rise=Q_(0, 'm')),		# H_rise is solved
			NPSH_a=Q_(52.3, 'm'),		# the pump's NPSH_a (Pump.NPSH_a)
			baseline=ImpellerInlet(d_1=Q_(1, 'in'), d_hub=Q_(0.6, 'in'))		# from Main.py basic_impeller_geometry
		),
		'IPA': dict(
			design_point=InducerDesignPoint(Q=Q_(2.735, 'L/s'), N_shaft=Q_(35000, 'rpm'), H_rise=Q_(0, 'm')),
			NPSH_a=Q_(88.8, 'm'),
			baseline=ImpellerInlet(d_1=Q_(1.2, 'in'), d_hub=Q_(0.7, 'in'))
		),
	}
	tolerance = InletTolerance(d_1=0.20, d_hub=0.20)		# how far the impeller inlet may move from the baseline
	# ------------------------------------------------------------------------------------------

	USE_PUMP = True		# True: take the impeller inlet, flow, speed and NPSH_a from the pump built in Main.py instead of the cases above

	if USE_PUMP:
		import sys
		from pathlib import Path
		# Put Code/ first and drop this folder so the Inducer package (not Inducer.py) is imported by Pump.py
		_here = Path(__file__).resolve().parent
		sys.path[:] = [p for p in sys.path if Path(p or '.').resolve() != _here]
		sys.path.insert(0, str(_here.parent))
		import Main
		for name in (Main.Propellant.LOX, Main.Propellant.IPA):
			sizing = SizeInducerFromPump(Main.CreatePump(name), tolerance=tolerance)		# tolerance=None keeps the impeller as it is
			PrintInletSizing(sizing, name)
			PlotContour(sizing, name=name).show()
		raise SystemExit

	for name, case in cases.items():
		result = OptimizeImpellerInlet(
			case['design_point'], constants, case['baseline'], HeadRequirementModel(NPSH_a=case['NPSH_a']), tolerance
		)
		PrintOptimizedInlet(result, case['baseline'], name)
		PlotContour(result.Sizing, name=name).show()