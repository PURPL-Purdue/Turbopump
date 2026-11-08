# Maximum inducer blade thickness at the hub from two ways of getting the hub relative velocity ratio w1/w2 (NASA SP-8052 sec. 2.1.5, eq. 28):
#   t_max = s (sin(beta) - (w1/w2) sin(beta - alpha)),  s = 2 pi r_hub_1 / Z
#   Method 1: hub inlet and outlet triangles from the radial equilibrium, iterated because the blade blockage depends on t_max
#   Method 2: SP-8052 cavity relations (eqs. 5, 25-27) with K from the optimum flow coefficient of lambda_w, lambda_c
import sys
from dataclasses import dataclass
from pathlib import Path

import numpy as np
from pint import Quantity as Q_

if not __package__:
    # Run directly as a script: put Code/ first and drop this folder so the Inducer package (not Inducer.py) is imported
    _here = Path(__file__).resolve().parent
    sys.path[:] = [p for p in sys.path if Path(p or '.').resolve() != _here]
    sys.path.insert(0, str(_here.parent))

from Inducer.Inducer import Inducer, InducerInputGeometry
from Inducer.MeanlineTriangle import INCIDENCE
from Inducer.RadialEquilibrium import InducerRadialEquilibrium
from Inducer.Sizing import InletSizing

Z_BLADE = 2                 # same as Pump.AddInducer


@dataclass
class ThicknessResult:
    method: str
    w_ratio: float          # w1 / w2 at the hub
    beta: float             # hub blade angle at the inlet (rad)
    alpha: float            # hub incidence angle at the inlet (rad)
    pitch: Q_               # blade spacing at the hub inlet
    t_max: Q_
    iterations: int | None = None
    K: float | None = None  # cavitation number from phi_opt (method 2)
    K_used: float | None = None     # K, raised to K_min where the hub is supercavitating
    K_min: float | None = None      # eq. 29 at the hub
    wedge: float | None = None      # blade wedge angle at the hub inlet (rad), eq. 69


def WedgeAngle(inducer: Inducer, beta: float) -> float:
    """SP-8052 eqs. 69-70: blade wedge angle alpha_w = beta - arctan(1.10 phi_d), with phi_d the design flow coefficient. beta and the result in rad"""
    return beta - np.arctan(1.10 * inducer.FlowCoefficient)


def InducerFromSizing(sizing: InletSizing, Z_blade: int = Z_BLADE) -> Inducer:
    """Inducer with the hub geometry solved in Sizing.py"""
    geometry = InducerInputGeometry(
        d_tip=sizing.D_tip, d_hub_2=sizing.D_hub_2, d_hub_1=sizing.D_hub_1, Z_blade=Z_blade,
        axial_length_ratio=float((sizing.L_ax / sizing.D_tip).to('').m),
    )
    return Inducer(sizing.design_point, geometry, lambda_c=sizing.constants.lambda_c, lambda_w=sizing.lambda_w)


def _MaxThickness(pitch: float, beta: float, alpha: float, w_ratio: float) -> float:
    """SP-8052 eq. 28: h_c / s = sin(beta) - (w1/w2) sin(beta - alpha), the thickest blade the cavity can hold"""
    return pitch * (np.sin(beta) - w_ratio * np.sin(beta - alpha))


def Method1RadialEquilibrium(inducer: Inducer, incidence: Q_ = Q_(INCIDENCE, 'deg'), tol: float = 1e-6,
                             max_iter: int = 100, relax: float = 0.7) -> ThicknessResult:
    """Iterate t -> radial equilibrium with blockage t -> hub w1/w2, beta, alpha -> t until t and w1/w2 stop changing"""
    r_h1, r_h2, L = inducer.R_hub_1, inducer.R_hub_2, inducer.L_ax
    z_in = Q_(0, 'm')
    pitch = 2 * np.pi * r_h1.to('m').m / inducer.Z_blade
    t, w_ratio = 0.0, np.nan
    for i in range(1, max_iter + 1):
        eq = InducerRadialEquilibrium(inducer, incidence=incidence, blade_thickness=Q_(t, 'm'))
        w_1 = eq.RelativeVelocity(r_h1, z_in)[2].to('m/s').m
        w_2 = eq.RelativeVelocity(r_h2, L)[2].to('m/s').m
        beta = eq.BladeAngle(r_h1, z_in).to('rad').m
        alpha = beta - eq.FlowAngle(r_h1, z_in).to('rad').m
        w_new = w_1 / w_2
        t_new = max(_MaxThickness(pitch, beta, alpha, w_new), 0.0)
        converged = abs(t_new - t) < tol * pitch and abs(w_new - w_ratio) < tol
        t, w_ratio = t + relax * (t_new - t), w_new
        if converged:
            return ThicknessResult('1: radial equilibrium', w_ratio, beta, alpha, Q_(pitch, 'm'), Q_(t_new, 'm'), iterations=i,
                                 wedge=WedgeAngle(inducer, beta))
    raise RuntimeError(f"Method 1 did not converge in {max_iter} iterations (t = {t:.3e} m, w1/w2 = {w_ratio:.4f})")


def OptimumCavitationNumber(phi: float) -> float:
    """SP-8052 eq. 5, K = phi_opt^2 / (1 - 2 phi_opt^2), i.e. eq. 6 solved for K"""
    return phi**2 / (1 - 2 * phi**2)


def Method2CavitationRelations(inducer: Inducer, incidence: Q_ = Q_(INCIDENCE, 'deg')) -> ThicknessResult:
    """w1/w2 from the cavity (wake) model, SP-8052 eqs. 25-27, with K from phi(lambda_w, lambda_c) as in 1DPerformance.py"""
    K = OptimumCavitationNumber(inducer.FlowCoefficient)
    r_h1 = inducer.R_hub_1.to('m').m
    omega = inducer.DP.N_shaft.to('rad/s').m
    c_m = inducer.DP.Q.to('m^3/s').m / (np.pi * (inducer.R_tip.to('m').m**2 - r_h1**2))   # eq. 23, zero blockage
    gamma = np.arctan(c_m / (omega * r_h1))                       # flow angle at the hub, zero prewhirl
    alpha = incidence.to('rad').m
    beta = gamma + alpha                                          # blade angle = flow angle + incidence

    # eq. 27 (the sign between the two terms was checked against eq. 29: F = 1 exactly at K = K_min)
    F_of = lambda K: ((1 + K)**0.5 * np.sin(beta) + (1 + K)**-0.5 * np.sin(beta - 2 * alpha)) / (2 * np.sin(beta - alpha))
    # eq. 29: supercavitation limit. The tip K can be below it at the hub (larger blade angle), where eq. 26 has no solution
    K_min = 2 * np.sin(alpha) * np.sin(beta - alpha) / (1 + np.cos(beta))
    K_used = max(K, K_min)
    F = F_of(K_used)
    w_c_over_w2 = F + np.sqrt(max(F**2 - 1, 0.0))   # eq. 26
    w_ratio = w_c_over_w2 / np.sqrt(1 + K_used)     # eq. 25, w_c = w_1 sqrt(1 + K)
    pitch = 2 * np.pi * r_h1 / inducer.Z_blade
    t = max(_MaxThickness(pitch, beta, alpha, w_ratio), 0.0)
    return ThicknessResult('2: SP-8052 cavitation', w_ratio, beta, alpha, Q_(pitch, 'm'), Q_(t, 'm'), K=K, K_used=K_used, K_min=K_min,
                           wedge=WedgeAngle(inducer, beta))


def CompareInducerBladeThickness(inducer: Inducer, name: str = 'Inducer',
                                 incidence: Q_ = Q_(INCIDENCE, 'deg')) -> tuple[ThicknessResult, ThicknessResult]:
    m1 = Method1RadialEquilibrium(inducer, incidence)
    m2 = Method2CavitationRelations(inducer, incidence)

    print(f"\n===== {name} inducer hub blade thickness =====")
    print(f"Hub radius inlet / outlet : {inducer.R_hub_1.to('mm'):.3f} / {inducer.R_hub_2.to('mm'):.3f}, tip {inducer.R_tip.to('mm'):.3f}, Z = {inducer.Z_blade}")
    print(f"Incidence                 : {incidence.to('deg'):.2f}")
    print(f"phi, K, K_min             : {inducer.FlowCoefficient:.4f}, {m2.K:.4f}, {m2.K_min:.4f}")
    print(f"{'Method':<26}{'w1/w2':>9}{'beta':>10}{'alpha':>10}{'wedge':>10}{'pitch':>12}{'t_max':>12}")
    for m in (m1, m2):
        print(f"{m.method:<26}{m.w_ratio:>9.4f}{np.degrees(m.beta):>8.2f}deg{np.degrees(m.alpha):>7.2f}deg"
                  f"{np.degrees(m.wedge):>7.2f}deg{m.pitch.to('mm').m:>10.3f}mm{m.t_max.to('mm').m:>10.4f}mm")
    print(f"Method 1 converged in {m1.iterations} iterations; t_max difference (1 - 2) = {(m1.t_max - m2.t_max).to('mm'):.4f}")
    if m2.K_used > m2.K:
        print(f"Note: K = {m2.K:.4f} is below the hub supercavitation limit K_min = {m2.K_min:.4f} (eq. 29), method 2 uses K = K_min at the hub")
    return m1, m2


def CompareBladeThickness(sizing: InletSizing, name: str = 'Inducer', Z_blade: int = Z_BLADE) -> tuple[ThicknessResult, ThicknessResult]:
    return CompareInducerBladeThickness(InducerFromSizing(sizing, Z_blade), name)


if __name__ == "__main__":
    import Main
    import Pump1D

    incidence = Q_(INCIDENCE, 'deg')
    for propellant in (Main.Propellant.LOX, Main.Propellant.IPA):
        pump = Main.CreatePump(propellant)
        inducer_head = Pump1D.GetInducerlessNPSH_r(pump)
        pump.AddInducer(inducer_head, geometry=Main.basic_inducer_geometry[propellant])
        CompareInducerBladeThickness(pump.Inducer, propellant.name, incidence)
