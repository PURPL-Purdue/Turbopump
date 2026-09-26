"""Preliminary inducer blade sizing (Huzel & Huang, "Modern Engineering for
Design of Liquid-Propellant Rocket Engines", inducer design section).

All calculations use US customary units, matching the empirical constants in
Huzel & Huang:
    flow        gal/min (gpm)
    head, NPSH  ft
    speed       rpm
    lengths     ft
    velocities  ft/s
    angles      rad internally, printed in deg
"""

import math
import warnings
from VelocityTriangle import VelocityTriangle
from dataclasses import dataclass

# --- Unit conversions -------------------------------------------------------
M_TO_FT = 3.28084
LPS_TO_GPM = 15.85032
GPM_TO_CFS = 1.0 / 448.83
G = 32.174  # gravitational acceleration, ft/s^2

# --- Empirical constants (Huzel & Huang) ------------------------------------
# D_t = K_DT * (Q / ((1 - rho^2) * n * phi))^(1/3)  ->  D_t in ft for Q in gpm.
# K_DT = (4 * 60 / (448.83 * pi^2))^(1/3), where 448.83 gpm = 1 ft^3/s.
K_DT = 0.37843
# Constant in the corrected suction specific speed vs. flow coefficient
# relation: N_ss' = K_NSS * (1 - 2 phi^2)^(3/4) / phi
K_NSS = 3574.0
# Flow angle / blade angle ratio (attack angle is ~42.5% of blade angle).
FLOW_TO_BLADE_ANGLE = 0.575


@dataclass
class InducerInputs:
    flow_gpm: float           # Q, volumetric flow
    speed_rpm: float          # n, shaft speed
    npsh_ft: float            # available NPSH at inducer inlet
    inducer_head_ft: float    # head rise across the inducer
    hub_tip_ratio: float      # rho = D_h / D_t
    blade_efficiency: float   # eta_bl, inducer hydraulic efficiency
    num_blades: int           # N
    solidity: float           # sigma = chord / blade spacing


@dataclass
class InducerResults:
    suction_specific_speed: float
    corrected_suction_specific_speed: float
    specific_speed: float
    flow_coefficient: float
    tip_diameter_ft: float
    hub_diameter_ft: float
    tip_speed_fps: float
    meridional_velocity_fps: float
    tangential_velocity_fps: float
    inlet_flow_angle_rad: float
    inlet_blade_angle_rad: float
    outlet_flow_angle_rad: float
    blade_spacing_ft: float
    chord_ft: float


# --- Performance parameters -------------------------------------------------
def specific_speed(n: float, q: float, head: float) -> float:
    """Pump specific speed: N_s = n * Q^0.5 / H^0.75."""
    return n * math.sqrt(q) / head**0.75


def suction_specific_speed(n: float, q: float, npsh: float) -> float:
    """Suction specific speed: N_ss = n * Q^0.5 / NPSH^0.75."""
    return n * math.sqrt(q) / npsh**0.75
    # TODO this uses critical NPSH, not avaliable NPSH. Which one do we need here?  


def npsh_for_suction_specific_speed(n: float, q: float, n_ss: float) -> float:
    """NPSH that gives a target N_ss: inverse of suction_specific_speed()."""
    return (n * math.sqrt(q) / n_ss) ** (4.0 / 3.0)


def npsh_required(
    d_tip: float, q: float, n: float, rho: float, sigma_b: float,
    lambda_c: float = 1.2,
) -> float:
    """Required NPSH (ft) at the inducer tip from the blade cavitation number.

    NPSH_r = lambda_c * C_m^2 / 2g + sigma_b * W_t^2 / 2g, with
    W_t^2 = C_m^2 + u_t^2 (no inlet prerotation).

    d_tip in ft, q in gpm, n in rpm. lambda_c covers inlet acceleration and
    losses on the absolute velocity (Gulich gives ~1.1 to 1.35).
    """
    c_m = meridional_velocity(q, d_tip, rho)
    u = tip_speed(d_tip, n)
    return (lambda_c * c_m**2 + sigma_b * (c_m**2 + u**2)) / (2.0 * G)


def corrected_suction_specific_speed(n_ss: float, rho: float) -> float:
    """N_ss based on the net (annular) inlet flow area: N_ss / sqrt(1 - rho^2)."""
    return n_ss / math.sqrt(1.0 - rho**2)


def flow_coefficient(n_ss_corrected: float, exact: bool = False) -> float:
    """Inlet flow coefficient phi = C_m / u_t for a given corrected N_ss.

    The underlying relation is N_ss' = K_NSS * (1 - 2 phi^2)^(3/4) / phi.

    exact=False: closed-form solution of the small-phi linearization
        (1 - 2 phi^2)^(3/4) ~= 1 - 1.5 phi^2, accurate to <0.1% for phi < 0.25.
    exact=True: solve the full relation numerically by bisection.
    """
    if not exact:
        k = K_NSS / n_ss_corrected
        return k / ((1.0 + math.sqrt(1.0 + 6.0 * k**2)) / 2.0)

    def residual(phi: float) -> float:
        return K_NSS * (1.0 - 2.0 * phi**2) ** 0.75 / phi - n_ss_corrected

    # residual is monotonically decreasing on (0, 1/sqrt(2)).
    lo, hi = 1e-9, 1.0 / math.sqrt(2.0) - 1e-9
    for _ in range(200):
        mid = 0.5 * (lo + hi)
        if residual(mid) > 0.0:
            lo = mid
        else:
            hi = mid
    return 0.5 * (lo + hi)


# --- Geometry ---------------------------------------------------------------
def tip_diameter(q: float, n: float, rho: float, phi: float) -> float:
    """Inducer tip diameter in ft (Q in gpm, n in rpm)."""
    return K_DT * (q / ((1.0 - rho**2) * n * phi)) ** (1.0 / 3.0)


def tip_speed(d_tip: float, n: float) -> float:
    """Blade tip speed u_t = pi * D_t * n / 60, in ft/s."""
    return math.pi * d_tip * n / 60.0


def meridional_velocity(q: float, d_tip: float, rho: float) -> float:
    """Axial inlet velocity C_m through the annulus, in ft/s (Q in gpm)."""
    area = math.pi / 4.0 * d_tip**2 * (1.0 - rho**2)
    return q * GPM_TO_CFS / area


def blade_spacing(diameter: float, num_blades: int) -> float:
    """Circumferential blade spacing S = pi * D / N."""
    return math.pi * diameter / num_blades


# --- Velocity triangles (at tip, no inlet prerotation) ----------------------
def tangential_velocity(head: float, u: float, eta: float) -> float:
    """Outlet absolute swirl from Euler's equation: C_u = g * H / (eta * u)."""
    return G * head / (eta * u)


def inlet_flow_angle(c_m: float, u: float) -> float:
    """Inlet relative flow angle measured from the tangential direction."""
    return math.atan2(c_m, u)


def outlet_flow_angle(c_m: float, u: float, c_u: float) -> float:
    """Outlet relative flow angle measured from the tangential direction."""
    if c_u >= u:
        warnings.warn(
            f"C_u ({c_u:.1f} ft/s) >= u ({u:.1f} ft/s): the inducer head is "
            "too high for this tip speed; the outlet velocity triangle is "
            "not physical."
        )
    return math.atan2(c_m, u - c_u)


# --- Top-level sizing -------------------------------------------------------
def size_inducer(inp: InducerInputs, exact_phi: bool = False) -> InducerResults:
    n, q, rho = inp.speed_rpm, inp.flow_gpm, inp.hub_tip_ratio

    n_ss = suction_specific_speed(n, q, inp.npsh_ft)
    n_ss_corr = corrected_suction_specific_speed(n_ss, rho)
    phi = flow_coefficient(n_ss_corr, exact=exact_phi)

    d_tip = tip_diameter(q, n, rho, phi)
    d_hub = rho * d_tip

    u = tip_speed(d_tip, n)
    c_m = phi * u
    c_u = tangential_velocity(inp.inducer_head_ft, u, inp.blade_efficiency)

    gamma_in = inlet_flow_angle(c_m, u)
    beta_in = gamma_in / FLOW_TO_BLADE_ANGLE
    gamma_out = outlet_flow_angle(c_m, u, c_u)

    # Solidity is conventionally defined at the tip for inducers.
    spacing = blade_spacing(d_tip, inp.num_blades)
    chord = inp.solidity * spacing

    return InducerResults(
        suction_specific_speed=n_ss,
        corrected_suction_specific_speed=n_ss_corr,
        specific_speed=specific_speed(n, q, inp.inducer_head_ft),
        flow_coefficient=phi,
        tip_diameter_ft=d_tip,
        hub_diameter_ft=d_hub,
        tip_speed_fps=u,
        meridional_velocity_fps=c_m,
        tangential_velocity_fps=c_u,
        inlet_flow_angle_rad=gamma_in,
        inlet_blade_angle_rad=beta_in,
        outlet_flow_angle_rad=gamma_out,
        blade_spacing_ft=spacing,
        chord_ft=chord,
    )


def print_results(res: InducerResults) -> None:
    rows = [
        ("Suction specific speed N_ss", res.suction_specific_speed, ""),
        ("Corrected N_ss'", res.corrected_suction_specific_speed, ""),
        ("Specific speed N_s", res.specific_speed, ""),
        ("Flow coefficient phi", res.flow_coefficient, ""),
        ("Tip diameter D_t", res.tip_diameter_ft * 12.0, "in"),
        ("Hub diameter D_h", res.hub_diameter_ft * 12.0, "in"),
        ("Tip speed u", res.tip_speed_fps, "ft/s"),
        ("Meridional velocity C_m", res.meridional_velocity_fps, "ft/s"),
        ("Tangential velocity C_u", res.tangential_velocity_fps, "ft/s"),
        ("Inlet flow angle gamma_1", math.degrees(res.inlet_flow_angle_rad), "deg"),
        ("Inlet blade angle beta_1", math.degrees(res.inlet_blade_angle_rad), "deg"),
        ("Outlet flow angle gamma_2", math.degrees(res.outlet_flow_angle_rad), "deg"),
        ("Blade spacing S (tip)", res.blade_spacing_ft * 12.0, "in"),
        ("Chord C", res.chord_ft * 12.0, "in"),
    ]
    for label, value, unit in rows:
        print(f"{label:<28} {value:>12.4f} {unit}")


# TODO: blade lead/pitch, delta = 2 * pi * r * tan(gamma), at hub and tip
# TODO: blade angle change, beta_2 - beta_1

DESIGN_INPUTS = InducerInputs(
    flow_gpm=1.69 * LPS_TO_GPM,
    speed_rpm=40_000,
    npsh_ft=10.0 * M_TO_FT,
    # NOTE: this is the full pump head (140 m). The inducer normally
    # supplies only a fraction of this; update once the inducer head
    # rise is chosen.
    inducer_head_ft=140.0 * M_TO_FT,
    hub_tip_ratio=0.3,
    blade_efficiency=0.85,
    num_blades=3,
    solidity=2.5,
)

if __name__ == "__main__":
    print_results(size_inducer(DESIGN_INPUTS))
