# ===========================
# Manifold + Orifice Utilities
# (Fuel and Ox can have different heights / widths / radii)
# ===========================

import math
import numpy as np


# ---------------------------
# Geometry helpers
# ---------------------------

def rounded_rect_section(w_m: float, h_m: float, r_m: float):
    """
    Cross-section area A, wetted perimeter P, and hydraulic diameter Dh
    for a rectangle with 4 rounded corners of radius r.

    Units: meters.
    """
    w = float(w_m)
    h = float(h_m)
    r = float(r_m)

    if w <= 0 or h <= 0:
        raise ValueError(f"w and h must be > 0. Got w={w}, h={h}")
    if r < 0:
        raise ValueError(f"r must be >= 0. Got r={r}")

    rmax = 0.5 * min(w, h)
    if r > rmax:
        raise ValueError(f"Corner radius too large: r={r} > {rmax} (half of min(w,h)).")

    # Area = full rectangle - 4 corner squares + 4 quarter circles = w*h - 4 r^2 + pi r^2
    A = (w * h) - (4.0 * r * r) + (math.pi * r * r)

    # Wetted perimeter = straight segments + arcs
    # Straight lengths: 2*(w - 2r) + 2*(h - 2r) = 2*(w + h - 4r)
    # Arcs: 4 quarter circles = 2*pi*r
    P = 2.0 * (w + h - 4.0 * r) + (2.0 * math.pi * r)

    Dh = 4.0 * A / P
    return A, P, Dh


def design_two_manifolds(
    *,
    # Fuel section
    w_fuel_m: float,
    h_fuel_m: float,
    r_fuel_m: float,
    R_fuel_m: float,
    # Ox section
    w_ox_m: float,
    h_ox_m: float,
    r_ox_m: float,
    R_ox_m: float,
    # Ring effective travel factor
    # 1.0 -> use full circumference (2*pi*R)
    # 0.5 -> use half circumference (pi*R) for "average path" ring model
    travel_factor: float = 1.0,
):
    """
    Returns geometry dict for two separate annular manifolds (fuel and ox).

    travel_factor lets you tune effective path length:
      L = travel_factor * (2*pi*R)
    """
    tf = float(travel_factor)
    if tf <= 0:
        raise ValueError("travel_factor must be > 0")

    Rf = float(R_fuel_m)
    Ro = float(R_ox_m)
    if Rf <= 0 or Ro <= 0:
        raise ValueError("R_fuel_m and R_ox_m must be > 0")

    A_f, P_f, Dh_f = rounded_rect_section(w_fuel_m, h_fuel_m, r_fuel_m)
    A_o, P_o, Dh_o = rounded_rect_section(w_ox_m, h_ox_m, r_ox_m)

    L_f = tf * (2.0 * math.pi * Rf)
    L_o = tf * (2.0 * math.pi * Ro)

    return {
        "fuel": {"A": A_f, "P": P_f, "Dh": Dh_f, "R": Rf, "L": L_f, "w": w_fuel_m, "h": h_fuel_m, "r": r_fuel_m},
        "ox":   {"A": A_o, "P": P_o, "Dh": Dh_o, "R": Ro, "L": L_o, "w": w_ox_m,   "h": h_ox_m,   "r": r_ox_m},
    }


# ---------------------------
# Fluid / pressure drop
# ---------------------------

def friction_factor_darcy(Re: float, rel_rough: float):
    """
    Darcy friction factor:
      - Laminar: 64/Re
      - Turbulent: Haaland explicit approximation

    Re: Reynolds number
    rel_rough: epsilon/Dh
    """
    Re = float(Re)
    rr = float(rel_rough)

    if Re <= 0:
        return float("nan")
    if rr < 0:
        raise ValueError("rel_rough must be >= 0")

    if Re < 2300.0:
        return 64.0 / Re

    # Haaland (Darcy)
    # 1/sqrt(f) = -1.8 log10( ( (rr/3.7)^1.11 ) + 6.9/Re )
    inv_sqrt_f = -1.8 * math.log10((rr / 3.7) ** 1.11 + 6.9 / Re)
    return 1.0 / (inv_sqrt_f ** 2)


def manifold_dp_darcy(
    *,
    mdot_kg_s: float,
    rho_kg_m3: float,
    mu_Pa_s: float,
    L_m: float,
    A_m2: float,
    Dh_m: float,
    rough_m: float = 1.5e-5,
    K_total: float = 0.0,
):
    """
    Pressure drop (Pa) using Darcy–Weisbach + minor losses:
      dp = f (L/Dh) (rho v^2 / 2) + K_total (rho v^2 / 2)

    All SI units.
    """
    mdot = float(mdot_kg_s)
    rho = float(rho_kg_m3)
    mu = float(mu_Pa_s)
    L = float(L_m)
    A = float(A_m2)
    Dh = float(Dh_m)
    eps = float(rough_m)
    K = float(K_total)

    if mdot <= 0:
        return 0.0
    if rho <= 0 or mu <= 0:
        raise ValueError(f"rho and mu must be > 0. Got rho={rho}, mu={mu}")
    if L <= 0 or A <= 0 or Dh <= 0:
        raise ValueError(f"L, A, Dh must be > 0. Got L={L}, A={A}, Dh={Dh}")
    if eps < 0 or K < 0:
        raise ValueError(f"rough_m and K_total must be >= 0. Got rough_m={eps}, K_total={K}")

    v = mdot / (rho * A)  # m/s
    Re = rho * v * Dh / mu
    rr = eps / Dh

    f = friction_factor_darcy(Re, rr)
    if not np.isfinite(f):
        raise ValueError(f"Bad friction factor. Re={Re}, rr={rr}, v={v}, Dh={Dh}")

    q = 0.5 * rho * v * v
    dp_fric = f * (L / Dh) * q
    dp_minor = K * q
    return dp_fric + dp_minor


# ---------------------------
# Orifice / injector helpers
# ---------------------------

def orifice_diameter_from_mdot(
    *,
    mdot_kg_s: float,
    rho_kg_m3: float,
    Cd: float,
    deltaP_Pa: float,
    n_holes: int,
):
    """
    Solves for hole diameter (m) for a set of identical holes:
      mdot_total = n * Cd * A * sqrt(2 rho dP)

    Returns diameter in meters.
    """
    mdot = float(mdot_kg_s)
    rho = float(rho_kg_m3)
    Cd = float(Cd)
    dP = float(deltaP_Pa)
    n = int(n_holes)

    if n <= 0:
        raise ValueError("n_holes must be >= 1")
    if mdot <= 0:
        return 0.0
    if rho <= 0 or Cd <= 0 or dP <= 0:
        raise ValueError(f"rho, Cd, deltaP must be > 0. Got rho={rho}, Cd={Cd}, dP={dP}")

    A_one = (mdot / n) / (Cd * math.sqrt(2.0 * rho * dP))
    if A_one <= 0:
        raise ValueError("Computed non-positive area for orifice.")

    d = math.sqrt(4.0 * A_one / math.pi)
    return d


def mdot_through_orifices(
    *,
    d_m: float,
    rho_kg_m3: float,
    Cd: float,
    deltaP_Pa: float,
    n_holes: int,
):
    """
    Forward model:
      mdot_total = n * Cd * A * sqrt(2 rho dP)
    """
    d = float(d_m)
    rho = float(rho_kg_m3)
    Cd = float(Cd)
    dP = float(deltaP_Pa)
    n = int(n_holes)

    if n <= 0:
        raise ValueError("n_holes must be >= 1")
    if d <= 0:
        return 0.0
    if rho <= 0 or Cd <= 0 or dP <= 0:
        raise ValueError(f"rho, Cd, deltaP must be > 0. Got rho={rho}, Cd={Cd}, dP={dP}")

    A_one = math.pi * (d * d) / 4.0
    return n * Cd * A_one * math.sqrt(2.0 * rho * dP)


# ---------------------------
# Convenience wrapper:
# compute fuel + ox DP given geometry
# ---------------------------

def compute_two_manifold_dps(
    *,
    geo: dict,
    # Fuel fluid props
    mdot_fuel: float,
    rho_fuel: float,
    mu_fuel: float,
    rough_fuel_m: float = 1.5e-5,
    K_fuel: float = 0.0,
    # Ox fluid props
    mdot_ox: float,
    rho_ox: float,
    mu_ox: float,
    rough_ox_m: float = 1.5e-5,
    K_ox: float = 0.0,
):
    """
    Returns (dp_fuel_Pa, dp_ox_Pa)
    """
    gF = geo["fuel"]
    gO = geo["ox"]

    dp_fuel = manifold_dp_darcy(
        mdot_kg_s=mdot_fuel,
        rho_kg_m3=rho_fuel,
        mu_Pa_s=mu_fuel,
        L_m=gF["L"],
        A_m2=gF["A"],
        Dh_m=gF["Dh"],
        rough_m=rough_fuel_m,
        K_total=K_fuel,
    )

    dp_ox = manifold_dp_darcy(
        mdot_kg_s=mdot_ox,
        rho_kg_m3=rho_ox,
        mu_Pa_s=mu_ox,
        L_m=gO["L"],
        A_m2=gO["A"],
        Dh_m=gO["Dh"],
        rough_m=rough_ox_m,
        K_total=K_ox,
    )

    return dp_fuel, dp_ox


# ---------------------------
# Example usage (edit these values)
# ---------------------------
if __name__ == "__main__":
    # Example ring radii (meters)
    R_fuel_outer = 0.055  # 55 mm
    R_ox_outer   = 0.052  # 52 mm

    # Example: different heights for fuel and ox
    geo = design_two_manifolds(
        w_fuel_m=0.006,  # 6 mm width
        h_fuel_m=0.004,  # 4 mm height (fuel)
        r_fuel_m=0.0008, # 0.8 mm corner radius
        R_fuel_m=R_fuel_outer,

        w_ox_m=0.006,    # 6 mm width
        h_ox_m=0.003,    # 3 mm height (ox)
        r_ox_m=0.0006,   # 0.6 mm corner radius
        R_ox_m=R_ox_outer,

        travel_factor=0.5,  # use pi*R as "average path length" model
    )

    # Example fluid properties (NOT authoritative; plug your real values)
    # RP-1-ish at room temp
    mdot_fuel = 0.9
    rho_fuel  = 810.0
    mu_fuel   = 1.8e-3

    # LOX-ish
    mdot_ox = 1.9
    rho_ox  = 1140.0
    mu_ox   = 2.0e-4

    dp_f, dp_o = compute_two_manifold_dps(
        geo=geo,
        mdot_fuel=mdot_fuel, rho_fuel=rho_fuel, mu_fuel=mu_fuel, K_fuel=1.0,
        mdot_ox=mdot_ox,     rho_ox=rho_ox,     mu_ox=mu_ox,     K_ox=1.0,
    )

    print("=== Geometry ===")
    print("Fuel:", geo["fuel"])
    print("Ox  :", geo["ox"])
    print("\n=== DP ===")
    print(f"dp_fuel = {dp_f/1e5:.4f} bar ({dp_f:.1f} Pa)")
    print(f"dp_ox   = {dp_o/1e5:.4f} bar ({dp_o:.1f} Pa)")

    # Example orifice sizing
    d_fuel = orifice_diameter_from_mdot(
        mdot_kg_s=mdot_fuel,
        rho_kg_m3=rho_fuel,
        Cd=0.75,
        deltaP_Pa=2.0e6,  # 20 bar
        n_holes=120
    )
    print(f"\nFuel orifice diameter for 120 holes @20 bar: {d_fuel*1e3:.3f} mm")
