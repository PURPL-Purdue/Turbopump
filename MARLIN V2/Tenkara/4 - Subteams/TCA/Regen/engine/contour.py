"""
Chamber / nozzle wall contour from the sized engine: geometry only.
Printing is in util/report.py, the DXF export in util/export.py.

The divergent section is a Rao bell (thrust-optimised parabola). Its wall
angles come from the team's contour_script.py, vendored as engine/rao_bell.py.
"""
import math
import numpy as np

from engine import engine_inputs as ei
from engine.rao_bell import find_wall_angles


def _bezier_t(xq, xN, xQ, xE):
    """Parameter t in [0, 1] where the bell Bezier's x(t) = xq (x(t) is monotonic)."""
    a = xN - 2 * xQ + xE
    b = 2 * (xQ - xN)
    c = xN - np.asarray(xq, dtype=float)
    if abs(a) < 1e-15:
        t = -c / b
    else:
        t = (-b + np.sqrt(np.clip(b * b - 4 * a * c, 0, None))) / (2 * a)
    return np.clip(t, 0, 1)


def wall_angle(xq, d):
    """
    Wall slope angle phi(x) in rad from the dims of build_contour, exact per
    segment: positive where r grows with x. The outward wall normal (away
    from the axis) is (-sin phi, cos phi).
    """
    xq = np.asarray(xq, dtype=float)
    xt, th_c = d["x_t"], math.radians(d["theta_conv_deg"])
    (xN, rN), (xQ, rQ), (xE, rE) = d["N"], d["Q"], d["E"]

    def bell(x):
        t = _bezier_t(x, xN, xQ, xE)
        dx = 2 * (1 - t) * (xQ - xN) + 2 * t * (xE - xQ)
        dr = 2 * (1 - t) * (rQ - rN) + 2 * t * (rE - rQ)
        return np.arctan2(dr, dx)

    return np.piecewise(xq,
        [xq <= d["x0"],
         (xq > d["x0"]) & (xq <= d["x1"]),
         (xq > d["x1"]) & (xq <= xt),
         (xq > xt) & (xq <= xN),
         xq > xN],
        [0.0,
         -th_c,
         lambda x: np.arcsin((x - xt) / d["R_up"]),
         lambda x: np.arcsin((x - xt) / d["R_dn"]),
         bell])


def build_contour(eng, m, cfg=ei):
    """
    Wall profile, injector face -> nozzle exit. Six segments:

        1. cylinder      r = Rc,  length L_cyl
        2. cone          half-angle theta_conv, straight line down toward the throat
        3. upstream arc  radius R_up, centred above the throat, tangent
                         to the cone at one end and horizontal at the throat
        4. downstream arc radius R_down = 0.382 Rt, same centre line, horizontal
                         at the throat, ending at the bell's initial angle theta_n (point N)
        5. bell          quadratic Bezier N -> Q -> E. Leaves N at theta_n and
                         reaches the exit E at theta_e. Q is where the two tangent
                         lines cross.

    Bell (Rao 1958 approximation, as in engine/rao_bell.py):
        theta_n, theta_e   Rao angle tables at eps = (Re/Rt)^2, bell_percent
        E                  x = (bell_percent/100) (sqrt(eps) - 1) Rt / tan 15deg past the throat,
                           r = Re
        Q                  intersection of the line through N at theta_n and
                           the line through E at theta_e

    Coordinates: x = 0 at injector face, x_t at throat, x[-1] at exit plane.
    The cylinder length is chosen so that the volume from injector to throat
    equals L* x At (that's the definition of L*), i.e. the convergent volume is
    subtracted, not ignored.

    Returns x, r (m,), and a dict of the sketch dimensions (also used by
    util/export.write_contour_dxf).
    """
    Rt, Rc, Re = eng["Rt"], eng["Rc"], eng["Re"]
    R_up, R_dn = cfg.R_up_over_Rt * Rt, cfg.R_down_over_Rt * Rt
    th_c = cfg.theta_conv
    eps = (Re / Rt)**2

    # --- bell wall angles from the Rao tables (team code, vendored) ---
    _, th_n, th_e = find_wall_angles(eps, Rt * 1e3, cfg.bell_percent)

    # --- geometry relative to the throat (x_t = 0 for now) ---
    # Upstream arc: centre at (0, Rt + R_up). Point where its tangent is at th_c:
    x1, r1 = -R_up * math.sin(th_c), Rt + R_up * (1 - math.cos(th_c))
    # Cone from (x1, r1) back up to r = Rc at angle th_c:
    L_cone_c = (Rc - r1) / math.tan(th_c)
    x0 = x1 - L_cone_c                       # end of cylinder / start of cone
    # Downstream arc: centre at (0, Rt + R_dn). Ends at N, where its tangent is at th_n:
    xN, rN = R_dn * math.sin(th_n), Rt + R_dn * (1 - math.cos(th_n))
    # Exit point E:
    xE = cfg.bell_percent / 100 * (math.sqrt(eps) - 1) * Rt / math.tan(math.radians(15))
    rE = Re
    # Control point Q: tangent at N meets tangent at E
    m1, m2 = math.tan(th_n), math.tan(th_e)
    C1, C2 = rN - m1 * xN, rE - m2 * xE
    xQ = (C2 - C1) / (m1 - m2)
    rQ = (m1 * C2 - m2 * C1) / (m1 - m2)
    if not xN < xQ < xE:
        raise ValueError(f"Bell control point out of order: N {xN*1e3:.2f}, "
                         f"Q {xQ*1e3:.2f}, E {xE*1e3:.2f} mm from throat")

    # --- convergent volume (cone + upstream arc) by numerical integration ---
    xs = np.linspace(x0, 0, 2000)
    rs = np.where(xs < x1,
                  Rc - (xs - x0) * math.tan(th_c),
                  (Rt + R_up) - np.sqrt(np.clip(R_up**2 - xs**2, 0, None)))
    V_conv = np.trapezoid(math.pi * rs**2, xs)
    L_cyl = (eng["Vc"] - V_conv) / eng["Ac"]
    if L_cyl <= 0:
        raise ValueError(f"L* too small: convergent section alone exceeds L*·At "
                         f"(L_cyl = {L_cyl*1e3:.1f} mm)")

    # --- shift so x = 0 at injector face ---
    x_t = L_cyl - x0
    x0, x1, xN, xQ, xE = x0 + x_t, x1 + x_t, xN + x_t, xQ + x_t, xE + x_t

    def bell_r(xq):
        """r on the Bezier at axial position xq."""
        t = _bezier_t(xq, xN, xQ, xE)
        return (1 - t)**2 * rN + 2 * (1 - t) * t * rQ + t**2 * rE

    # --- sample: uniform in x, nearest station snapped onto the throat exactly ---
    x = np.linspace(0, xE, m)
    x[np.argmin(np.abs(x - x_t))] = x_t
    r = np.piecewise(x,
        [x <= x0,
         (x > x0) & (x <= x1),
         (x > x1) & (x <= x_t),
         (x > x_t) & (x <= xN),
         x > xN],
        [lambda x: Rc,
         lambda x: Rc - (x - x0) * math.tan(th_c),
         lambda x: (Rt + R_up) - np.sqrt(R_up**2 - (x - x_t)**2),
         lambda x: (Rt + R_dn) - np.sqrt(R_dn**2 - (x - x_t)**2),
         bell_r])

    dims = dict(Rc=Rc, Rt=Rt, Re=Re, eps=eps, L_cyl=L_cyl, L_cone_c=L_cone_c,
                R_up=R_up, R_dn=R_dn, x_t=x_t, L_total=xE, L_bell=xE - x_t,
                theta_conv_deg=cfg.theta_conv_deg, bell_percent=cfg.bell_percent,
                theta_n_deg=math.degrees(th_n), theta_e_deg=math.degrees(th_e),
                x0=x0, x1=x1, r1=r1, N=(xN, rN), Q=(xQ, rQ), E=(xE, rE))
    return x, r, dims
