"""
Centrifugal impeller generator (CadQuery) - shrouded / unshrouded, splitter blades.

UNITS: geometry in mm, angles in degrees on the user side (radians internally),
       hydraulics in SI (m, m/s, m^3/s, rpm).  Velocities: v = absolute, w = relative, u = blade speed.
Angles: beta is measured FROM THE TANGENTIAL direction in the meridional (m-theta) plane.

Meridional frame: (z, r).  Flow enters axially at z = 0 and leaves radially.
    hub    : (0, r1h)        -> (L,     r2)   hub plane at z = L
    shroud : (0, r1s)        -> (L-b2,  r2)   b2 = outlet passage width
Bezier handles (hub_handles / shroud_handles = (s1, s2), each 0..1):
    P1 = P0 + s1 * (axial extent) along +z        (keeps axial inlet tangent)
    P2 = P3 - s2 * (radial extent) along r         (keeps radial outlet tangent)
Loading: the angular momentum r*v_theta is distributed along the meridional length m in [0,1]
    r*v_theta(m) = rVt1 + f(m) * (rVt2 - rVt1),  f = cumulative of the loading density l(m).
Blade loading dW = W_ss - W_ps = (2*pi/Z) * (v_m / w) * d(r v_theta)/ds   (from the momentum eq.)
"""
from __future__ import annotations
from dataclasses import dataclass, field
from typing import Optional, Callable
import math
import numpy as np

G = 9.80665


# --------------------------------------------------------------------------------------
# Parameters
# --------------------------------------------------------------------------------------
@dataclass
class ImpellerParams:
    # ---- duty ----
    rpm: float
    flow_rate: float                     # [m^3/s] delivered flow
    desired_head: float                  # [m]
    head_is_static: bool = False         # False: head = total head; True: static -> adds (v2^2-v1^2)/2
    rho: float = 998.0                   # [kg/m^3] only for power estimate

    # ---- efficiencies ----
    n_h: float = 0.85                    # hydraulic (Euler work = g*H / n_h)
    n_v: float = 1.0                     # volumetric (impeller flow = Q / n_v; leakage)
    n_m: float = 0.95                    # mechanical (bearings/seals)
    n_p: Optional[float] = None          # overall; default n_h*n_v*n_m

    # ---- inlet velocity triangle ----
    v1_m: Optional[float] = None         # meridional inlet velocity [m/s]; if d1_hub is None -> sets hub dia
    v1_t: float = 0.0                    # swirl at mean (area-weighted) radius [m/s] (from inducer)
    inlet_swirl_model: str = "free_vortex"   # see SWIRL_MODELS
    inlet_swirl_exp: float = -1.0        # exponent for "power_law"
    incidence_deg: float = 0.0           # blade angle - flow angle at LE (decays over first ~30 % of m)
    v2_m: Optional[float] = None         # outlet meridional velocity; default = v1_m
    blockage: float = 0.90               # blade blockage used in continuity

    # ---- slip ----
    slip_factor: Optional[float] = None  # fixed value overrides the model
    slip_model: str = "gulich_wiesner"   # "wiesner" | "gulich_wiesner" | "stodola"
    count_splitters_in_slip: bool = True # Z_eff = Z*(1+n_split) at the outlet

    # ---- geometry [mm] ----
    d1: float = 50.0                     # inlet (eye / shroud) diameter
    d2: float = 100.0                    # outlet diameter
    d1_hub: Optional[float] = None       # hub diameter at inlet (min = bore). None -> from continuity
    axial_length: Optional[float] = None # inlet plane -> outlet hub plane
    blade_thickness: float = 2.0
    n_blades: int = 8
    hub_handles: tuple = (0.5, 0.5)
    shroud_handles: tuple = (0.5, 0.5)

    shrouded: bool = True
    shroud_thickness: float = 2.0
    hub_thickness: float = 4.0           # solid hub extends this far behind the outlet hub plane (z = L .. L+t)
    bore_diameter: float = 12.0
    inlet_lip: float = 0.0               # fillet radius on shroud inlet edge (0 = none)

    # ---- loading ----
    loading_fn: Optional[Callable[[float], float]] = None   # density l(m), m in [0,1]; overrides loading_type
    loading_type: str = "double_linear"  # "uniform" | "linear_aft" | "linear_fore" | "double_linear"
    loading_peak: float = 0.4            # peak location for double_linear
    loading_radius_exp: float = 1.5      # per-streamline density *= (r/r2)^exp: shifts work to larger radius
                                         # (0 = identical r*vt(m) on every streamline -> very high hub angles)

    # ---- splitters ----
    splitters: bool = True
    splitter_start: float = 0.35         # LE position as fraction of meridional length
    num_splitters_per_blade: int = 1
    splitter_pitch_offset: Optional[float] = None   # fraction of pitch of 1st splitter; None = equal spacing

    # ---- other ----
    rotation: str = "ccw"                # "ccw" | "cw"
    n_span: int = 5                      # hub, shroud + (n_span-2) intermediate
    n_chord: int = 36                    # points per blade side per section
    n_fine: int = 400                    # integration grid
    blade_le_m: float = 0.01             # main-blade LE start (m fraction); >0 avoids coplanar faces with the z=0 inlet plane
    trim_eps: float = 0.3                # [mm] pull TE back from the outlet radius
    root_penetration: float = 0.5        # [mm] push blade root into hub (and tip into shroud if shrouded)
    perform_union: bool = True

    # ---- derived (filled by solve) ----
    alpha1: Optional[float] = None
    beta1: Optional[float] = None
    w1: Optional[float] = None
    w1_m: Optional[float] = None
    w1_t: Optional[float] = None
    v1: Optional[float] = None
    v1_components: Optional[tuple] = None   # (v_t, v_ax, v_r) at mean radius
    v2: Optional[float] = None
    v2_u: Optional[float] = None
    alpha2: Optional[float] = None
    beta2: Optional[float] = None
    w2: Optional[float] = None
    w2_m: Optional[float] = None
    w2_t: Optional[float] = None
    v2_components: Optional[tuple] = None


# --------------------------------------------------------------------------------------
# Bezier helpers (P: (4,2) array of (z, r))
# --------------------------------------------------------------------------------------
def _bez(P, t):
    t = np.asarray(t, float)[:, None]
    return ((1 - t) ** 3) * P[0] + 3 * ((1 - t) ** 2) * t * P[1] \
        + 3 * (1 - t) * t ** 2 * P[2] + t ** 3 * P[3]


def _dbez(P, t):
    t = np.asarray(t, float)[:, None]
    return 3 * ((1 - t) ** 2 * (P[1] - P[0]) + 2 * (1 - t) * t * (P[2] - P[1])
                + t ** 2 * (P[3] - P[2]))


def _cumtrapz(y, x):
    return np.concatenate([[0.0], np.cumsum(0.5 * (y[1:] + y[:-1]) * np.diff(x))])


def _resample_by_arclength(P, n):
    """Return points, unit tangents (z,r) and total length at n points uniform in arclength."""
    t = np.linspace(0, 1, 4000)
    pts = _bez(P, t)
    s = np.concatenate([[0], np.cumsum(np.hypot(*np.diff(pts, axis=0).T))])
    tm = np.interp(np.linspace(0, 1, n) * s[-1], s, t)
    d = _dbez(P, tm)
    d /= np.linalg.norm(d, axis=1, keepdims=True)
    return _bez(P, tm), d, s[-1]


# --------------------------------------------------------------------------------------
# Inlet radial equilibrium models:  v_theta(r) given value at mean radius
# --------------------------------------------------------------------------------------
def _swirl_models():
    return {
        "free_vortex": lambda x, vm, n: vm / x,                    # r*vt = const
        "solid_body":  lambda x, vm, n: vm * x,                    # vt ~ r
        "constant":    lambda x, vm, n: vm + 0 * x,                # vt = const (= const-alpha for uniform vm)
        "compound":    lambda x, vm, n: vm * 0.5 * (x + 1 / x),    # 50/50 free + forced vortex
        "power_law":   lambda x, vm, n: vm * x ** n,               # vt ~ r^n (n=-1 free, +1 solid, 0 const)
    }


SWIRL_MODELS = _swirl_models()   # add your own: SWIRL_MODELS["mine"] = lambda x, vm, n: ...


def inlet_swirl(r, rm, vu_m, model, n=-1.0):
    x = np.maximum(np.asarray(r, float), 1e-6) / rm
    return SWIRL_MODELS[model](x, vu_m, n)


# --------------------------------------------------------------------------------------
# Slip
# --------------------------------------------------------------------------------------
def slip_factor_model(beta2, z, d1m, d2, nq, model="gulich_wiesner"):
    sb = math.sin(beta2)
    wies = 1.0 - math.sqrt(sb) / z ** 0.7
    if model == "wiesner":
        return wies
    if model == "stodola":
        return 1.0 - math.pi * sb / z
    f1 = max(0.98, 1.02 + 1.02e-3 * (nq - 50.0))
    eps_lim = math.exp(-8.16 * sb / z)
    ratio = d1m / d2
    kw = 1.0 if ratio <= eps_lim else 1.0 - ((ratio - eps_lim) / (1.0 - eps_lim)) ** 3
    return f1 * wies * kw


# --------------------------------------------------------------------------------------
# Loading
# --------------------------------------------------------------------------------------
def loading_distribution(p: ImpellerParams, m):
    """Returns (density l(m) normalised to unit integral, cumulative f(m) 0..1)."""
    if p.loading_fn is not None:
        l = np.array([p.loading_fn(float(x)) for x in m], float)
    elif p.loading_type == "uniform":
        l = np.ones_like(m)
    elif p.loading_type == "linear_aft":
        l = m.copy()
    elif p.loading_type == "linear_fore":
        l = 1 - m
    elif p.loading_type == "double_linear":
        pk = min(max(p.loading_peak, 1e-3), 1 - 1e-3)
        l = np.where(m < pk, m / pk, (1 - m) / (1 - pk))
    else:
        raise ValueError(f"unknown loading_type {p.loading_type}")
    l = np.clip(l, 0, None) + 1e-4              # tiny floor keeps f strictly monotone
    f = _cumtrapz(l, m)
    l = l / f[-1]
    f = f / f[-1]
    return l, f


# --------------------------------------------------------------------------------------
# Design result
# --------------------------------------------------------------------------------------
@dataclass
class Design:
    p: ImpellerParams
    m: np.ndarray = None
    lam: np.ndarray = None
    hub: np.ndarray = None            # (n,2) z,r
    shroud: np.ndarray = None
    hub_t: np.ndarray = None
    shroud_t: np.ndarray = None
    stream: list = field(default_factory=list)       # per span: dict of arrays
    load_density: np.ndarray = None
    load_cum: np.ndarray = None
    info: dict = field(default_factory=dict)


def solve_design(p: ImpellerParams) -> Design:
    n = p.n_fine
    omega = p.rpm * 2 * math.pi / 60
    Q = p.flow_rate / p.n_v
    r1s = p.d1 / 2e3                                        # m
    r2 = p.d2 / 2e3
    u2 = omega * r2
    B = p.blockage

    # ---------- inlet hub from continuity ----------
    if p.d1_hub is None:
        if p.v1_m is not None:
            dh2 = (p.d1 / 1e3) ** 2 - 4 * Q / (math.pi * p.v1_m * B)
            if dh2 < 0:
                raise ValueError("v1_m too low for this d1 (hub dia^2 < 0): increase d1 or v1_m")
            d1h = math.sqrt(dh2) * 1e3
        else:
            d1h = 0.3 * p.d1
    else:
        d1h = p.d1_hub
    d1h = max(d1h, p.bore_diameter)                         # must fit shaft
    if d1h >= p.d1:
        raise ValueError("hub diameter >= eye diameter; increase d1")
    p.d1_hub = d1h
    r1h = d1h / 2e3
    A1 = math.pi / 4 * ((2 * r1s) ** 2 - d1h ** 2 / 1e6) * B
    v1_m = Q / A1
    p.v1_m = v1_m
    rm1 = math.sqrt(0.5 * (r1h ** 2 + r1s ** 2))            # area-weighted mean radius
    u1m = omega * rm1
    vt1m = p.v1_t

    # ---------- outlet width from continuity ----------
    v2_m = p.v2_m if p.v2_m is not None else v1_m
    b2 = Q / (2 * math.pi * r2 * v2_m * B) * 1e3            # mm
    p.v2_m = v2_m

    # ---------- Euler + slip -> outlet triangle ----------
    nq = p.rpm * math.sqrt(p.flow_rate) / p.desired_head ** 0.75
    n_split = p.num_splitters_per_blade if p.splitters else 0
    z_eff = p.n_blades * (1 + n_split) if p.count_splitters_in_slip else p.n_blades
    v1_abs = math.hypot(v1_m, vt1m)
    beta2 = math.radians(25.0)
    v2_abs = v2_m
    gamma = 0.9
    for _ in range(200):
        gH = G * p.desired_head
        if p.head_is_static:
            gH += 0.5 * (v2_abs ** 2 - v1_abs ** 2)
        work = gH / p.n_h                                   # Euler work u2 vt2 - u1 vt1
        vt2 = (work + u1m * vt1m) / u2
        if p.slip_factor is not None:
            gamma_new = p.slip_factor
            gh = gs = gamma_new
        else:
            gh = slip_factor_model(beta2, z_eff, 2 * r1h, 2 * r2, nq, p.slip_model)
            gs = slip_factor_model(beta2, z_eff, 2 * r1s, 2 * r2, nq, p.slip_model)
            gamma_new = 0.5 * (gh + gs)
        arg = gamma_new * u2 - vt2                          # = v_m2 * cot(beta2)
        if arg <= 0:
            raise ValueError("Head not reachable with this slip/rpm/d2 (beta2 >= 90 deg). "
                             "Increase rpm or d2, or lower the head.")
        beta2_new = math.atan2(v2_m, arg)
        v2_new = math.hypot(v2_m, vt2)
        done = abs(beta2_new - beta2) < 1e-9 and abs(gamma_new - gamma) < 1e-9 and abs(v2_new - v2_abs) < 1e-9
        beta2 = 0.5 * beta2 + 0.5 * beta2_new
        gamma, v2_abs = gamma_new, v2_new
        if done:
            break
    beta2 = beta2_new

    n_p = p.n_p if p.n_p is not None else p.n_h * p.n_v * p.n_m
    p.n_p = n_p
    shaft_power = p.rho * G * p.flow_rate * p.desired_head / n_p

    # ---------- fill triangle fields ----------
    p.v1, p.alpha1 = v1_abs, math.degrees(math.atan2(v1_m, vt1m))
    p.w1_m, p.w1_t = v1_m, vt1m - u1m
    p.w1 = math.hypot(p.w1_m, p.w1_t)
    p.beta1 = math.degrees(math.atan2(v1_m, u1m - vt1m))
    p.v1_components = (vt1m, v1_m, 0.0)
    p.v2, p.v2_u = v2_abs, vt2
    p.alpha2 = math.degrees(math.atan2(v2_m, vt2))
    p.beta2 = math.degrees(beta2)
    p.w2_m, p.w2_t = v2_m, vt2 - u2
    p.w2 = math.hypot(p.w2_m, p.w2_t)
    p.v2_components = (vt2, 0.0, v2_m)

    # ---------- meridional geometry ----------
    r1h_mm, r1s_mm, r2_mm = r1h * 1e3, r1s * 1e3, r2 * 1e3
    L = p.axial_length if p.axial_length is not None else b2 + 0.6 * (r2_mm - r1s_mm)
    p.axial_length = L
    zs_end = L - b2
    if zs_end <= 0:
        raise ValueError("axial_length must exceed outlet width b2")
    hs, ss = p.hub_handles, p.shroud_handles
    Ph = np.array([[0, r1h_mm], [hs[0] * L, r1h_mm], [L, r2_mm - hs[1] * (r2_mm - r1h_mm)], [L, r2_mm]], float)
    Ps = np.array([[0, r1s_mm], [ss[0] * zs_end, r1s_mm], [zs_end, r2_mm - ss[1] * (r2_mm - r1s_mm)],
                   [zs_end, r2_mm]], float)
    hub, hub_t, _ = _resample_by_arclength(Ph, n)
    shr, shr_t, _ = _resample_by_arclength(Ps, n)
    m = np.linspace(0, 1, n)

    # passage area & meridional velocity (continuity with blockage)
    b_loc = np.linalg.norm(shr - hub, axis=1)
    r_mean = 0.5 * (hub[:, 1] + shr[:, 1])
    A = 2 * math.pi * r_mean * b_loc * B / 1e6              # m^2
    vm = Q / A

    # ---------- loading ----------
    dens, f = loading_distribution(p, m)
    lam = np.linspace(0, 1, p.n_span)
    rVt2 = r2 * vt2
    streams = []
    for lk in lam:
        pts = hub + lk * (shr - hub)
        zz, rr = pts[:, 0], pts[:, 1]
        s = np.concatenate([[0], np.cumsum(np.hypot(np.diff(zz), np.diff(rr)))])      # mm
        r_si = rr / 1e3
        vt_in = float(inlet_swirl(rr[0] / 1e3, rm1, vt1m, p.inlet_swirl_model, p.inlet_swirl_exp))
        rVt1 = rr[0] / 1e3 * vt_in
        dens_k = dens * (rr / r2_mm) ** p.loading_radius_exp
        f_k = _cumtrapz(dens_k, m)
        norm_k = f_k[-1]
        f_k = f_k / norm_k
        rVt = rVt1 + f_k * (rVt2 - rVt1)
        vt = rVt / np.maximum(r_si, 1e-6)
        u = omega * r_si
        beta_flow = np.arctan2(vm, u - vt)
        inc = np.radians(p.incidence_deg) * np.clip(1 - m / 0.3, 0, 1) ** 2
        beta_b = beta_flow + inc
        w = np.hypot(vm, u - vt)
        drVt_ds = np.gradient(rVt, s / 1e3)
        dW = (2 * math.pi / z_eff) * (vm / w) * drVt_ds
        dW_full = (2 * math.pi / p.n_blades) * (vm / w) * drVt_ds     # loading of the main blades only
        sgn = 1.0 if p.rotation.lower() == "ccw" else -1.0
        # wrap angle: r dtheta = (w_theta / w_m) ds , blade lags the rotation
        theta = -sgn * _cumtrapz(1.0 / np.tan(beta_b) / rr, s)
        streams.append(dict(lam=lk, z=zz, r=rr, s=s, vm=vm, vt=vt, u=u, w=w, beta_flow=beta_flow,
                            beta_b=beta_b, theta=theta, dW=dW, dW_full=dW_full, rVt=rVt, f=f_k, dens=dens_k / norm_k))

    d = Design(p=p, m=m, lam=lam, hub=hub, shroud=shr, hub_t=hub_t, shroud_t=shr_t, stream=streams,
               load_density=dens, load_cum=f)
    d.info = dict(omega=omega, u2=u2, nq=nq, b2_mm=b2, L_mm=L, d1_hub_mm=d1h, v1_m=v1_m, v2_m=v2_m,
                  gamma=gamma, gamma_hub=gh, gamma_shroud=gs, euler_work=work, vt2=vt2,
                  beta1_deg=p.beta1, beta2_deg=p.beta2, shaft_power_W=shaft_power,
                  rm1_mm=rm1 * 1e3, z_eff=z_eff, n_p=n_p)
    return d


# --------------------------------------------------------------------------------------
# Blade section sampling
# --------------------------------------------------------------------------------------
def _blade_sections(d: Design, m_start: float, m_end: float):
    """Per span station: dict(ps, ss, camber) arrays of (x,y,z) [mm], LE->TE."""
    p = d.p
    u = np.linspace(0, 1, p.n_chord)
    ms = m_start + (m_end - m_start) * 0.5 * (1 - np.cos(math.pi * u))     # LE/TE clustering
    out = []
    for k, st in enumerate(d.stream):
        z = np.interp(ms, d.m, st["z"])
        r = np.interp(ms, d.m, st["r"])
        th = np.interp(ms, d.m, st["theta"])
        bb = np.interp(ms, d.m, st["beta_b"])
        s = np.interp(ms, d.m, st["s"])
        s_loc = s - s[0]
        le_len = max(0.04 * s_loc[-1], 1e-6)
        x = np.clip(s_loc / le_len, 0, 1)
        tfac = np.maximum(np.sqrt(1 - (1 - x) ** 2), 0.12)                 # elliptic LE, blunt TE
        t_tan = p.blade_thickness * tfac / np.maximum(np.abs(np.sin(bb)), 0.25)
        dth = 0.5 * t_tan / r
        sgn = 1.0 if p.rotation.lower() == "ccw" else -1.0
        th_ps, th_ss = th + sgn * dth, th - sgn * dth                       # PS leads in the rotation direction

        # push root / tip into hub / shroud (meridional normal)
        zr = np.stack([z, r], 1)
        if k == 0 or (k == len(d.stream) - 1 and p.shrouded):
            if k == 0:
                tan = d.hub_t
                nrm = np.stack([tan[:, 1], -tan[:, 0]], 1)                  # away from passage (hub)
                src = d.hub
            else:
                tan = d.shroud_t
                nrm = np.stack([-tan[:, 1], tan[:, 0]], 1)                  # away from passage (shroud)
                src = d.shroud
            nn = np.stack([np.interp(ms, d.m, nrm[:, 0]), np.interp(ms, d.m, nrm[:, 1])], 1)
            zr = zr + nn * p.root_penetration
        z2, r2_ = zr[:, 0], zr[:, 1]

        def xyz(th_):
            return np.stack([r2_ * np.cos(th_), r2_ * np.sin(th_), z2], 1)
        out.append(dict(ps=xyz(th_ps), ss=xyz(th_ss), camber=xyz(th), z=z2, r=r2_, theta=th))
    return out


def _trim_end(d: Design):
    p = d.p
    rmean = 0.5 * (d.hub[:, 1] + d.shroud[:, 1])
    r_target = p.d2 / 2 - p.trim_eps
    return float(np.interp(r_target, rmean, d.m))


# --------------------------------------------------------------------------------------
# CAD
# --------------------------------------------------------------------------------------
def build_cad(d: Design):
    import cadquery as cq
    p = d.p
    L = p.axial_length
    r_b = p.bore_diameter / 2

    # ---------- hub (solid of revolution; flat back at z = L + hub_thickness) ----------
    idx = np.linspace(0, len(d.hub) - 1, 80).astype(int)
    hp = d.hub[idx]
    prof = (cq.Workplane("XZ").moveTo(r_b, 0).lineTo(hp[0, 1], 0)
            .spline([(r, z) for z, r in hp[1:]], includeCurrent=True)
            .lineTo(hp[-1, 1], L + p.hub_thickness).lineTo(r_b, L + p.hub_thickness).close())
    hub = prof.revolve(360, (0, 0, 0), (0, 1, 0))

    # ---------- shroud (shell by normal offset) ----------
    shroud = None
    if p.shrouded:
        sp = d.shroud[idx]
        tg = d.shroud_t[idx]
        nrm = np.stack([-tg[:, 1], tg[:, 0]], 1)
        so = sp + nrm * p.shroud_thickness
        w = (cq.Workplane("XZ").moveTo(sp[0, 1], sp[0, 0])
             .spline([(r, z) for z, r in sp[1:]], includeCurrent=True)
             .lineTo(so[-1, 1], so[-1, 0])
             .spline([(r, z) for z, r in so[::-1][1:]], includeCurrent=True).close())
        shroud = w.revolve(360, (0, 0, 0), (0, 1, 0))
        if p.inlet_lip > 0:
            try:
                shroud = shroud.faces("<Z").edges().fillet(min(p.inlet_lip, 0.45 * p.shroud_thickness))
            except Exception as e:                     # fillet is cosmetic - never fail the build
                print("inlet_lip fillet skipped:", e)

    # ---------- blades ----------
    def loft(secs):
        """Loft the sections, cap both ends (non-planar -> ruled PS/SS surface), sew into a closed solid."""
        from OCP.BRepBuilderAPI import BRepBuilderAPI_Sewing, BRepBuilderAPI_MakeSolid
        from OCP.TopoDS import TopoDS
        wires, ends = [], []
        for s in secs:
            ps = [cq.Vector(*q) for q in s["ps"]]
            ss = [cq.Vector(*q) for q in s["ss"]]
            pe, se = cq.Edge.makeSpline(ps), cq.Edge.makeSpline(ss)
            ends.append((pe, se))
            e = [pe, cq.Edge.makeLine(ps[-1], ss[-1]), cq.Edge.makeSpline(ss[::-1]), cq.Edge.makeLine(ss[0], ps[0])]
            wires.append(cq.Wire.assembleEdges(e))
        sides = cq.Solid.makeLoft(wires, False)
        caps = [cq.Face.makeRuledSurface(*ends[0]), cq.Face.makeRuledSurface(*ends[-1])]
        sew = BRepBuilderAPI_Sewing(1e-3)
        for f in list(sides.Faces()) + caps:
            sew.Add(f.wrapped)
        sew.Perform()
        shell = cq.Shape.cast(sew.SewedShape()).Shells()[0]
        sol = cq.Solid(BRepBuilderAPI_MakeSolid(TopoDS.Shell_s(shell.wrapped)).Solid())
        if sol.Volume() < 0:
            sol = sol.fix()
        return sol

    m_end = _trim_end(d)
    pitch = 360.0 / p.n_blades
    sgn = 1.0 if p.rotation.lower() == "ccw" else -1.0
    main = loft(_blade_sections(d, p.blade_le_m, m_end))
    blades = [main.rotate(cq.Vector(0, 0, 0), cq.Vector(0, 0, 1), i * pitch) for i in range(p.n_blades)]

    split = []
    if p.splitters and p.num_splitters_per_blade > 0:
        spl0 = loft(_blade_sections(d, p.splitter_start, m_end))
        N = p.num_splitters_per_blade
        off0 = p.splitter_pitch_offset if p.splitter_pitch_offset is not None else 1.0 / (N + 1)
        for i in range(p.n_blades):
            for j in range(N):
                frac = off0 + j * (1 - off0) / N
                ang = (i + frac) * pitch
                split.append(spl0.rotate(cq.Vector(0, 0, 0), cq.Vector(0, 0, 1), ang))

    parts = dict(hub=hub, shroud=shroud, blades=blades, splitters=split)
    all_blades = cq.Compound.makeCompound(blades + split)
    if p.perform_union:
        res = hub.union(all_blades)
        if shroud is not None:
            res = res.union(shroud)
        parts["impeller"] = res
    else:
        objs = [hub.val(), all_blades] + ([shroud.val()] if shroud is not None else [])
        parts["impeller"] = cq.Workplane("XY").newObject([cq.Compound.makeCompound(objs)])
    return parts


# --------------------------------------------------------------------------------------
# Plots
# --------------------------------------------------------------------------------------
def make_plots(d: Design, show=True, save_prefix: Optional[str] = None):
    import matplotlib.pyplot as plt
    p = d.p
    m = d.m
    cols = plt.cm.viridis(np.linspace(0, 0.9, len(d.stream)))
    fig, ax = plt.subplots(2, 3, figsize=(17, 9))

    a = ax[0, 0]                                           # meridional contours
    a.plot(d.hub[:, 0], d.hub[:, 1], "k", lw=2, label="hub")
    a.plot(d.shroud[:, 0], d.shroud[:, 1], "r", lw=2, label="shroud")
    for st, c in zip(d.stream[1:-1], cols[1:-1]):
        a.plot(st["z"], st["r"], color=c, lw=1, ls="--")
    a.set(title="Meridional contours", xlabel="z [mm]", ylabel="r [mm]")
    a.set_aspect("equal"); a.legend(); a.grid(alpha=.3)

    a = ax[0, 1]                                           # loading
    a.plot(m, d.load_density, "k", label="base density l(m)")
    a.plot(m, d.load_cum, "k--", label="base cumulative")
    for st, c in zip(d.stream, cols):
        a.plot(m, st["dens"], color=c, lw=0.8, ls=":")
    a.set(title="Loading distribution", xlabel="m (meridional fraction)")
    a.grid(alpha=.3)
    a2 = a.twinx()
    for st, c in zip(d.stream, cols):
        a2.plot(m, st["dW_full"], color=c, lw=1)
    a2.set_ylabel("W_ss - W_ps [m/s] (hub->shroud)")
    a.legend(loc="upper right")

    a = ax[0, 2]                                           # blade angle
    for st, c in zip(d.stream, cols):
        a.plot(m, np.degrees(st["beta_b"]), color=c, label=f"span {st['lam']:.2f}")
    a.set(title="Blade angle from tangential", xlabel="m", ylabel="beta_b [deg]")
    a.legend(); a.grid(alpha=.3)

    a = ax[1, 0]                                           # relative velocity
    for st, c in zip(d.stream, cols):
        a.plot(m, st["w"], color=c, label=f"span {st['lam']:.2f}")
    mid = d.stream[len(d.stream) // 2]
    a.plot(m, mid["w"] + 0.5 * mid["dW_full"], "k:", lw=1, label="mid: suction")
    a.plot(m, mid["w"] - 0.5 * mid["dW_full"], "k-.", lw=1, label="mid: pressure")
    a.set(title="Relative velocity", xlabel="m", ylabel="w [m/s]")
    a.legend(fontsize=7); a.grid(alpha=.3)

    a = ax[1, 1]                                           # meridional velocity & rVt
    a.plot(m, d.stream[0]["vm"], "b", label="v_m")
    a.set_xlabel("m"); a.set_ylabel("v_m [m/s]"); a.grid(alpha=.3)
    a3 = a.twinx()
    for st, c in zip(d.stream, cols):
        a3.plot(m, st["rVt"], color=c, lw=1)
    a3.set_ylabel("r*v_theta [m^2/s]")
    a.set_title("Continuity v_m and r*v_theta"); a.legend(loc="upper left")

    a = ax[1, 2]                                           # plan projection
    th_ring = np.linspace(0, 2 * math.pi, 400)
    for rr in (p.d1_hub / 2, p.d1 / 2, p.d2 / 2):
        a.plot(rr * np.cos(th_ring), rr * np.sin(th_ring), color="gray", lw=.6)
    m_end = _trim_end(d)
    pitch = 2 * math.pi / p.n_blades
    sgn = 1.0 if p.rotation.lower() == "ccw" else -1.0
    sel = [(0, "k", "hub"), (len(d.stream) - 1, "r", "shroud")]
    for k, c, nm in sel:
        for i in range(p.n_blades):
            sec = _blade_sections(d, p.blade_le_m, m_end)[k] if i == 0 else sec
            x, y = sec["camber"][:, 0], sec["camber"][:, 1]
            cs, sn = math.cos(i * pitch), math.sin(i * pitch)
            a.plot(x * cs - y * sn, x * sn + y * cs, color=c, lw=1.3, label=nm if i == 0 else None)
        if p.splitters and p.num_splitters_per_blade > 0:
            sps = _blade_sections(d, p.splitter_start, m_end)[k]
            N = p.num_splitters_per_blade
            off0 = p.splitter_pitch_offset if p.splitter_pitch_offset is not None else 1.0 / (N + 1)
            for i in range(p.n_blades):
                for j in range(N):
                    ang = (i + off0 + j * (1 - off0) / N) * pitch
                    cs, sn = math.cos(ang), math.sin(ang)
                    x, y = sps["camber"][:, 0], sps["camber"][:, 1]
                    a.plot(x * cs - y * sn, x * sn + y * cs, color=c, lw=1, ls="--")
    a.set_aspect("equal"); a.set(title="Plan view: hub (k) / shroud (r) camber lines, splitters dashed",
                                 xlabel="x [mm]", ylabel="y [mm]")
    a.legend(); a.grid(alpha=.3)
    fig.tight_layout()
    if save_prefix:
        fig.savefig(f"{save_prefix}_plots.png", dpi=150)
    if show:
        plt.show()
    return fig


def print_summary(d: Design):
    p, i = d.p, d.info
    print("---- design summary ----")
    print(f"n_q={i['nq']:.1f}  u2={i['u2']:.1f} m/s  Euler work={i['euler_work']:.0f} J/kg  eff_overall={i['n_p']:.3f}")
    print(f"hub dia (inlet)={i['d1_hub_mm']:.2f} mm  b2={i['b2_mm']:.2f} mm  L={i['L_mm']:.2f} mm  rm1={i['rm1_mm']:.2f} mm")
    print(f"v1_m={i['v1_m']:.2f}  v2_m={i['v2_m']:.2f}  vt2={i['vt2']:.2f} m/s")
    print(f"slip: gamma={i['gamma']:.3f} (hub {i['gamma_hub']:.3f}, shroud {i['gamma_shroud']:.3f}), Z_eff={i['z_eff']}")
    print(f"beta1(flow, mean)={p.beta1:.1f} deg  beta2={p.beta2:.1f} deg  alpha1={p.alpha1:.1f}  alpha2={p.alpha2:.1f}")
    print(f"shaft power ~ {i['shaft_power_W'] / 1e3:.2f} kW")


# --------------------------------------------------------------------------------------
if __name__ == "__main__":
    params = ImpellerParams(
        rpm=10000, flow_rate=0.008, desired_head=100.0,
        d1=50.0, d2=100.0, n_blades=6, splitters=True, num_splitters_per_blade=1,
        shrouded=True, hub_handles=(0.5, 0.5), shroud_handles=(0.5, 0.5),
        loading_type="double_linear", loading_peak=0.4, v1_t=0.0, perform_union=True,
    )
    design = solve_design(params)
    print_summary(design)
    make_plots(design, show=False, save_prefix="impeller")
    try:
        import cadquery as cq
        parts = build_cad(design)
        cq.exporters.export(parts["impeller"], "impeller.step")
        print("wrote impeller.step")
    except ImportError:
        print("cadquery not installed - skipped CAD")