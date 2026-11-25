"""
CEA interface for the regen analysis. Two jobs:

    size_engine()      - throat / chamber / exit dimensions from c* and Cf
    gas_side_states()  - hot-gas state at an arbitrary list of area ratios,
                         ONE CEA call, returned as arrays for the regen march

Both take the engine definition as `cfg` (default engine/engine_inputs.py). The CEA
mixture/solver is built once per propellant pair and cached.

Sanity check: scripts/engine_sizing.py.

CEA output conventions (checked against the library):
    station order   [chamber, throat, subar..., supar...]  (pi_p not used here)
    P               bar
    T               K
    viscosity       millipoise      -> x1e-4  = Pa-s
    conductivity    mW/cm-K         -> x0.1   = W/m-K
    cp              kJ/kg-K         -> x1e3   = J/kg-K
"""
import math
import numpy as np
import cea

from engine import engine_inputs as ei

# ============================== CEA SETUP ==================================
# Mixture + solver built once per propellant pair; reusable across solve() calls.
_cache = {}


def _mixture(cfg):
    key = tuple(cfg.reac_names)
    if key not in _cache:
        reac = cea.Mixture(list(cfg.reac_names))
        prod = cea.Mixture(list(cfg.reac_names), products_from_reactants=True)
        _cache[key] = (reac, cea.RocketSolver(prod, reactants=reac, transport=True))
    return _cache[key]


def _solve(cfg, **stations):
    """Run CEA at chamber pressure and cfg.OF with the requested station lists."""
    reac, solver = _mixture(cfg)
    weights = reac.of_ratio_to_weights(cfg.oxidant_weights, cfg.fuel_weights, cfg.OF)
    h = reac.calc_property(cea.ENTHALPY, weights, cfg.T_reactant)   # J/kg of mixture
    h -= cfg.fuel_dH_vap / (1 + cfg.OF)                              # liquid fuel (see engine_inputs)
    hc = h / cea.R
    sol = cea.RocketSolution(solver)
    solver.solve(sol, weights, cfg.p0_bar, hc=hc, **stations)
    if not sol.converged:
        raise RuntimeError("CEA did not converge")
    return sol


# ============================== ENGINE SIZING ==============================

def size_engine(cfg=ei):
    """
    Mass flow (pump-limited) -> throat from c* -> exit from Pe = Pa -> chamber from
    contraction ratio and L*. Returns a dict of SI values.
    """
    # Pump limit: ox = OF * fuel, so ox limits when OF > 1
    if cfg.OF >= 1:
        m_dot_ox = cfg.m_dot_pump_max
        m_dot_fuel = m_dot_ox / cfg.OF
    else:
        m_dot_fuel = cfg.m_dot_pump_max
        m_dot_ox = m_dot_fuel * cfg.OF
    m_dot = m_dot_ox + m_dot_fuel

    # Stations: [chamber, throat, exit (pi_p), chamber section (subar)]
    sol = _solve(cfg, pi_p=[cfg.p0_bar / cfg.pa_bar], subar=[cfg.contraction_ratio])
    CH, TH, EX = 0, 1, 2

    c_star = sol.c_star[TH] * cfg.c_star_eff
    Cf = sol.coefficient_of_thrust[EX] * cfg.c_thrust_eff
    eps = sol.ae_at[EX]

    At = m_dot * c_star / cfg.p0_Pa
    Ae = eps * At
    Ac = cfg.contraction_ratio * At
    Vc = cfg.lstar * At                     # volume injector -> throat

    F_N = Cf * cfg.p0_Pa * At
    return dict(
        m_dot_ox=m_dot_ox, m_dot_fuel=m_dot_fuel, m_dot=m_dot,
        c_star=c_star, Cf=Cf, eps=eps,
        At=At, Rt=math.sqrt(At / math.pi),
        Ae=Ae, Re=math.sqrt(Ae / math.pi),
        Ac=Ac, Rc=math.sqrt(Ac / math.pi),
        Vc=Vc,
        F_N=F_N, F_lbf=F_N * cfg.N_TO_LBF, Isp=F_N / (m_dot * cfg.G0),
        T0=sol.T[CH],
    )


# ============================== GAS-SIDE STATES ============================

def gas_side_states(area_ratio, supersonic, cfg=ei):
    """
    Hot-gas state at every station of the contour.

    Parameters
        area_ratio : (m,) A/At at each station, indexed injector -> exit
        supersonic : (m,) bool, True downstream of the throat

    CEA needs subsonic (subar) and supersonic (supar) area ratios as separate
    lists, and returns them in the order given after [chamber, throat]. We send
    them in contour order, then scatter back into contour-indexed arrays. The
    throat station (A/At == 1) is taken from CEA's own throat result.

    Returns dict of (m,) arrays in SI:
        T, p, M, gamma, mu, cp, Pr
    plus scalars at stagnation for Bartz: T0, p0, mu0, cp0, Pr0
    """
    area_ratio = np.asarray(area_ratio, dtype=float)
    supersonic = np.asarray(supersonic, dtype=bool)
    m = len(area_ratio)

    is_throat = np.isclose(area_ratio, 1.0)
    sub_idx = np.where(~supersonic & ~is_throat)[0]
    sup_idx = np.where(supersonic & ~is_throat)[0]

    sol = _solve(cfg, subar=list(area_ratio[sub_idx]), supar=list(area_ratio[sup_idx]))

    # Map CEA station -> contour index. CEA order: [0 chamber, 1 throat, sub..., sup...]
    cea_idx = np.empty(m, dtype=int)
    cea_idx[is_throat] = 1
    cea_idx[sub_idx] = 2 + np.arange(len(sub_idx))
    cea_idx[sup_idx] = 2 + len(sub_idx) + np.arange(len(sup_idx))

    def pull(arr, scale=1.0):
        return np.asarray(arr[:sol.num_pts], dtype=float)[cea_idx] * scale

    # Frozen transport properties: Bartz was derived for a non-reacting boundary
    # layer, and frozen cp/Pr are what Bartz (1957) and Huzel & Huang use.
    return dict(
        T=pull(sol.T),
        p=pull(sol.P, 1e5),
        M=pull(sol.Mach),
        gamma=pull(sol.gamma_s),
        mu=pull(sol.viscosity, 1e-4),
        cp=pull(sol.cp_fr, 1e3),
        Pr=pull(sol.Pr_fr),
        T0=sol.T[0],
        p0=sol.P[0] * 1e5,
        mu0=sol.viscosity[0] * 1e-4,
        cp0=sol.cp_fr[0] * 1e3,
        Pr0=sol.Pr_fr[0]
    )
