"""
Isopropanol liquid properties vs temperature, for the coolant-side correlations.

Source: `thermo` package (Caleb Bell), Chemical('isopropanol', T, P). CoolProp and
NIST don't carry IPA transport properties, so this is the practical option.

Chemical() costs ~17 ms per call, and the march evaluates properties inside a
root-find at every station. So we tabulate once at the design pressure and
np.interp thereafter (microseconds). Pressure dependence of a liquid's rho, mu,
k, cp over the jacket's few-hundred-psi drop is negligible next to the T effect.

Checked against handbook values at 20 C: rho 786, mu 2.4 mPa-s, k 0.14, cp 2.6 kJ/kg-K.

Valid range: liquid only. Tc = 508 K, so the table stops at 480 K. Lookups
outside [T_MIN, T_MAX] are CLAMPED to the table ends, not raised: the T_wg
root-find in the march probes the ends of its bracket (T_wg = T_l gives T_wl
below the coolant; T_wg = T_r gives T_wl ~ 3000 K), and it only needs the sign
of the residual there. The CONVERGED T_wl is what matters - report() checks
that against T_MAX and warns, since Sieder-Tate isn't valid near/above Tc.
"""
import warnings
import numpy as np

NAME = "IPA"
T_MIN, T_MAX, N_PTS = 270.0, 480.0, 106       # 2 K spacing
LIMIT_NOTE = "Tc = 508 K"


def _build_table(p_Pa):
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        from thermo import Chemical
        T = np.linspace(T_MIN, T_MAX, N_PTS)
        rows = [Chemical("isopropanol", T=Ti, P=p_Pa) for Ti in T]
    return dict(T=T,
                rho=np.array([c.rho for c in rows]),
                mu=np.array([c.mu for c in rows]),
                k=np.array([c.k for c in rows]),
                cp=np.array([c.Cp for c in rows]))


_table = None


def init(p_Pa):
    """Build the property table at the design coolant pressure. Call once before the march."""
    global _table
    _table = _build_table(p_Pa)


def _lookup(name, T):
    if _table is None:
        raise RuntimeError("ipa_props.init(p_Pa) must be called before property lookups")
    # np.interp clamps to the end values outside the table - intentional, see docstring.
    return np.interp(T, _table["T"], _table[name])


def rho(T): return _lookup("rho", T)   # kg/m^3
def mu(T):  return _lookup("mu", T)    # Pa-s
def k(T):   return _lookup("k", T)     # W/m-K
def cp(T):  return _lookup("cp", T)    # J/kg-K


if __name__ == "__main__":
    init(850 * 6894.76)
    print("   T[K]   rho[kg/m3]   mu[mPa-s]   k[W/mK]   cp[J/kgK]")
    for T in (293, 320, 350, 380, 410, 440, 470):
        print(f"{T:6.0f}   {rho(T):9.1f}   {mu(T)*1e3:9.3f}   {k(T):8.4f}   {cp(T):9.1f}")
