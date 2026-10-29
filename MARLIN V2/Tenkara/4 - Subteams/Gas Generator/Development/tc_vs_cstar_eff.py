"""
Combustion temperature vs. c* efficiency for LOX/IPA using RocketCEA.

Model: "reduced energy release" is imposed as an enthalpy deficit on the reactants.
CEA is then run as a normal HP (constant-enthalpy, constant-pressure) problem, so
the product composition, MW, gamma and c* are all re-equilibrated at the lower
temperature (self-consistent, no constant-R / constant-gamma assumption).

    eta_c* = c*_CEA(deficit) / c*_CEA(0)

The curve is built parametrically (sweep the deficit, compute eta_c*), so no root
finding or registration of hundreds of fuels is needed.

Units used by RocketCEA: Pc [psia], T [Rankine], c* [ft/s], MW [lbm/lbmol].
Everything is converted to SI inside this script.
"""
import itertools

import matplotlib.pyplot as plt
import numpy as np
from rocketcea.cea_obj import CEA_Obj, add_new_fuel

# --------------------------------------------------------------------------
# USER INPUTS
# --------------------------------------------------------------------------
PC_PSIA = 330.0        # chamber pressure [psia]
OF = 0.5               # oxidizer/fuel mass ratio
EPS = 13.0             # nozzle area ratio (only needed as an argument; chamber
                       # properties do not depend on it)

# Isopropanol (liquid) data. CHECK these against your own source / thermo.inp.
IPA_FORMULA = "C 3 H 8 O 1"
IPA_MW_KG_MOL = 0.060095          # [kg/mol]
IPA_HF_KJ_MOL = -318.1            # standard enthalpy of formation, liquid [kJ/mol]
IPA_T_K = 298.15                  # fuel inlet temperature [K]

D_MAX_KJ_KG = 4200.0   # largest enthalpy deficit swept [kJ per kg of propellant]
N_POINTS = 60
ETA_MIN_PLOT = 0.4    # drop points below this c* efficiency

# Optional: if you define an "ideal heat release" Q_ideal [kJ/kg propellant]
# you get eta_Q = 1 - deficit/Q_ideal in the output table. Left as None because
# Q_ideal has no unique definition (reference temperature, product set...).
Q_IDEAL_KJ_KG = None

# --------------------------------------------------------------------------
# CONSTANTS
# --------------------------------------------------------------------------
R_UNIV = 8314.462          # [J/(kmol K)]
FT_TO_M = 0.3048
RANKINE_TO_K = 1.0 / 1.8
KJ_TO_CAL = 1000.0 / 4.184

_counter = itertools.count()


# --------------------------------------------------------------------------
# CEA HELPERS
# --------------------------------------------------------------------------
def make_cea_obj(deficit_kj_per_kg_prop: float) -> CEA_Obj:
    """LOX + IPA whose enthalpy is lowered so the mixture releases less energy.

    The deficit is given per kg of total propellant. Because CEA takes the
    enthalpy per mole of FUEL, it is rescaled with the O/F ratio:
        kJ/kg_fuel = deficit * (1 + O/F)
        kJ/mol     = kJ/kg_fuel * M_IPA
    """
    deficit_kj_per_mol_fuel = deficit_kj_per_kg_prop * (1.0 + OF) * IPA_MW_KG_MOL
    h_cal_per_mol = (IPA_HF_KJ_MOL - deficit_kj_per_mol_fuel) * KJ_TO_CAL

    name = f"IPAd{next(_counter)}"          # unique name, <= 15 characters
    card = f"""
fuel {name}  {IPA_FORMULA}   wt%=100.
h,cal={h_cal_per_mol:.3f}  t(k)={IPA_T_K}
"""
    add_new_fuel(name, card)
    return CEA_Obj(oxName="LOX", fuelName=name)


def chamber_state(deficit_kj_per_kg_prop: float) -> dict:
    cea = make_cea_obj(deficit_kj_per_kg_prop)
    t_k = float(cea.get_Tcomb(PC_PSIA, OF)) * RANKINE_TO_K
    cstar = float(cea.get_Cstar(PC_PSIA, OF)) * FT_TO_M
    mw, gam = cea.get_Chamber_MolWt_gamma(PC_PSIA, OF, EPS)
    mw, gam = float(mw), float(gam)
    r_gas = R_UNIV / mw                                     # [J/(kg K)]

    # Closed-form c* from chamber R, T, gamma, for cross-checking CEA.
    # (CEA's c* uses throat conditions, so expect a small mismatch.)
    f_gamma = gam * (2.0 / (gam + 1.0)) ** ((gam + 1.0) / (gam - 1.0))
    cstar_formula = np.sqrt(r_gas * t_k / f_gamma)

    return dict(deficit=deficit_kj_per_kg_prop, T=t_k, cstar=cstar, MW=mw,
                gamma=gam, R=r_gas, cstar_formula=cstar_formula)


# --------------------------------------------------------------------------
# SWEEP
# --------------------------------------------------------------------------
def run_sweep() -> dict:
    rows = [chamber_state(d) for d in np.linspace(0.0, D_MAX_KJ_KG, N_POINTS)]
    data = {k: np.array([r[k] for r in rows]) for k in rows[0]}
    data["eta"] = data["cstar"] / data["cstar"][0]

    keep = data["eta"] >= ETA_MIN_PLOT
    data = {k: v[keep] for k, v in data.items()}
    return data


def t_at_eta(data: dict, eta_target: float) -> float:
    """Interpolate chamber temperature [K] at a requested c* efficiency."""
    order = np.argsort(data["eta"])
    return float(np.interp(eta_target, data["eta"][order], data["T"][order]))


# --------------------------------------------------------------------------
# OUTPUT
# --------------------------------------------------------------------------
def print_table(d: dict) -> None:
    hdr = (f"{'deficit':>9} {'eta_c*':>7} {'Tc [K]':>8} {'MW':>7} {'gamma':>7} "
           f"{'R':>7} {'c*_CEA':>8} {'c*_form':>8}")
    if Q_IDEAL_KJ_KG:
        hdr += f" {'eta_Q':>7}"
    print(hdr)
    print(f"{'kJ/kg':>9} {'-':>7} {'K':>8} {'g/mol':>7} {'-':>7} "
          f"{'J/kgK':>7} {'m/s':>8} {'m/s':>8}")
    for i in range(len(d["eta"])):
        line = (f"{d['deficit'][i]:9.0f} {d['eta'][i]:7.4f} {d['T'][i]:8.1f} "
                f"{d['MW'][i]:7.3f} {d['gamma'][i]:7.4f} {d['R'][i]:7.1f} "
                f"{d['cstar'][i]:8.1f} {d['cstar_formula'][i]:8.1f}")
        if Q_IDEAL_KJ_KG:
            line += f" {1.0 - d['deficit'][i] / Q_IDEAL_KJ_KG:7.4f}"
        print(line)


def make_plots(d: dict, fname: str = "tc_vs_cstar_efficiency.png") -> None:
    eta = d["eta"]
    t0 = d["T"][0]

    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12, 4.8))

    ax1.plot(eta * 100, d["T"], "o-", label="CEA, reduced energy release (R, gamma updated)")
    ax1.plot(eta * 100, t0 * eta**2, "--",
             label=r"constant R, gamma:  T = eta$^2$ T$_{CEA}$")
    ax1.set_xlabel("c* efficiency [%]")
    ax1.set_ylabel("Chamber temperature [K]")
    ax1.set_title(f"LOX/IPA, Pc = {PC_PSIA:.0f} psia, O/F = {OF:.2f}")
    ax1.grid(True, alpha=0.3)
    ax1.legend(fontsize=8)

    ax2.plot(eta * 100, d["MW"], "o-", color="tab:green", label="MW [g/mol]")
    ax2.set_xlabel("c* efficiency [%]")
    ax2.set_ylabel("Chamber MW [g/mol]", color="tab:green")
    ax2.grid(True, alpha=0.3)
    ax3 = ax2.twinx()
    ax3.plot(eta * 100, d["gamma"], "s-", color="tab:red", label="gamma")
    ax3.set_ylabel("Chamber gamma [-]", color="tab:red")
    ax2.set_title("Gas properties along the same curve")

    fig.tight_layout()
    fig.savefig(fname, dpi=200)
    print(f"\nSaved {fname}")


if __name__ == "__main__":
    data = run_sweep()
    print_table(data)
    for e in (0.95, 0.90, 0.85, 0.80, 0.75, 0.70, 0.65, 0.60):
        print(f"eta_c* = {e:.2f}  ->  Tc = {t_at_eta(data, e):7.1f} K")
    make_plots(data)
    plt.show()
