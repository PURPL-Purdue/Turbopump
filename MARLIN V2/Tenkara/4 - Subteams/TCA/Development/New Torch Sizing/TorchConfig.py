"""Torch igniter chamber / injector pressure analysis (CH4 / O2).

Changes vs. the original:
  * Chamber temperature, gamma and R all come from one CEA solve
  * Solver converges on chamber pressure and raises if it does not converge
  * Injector pressures are treated as the source of truth; stiffness is
    reported as a check rather than used to override them
  * Discharge coefficients and a c* efficiency are explicit parameters
  * CEA objects are built once and reused
  * Results are returned in a dataclass instead of only printed
  * c* is computed in `analyze()` (Pc * A_throat / mdot) and carried in
    TorchResults instead of being referenced as an undefined name in
    print_results
"""

from dataclasses import dataclass

import numpy as np
import cea

# ---------------------------------------------------------------- constants
PSI_TO_PA = 6894.757293
PA_TO_PSI = 1 / PSI_TO_PA
BAR_TO_PA = 1e5

R_UNIVERSAL = 8314.462618   # [J/(kmol*K)]
MW_O2 = 31.998              # [kg/kmol]
MW_CH4 = 16.043             # [kg/kmol]
R_OX = R_UNIVERSAL / MW_O2      # [J/(kg*K)]
R_FUEL = R_UNIVERSAL / MW_CH4   # [J/(kg*K)]


# ------------------------------------------------------------------- config
@dataclass
class TorchConfig:
    # Geometry [m]
    D_throat: float = 0.004781
    D_fuel: float = 0.000849
    D_ox: float = 0.001135

    # Line temperatures [K]. Also used as the CEA reactant temperatures.
    T_ox: float = 22 + 273.15
    T_fuel: float = 22 + 273.15

    # Injector gas properties (ideal-gas approximation)
    k_ox: float = 1.40
    k_fuel: float = 1.31

    # Operating point
    mdot: float = 0.01765       # total mass flow [kg/s]
    of_ratio: float = 2.5       # O/F [--]
    stiffness_min: float = 0.3  # required (P_line - Pc) / Pc

    # Loss factors (1.0 = ideal; real injectors are typically 0.6-0.9)
    Cd_ox: float = 0.9
    Cd_fuel: float = 0.9
    Cd_throat: float = 0.9
    cstar_eff: float = 1.0      # applied as Tc_eff = Tc * eff^2 (c* ~ sqrt(Tc))

    # Solver
    p_guess: float = 200 * PSI_TO_PA
    tol: float = 1e-6
    max_iter: int = 100


@dataclass
class TorchResults:
    Pc: float
    k_c: float
    R_c: float
    Tc: float
    cstar: float
    mdot_ox: float
    mdot_f: float
    P_ox_line: float
    P_f_line: float
    stiffness_ox: float
    stiffness_f: float
    ox_choked: bool
    fuel_choked: bool
    stiffness_ok: bool


# ------------------------------------------------------------- basic helpers
def circle_area(diameter):
    """Area [m^2] of a circle from diameter [m]."""
    return np.pi * diameter ** 2 / 4


def choked_flow_factor(k):
    """(2/(k+1))^((k+1)/(2(k-1))) used in the choked mass flow equation."""
    return (2 / (k + 1)) ** ((k + 1) / (2 * (k - 1)))


def choked_pressure(mdot, area, k, T, R, Cd=1.0):
    """Upstream stagnation pressure [Pa] that chokes `mdot` through an orifice.

    mdot = Cd * A * P * sqrt(k / (R T)) * (2/(k+1))^((k+1)/(2(k-1)))
    """
    return mdot / (Cd * area * np.sqrt(k / (R * T)) * choked_flow_factor(k))


def critical_pressure_ratio(k):
    """Minimum P_upstream / P_downstream for choked flow."""
    return ((k + 1) / 2) ** (k / (k - 1))


# --------------------------------------------------------------- CEA wrapper
class CEAModel:
    """Builds the CEA objects once; call `props(p_c, OF)` inside the loop."""

    def __init__(self, T_fuel, T_ox):
        names = ["CH4", "O2"]
        self.T_reactant = np.array([T_fuel, T_ox])
        self.fuel_weights = np.array([1.0, 0.0])
        self.ox_weights = np.array([0.0, 1.0])

        self.reac = cea.Mixture(names)
        self.prod = cea.Mixture(names, products_from_reactants=True)
        self.solver = cea.RocketSolver(self.prod, reactants=self.reac)
        self.solution = cea.RocketSolution(self.solver)

    def props(self, p_c, of_ratio):
        """Return (gamma, R [J/kg/K], Tc [K]) at chamber pressure p_c [Pa]."""
        weights = self.reac.of_ratio_to_weights(
            self.ox_weights, self.fuel_weights, of_ratio=of_ratio
        )
        hc = self.reac.calc_property(
            cea.ENTHALPY, weights, self.T_reactant
        ) / cea.R

        self.solver.solve(
            self.solution, weights, p_c / BAR_TO_PA, hc=hc, iac=True
        )

        k = self.solution.gamma_s[0]
        R = cea.R / self.solution.MW[0]
        # NOTE: confirm `solution.T` exists in your CEA version. [0] is the
        # chamber (combustion) station, matching gamma_s[0] / MW[0] above.
        Tc = self.solution.T[0]
        return k, R, Tc


# ------------------------------------------------------------------- solver
def solve_chamber_pressure(cfg, cea_model):
    """Fixed-point iteration on chamber pressure.

    Returns (Pc [Pa], gamma, R, Tc_effective). Raises if not converged.
    """
    A_t = circle_area(cfg.D_throat)
    p_c = cfg.p_guess

    for _ in range(cfg.max_iter):
        k, R, Tc = cea_model.props(p_c, cfg.of_ratio)
        Tc_eff = Tc * cfg.cstar_eff ** 2

        # Choked throat: Pc = mdot * sqrt(R T) / (Cd * A * sqrt(k) * factor)
        p_new = choked_pressure(cfg.mdot, A_t, k, Tc_eff, R, cfg.Cd_throat)

        if abs(p_new - p_c) / p_c < cfg.tol:
            return p_new, k, R, Tc_eff
        p_c = p_new

    raise RuntimeError(
        f"Chamber pressure did not converge in {cfg.max_iter} iterations "
        f"(last Pc = {p_c * PA_TO_PSI:.1f} psi)"
    )


# --------------------------------------------------------------------- main
def analyze(cfg: TorchConfig) -> TorchResults:
    # Split total mass flow
    mdot_ox = cfg.mdot * cfg.of_ratio / (cfg.of_ratio + 1)
    mdot_f = cfg.mdot - mdot_ox

    # Chamber state (Pc, gamma, R, Tc all from the same CEA solution)
    cea_model = CEAModel(cfg.T_fuel, cfg.T_ox)
    Pc, k_c, R_c, Tc = solve_chamber_pressure(cfg, cea_model)

    # Characteristic velocity, c* = Pc * A_throat / mdot [m/s]
    A_t = circle_area(cfg.D_throat)
    cstar = Pc * A_t / cfg.mdot

    # Line pressures required to push the target flow through the injectors
    P_ox = choked_pressure(
        mdot_ox, circle_area(cfg.D_ox), cfg.k_ox, cfg.T_ox, R_OX, cfg.Cd_ox
    )
    P_f = choked_pressure(
        mdot_f, circle_area(cfg.D_fuel), cfg.k_fuel, cfg.T_fuel, R_FUEL, cfg.Cd_fuel
    )

    # Checks on the final pressures
    ox_choked = P_ox / Pc >= critical_pressure_ratio(cfg.k_ox)
    fuel_choked = P_f / Pc >= critical_pressure_ratio(cfg.k_fuel)

    stiff_ox = P_ox / Pc - 1
    stiff_f = P_f / Pc - 1
    stiffness_ok = min(stiff_ox, stiff_f) >= cfg.stiffness_min

    return TorchResults(
        Pc=Pc, k_c=k_c, R_c=R_c, Tc=Tc, cstar=cstar,
        mdot_ox=mdot_ox, mdot_f=mdot_f,
        P_ox_line=P_ox, P_f_line=P_f,
        stiffness_ox=stiff_ox, stiffness_f=stiff_f,
        ox_choked=ox_choked, fuel_choked=fuel_choked,
        stiffness_ok=stiffness_ok,
    )


def print_results(r: TorchResults, cfg: TorchConfig):
    print(f"Chamber Pressure:        {r.Pc * PA_TO_PSI:8.1f} psi")
    print(f"Chamber Temperature:     {r.Tc:8.1f} K")
    print(f"Chamber Gamma:           {r.k_c:8.4f}")
    print(f"Chamber R:               {r.R_c:8.2f} J/(kg*K)")
    print(f"Oxidizer Line Pressure:  {r.P_ox_line * PA_TO_PSI:8.1f} psi "
          f"(stiffness {r.stiffness_ox:.2f})")
    print(f"Fuel Line Pressure:      {r.P_f_line * PA_TO_PSI:8.1f} psi "
          f"(stiffness {r.stiffness_f:.2f})")
    print(f"Oxidizer injector choked: {r.ox_choked}")
    print(f"Fuel injector choked:     {r.fuel_choked}")
    print(f"Cstar: {r.cstar:0.1f} m/s")
    if not r.ox_choked or not r.fuel_choked:
        print("WARNING: an injector is not choked; the choked-flow sizing "
              "does not apply. Resize the injector.")
    if not r.stiffness_ok:
        print(f"WARNING: stiffness below required {cfg.stiffness_min:.2f}. "
              "Reduce injector diameter(s).")


def main():
    print("\nActual Hardware")
    cfg = TorchConfig(
    D_throat=6.8072e-3,          # 0.268 in -> m
    D_ox=1.6e-3,         # m
    D_fuel=1.5e-3,       # m
    mdot=0.0290,          # kg/s
    of_ratio=1.5,
    Cd_ox=0.8, Cd_fuel=0.8, Cd_throat=0.75,
    cstar_eff=1,
    stiffness_min=0.30,
)
    print_results(analyze(cfg), cfg)
    return


if __name__ == "__main__":
    main()