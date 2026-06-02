import numpy as np
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
import pint
import cea


# ============================================================
# UNIT SETUP
# ============================================================

ureg = pint.UnitRegistry()

g0 = 9.81 * ureg.m / ureg.s**2


# ============================================================
# ENGINE INPUTS
# ============================================================

of_range = (0.55, 6.0)

Pc_range = (1, 60) * ureg.bar

Pamb = 1.01325 * ureg.bar
Pexit = 18 * ureg.psi

# Constant thrust contour levels
thrust_levels = np.arange(3000,5500,500) * ureg.lbf


# ============================================================
# PROPELLANTS
# ============================================================

ipa_reac_names = [
    "C3H8O,2propanol",
    "O2(L)"
]

eth_reac_names = [
    "C2H5OH(L)",
    "O2(L)"
]

reac_temps = np.array(
    [298.15, 90.17]
) * ureg.K

fuel_weights = np.array([1.0, 0.0])

oxidant_weights = np.array([0.0, 1.0])

cont_ratio = 1.8


# ============================================================
# CEA SETUP FUNCTION
# ============================================================

def cea_setup(
    reac_names,
    oxidant_weights,
    fuel_weights,
    of_ratio,
    reac_temps,
    Pc,
    Prat,
    ac_at
):

    reac = cea.Mixture(reac_names)

    prod = cea.Mixture(
        reac_names,
        products_from_reactants=True
    )

    solver = cea.RocketSolver(
        prod,
        reactants=reac
    )

    solution = cea.RocketSolution(solver)

    weights = reac.of_ratio_to_weights(
        oxidant_weights,
        fuel_weights,
        of_ratio
    )

    hc = (
        reac.calc_property(
            cea.ENTHALPY,
            weights,
            reac_temps
        ) / cea.R
    )

    solver.solve(
        solution,
        weights,
        float(Pc),
        float(Prat),
        ac_at=float(ac_at),
        iac=True,
        hc=hc
    )

    return solution


# ============================================================
# BUILD GRIDS
# ============================================================

def build_grids(
    of_range,
    Pc_range,
    reac_names,
    oxidant_weights,
    fuel_weights,
    reac_temps,
    Pexit,
    ac_at
):

    of_set = np.arange(
        of_range[0],
        of_range[1],
        0.1
    )

    Pc_set = (
        np.arange(
            Pc_range[0].magnitude,
            Pc_range[1].magnitude + 1,
            1
        ) * Pc_range[0].units
    )

    Tc_grid = np.zeros(
        (len(of_set), len(Pc_set))
    )

    ve_grid = np.zeros(
        (len(of_set), len(Pc_set))
    )

    for i, of in enumerate(of_set):

        for j, Pc in enumerate(Pc_set):

            Prat = float(
                (Pc / Pexit)
                .to_base_units()
                .magnitude
            )

            solution = cea_setup(
                reac_names=reac_names,
                oxidant_weights=oxidant_weights,
                fuel_weights=fuel_weights,
                of_ratio=float(of),
                reac_temps=reac_temps
                .to(ureg.K)
                .magnitude,
                Pc=float(
                    Pc.to(ureg.bar)
                    .magnitude
                ),
                Prat=Prat,
                ac_at=ac_at
            )

            # ------------------------------------------------
            # CHAMBER TEMPERATURE
            # ------------------------------------------------

            Tc = solution.T[0] * ureg.K

            Tc_grid[i, j] = (
                Tc.to(ureg.K)
                .magnitude
            )

            # ------------------------------------------------
            # EFFECTIVE EXHAUST VELOCITY
            # ------------------------------------------------

            ve = (
                solution.Isp[2]
                * ureg.m
                / ureg.s
            )

            ve_grid[i, j] = (
                ve.to(ureg.m / ureg.s)
                .magnitude
            )

    return (
        of_set,
        Pc_set,
        Tc_grid,
        ve_grid
    )


# ============================================================
# PLOT FUNCTION
# ============================================================

def plot_temperature_map(
    of_set,
    Pc_set,
    Tc_grid,
    ve_grid,
    thrust_levels
):

    X, Y = np.meshgrid(
        Pc_set.magnitude,
        of_set
    )

    plt.figure(figsize=(12, 8))

    # --------------------------------------------------------
    # TEMPERATURE HEATMAP
    # --------------------------------------------------------

    heatmap = plt.contourf(
        X,
        Y,
        Tc_grid,
        levels=50,
    )

    cbar = plt.colorbar(heatmap)

    cbar.set_label(
        "Chamber Temperature [K]"
    )

    # --------------------------------------------------------
    # MDOT VALUES
    # --------------------------------------------------------

    mdot_values = [10] * ureg.kg / ureg.s


    legend_handles = []

    # --------------------------------------------------------
    # CONSTANT THRUST CONTOURS
    # --------------------------------------------------------

    for mdot in mdot_values:

        thrust_grid = (
            ve_grid
            * mdot.to(ureg.kg / ureg.s)
            .magnitude
        )

        contours = plt.contour(
            X,
            Y,
            thrust_grid,
            levels=thrust_levels
            .to(ureg.N)
            .magnitude,
            colors="black",
            linewidths=1.5
        )

        plt.clabel(
            contours,
            inline=True,
            fontsize=8,
            colors="black",
            fmt=lambda x:
            f"{x/1000:.1f} kN"
        )

    # --------------------------------------------------------
    # LABELS
    # --------------------------------------------------------

    plt.xlabel(
        "Chamber Pressure [bar]"
    )

    plt.ylabel(
        "O/F Ratio"
    )

    plt.title(
        "Chamber Temperature Heatmap\n"
        "with Constant Thrust Contours\n"
        "mdot=10kg/s"
    )

    plt.tight_layout()

    plt.show()


# ============================================================
# RUN
# ============================================================

of_set, Pc_set, Tc_grid, ve_grid = build_grids(
    of_range=of_range,
    Pc_range=Pc_range,
    reac_names=ipa_reac_names,
    oxidant_weights=oxidant_weights,
    fuel_weights=fuel_weights,
    reac_temps=reac_temps,
    Pexit=Pexit,
    ac_at=cont_ratio
)

plot_temperature_map(
    of_set,
    Pc_set,
    Tc_grid,
    ve_grid,
    thrust_levels
)