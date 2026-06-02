import numpy as np
import matplotlib.pyplot as plt
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

Pc_range = (15, 60) * ureg.bar

Pamb = 1.01325 * ureg.bar
Pexit = 18 * ureg.psi

thrust_levels = [
    10000,
    20000,
    30000,
    40000,
    50000
] * ureg.N


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
# BUILD EXHAUST VELOCITY GRID
# ============================================================

def get_ve_grid(
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
        ve_grid
    )


# ============================================================
# PLOT MDOT MAP
# ============================================================

def plot_mdot_map(
    of_set,
    Pc_set,
    ve_grid,
    thrust_levels
):

    X, Y = np.meshgrid(
        Pc_set.magnitude,
        of_set
    )

    # --------------------------------------------------------
    # CHOOSE REFERENCE THRUST FOR HEATMAP
    # --------------------------------------------------------

    reference_thrust = 30000 * ureg.N

    mdot_grid = (
        reference_thrust.to(ureg.N).magnitude
        / ve_grid
    )

    plt.figure(figsize=(12, 8))

    # --------------------------------------------------------
    # MDOT HEATMAP
    # --------------------------------------------------------

    heatmap = plt.contourf(
        X,
        Y,
        mdot_grid,
        levels=50
    )

    cbar = plt.colorbar(heatmap)

    cbar.set_label(
        "Required Mass Flow Rate [kg/s]"
    )

    # --------------------------------------------------------
    # CONSTANT THRUST CONTOURS
    # --------------------------------------------------------

    for thrust in thrust_levels:

        thrust_grid = (
            thrust.to(ureg.N).magnitude
            / ve_grid
        )

        contours = plt.contour(
            X,
            Y,
            thrust_grid,
            levels=[5, 10, 15, 20],
            linewidths=1.5
        )

        plt.clabel(
            contours,
            inline=True,
            fontsize=8,
            fmt=lambda x:
            f"{x:.0f} kg/s"
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
        "Mass Flow Rate Map\n"
        "with Constant-Thrust Contours"
    )

    plt.tight_layout()

    plt.show()


# ============================================================
# RUN
# ============================================================

of_set, Pc_set, ve_grid = get_ve_grid(
    of_range=of_range,
    Pc_range=Pc_range,
    reac_names=ipa_reac_names,
    oxidant_weights=oxidant_weights,
    fuel_weights=fuel_weights,
    reac_temps=reac_temps,
    Pexit=Pexit,
    ac_at=cont_ratio
)

plot_mdot_map(
    of_set,
    Pc_set,
    ve_grid,
    thrust_levels
)