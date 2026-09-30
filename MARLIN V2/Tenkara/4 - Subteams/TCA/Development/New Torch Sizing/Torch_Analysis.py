# Torch_Analysis
import numpy as np
from pyfluids import Fluid, FluidsList, Input
from New_TCA_Igniter_Sizing import diameter_to_area
from New_TCA_Igniter_Sizing import area_to_diameter
from New_TCA_Igniter_Sizing import choked_backpressure
import cea as cea


CM_TO_IN = 1 / 2.54          # Centimeters to Inches
IN_TO_CM = 2.54              # Inches to Centimeters
BAR_TO_PA = 1e5              # Bar to Pascals
PSI_TO_PA = 6894.757293      # PSI to Pascals
PA_TO_PSI = 1 / PSI_TO_PA    # Pascals to PSI

T_AMB_CELSIUS = 20           # [deg C]
T_AMB_KELVIN = 293.15        # [K]

R_UNIVERSAL = 8314.462618    # [J/(kmol*K)]

MW_O2 = 31.998               # [kg/kmol]
MW_CH4 = 16.043              # [kg/kmol]

R_ox = R_UNIVERSAL / MW_O2       # [J/(kg*K)]
R_fuel = R_UNIVERSAL / MW_CH4    # [J/(kg*K)]


def chamber_pressure(throatArea, k, massflow, chamberTemp, R):
#   Input:
    # Throat area [m^2]
    # Heat capacity ratio [--]
    # Mass flow rate [kg/s]
    # Chamber temperature [K]
    # Specific gas constant [J/(kg*K)]
#   Output:
    # Chamber Pressure [Pa]

    return (massflow / (throatArea * k)) * (np.sqrt(k * R * chamberTemp) / np.sqrt((2 / (k + 1)) ** ((k + 1) / (k - 1))))


def get_k(p_c, OF):

    reac_names = ["CH4", "O2"]
    T_reactant = np.array([T_AMB_KELVIN, T_AMB_KELVIN])
    fuel_weights = np.array([1.0, 0.0])
    ox_weights = np.array([0.0, 1.0])

    # Convert pressure from Pa to bar for CEA
    p_c = p_c / BAR_TO_PA

    reac = cea.Mixture(reac_names)
    prod = cea.Mixture(reac_names, products_from_reactants=True)
    solver = cea.RocketSolver(prod, reactants=reac)
    solution = cea.RocketSolution(solver)

    weights = reac.of_ratio_to_weights(ox_weights, fuel_weights, of_ratio=OF)
    hc = reac.calc_property(cea.ENTHALPY, weights, T_reactant) / cea.R

    solver.solve(solution, weights, p_c, hc=hc, iac=True)

    k_c = solution.gamma_s[0]

    # Specific gas constant [J/(kg*K)]
    R = cea.R / solution.MW[0]

    return k_c, R


def critical_pressure(k):
#   Input:
    # Gamma [--]
#   Output:
    # Critical pressure ratio [--]

    return ((k + 1) / 2) ** (k / (k - 1))


def solve_chamber_pressure(throatArea, massflow, chamberTemp, OF, k, R):

    for i in range(100):

        p_c = chamber_pressure(throatArea, k, massflow, chamberTemp, R)
        k_new, R_new = get_k(p_c, OF)

        if abs(k_new - k) < 1e-5:

            k = k_new
            R = R_new
            p_c = chamber_pressure(throatArea, k, massflow, chamberTemp, R)

            break

        k = k_new
        R = R_new

    return k, R, p_c


def main():

    ## Fixed torch parameters

    D_t = 0.004781        # Throat diameter [m]
    Df = 0.000849         # Fuel injector diameter [m]
    Dox = 0.001135        # Oxidizer injector diameter [m]

    Tc = 3205             # Chamber temperature [K]
    Tox = 22 + 273.15     # Oxidizer line temperature [K]
    Tf = 22 + 273.15      # Fuel line temperature [K]

    k_c = 1.2             # Chamber gamma [--]
    k_ox = 1.40           # Oxidizer gamma [--]
    k_fuel = 1.31         # Fuel gamma [--]

    mdot = 0.01765        # Total mass flow rate [kg/s]
    of_ratio = 2.5        # O/F Ratio [--]
    stiffness = 1.65      # System stiffness [--]


    # Split total mass flow
    mdot_ox = mdot * of_ratio / (of_ratio + 1)
    mdot_f = mdot - mdot_ox


    # Initial CEA estimate at 200 psi
    k_c, R = get_k(200 * PSI_TO_PA, of_ratio)


    # Iterate chamber pressure, gamma, and R
    k_c, R, Pc = solve_chamber_pressure(diameter_to_area(D_t), mdot, Tc, of_ratio, k_c, R)


    # Injector pressures required for given mass flow
    Pox = chamber_pressure(diameter_to_area(Dox), k_ox, mdot_ox, Tox, R_ox)
    Pf = chamber_pressure(diameter_to_area(Df), k_fuel, mdot_f, Tf, R_fuel)


    # Critical pressure ratios
    Pcrit_ox = critical_pressure(k_ox)
    Pcrit_f = critical_pressure(k_fuel)


    # Minimum line pressures required for choking
    Pcrit_line_ox = Pc * Pcrit_ox
    Pcrit_line_f = Pc * Pcrit_f


    # Line pressure from stiffness
    pox_line_stiffness = Pc * (1 + stiffness)
    pf_line_stiffness = Pc * (1 + stiffness)


    if (pox_line_stiffness / Pc) >= Pcrit_ox:

        print("Oxidizer is choked")
        pox_line = max(Pox, Pcrit_line_ox, pox_line_stiffness)

    else:

        print("Oxidizer is not choked")
        pox_line = pox_line_stiffness


    if (pf_line_stiffness / Pc) >= Pcrit_f:

        print("Fuel is choked")
        pf_line = max(Pf, Pcrit_line_f, pf_line_stiffness)

    else:

        print("Fuel is not choked")
        pf_line = pf_line_stiffness


    # Final results only
    print("Chamber Pressure:", Pc * PA_TO_PSI, "psi")
    print("Chamber Gamma:", k_c)
    print("Chamber R:", R, "J/(kg*K)")
    print("Oxidizer Line Pressure:", pox_line * PA_TO_PSI, "psi")
    print("Fuel Line Pressure:", pf_line * PA_TO_PSI, "psi")


if __name__ == "__main__":
    main()