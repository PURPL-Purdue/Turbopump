import numpy as np
import pint 
import matplotlib.pyplot as plt
from cea_setup import cea_setup

def format_values(values, skip_index=1, width=10, precision=3):
    return " ".join(
        f"{float(values[i]):{width}.{precision}f}"
        for i in range(len(values))
        if i != skip_index
    )


def get_isp_heatmap(ureg, of_range, Pc_range, reac_names, oxidant_weights, fuel_weights, reac_temps, Pexit, ac_at):
    g0 = 9.81 * ureg.m / (ureg.s)**2
    Pamb = 1 * ureg.bar

    of_set = np.arange(of_range[0], of_range[1], 0.1) 
    Pc_set = np.arange(Pc_range[0].magnitude, Pc_range[1].magnitude, 1) * Pc_range[0].units
    isp_grid = np.zeros((len(of_set), len(Pc_set)))


    for i, of in enumerate(of_set):
        for j, Pc in enumerate(Pc_set):            
            solution = cea_setup(
                reac_names=reac_names, 
                oxidant_weights=oxidant_weights, 
                fuel_weights=fuel_weights, 
                of_ratio=of, 
                reac_temps=reac_temps.to(ureg.K).magnitude, 
                Pc=Pc.to(ureg.bar).magnitude, 
                Prat = float((Pc / Pexit).to_base_units().magnitude),
                ac_at=ac_at)

            # Specific impulse at sea level (s)
            isp = solution.Isp[2] * (ureg.m / ureg.s)
            isp = isp / g0
            
            isp_grid[i, j] = isp.magnitude

    return of_set, Pc_set, isp_grid

def get_Tc_heatmap(ureg, of_range, Pc_range, reac_names, oxidant_weights, fuel_weights, reac_temps, Pexit, ac_at):
    Pamb = 1 * ureg.bar

    of_set = np.arange(of_range[0], of_range[1], 0.1) 
    Pc_set = np.arange(Pc_range[0].magnitude, Pc_range[1].magnitude, 1) * Pc_range[0].units
    Tc_grid = np.zeros((len(of_set), len(Pc_set)))

    for i, of in enumerate(of_set):
        for j, Pc in enumerate(Pc_set):            
            solution = cea_setup(
                reac_names=reac_names, 
                oxidant_weights=oxidant_weights, 
                fuel_weights=fuel_weights, 
                of_ratio=of, 
                reac_temps=reac_temps.to(ureg.K).magnitude, 
                Pc=Pc.to(ureg.bar).magnitude, 
                Prat = float((Pc / Pexit).to_base_units().magnitude),
                ac_at=ac_at)

            # Specific impulse at sea level (s)
            Tc = solution.T[0] * ureg.K
            Tc_grid[i, j] = Tc.magnitude

    return of_set, Pc_set, Tc_grid

def get_isp_density_heatmap(ureg, rho_ox, rho_f, of_range, Pc_range, reac_names, oxidant_weights, fuel_weights, reac_temps, Pexit, ac_at):
    g0 = 9.81 * ureg.m / (ureg.s)**2
    Pamb = 1 * ureg.bar

    of_set = np.arange(of_range[0], of_range[1], 0.1) 
    Pc_set = np.arange(Pc_range[0].magnitude, Pc_range[1].magnitude, 1) * Pc_range[0].units
    isp_density_grid = np.zeros((len(of_set), len(Pc_set)))


    for i, of in enumerate(of_set):
        for j, Pc in enumerate(Pc_set):            
            solution = cea_setup(
                reac_names=reac_names, 
                oxidant_weights=oxidant_weights, 
                fuel_weights=fuel_weights, 
                of_ratio=of, 
                reac_temps=reac_temps.to(ureg.K).magnitude, 
                Pc=Pc.to(ureg.bar).magnitude, 
                Prat = float((Pc / Pexit).to_base_units().magnitude),
                ac_at=ac_at)

            # Specific impulse at sea level (s)
            isp = solution.Isp[2] * (ureg.m / ureg.s)
            isp = isp / g0

            mass_frac_ox = of / (1 + of)
            mass_frac_fuel = 1 / (1 + of)

            isp_denisty = isp * 1 / (mass_frac_ox/rho_ox + mass_frac_fuel/rho_f)
            isp_density_grid[i, j] = isp_denisty.magnitude

    return of_set, Pc_set, isp_density_grid

#thrust plots with lines of constant mdot
def get_thrust_heatmap(ureg, of_range, Pc_range, reac_names, oxidant_weights, fuel_weights, reac_temps, Pexit, ac_at):
    g0 = 9.81 * ureg.m / (ureg.s)**2

    of_set = np.arange(of_range[0], of_range[1], 0.1) 
    Pc_set = np.arange(Pc_range[0].magnitude, Pc_range[1].magnitude, 1) * Pc_range[0].units
    thrust_grid = np.zeros((len(of_set), len(Pc_set)))


    for i, of in enumerate(of_set):
        for j, Pc in enumerate(Pc_set):            
            solution = cea_setup(
                reac_names=reac_names, 
                oxidant_weights=oxidant_weights, 
                fuel_weights=fuel_weights, 
                of_ratio=of, 
                reac_temps=reac_temps.to(ureg.K).magnitude, 
                Pc=Pc.to(ureg.bar).magnitude, 
                Prat = float((Pc / Pexit).to_base_units().magnitude),
                ac_at=ac_at)

            # Specific impulse at sea level (m/s)
            isp = solution.Isp[2] * (ureg.m / ureg.s)
            thrust = isp * mdot 
            
            thrust_grid[i, j] = thrust.magnitude


    return of_set, Pc_set, thrust_grid
