"""
Gas-side heat transfer: recovery temperature and the Bartz coefficient.
Heister et al., Rocket Propulsion, 6.3 (Eq. 6.13, 6.15).

Pure physics: no config reads. The caller passes everything in.
"""
import math
import numpy as np

from engine.cea_interface import gas_side_states


def recovery_temperature(T, M, gamma, Pr):
    """Eq. 6.15: T_r = T (1 + (gamma-1)/2 r M^2), recovery factor r = Pr^(1/3) (turbulent)."""
    return T * (1 + (gamma - 1) / 2 * Pr**(1 / 3) * M**2)


def gas_environment(eng, x, r, dims, cfg):
    """
    Hot-gas conditions at every station, the inputs to step 5 of 6.5.2.
    CEA gives T, M, gamma, mu, cp, Pr along the contour; added here:
        r      wall radius (Bartz uses D = 2r)
        G      rho V = m_dot / A, by continuity
        omega  mu ~ T^omega fit to the CEA data (for Bartz's mu_am/mu_0)
        T_r    recovery temperature, Eq. 6.15
    """
    gas = gas_side_states((r / eng["Rt"])**2, supersonic=(x > dims["x_t"]), cfg=cfg)
    gas["r"] = r
    gas["G"] = eng["m_dot"] / (math.pi * r**2)
    gas["omega"] = np.polyfit(np.log(gas["T"]), np.log(gas["mu"]), 1)[0]  # check this
    gas["T_r"] = recovery_temperature(gas["T"], gas["M"], gas["gamma"], gas["Pr"])
    return gas


def bartz_hg(D, T_wg, T_inf, G, T0, mu0, cp0, Pr0, omega):
    """
    Bartz (1957), Eq. 6.13:

        h_g = 0.026 / D^0.2 * (mu^0.2 cp / Pr^0.6)_0 * (rho V)^0.8
              * (rho_am / rho)^0.8 * (mu_am / mu_0)^0.2

    ()_0 stagnation, ()_am at T_am = (T + T_wg)/2. Closures:
        rho_am/rho = T/T_am        ideal gas, p constant across the boundary layer
        mu_am/mu_0 = (T_am/T0)^omega
    """
    T_am = 0.5 * (T_inf + T_wg)
    stag = mu0**0.2 * cp0 / Pr0**0.6
    return 0.026 / D**0.2 * stag * G**0.8 * (T_inf / T_am)**0.8 * ((T_am / T0)**omega)**0.2  # need to validate how the am properties are calculated


def bartz_at_station(gas, i, T_wg):
    """Bartz h_g at station i of a gas_environment() dict, for a wall temperature T_wg."""
    return bartz_hg(2 * gas["r"][i], T_wg, gas["T"][i], gas["G"][i],
                    gas["T0"], gas["mu0"], gas["cp0"], gas["Pr0"], gas["omega"])
