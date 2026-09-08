# ============================================================
# 0-D ADIABATIC WELL-STIRRED REACTOR (WSR)
# MINIMUM IGNITION ENERGY SOLVER
# Paper-faithful implementation
# ============================================================

import cantera as ct
import numpy as np

# -------------------------------
# USER INPUTS
# -------------------------------

mech = "A2NTC_skeletal.yaml"     # reaction mechanism
P = 20e5                        # chamber pressure [Pa]
T_in = 300.0                    # inlet temperature [K]
OF = 2.0                        # oxidizer-to-fuel ratio (mass)
tau = 5e-3                      # residence time [s]
mdot_in = 1.0                   # inlet mass flow rate [kg/s]

T_ign_crit = 2000.0             # ignition temperature threshold [K]

# -------------------------------
# MIXTURE (MASS FRACTIONS)
# -------------------------------

Y = {
    "O2": OF / (OF + 1.0),
    "POSF10325": 1.0 / (OF + 1.0)
}

# -------------------------------
# IGNITION ENERGY INJECTION
# -------------------------------

def apply_ignition_energy(reactor, E_ign):
    """
    Apply a delta-function ignition energy at t = 0
    """
    gas = reactor.thermo
    dh = E_ign / reactor.mass
    gas.HP = gas.h + dh, gas.P
    

# -------------------------------
# IGNITION TEST
# -------------------------------

def ignites(E_ign):
    """
    Returns True if ignition occurs within residence time
    """

    # Fresh gas
    gas = ct.Solution(mech)
    gas.TPY = T_in, P, Y

    # Reactor
    r = ct.IdealGasReactor(gas, energy='on')
    r.volume = mdot_in * tau / gas.density

    # Reservoirs
    res_in = ct.Reservoir(gas)
    res_out = ct.Reservoir(gas)

    # Fixed mass flows (choked outlet assumed)
    ct.MassFlowController(res_in, r, mdot=mdot_in)
    ct.MassFlowController(r, res_out, mdot=mdot_in)

    sim = ct.ReactorNet([r])

    # Delta-function ignition energy
    apply_ignition_energy(r, E_ign)

    Tmax = gas.T
    t = 0.0

    while t < tau:
        t = sim.step()
        Tmax = max(Tmax, gas.T)

        if gas.T > T_ign_crit:
            return True

    return False


# -------------------------------
# BISECTION SEARCH FOR E_min
# -------------------------------

E_low = 0.0
E_high = 1e5   # J (increase if no ignition)

# Ensure upper bound ignites
if not ignites(E_high):
    raise RuntimeError("Upper bound too low — increase E_high")

for _ in range(30):
    E_mid = 0.5 * (E_low + E_high)

    if ignites(E_mid):
        E_high = E_mid
    else:
        E_low = E_mid

E_min = E_high

# -------------------------------
# TORCH MASS FLOW ESTIMATE
# -------------------------------

Δh_rxn = 1.2e7   # J/kg (torch reaction energy — OF dependent)

mdot_torch = mdot_in * E_min / Δh_rxn

# -------------------------------
# RESULTS
# -------------------------------

print("\n================ RESULTS ================")
print(f"Minimum ignition energy : {E_min:.3e} J")
print(f"Torch mass flow rate    : {mdot_torch:.3e} kg/s")
print("========================================\n")
