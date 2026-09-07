import cantera as ct
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import time

PSI2PA = 6894.76

## Experiment Set Up ##
OF = 1
residence_time = 100 # [s]

# Propellant Temperatures [K]
T_comb = 1000
T_amb = 298
T_ox = 298
T_fuel = 298
T_mix = T_ox *(OF/(OF+1)) + T_fuel *(1/(OF+1))

# Propellant Pressures [Pa]
P_wsr = 500 * PSI2PA
P_amb = 14.7 * PSI2PA

# Mass Flow Rates [kg/s]
mdot_total = 0.12
mdot_ox = mdot_total*(OF/(OF+1))
mdot_fuel = mdot_total*(1/(OF+1))

# Jet-A reaction mechanism
mech = "A2NTC_skeletal.yaml"

# Gas Compositions #
# Pre-Mixed Inlet
Y_mix = {"POSF10325": (1/(OF+1)), "O2": (OF/(OF+1))}
sol_inlet = ct.Solution(mech)
sol_inlet.TPY = T_mix, P_wsr, Y_mix

# Arbitrary Exhaust 
sol_exhaust = ct.Solution(mech)
sol_exhaust.TP = T_amb, P_amb

def define_network(inlet):
## Define Reactor Network ##

    # Reservoirs
    res_inlet = ct.Reservoir(inlet, clone=True)
    res_exhaust = ct.Reservoir(sol_exhaust, clone=True)

    # Well Stirred Reactor & set residence time to be internally conistent with volume
    rctr_wsr = ct.ConstPressureReactor(inlet, energy='on')
    rho_inlet = rctr_wsr.density 
    rctr_wsr.volume = mdot_total * residence_time / rho_inlet
    
    # Mass Flow Controllers
    mfc_inlet = ct.MassFlowController(upstream=res_inlet, downstream=rctr_wsr, mdot=mdot_total)
    mfc_outlet = ct.MassFlowController(upstream=rctr_wsr, downstream=res_exhaust, mdot=mdot_total)
    #pc_outlet = ct.PressureController(upstream=rctr_wsr, downstream=res_exhaust, primary=mfc_inlet, K=1e-5)

    sim = ct.ReactorNet([rctr_wsr])
    
    return sim, rctr_wsr

def add_energy(E_ign):
    # Returns Clone with Added Enthalpy to Inlet Mixture, keep other properties the same
    inlet_clone = ct.Solution(mech)
    inlet_clone.TPY = T_mix, P_wsr, Y_mix
    dh = E_ign / mdot_total
    inlet_clone.HPY = inlet_clone.h + dh , inlet_clone.P, inlet_clone.Y

    return inlet_clone

def ignites(E_ign, plotFlag=False):
    inlet_clone = add_energy(E_ign)
    (sim, rctr_wsr) = define_network(inlet_clone)
    max_sim_time = 5000 #[s]
    if plotFlag:
        time_history = ct.SolutionArray(sol_inlet, extra=["t"])
        step_solution(time_history, max_sim_time, rctr_wsr, sim)
        print(time_history)
        return True

    else:
        t = 0.0
        sim.rtol = 1e-12
        sim.atol = 1e-12
        while t < 500:
            t = sim.step()
            if rctr_wsr.phase.T > T_comb:
                return True

def step_solution(time_history, max_sim_time, reactor, sim):
    # Start the stopwatch
    tic = time.time()

    # Set simulation start time to zero
    t = 0
    counter = 1
    while t < max_sim_time:
        t = sim.step()

        # We will store only every 10th value. Remember, we have 1200+ species, so there
        # will be 1200+ columns for us to work with
        if counter % 10 == 0:
            # Extract the state of the reactor
            time_history.append(reactor.phase.state, t=t)

        counter += 1

    # Stop the stopwatch
    toc = time.time()
    print(f"Simulation Took {toc-tic:3.2f}s to compute, with {counter} steps")
    
    #plot results
    fig, axes = plt.subplots(nrows=2, ncols=2, figsize=(10, 8))
    axes[0,0].semilogx(time_history.t, time_history.T,"-o", label="Temperature")
    axes[0,0].set_xlabel("Time [s]")
    axes[0,0].set_ylabel("Mixture Temperature [K]")
    axes[0,0].set_title("Mixture Temperature vs. Time")
    axes[0,0].legend()
    axes[0,0].grid(True)
    
    axes[1,1].semilogx(time_history.t, time_history.P / PSI2PA,"-o", label="Pressure")
    axes[1,1].set_xlabel("Time [s]")
    axes[1,1].set_ylabel("Mixture Pressure [psi]")
    axes[1,1].set_title("Mixture Pressure vs. Time")
    axes[1,1].legend()
    axes[1,1].grid(True)

    axes[0,1].semilogx(time_history.t, time_history.O2, "-o", label="O2")
    axes[0,1].semilogx(time_history.t, time_history.POSF10325, "-o", label="Jet-A")
    axes[0,1].set_xlabel("Time [s]")
    axes[0,1].set_ylabel("Mass Fraction")
    axes[0,1].set_title("Mixture Composition vs. Time")  
    axes[0,1].legend()
    axes[0,1].grid(True)
    
    axes[1,0].semilogx(time_history.t, time_history.u, "-o", label="Internal Energy")
    axes[1,0].semilogx(time_history.t, time_history.h, "-o", label='Enthalpy')
    axes[1,0].set_xlabel("Time [s]")
    axes[1,0].set_ylabel("Enthalpy, Internal Energy [J/kg]")
    axes[1,0].set_title("Mixture Energy vs. Time")  
    axes[1,0].legend()
    axes[1,0].grid(True)
    plt.show()

    return time_history

def steady_solve(sim, reactor):
    # Solve Sim to Steady State
    tic = time.time()
    sim.solve_steady()
    toc = time.time()
    print(f"Simulation Completed in {toc-tic:3.2f} seconds.")
    print(f"Final State: {reactor.T} K, {reactor.thermo.P / PSI2PA} psi")
    
    return reactor.T, reactor.thermo.P

def main():
    # Find Minimum Ignition Energy Case

    # MIE upper/lower bound guesses [J]
    E_low = 0.0
    E_high = 1e5

    # Bisecton Search to find MIE #
    if not ignites(E_high):
        raise RuntimeError("Upper bound too low — increase E_high")
    for xx in range(30):
        E_mid = 0.5 * (E_low + E_high)

        if ignites(E_mid):
            E_high = E_mid
        else:
            E_low = E_mid

        E_min = E_high
    print(f"Minimum Ignition Energy: {E_min:0.3f} J")

    # Once Minimum Ignition Energy Case is Found
    ## Step Through Solution to Get Temporal Data ## 
    ignites(E_min, plotFlag=True)
    
    
    return

if __name__ == "__main__":
    main()