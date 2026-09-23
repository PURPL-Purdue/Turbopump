import cantera as ct
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import time
from pyfluids import Fluid, FluidsList, Input
from rocketcea.cea_obj import CEA_Obj
import warnings
import os
warnings.filterwarnings("ignore")

# Helper Functions:
# function to get mass averaged chemical potential of mixture 
def mass_avg_chem_pot(mixture):
    mu = mixture.chemical_potentials      # J/kmol
    W  = mixture.molecular_weights        # kg/kmol
    Y  = mixture.Y                        # mass fractions
    return np.dot(Y, mu / W)              # J/kg mixture

def define_network(inlet, mdot_total, residence_time, sol_exhaust):
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

def add_energy(e_ign, mech, T_mix, P_wsr, Y_mix):
    # Returns Clone with Added Enthalpy to Inlet Mixture, keep other properties the same
    inlet_clone = ct.Solution(mech)
    inlet_clone.TPY = T_mix, P_wsr, Y_mix
    dh = e_ign # [J/kg]
    inlet_clone.HP = inlet_clone.h + dh, inlet_clone.P

    return inlet_clone

def ignites(e_ign, residence_time, sol_inlet, T_comb, mdot_total, sol_exhaust, mech, prop_names, save_path, plotFlag=False):
    T_mix, P_wsr, Y_mix = sol_inlet.TPY 
    inlet_clone = add_energy(e_ign, mech, T_mix, P_wsr, Y_mix)
    (sim, rctr_wsr) = define_network(inlet_clone, mdot_total, residence_time, sol_exhaust)
    max_sim_time = 500 #[s]
    if plotFlag:
        time_history = ct.SolutionArray(sol_inlet, extra=["t"])
        step_solution(time_history, max_sim_time, rctr_wsr, sim, prop_names, save_path)
        #print(time_history)
        return True

    else:
        t = 0.0
        sim.rtol = 1e-12
        sim.atol = 1e-12
        while t < residence_time:
            t = sim.step()
            if rctr_wsr.phase.T > T_comb:
                return True

def step_solution(time_history, max_sim_time, reactor, sim, prop_names, save_path):
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
    
    makeplot(time_history, prop_names, save_path)

    return time_history

def steady_solve(sim, reactor):
    PSI2PA = 6894.76
    # Solve Sim to Steady State
    tic = time.time()
    sim.solve_steady()
    toc = time.time()
    print(f"Simulation Completed in {toc-tic:3.2f} seconds.")
    print(f"Final State: {reactor.T:0.2f} K, {reactor.thermo.P / PSI2PA :0.2f} psi")
    
    return reactor.T, reactor.thermo.P

def makeplot(time_history, prop_names, save_path):
    PSI2PA = 6894.76 # psi to pa
    R2K = 0.555556 # Rankine to Kelvin
    BTU2J = 1055.06 # BTU to J
    LBM2KG = 0.453592 # lbm to kg
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

    ox = getattr(time_history, prop_names['ox'])
    fuel = getattr(time_history, prop_names['fuel'])
    axes[0,1].semilogx(time_history.t, ox, "-o", label=prop_names['ox'])
    axes[0,1].semilogx(time_history.t, fuel, "-o", label=prop_names['fuel'])
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

    fig_path = os.path.join(save_path, "reaction.png")
    print(fig_path)
    plt.savefig(fig_path)
    plt.close()
    return

def Auto_Igntition_Temp(bisect_count, mdot, OF, residence_time, Pc, mech, prop_temps, save_path):
# Arguments: 
    # mdot_total: mass flow rate to be ignited [lbm/s]
    # OF: oxidizer-fuel ratio
    # residence_time: Combustion Residence Time [s]
    # P_wsr: Combustor Pressure [psi]
    # mech: chemical kinetics reaction mechanism
    # prop_strings: Dictionary of Form fuel: {"fuel_name", T_fuel_init, "ox_name": T_ox_init} [K] 
    
    PSI2PA = 6894.76 # psi to pa
    R2K = 0.555556 # Rankine to Kelvin
    BTU2J = 1055.06 # BTU to J
    LBM2KG = 0.453592 # lbm to kg

    # Parse Propellant Names/Temperatures
    fuel,ox = prop_temps.keys()
    T_ox = prop_temps[ox]
    T_fuel = prop_temps[fuel]
    prop_names = {"fuel":fuel, "ox":ox}

    # WSR Temperatures [K]
    T_comb = 1000
    T_amb = 298
    T_mix = T_ox *(OF/(OF+1)) + T_fuel *(1/(OF+1))

    # Propellant Pressures [Pa]
    P_wsr = Pc * PSI2PA
    P_amb = 14.7 * PSI2PA

    # Mass Flow Rates [kg/s]
    mdot_total = mdot * LBM2KG

    # Gas Compositions #
    # Pre-Mixed Inlet
    Y_mix = {fuel: (1/(OF+1)), ox: (OF/(OF+1))}
    sol_inlet = ct.Solution(mech)
    sol_inlet.TPY = T_mix, P_wsr, Y_mix

    # Arbitrary Exhaust 
    sol_exhaust = ct.Solution(mech)
    sol_exhaust.TP = T_amb, P_amb

############### MAIN LOOP ###############
# Find Minimum Ignition Energy Case

    # MIE upper/lower bound guesses [J/kg]
    e_low = 0
    e_high = 8e6

    # Bisecton Search to find MIE #
    tic = time.time()
    if not ignites(e_high, residence_time, sol_inlet, T_comb, mdot_total, sol_exhaust, mech, prop_names, save_path, plotFlag=False):
        raise RuntimeError("Upper bound too low — increase e_high")
    for xx in range(bisect_count):
        print(f"Bisection Progress: {(xx+1)/bisect_count * 100:0.2f}%")
        e_mid = 0.5 * (e_low + e_high)

        if ignites(e_mid, residence_time, sol_inlet, T_comb, mdot_total, sol_exhaust, mech, prop_names, save_path, plotFlag=False):
            e_high = e_mid
        else:
            e_low = e_mid

        e_min = e_high
    toc = time.time()
    print(f"Bisection Search Completed in {toc-tic:0.2f}s\n")
    
    # Once Minimum Ignition Energy Case is Found
    ## Step Through Solution to Get Temporal Data and Plots ## 
    ignites(e_min, residence_time, sol_inlet, T_comb, mdot_total, sol_exhaust, mech, prop_names, save_path, plotFlag=True)
    
    ## Return Auto Ignition Temperature of WSR 
    inlet_ignition = add_energy(e_min, mech, T_mix, P_wsr, Y_mix)
    T_ignition = inlet_ignition.T
    
    print(f'Inlet State: \n T [K]  {T_mix:0.3f}, P [psi]  {Pc:0.3f}')
    print(f'MIE [kJ/kg]  {e_min/1000:0.3f}')
    print(f'Autoignition Mixture State: \n T [K] {T_ignition:0.3f},  P [psi]  {Pc:0.0f}')
    return T_ignition

    """ NEED TO ADD NEW FUNCTION CALLS IN THE MAIN FUNCTION"""
###################################
    