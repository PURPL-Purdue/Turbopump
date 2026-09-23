import cantera as ct
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import time
from pyfluids import Fluid, FluidsList, Input
from rocketcea.cea_obj import CEA_Obj
import warnings

warnings.filterwarnings("ignore")

PSI2PA = 6894.76 # psi to pa
R2K = 0.555556 # Rankine to Kelvin
BTU2J = 1055.06 # BTU to J
LBM2KG = 0.453592 # lbm to kg

## Experiment Set Up ##
OF = 1.5
residence_time = 0.6 # [s]

# Propellant Temperatures [K]
T_comb = 1000
T_amb = 298
T_ox = 90
T_fuel = 298
T_mix = T_ox *(OF/(OF+1)) + T_fuel *(1/(OF+1))

# Propellant Pressures [Pa]
P_wsr = 500 * PSI2PA
P_amb = 14.7 * PSI2PA

# Mass Flow Rates [kg/s]
mdot_total = 20 * LBM2KG * (1.420/100)
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

# function to get masss averaged chemical potential of mixture 
def mass_avg_chem_pot(mixture):
    mu = mixture.chemical_potentials      # J/kmol
    W  = mixture.molecular_weights        # kg/kmol
    Y  = mixture.Y                        # mass fractions
    return np.dot(Y, mu / W)              # J/kg mixture

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

def add_energy(e_ign):
    # Returns Clone with Added Enthalpy to Inlet Mixture, keep other properties the same
    inlet_clone = ct.Solution(mech)
    inlet_clone.TPY = T_mix, P_wsr, Y_mix
    dh = e_ign # [J/kg]
    inlet_clone.HP = inlet_clone.h + dh, inlet_clone.P

    return inlet_clone

def ignites(e_ign, plotFlag=False):
    inlet_clone = add_energy(e_ign)
    (sim, rctr_wsr) = define_network(inlet_clone)
    max_sim_time = 500 #[s]
    if plotFlag:
        time_history = ct.SolutionArray(sol_inlet, extra=["t"])
        step_solution(time_history, max_sim_time, rctr_wsr, sim)
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
        """(T,P) = steady_solve(sim,rctr_wsr)
        if T > T_comb:
            return True"""

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
    
    makeplot(time_history)

    return time_history

def steady_solve(sim, reactor):
    # Solve Sim to Steady State
    tic = time.time()
    sim.solve_steady()
    toc = time.time()
    print(f"Simulation Completed in {toc-tic:3.2f} seconds.")
    print(f"Final State: {reactor.T:0.2f} K, {reactor.thermo.P / PSI2PA :0.2f} psi")
    
    return reactor.T, reactor.thermo.P

def makeplot(time_history):
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
    return

def main():
    # Find Minimum Ignition Energy Case

    # MIE upper/lower bound guesses [J/kg]
    e_low = 0
    e_high = 8e6

    # Bisecton Search to find MIE #
    tic = time.time()
    if not ignites(e_high, plotFlag=False):
        raise RuntimeError("Upper bound too low — increase e_high")
    for xx in range(30):
        e_mid = 0.5 * (e_low + e_high)

        if ignites(e_mid):
            e_high = e_mid
        else:
            e_low = e_mid

        e_min = e_high
    toc = time.time()
    print(f"Bisection Search Completed in {toc-tic:0.2f}s\n")
    # Once Minimum Ignition Energy Case is Found
    ## Step Through Solution to Get Temporal Data and Plots ## 
    ignites(e_min, plotFlag=True)
    
    # Use MIE to compute Torch Mass Flow Rate (uses torch propellants)
    """ox = 'GOX'
    fuel = 'GH2'
    pc_torch = 300 * PSI2PA # Torch Chamber Pressure [Pa]
    OF_torch = 5
    cea = CEA_Obj(oxName=ox, fuelName=fuel)

    T_torch = cea.get_Tcomb(Pc=pc_torch/PSI2PA, MR=OF_torch) * R2K # Torch Combustion Temp [K]
    h_torch_set = cea.get_Enthalpies(Pc=pc_torch/PSI2PA, MR=OF_torch, eps=1)# get enthalpies (BTU/lbm) [chamber (reference), throat, exit(same as throat, eps=1)]
    h_rxn_torch = (h_torch_set[0] - h_torch_set[1]) * (BTU2J / LBM2KG) # get change in enthalpy due to reaction (BTU/lbm)-> [J/kg]"""
    
    h_rxn_torch = 21 * 1e6 #J/kg GOX/GH2
    mdot_torch = mdot_total * e_min / h_rxn_torch # mass flow of torch needed to ignite mixture [kg/s]

    ## Print Results
    pd.options.display.float_format = "{:.5f}".format
    df = pd.DataFrame(
    {
        "Simulation": [e_min/1000, mdot_torch] ,
    },
    index=["MIE [kJ/kg]", "mdot [kg/s]"]
)
    print()
    print(df)
    print()
    return

if __name__ == "__main__":
    main()