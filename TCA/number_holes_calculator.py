"""
This script calculates the number of holes depending on a target 
hole diameter for a bi-propellant injector plate. Change the initial
guesses for hole diameter as needed. 
Author: Joaquin Alarcon
"""

import numpy as np

psi_to_pa = 6894.76 #conversion factor from psi to Pa

Cd = 0.8 #discharge coefficient
delta_P = 350*psi_to_pa #pressure drop across injector plate [Pa]
OF_Ratio = 2.1

mdot = 9 #[kg/s]
rho_lox = 1141 #[kg/m3]
rho_rp1 = 810 #[kg/m3]
d_init_guess_lox = 0.001 #[m] (CHANGE AS NEEDED)
d_init_guess_rp1 = 0.001 #[m] (CHANGE AS NEEDED)

choice = int(input("Enter 1 for LOX calculation, 2 for RP-1 calculation: "))

mdot_kero = mdot/(1+OF_Ratio) #Fuel mass flow (RP-1), from global O/F
mdot_lox = mdot*OF_Ratio/(1+OF_Ratio) #Oxidizer mass flow (LOX), from global O/F

def mdot_per_hole(d, rho):
    A = np.pi*(d/2)**2
    return Cd*A*np.sqrt(2*rho*delta_P)

def d_for_Nholes(mdot_total, N, rho):
    mdot_hole = mdot_total / N
    A = mdot_hole / (Cd*np.sqrt(2*rho*delta_P))
    d = 2 * np.sqrt(A/np.pi)
    return d

def solve_N_and_d_near_target(mdot_total, rho, d_target, N_search = 200):
    #Estimate N from the target diameter
    mdot_hole_target = mdot_per_hole(d_target, rho)
    N0 = max(1, int(np.round(mdot_total/mdot_hole_target)))
    #search around N0 to find best closeness in diameter
    N_min = max(1, N0 - N_search)
    N_max = N0 + N_search

    best = None
    for N in range(N_min, N_max+1):
        d_req = d_for_Nholes(mdot_total, N, rho)
        error = abs(d_req - d_target)
        score = (error, N)
        if best is None or score < best["score"]:
            best = {
                "N": N,
                "d": d_req,
                "score": score
            }
    N_best = best["N"]
    d_best = best["d"]
    mdot_hole_best = mdot_total / N_best
    return N_best, d_best, mdot_hole_best

if choice == 2:
    N_rp1, d_rp1, mdot_hole_rp1 = solve_N_and_d_near_target(mdot_kero, rho_rp1, d_init_guess_rp1)
    print(f"RP-1 Holes: {N_rp1}, Diameter: {d_rp1*1000:.3f} mm")

if choice == 1:
    N_lox, d_lox, mdot_hole_lox = solve_N_and_d_near_target(mdot_lox, rho_lox, d_init_guess_lox)
    print(f"LOX Holes: {N_lox}, Diameter: {d_lox*1000:.3f} mm")