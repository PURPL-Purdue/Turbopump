import cea

def cea_setup(reac_names, oxidant_weights, fuel_weights, of_ratio, reac_temps, Pc, Prat, ac_at):
    reac = cea.Mixture(reac_names)
    prod = cea.Mixture(reac_names, products_from_reactants=True)

    solver = cea.RocketSolver(prod, reactants=reac)
    solution = cea.RocketSolution(solver)

    weights = reac.of_ratio_to_weights(oxidant_weights, fuel_weights, of_ratio) #Convert OF to weights.

    #Compute chamber enthalpy. Normalized.
    hc = reac.calc_property(cea.ENTHALPY, weights, reac_temps)/cea.R
    
    #Solve the rocket problem for given inputs. Normalized.
    solver.solve(solution, weights, Pc, Prat, ac_at=ac_at, iac=True, hc=hc)

    return solution