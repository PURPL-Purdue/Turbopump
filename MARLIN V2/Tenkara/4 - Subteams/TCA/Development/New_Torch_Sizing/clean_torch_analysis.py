import cea 
import numpy as np

IN2M = 0.0254
LBM2KG = 0.453592
BAR2PA = 100000
MAX_ITER = 50
PSI2PA = 6894.76

class CEAModel:
    """Builds the CEA object; call `solve(p_c, OF)` to retrieve solution."""
    def __init__(self, props):
        self.T_reactant = np.array([273, 273])
        self.fuel_weights = np.array([1.0, 0.0])
        self.ox_weights = np.array([0.0, 1.0])

        self.reac = cea.Mixture(props)
        self.prod = cea.Mixture(props, products_from_reactants=True)
        self.solver = cea.RocketSolver(self.prod, reactants=self.reac)
        self.solution = cea.RocketSolution(self.solver)

    def solve(self, p_c, OF):
        """Return solution element of cea solve given chamber pressure and OF"""
        weights = self.reac.of_ratio_to_weights(self.ox_weights, self.fuel_weights, of_ratio=OF)

        hc = self.reac.calc_property(cea.ENTHALPY, weights, self.T_reactant) / cea.R

        self.solver.solve(self.solution, weights, p_c / BAR2PA, hc=hc, iac=True)

        return self.solution

class TorchSizer:
    # takes all standard SI inputs
    def __init__(self,fuel, ox, d_f, d_ox, d_t, cd_f, cd_ox, cd_t, mdot, OF, Pc_guess, cstar_eff):
        self.props = [fuel, ox]
        self.d_f = d_f
        self.d_ox = d_ox
        self.d_t = d_t
        self.cd_f = cd_f
        self.cd_ox = cd_ox
        self.cd_t = cd_t
        self.mdot = mdot
        self.OF = OF
        self.Pc_guess = Pc_guess
        self.cstar_eff = cstar_eff

        
        #split mass flows
        self.mdot_ox = self.mdot * self.OF / (self.OF + 1)
        self.mdot_f = self.mdot * 1/(self.OF+1)

        # find areas
        self.A_ox = self.circle_area(d_ox)
        self.A_f = self.circle_area(d_f)
        self.A_t = self.circle_area(d_t)

        # set up CEA solve & converge on a nominal pressure
        self.cea_model = CEAModel(self.props)
        self.solve_pressure()
        
        return

    def solve_pressure(self):
        # converge on chamber pressure
        Pc_guess = self.Pc_guess 
        for _ in range(MAX_ITER):
            solution = self.cea_model.solve(Pc_guess, self.OF)
        
            #k_c = solution.gamma_s[0]
            #R_c = cea.R / solution.MW[0]
            #Tc = solution.T[0]
            #Pc_new = self.choked_pressure(self.mdot, self.A_t, k_c, Tc, R_c, Cd=self.cd_t)
            Pc_new = (solution.c_star[0] * self.cstar_eff * self.mdot) / (self.A_t * self.cd_t)
            print(Pc_new/PSI2PA, Pc_guess/PSI2PA)
            if abs((Pc_new - Pc_guess) / Pc_guess) < 1e-6:
                self.Pc = Pc_new
                print(f"Pressure Converged: {self.Pc/PSI2PA:0.3f} psi")
                return 
            Pc_guess = Pc_new
        raise(ValueError("Pressure did not converge."))
        # calculate upstream line pressures for choked flow
        #self.P_ox = (self.mdot_ox / (self.cd_ox * self.A_ox))**2 / (self.k_ox * self.rho_ox * (2 / (self.k_ox + 1))**((self.k_ox + 1)/(self.k_ox - 1))) 
        #self.P_f = (self.mdot_f / (self.cd_f * self.A_f))**2 / (self.k_f * self.rho_f * (2 / (self.k_f + 1))**((self.k_f + 1)/(self.k_f - 1))) 

    
    ## Helper functions ##
    
    def choked_pressure(self, mdot, area, k, T, R, Cd=1.0):
        """Upstream stagnation pressure [Pa] that chokes mdot through an orifice.

        mdot = Cd * A * P * sqrt(k / (R T)) * (2/(k+1))^((k+1)/(2(k-1)))
        """
        factor = (2 / (k + 1)) ** ((k + 1) / (2 * (k - 1)))
        return mdot / (Cd * area * np.sqrt(k / (R * T)) * factor)

    def circle_area(self, diameter):
        """Area [m^2] of a circle from diameter [m]."""
        return np.pi * diameter ** 2 / 4
    
def main():
    # define sizes and convert to metric
    
    d_f =  0.059 * IN2M
    d_ox =  0.063 * IN2M 
    d_t = 0.268 * IN2M
    mdot = 0.055 * LBM2KG

    torch = TorchSizer(
        fuel = "CH4",
        ox = "O2",
        d_f = d_f,
        d_ox = d_ox,
        d_t = d_t,
        cd_f = 0.8,
        cd_ox = 0.8,
        cd_t = 1,
        mdot = mdot,
        OF = 1.5,
        Pc_guess = 500 * PSI2PA,
        cstar_eff = 1
    )

    return

if __name__ == "__main__":
    main()