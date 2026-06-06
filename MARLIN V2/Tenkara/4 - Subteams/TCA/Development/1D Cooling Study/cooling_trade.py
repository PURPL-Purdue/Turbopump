
import pint 
import cea 
import numpy as np
import CoolProp as CoolProp

from cea_setup import cea_setup
import gamma_from_of

def steady_state_throat_temp(ureg, reac_names, reac_temps, of_ratio, Dt, Pc, L, k):

    fuel_weights = np.array([1.0, 0.0])
    oxidant_weights = np.array([0.0, 1.0])
    cont_ratio = 1.8
    Prat = (12 * ureg.psi / (1 * ureg.bar)).to(ureg.dimensionless).magnitude

    solution = cea_setup(reac_names, oxidant_weights, fuel_weights, of_ratio, reac_temps.magnitude, Pc.to(ureg.bar).magnitude, Prat, cont_ratio)
    num_pts = solution.num_pts
    T = solution.T * ureg.K
    P = solution.P * ureg.bar
    rho = solution.density * ureg.kg / ureg.m**3

    visc = solution.viscosity * ureg.millipoise
    Cp = solution.cp_eq * ureg.kJ / (ureg.kg * ureg.K)
    Pr = solution.Pr_eq
    Ma = solution.Mach
    sonVel = solution.sonic_velocity * ureg.m / ureg.s
    V = Ma * sonVel


    visc = visc.to(ureg.Pa * ureg.s)
    Cp = Cp.to(ureg.J / (ureg.kg * ureg.K))

    T0 = T[0]
    Tinf = T[1]
    Mu0 = visc[0]
    Muinf = visc[1]
    Cp0 = Cp[0]
    Cpinf = Cp[1]
    Rho0 = rho[0]
    Rhoinf = rho[1]
    Pr0 = Pr[0]
    PrInf = Pr[1]
    Vinf = V[1]


    # Rhoam / Rhoinf = Tinf / Tam | ideal gas relation
    Tw = 1300 * ureg.K # initial guess for wall temp at throat
    i = 1
    error = 1
    while (error > 1/1000):
        print(f"Iteration {i} | Wall Temp: {round(Tw.magnitude, 3)} K")
        Tam = (Tinf + Tw) / 2
        RhoamqRhoinf = Tinf / Tam


        # establish power law for viscosity Mu/Muref = (T/Tref)**w
        # use best fit for a range of inf/stagnation properties
        # cea with multiple contraction ratios and log fit is ideal

        # currently aproximating using frestream and stagnation for parity
        w = np.log(Muinf/Mu0) / np.log(Tinf / T0)
        # use this power law for Muam/Mu0
        MuamqMu0 = (Tam / T0)**w


        # bartz equation for hg
        hg = (0.026 / Dt**0.2) * ((Mu0**0.2 * Cp0) / Pr0**0.6) * (Rhoinf * Vinf) ** 0.8 *(RhoamqRhoinf)**0.8 * (MuamqMu0)**0.2
        hg = hg.to(ureg.W / (ureg.m**2 * ureg.K))


        # Internal flow of coolant in channels
        #channel geometry
        AR = 0.8
        height = 0.2 * ureg.inch
        width = AR*height
        perim = 2*width + 2*height
        Acool = width * height
        Dh = 4*Acool/perim

        MaCool = 0.3
        Tcool = 150 * ureg.K
        Pcool = 850 * ureg.psi
        
        RhoCool = CoolProp.CoolProp.PropsSI("D", "T", Tcool.magnitude, "P", Pcool.to(ureg.Pa).magnitude, "H2") * ureg.kg / ureg.m**3
        KCool = CoolProp.CoolProp.PropsSI("L", "T", Tcool.magnitude, "P", Pcool.to(ureg.Pa).magnitude, "H2") * ureg.W /(ureg.m * ureg.K)
        MuCool = CoolProp.CoolProp.PropsSI("V", "T", Tcool.magnitude, "P", Pcool.to(ureg.Pa).magnitude, "H2") * ureg.Pa * ureg.s
        sonVelCool = CoolProp.CoolProp.PropsSI("A", "T", Tcool.magnitude, "P", Pcool.to(ureg.Pa).magnitude, "H2") * ureg.m / ureg.s
        Vcool = sonVelCool * MaCool

        ReCool = (RhoCool*Vcool*Dh / MuCool).to_base_units()
        PrCool = CoolProp.CoolProp.PropsSI("PRANDTL", "T", Tcool.magnitude, "P", Pcool.to(ureg.Pa).magnitude, "H2")
        hl = ((KCool/Dh) * (0.026*ReCool**0.8 * PrCool**0.4)).to(ureg.W / (ureg.m**2 * ureg.K)) # McAdams convection coorelation


        gam = gamma_from_of.gamma_lookup("1D Cooling Study\Gamma_Lookup.xlsx", of_ratio.round(1), reac_names)
        r = PrInf**(1/3)
        Tr = (1 + (gam-1)/2*Ma[1]**2*r)*Tinf
        qdot = (Tr - Tcool) / (1/hg + (L/k).to_base_units() + 1/hl)
        qdot.to(ureg.MW / ureg.m**2)
        Tw_new = Tr - qdot/hg
        
        error = np.abs((Tw_new - Tw)/Tw)
        Tw = Tw_new
        i += 1


    mdotCool = RhoCool*Vcool*Acool
    print(f"Dh: {round(Dh.magnitude, 2)} {Dh.units}")
    print(f"Coolant mdot: {round(mdotCool.to(ureg.kg / ureg.s).magnitude, 3)} kg/s")
    return Tw

## channel size for a temp, just sweep some params, maybe heat map. 
# stay time code for L*
# 
def combustion_temp_sweep(ureg, reac_names, reac_temps, Pc, of_range):

    fuel_weights = np.array([1.0, 0.0])
    oxidant_weights = np.array([0.0, 1.0])
    cont_ratio = 1.8
    Prat = (12 * ureg.psi / (1 * ureg.bar)).to(ureg.dimensionless).magnitude
    ofs = []
    temps = []


    for of_ratio in of_range:

        solution = cea_setup(reac_names, oxidant_weights, fuel_weights, of_ratio, reac_temps.magnitude, Pc.to(ureg.bar).magnitude, Prat, cont_ratio)
        Tc = solution.T[0]

        ofs.append(of_ratio)
        temps.append(Tc)

    tempMax = np.max(temps)
    ofMax = ofs[np.argmax(temps)]
    print(f"Max Tc: {tempMax:0.3f} at OF: {ofMax:0.2f}")
    return temps, ofs, tempMax, ofMax

def main():
    ureg = pint.UnitRegistry()

    L = 3 * ureg.mm # chamber wall thickness
    k = 30 * ureg.W / (ureg.cm * ureg.K) #chamber wall conductivity

    Dt = 2.3838 * ureg.inch
    Pc = 800 * ureg.psi

    reac_temps = np.array([298.0, 90.17]) * ureg.K

    of_range = np.arange(0.5,8,0.01)

    print("AT MAX TEMP")
    print("IPA")
    ipa_reac_names = ["C3H8O,2propanol", "O2(L)"]
    ipa_Tcs, ipa_ofs, ipa_tempMax, ipa_ofMax = combustion_temp_sweep(ureg, ipa_reac_names, reac_temps, Pc, of_range)
    steady_state_throat_temp(ureg, ipa_reac_names, reac_temps, ipa_ofMax, Dt, Pc, L, k)

    print("Ethanol")
    eth_reac_names = ["C2H5OH(L)","O2(L)"] 
    eth_Tcs, eth_ofs, eth_tempMax, eth_ofMax = combustion_temp_sweep(ureg, eth_reac_names, reac_temps, Pc, of_range)
    steady_state_throat_temp(ureg, eth_reac_names, reac_temps, eth_ofMax, Dt, Pc, L, k)

    print("AT SAME TEMP")
    eth_idx = np.argmin(np.abs(np.array(eth_Tcs) - 3000))
    ipa_idx = np.argmin(np.abs(np.array(ipa_Tcs) - 3000))

    eth_temp = eth_Tcs[eth_idx]
    eth_of = eth_ofs[eth_idx]

    ipa_temp = ipa_Tcs[ipa_idx]
    ipa_of = ipa_ofs[ipa_idx]

    

    print(f"IPA:     T = {ipa_temp:.2f} K at O/F = {ipa_of:.2f}")
    steady_state_throat_temp(ureg, ipa_reac_names, reac_temps, ipa_of, Dt, Pc, L, k)

    print(f"Ethanol: T = {eth_temp:.2f} K at O/F = {eth_of:.2f}")
    steady_state_throat_temp(ureg, eth_reac_names, reac_temps, eth_of, Dt, Pc, L, k)

    return

if __name__ == "__main__":
    main()