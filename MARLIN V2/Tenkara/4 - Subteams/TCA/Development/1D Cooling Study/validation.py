
import pint 
import cea 
import numpy as np
import CoolProp as CoolProp

from cea_setup import cea_setup

ureg = pint.UnitRegistry()

L = 3 * ureg.mm
k = 30 *ureg.W / (ureg.cm * ureg.K)

Dt = 12 * ureg.cm
Pc = 65 * ureg.bar

reac_names = ["H2", "O2(L)"]
reac_temps = np.array([298, 90.17]) * ureg.K
fuel_weights = np.array([1.0, 0.0])
oxidant_weights = np.array([0.0, 1.0])
cont_ratio = 1.8
of_ratio = 5.5
Prat = 1.7384

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
    print(f"Iteration {i} | Wall Temp: {Tw}")
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


    # Internal flow of LH2 in channels
    MaCool = 0.1
    Tcool = 100 * ureg.K
    Pcool = 10 * ureg.MPa
    Dh = 4 * ureg.mm

    RhoCool = CoolProp.CoolProp.PropsSI("D", "T", Tcool.magnitude, "P", Pcool.to(ureg.Pa).magnitude, "H2") * ureg.kg / ureg.m**3
    KCool = CoolProp.CoolProp.PropsSI("L", "T", Tcool.magnitude, "P", Pcool.to(ureg.Pa).magnitude, "H2") * ureg.W /(ureg.m * ureg.K)
    MuCool = CoolProp.CoolProp.PropsSI("V", "T", Tcool.magnitude, "P", Pcool.to(ureg.Pa).magnitude, "H2") * ureg.Pa * ureg.s
    sonVelCool = CoolProp.CoolProp.PropsSI("A", "T", Tcool.magnitude, "P", Pcool.to(ureg.Pa).magnitude, "H2") * ureg.m / ureg.s
    Vcool = sonVelCool * MaCool

    ReCool = (RhoCool*Vcool*Dh / MuCool).to_base_units()
    PrCool = CoolProp.CoolProp.PropsSI("PRANDTL", "T", Tcool.magnitude, "P", Pcool.to(ureg.Pa).magnitude, "H2")
    hl = ((KCool/Dh) * (0.026*ReCool**0.8 * PrCool**0.4)).to(ureg.W / (ureg.m**2 * ureg.K)) # McAdams convection coorelation


    gam = 1.1655
    r = PrInf**(1/3)
    Tr = (1 + (gam-1)/2*Ma[1]**2*r)*Tinf
    qdot = (Tr - Tcool) / (1/hg + (L/k).to_base_units() + 1/hl)
    qdot.to(ureg.MW / ureg.m**2)
    Tw_new = Tr - qdot/hg
    
    error = np.abs((Tw_new - Tw)/Tw)
    Tw = Tw_new
    i += 1
