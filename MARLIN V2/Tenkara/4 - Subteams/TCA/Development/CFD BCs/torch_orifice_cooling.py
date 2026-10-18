"""
Gnielinski HTC for GOX and GCH4 through an orifice (wall cooling).
All inputs are set in main(). No function has default values, so nothing is silently overwritten.
Not run. Requires CoolProp: pip install CoolProp
"""
import numpy as np
from CoolProp.CoolProp import PropsSI

IN2M   = 0.0254
LBM2KG = 0.453592
PSI2PA = 6894.757
R_UNIV = 8.314462618

def gnielinski_nu(Re, Pr, L_over_D):
    """Gnielinski (1976) with Petukhov friction factor.
    Valid ~ 3000 < Re < 5e6, 0.5 < Pr < 2000.
    L_over_D: float applies entrance correction (1 + (D/L)^(2/3)); None skips it."""
    f = (0.79 * np.log(Re) - 1.64) ** -2              # Darcy friction factor
    Nu = (f / 8) * (Re - 1000) * Pr / (1 + 12.7 * np.sqrt(f / 8) * (Pr ** (2 / 3) - 1))
    if L_over_D is not None:
        Nu *= 1 + (1 / L_over_D) ** (2 / 3)
    return Nu, f

def choked_mdot(P0, T0, gamma, R, A, Cd):
    """Choked mass flow [kg/s] through area A with discharge coefficient Cd."""
    return Cd * A * P0 * np.sqrt(gamma / (R * T0)) * (2 / (gamma + 1)) ** ((gamma + 1) / (2 * (gamma - 1)))

def htc(fluid, P0_psi, Pc_psi, mdot_lbm, d_in, n_orifices, T0, Cd, L_over_D, T_wall):
    """
    fluid       : CoolProp fluid name
    P0_psi      : upstream pressure [psia]
    Pc_psi      : chamber (downstream) pressure [psia]
    mdot_lbm    : total propellant flow [lbm/s]
    d_in        : orifice diameter [in]
    n_orifices  : number of orifices sharing mdot
    T0          : upstream temperature [K]
    Cd          : discharge coefficient (used only for the choked-flow check)
    L_over_D    : orifice length/diameter for entrance correction, or None
    T_wall      : wall temperature [K] for (T_b/T_w)^0.45 correction, or None
    """
    D    = d_in * IN2M
    A    = np.pi * D ** 2 / 4
    P0   = P0_psi * PSI2PA          # inputs are treated as psia
    Pc   = Pc_psi * PSI2PA
    mdot = mdot_lbm * LBM2KG / n_orifices      # per orifice [kg/s]

    # upstream state -> gamma, R
    cp = PropsSI("Cpmass", "T", T0, "P", P0, fluid)
    cv = PropsSI("Cvmass", "T", T0, "P", P0, fluid)
    gamma = cp / cv
    R = R_UNIV / PropsSI("molar_mass", fluid)

    # orifice-state estimate (isentropic): choked if P ratio below critical
    Pcrit = P0 * (2 / (gamma + 1)) ** (gamma / (gamma - 1))
    P_th  = max(Pc, Pcrit)
    T_th  = T0 * (P_th / P0) ** ((gamma - 1) / gamma)

    # properties at that state
    mu = PropsSI("viscosity",    "T", T_th, "P", P_th, fluid)
    k  = PropsSI("conductivity", "T", T_th, "P", P_th, fluid)
    cp = PropsSI("Cpmass",       "T", T_th, "P", P_th, fluid)
    Pr = cp * mu / k

    Re = 4 * mdot / (np.pi * D * mu)
    Nu, f = gnielinski_nu(Re, Pr, L_over_D)
    if T_wall is not None:
        Nu *= (T_th / T_wall) ** 0.45
    h = Nu * k / D

    mdot_choke = choked_mdot(P0, T0, gamma, R, A, Cd)
    print(f"--- {fluid} ---")
    print(f"gamma={gamma:.3f}  choked={Pc < Pcrit}  P_th={P_th/PSI2PA:.1f} psi  T_th={T_th:.1f} K")
    print(f"mdot/orifice={mdot:.5f} kg/s  (choked limit @Cd={Cd}: {mdot_choke:.5f} kg/s per orifice)")
    if mdot > mdot_choke:
        print("WARNING: flow per orifice exceeds choked limit; increase n_orifices")
    print(f"Re={Re:.0f}  Pr={Pr:.3f}  f={f:.4f}  Nu={Nu:.1f}")
    if not (3000 < Re < 5e6):
        print("WARNING: Re outside Gnielinski validity range")
    print(f"h = {h:.0f} W/m^2-K  ({h*0.17611:.1f} BTU/hr-ft^2-F)")
    return h

def main():
    # ---- common inputs ----
    d_in     = 0.063      # orifice diameter [in]
    Pc_psi   = 200.0      # chamber pressure [psia]
    T0       = 293.15     # upstream temperature [K] (room temp, edit if needed)
    Cd       = 0.9        # discharge coefficient
    L_over_D = None       # orifice L/D for entrance correction, or None
    T_wall   = None       # wall temp [K] for property correction, or None

    # ---- per-propellant inputs ----
    ox = dict(fluid="Oxygen",  P0_psi=505.1, mdot_lbm=0.0470, n_orifices=1)  # n_orifices: placeholder, set to actual
    fu = dict(fluid="Methane", P0_psi=490.8, mdot_lbm=0.0313, n_orifices=1)  # n_orifices: placeholder, set to actual

    for p in (ox, fu):
        htc(fluid=p["fluid"], P0_psi=p["P0_psi"], Pc_psi=Pc_psi,
            mdot_lbm=p["mdot_lbm"], d_in=d_in, n_orifices=p["n_orifices"],
            T0=T0, Cd=Cd, L_over_D=L_over_D, T_wall=T_wall)

if __name__ == "__main__":
    main()