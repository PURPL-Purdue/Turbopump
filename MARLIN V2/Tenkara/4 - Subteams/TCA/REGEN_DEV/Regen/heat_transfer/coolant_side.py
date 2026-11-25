"""
Coolant-side heat transfer coefficient h_l. Heister et al., 6.3 (Eq. 6.17).

Every correlation has the same signature:
    h_l = f(T_l, T_wl, mdot_chan, A_c, D_h, fluid, x_in=inf)
        T_l    coolant bulk temperature [K]
        T_wl   coolant-side wall temperature [K]
        mdot_chan  mass flow per channel [kg/s]
        A_c, D_h   channel flow area [m^2] and hydraulic diameter [m]
        fluid  coolant property module (mu, k, cp, rho as functions of T)
        x_in   distance along the channel from the coolant inlet [m]; inf =
               fully developed (entrance effects off). Correlations without an
               entrance term ignore it.

To add one: write the function, add it to CORRELATIONS, give it a row in
TRACEABILITY.md, then pick it in regen.yaml (jacket.coolant_htc).
"""
import math

from heat_transfer.friction import konakov_smooth_darcy


def sieder_tate(T_l, T_wl, mdot_chan, A_c, D_h, fluid, x_in=math.inf):
    """
    Sieder-Tate form, Eq. 6.17:  Nu = h_l D_h / k = a Re^m Pr^n (mu / mu_w)^b

    a, m, n from Table 6.3 - IPA is not tabulated; kerosene
    (0.023, 0.8, 0.4) is the closest. b = 0.114 (book p. 197). Re and Pr at the mean film temperature
    (T_l + T_wl)/2 (p. 198), mu_w at T_wl. Re = m_dot D_h / (A mu).
    Valid for Re > 10,000, 0.7 < Pr < 160. No entrance term: x_in is ignored.
    """
    a, m, n, b = 0.023, 0.8, 0.4, 0.114
    T_f = 0.5 * (T_l + T_wl)   # mean film temperature for Re and Pr
    mu, k, cp = fluid.mu(T_f), fluid.k(T_f), fluid.cp(T_f)
    Re = mdot_chan * D_h / (A_c * mu)
    Pr = mu * cp / k
    Nu = a * Re**m * Pr**n * (mu / fluid.mu(T_wl))**b
    return Nu * k / D_h

RE_LAM, RE_TURB = 2300.0, 1e4     # VDI G1 Eq. 30: the transition range


def gnielinski_nu(Re, Pr, Pr_w, d_over_x=0.0):
    """
    Turbulent local Nusselt number, VDI Heat Atlas (2010) chapter G1 (V. Gnielinski):

        Eq. 28  Nu = (xi/8) Re Pr / (1 + 12.7 sqrt(xi/8) (Pr^(2/3) - 1))
                     * [1 + (1/3) (d_i/x)^(2/3)]
        Eq. 27  xi = (1.8 log10 Re - 1.5)^-2      smooth pipe, Konakov (Darcy)
        Eq. 40  Nu_liquid = Nu (Pr / Pr_w)^0.11   liquids, 0.1 < Pr/Pr_w < 10

    d_over_x = D_h / x, x from the coolant inlet. 0 = fully developed: the
    bracket is 1 and this is Eq. 26 without its length factor.
    Valid 1e4 <= Re <= 1e6, 0.1 <= Pr <= 1000.
    """
    xi = konakov_smooth_darcy(Re)
    Nu = (xi / 8) * Re * Pr / (1 + 12.7 * math.sqrt(xi / 8) * (Pr**(2 / 3) - 1))
    Nu *= 1 + d_over_x**(2 / 3) / 3                                     # Eq. 28
    return Nu * (Pr / Pr_w)**0.11


def laminar_nu(Re, Pr, d_over_x=0.0):
    """
    Laminar local Nusselt number, constant heat flux, velocity and temperature
    developing together. VDI G1 3.2.2:

        Eq. 14  Nu_1 = 4.354                               fully developed
        Eq. 15  Nu_2 = 1.302 (Re Pr d_i/x)^(1/3)           thermally developing
        Eq. 22  Nu_3 = 0.462 Pr^(1/3) (Re d_i/x)^(1/2)     both developing
        Eq. 23  Nu   = (Nu_1^3 + (Nu_2 - 1)^3 + Nu_3^3)^(1/3)   (as printed)

    Valid 0.7 < Pr < 1000. d_over_x = 0 (fully developed) returns Eq. 14
    directly; Eq. 23 as printed gives 4.335 there (see TRACEABILITY).
    No wall-property correction: G1 3.2.3 gives none for constant heat flux.
    """
    if d_over_x == 0:
        return 4.354                                                    # Eq. 14
    nu2 = 1.302 * (Re * Pr * d_over_x)**(1 / 3)                         # Eq. 15
    nu3 = 0.462 * Pr**(1 / 3) * (Re * d_over_x)**0.5                    # Eq. 22
    return (4.354**3 + (nu2 - 1)**3 + nu3**3)**(1 / 3)                  # Eq. 23


def vdi_nu(Re, Pr, Pr_w, d_over_x=0.0):
    """
    Nusselt number for any Re, VDI G1 4.2:

        Re >= 1e4          turbulent: gnielinski_nu (Eqs. 27, 28, 40)
        Re <= 2300         laminar: laminar_nu at Re (Eq. 23)
        2300 < Re < 1e4    Eq. 29  Nu = (1 - g) Nu_lam,2300 + g Nu_turb,1e4
                           Eq. 30  g  = (Re - 2300) / (1e4 - 2300)

    G1 writes Eq. 29 with mean (whole-tube) values. Here both ends are local
    values at this station's d_i/x: our step, see TRACEABILITY.
    """
    if Re >= RE_TURB:
        return gnielinski_nu(Re, Pr, Pr_w, d_over_x)
    if Re <= RE_LAM:
        return laminar_nu(Re, Pr, d_over_x)
    g = (Re - RE_LAM) / (RE_TURB - RE_LAM)                              # Eq. 30
    return ((1 - g) * laminar_nu(RE_LAM, Pr, d_over_x)
            + g * gnielinski_nu(RE_TURB, Pr, Pr_w, d_over_x))           # Eq. 29


def gnielinski(T_l, T_wl, mdot_chan, A_c, D_h, fluid, x_in=math.inf):
    """
    Gnielinski h_l with VDI's laminar and transition ranges (vdi_nu). Re, Pr, k
    at the bulk coolant temperature T_l (G1 4.1: the mean fluid temperature; in
    a station march that's the local bulk temperature). Pr_w at the wall, T_wl.
    x_in = inf: fully developed everywhere (no entrance terms).
    """
    mu, k, cp = fluid.mu(T_l), fluid.k(T_l), fluid.cp(T_l)
    Re = mdot_chan * D_h / (A_c * mu)
    Pr = mu * cp / k
    Pr_w = fluid.mu(T_wl) * fluid.cp(T_wl) / fluid.k(T_wl)
    return vdi_nu(Re, Pr, Pr_w, D_h / x_in) * k / D_h


# name in regen.yaml -> (function, lowest Re it's valid for, for the report's check)
CORRELATIONS = {
    "sieder_tate": (sieder_tate, 1e4),
    "gnielinski": (gnielinski, 0.0),      # laminar and transition covered (vdi_nu)
}


def get(name):
    """The correlation function registered under `name`."""
    try:
        return CORRELATIONS[name][0]
    except KeyError:
        raise ValueError(f"Unknown coolant_htc '{name}'. Available: {', '.join(CORRELATIONS)}")


def Re_min(name):
    """Lowest Reynolds number the named correlation is valid for."""
    return CORRELATIONS[name][1]
