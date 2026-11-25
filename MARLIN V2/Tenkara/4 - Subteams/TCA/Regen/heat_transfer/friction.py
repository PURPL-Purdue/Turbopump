import math


def konakov_smooth_darcy(Re):
    """
    Darcy friction factor for turbulent flow in a smooth pipe, Konakov:
        xi = (1.8 log10(Re) - 1.5)^-2
    VDI Heat Atlas (2010) G1, Eq. 27. Used inside the Gnielinski correlation
    (heat_transfer/coolant_side.py), not for the pressure drop: the correlation
    calls for the smooth-pipe value.
    """
    return (1.8 * math.log10(Re) - 1.5)**-2


def friction_factor(Re, eD=0.0):
    """
    Fanning friction factor. Smooth wall (eD = 0), Eq. 6.37:
        16/Re                      Re < 2100
        0.046/Re^0.2               5000 < Re < 200,000
        0.0014 + 0.125/Re^0.32     3000 < Re < 3e6
    The turbulent fits overlap: first one up to 2e5, second above; the
    2100-3000 transition gap uses the turbulent fit.
    Rough wall (eD = eps/D_h > 0): Moody diagram via Colebrook, as the book suggests.
    """
    if eD > 0:
        from fluids import friction_factor as darcy
        return darcy(Re=Re, eD=eD) / 4
    if Re < 2100:
        return 16 / Re
    if Re < 2e5:
        return 0.046 / Re**0.2
    return 0.0014 + 0.125 / Re**0.32
