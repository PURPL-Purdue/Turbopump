"""
Fin conduction through the lands - Blong et al. (2024), "Investigation of Heat
Transfer in a Water-Cooled GOx/RP-1 RDRE", Eqs. 5 and 21-27 (straight fin,
Incropera fin theory).

The coolant touches the whole channel perimeter, not just the floor. Heat
reaching the side walls and the top (closeout) has to conduct up the land
first, so that area counts at the fin efficiency eta. Blong's fin is the land
plus half the closeout: up the land (h_c), then along the closeout to the
channel centreline (w_c / 2), a symmetry line, so the tip is adiabatic.

    Eq. 5   w_fin = 2 pi (r + t_w + h_c/2) / N - w_c     land width at mid-height
    Eq. 22  A_wl  = dx w_c                               channel floor
    Eq. 23  A_wg  = dx (w_c + w_fin)                     hot-wall area per channel
    Eq. 24  A_fin = 2 dx L_fin                           both fins
    Eq. 25  eta   = tanh(m L_fin) / (m L_fin)            efficiency of the fin compared to base surface
            m     = sqrt(2 h_l / (k_w w_fin))            NOT Blong Eq. 26, see below
    Eq. 27  L_fin = h_c + w_c / 2
    Eq. 21  q_l   = h_l (T_wl - T_l) (eta A_fin + A_wl) / A_wg

m: Blong's Eq. 26, m = sqrt(2 h_l (w_fin + dx) / (k_w dx w_fin)), uses the
perimeter 2 (w_fin + dx) of one dx-long slice of land, which counts its front
and back faces as wetted. They're cuts through a continuous land, not coolant
surface, and the term makes the result depend on the station spacing (eta -> 0
as dx -> 0; checked Nov 24: 100 -> 3000 stations took dp from 58 to 257 psi).
Only the two side faces are wetted: P = 2 dx, A_c = w_fin dx, so
m = sqrt(h P / (k A_c)) = sqrt(2 h_l / (k_w w_fin)), the standard straight
fin (Incropera / Bergman, the fin theory Blong cites).

fin_factor() returns phi = (eta A_fin + A_wl) / A_wg, so q_l = h_l phi (T_wl - T_l):
phi multiplies the correlation's h_l. dx cancels in phi.
Assumes (Blong): 1D conduction along the fin, constant k_w, adiabatic tip, the
same h_l on every wetted face.

Pure physics: no config reads. Switched on and off with jacket.fin_conduction
in regen.yaml (jacket/march.coolant_area_factor).
"""
import math


def fin_width(r, t_w, h, N, w):
    """Eq. 5: land width at the channel's mid-height [m]. r is the hot-gas wall radius."""
    return 2 * math.pi * (r + t_w + 0.5 * h) / N - w


def fin_efficiency(h_l, k_w, w_fin, L_fin):
    """Eq. 25: adiabatic-tip straight fin efficiency, 0 < eta <= 1. m without
    Blong's dx term (module docstring)."""
    m = math.sqrt(2 * h_l / (k_w * w_fin))
    mL = m * L_fin
    return math.tanh(mL) / mL


def fin_factor(h_l, w, h, r, t_w, N, k_w, dx):
    """
    Eqs. 5, 21-25, 27 at one station. Returns (phi, eta).
        h_l   coolant-side coefficient from the correlation [W/m^2-K]
        w, h  channel width and height [m]
        r     hot-gas wall radius at the station [m]
        t_w   hot wall thickness [m]
        N     number of channels
        k_w   wall conductivity [W/m-K]
        dx    segment length [m] (Blong's dx; the march's ds). Cancels in phi.
    """
    w_fin = fin_width(r, t_w, h, N, w)
    L_fin = h + 0.5 * w                                   # Eq. 27
    eta = fin_efficiency(h_l, k_w, w_fin, L_fin)          # Eq. 25, m without dx
    A_wl = dx * w                                         # Eq. 22
    A_wg = dx * (w + w_fin)                               # Eq. 23
    A_fin = 2 * dx * L_fin                                # Eq. 24
    return (eta * A_fin + A_wl) / A_wg, eta               # Eq. 21
