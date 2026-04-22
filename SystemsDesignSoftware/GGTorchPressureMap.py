import numpy as np
from rocketcea.cea_obj import CEA_Obj
from scipy.interpolate import LinearNDInterpolator

# Unit conversions
psi_to_pa = 6894.76
ft_to_m = 0.3048

# Oxygen gas constant
R_UNIV = 8.314462618
M_O2 = 0.031998 # kg/mol
M_H2 = 0.002016
R_O2 = R_UNIV / M_O2  # J/(kg*K)
R_H2 = R_UNIV / M_H2  # J/(kg*K)

def finite_minmax(a):
    a = np.asarray(a, dtype=float)
    m = np.isfinite(a)
    if not np.any(m):
        raise ValueError("Array has no finite values.")
    return float(np.min(a[m])), float(np.max(a[m]))

def required_gox_feed_pressure_orifice(
    m_dot, CdA, T, gamma, p_chamber, R=R_O2,
    max_iter=80, rtol=1e-8
):
    """
    Solve for upstream/feed pressure p0 (Pa) given:
      m_dot: kg/s
      CdA: m^2 (Cd*A lumped)
      T: K (stagnation temperature upstream)
      gamma: -
      p_chamber: Pa (downstream static pressure)
      R: J/(kg*K)

    Uses compressible orifice flow (ideal gas), automatically handling choked/unchoked.
    """

    if m_dot <= 0:
        return float(p_chamber)

    # Critical pressure ratio for choking (downstream/upstream)
    pr_crit = (2.0/(gamma + 1.0)) ** (gamma/(gamma - 1.0))

    # Choked mass flow coefficient
    Cg = np.sqrt(gamma/(R*T)) * (2.0/(gamma + 1.0)) ** ((gamma + 1.0)/(2.0*(gamma - 1.0)))

    # Mass flow as a function of upstream pressure p0
    def mdot_of_p0(p0):
        # enforce p0 > p_chamber
        if p0 <= p_chamber:
            return 0.0

        pr = p_chamber / p0  # downstream/upstream

        if pr <= pr_crit:
            # choked
            return CdA * p0 * Cg

        # unchoked (subsonic compressible)
        term = (pr ** (2.0/gamma) - pr ** ((gamma + 1.0)/gamma))
        if term <= 0.0:
            return 0.0

        return CdA * p0 * np.sqrt((2.0*gamma/(R*T*(gamma - 1.0))) * term)

    # --- Bracket the solution ---
    p_lo = p_chamber * 1.000001
    f_lo = mdot_of_p0(p_lo) - m_dot  # should be negative

    # start with a reasonable high guess: choked inversion
    p_hi = max(p_chamber * 1.1, m_dot/(CdA*Cg))
    f_hi = mdot_of_p0(p_hi) - m_dot

    # expand until we bracket (f_lo < 0 < f_hi)
    expand = 0
    while f_hi < 0 and expand < 60:
        p_hi *= 1.5
        f_hi = mdot_of_p0(p_hi) - m_dot
        expand += 1

    if f_hi < 0:
        raise RuntimeError("Could not bracket solution for feed pressure. Check CdA/T/gamma/m_dot.")

    # --- Bisection (monotonic, guaranteed) ---
    for _ in range(max_iter):
        p_mid = 0.5*(p_lo + p_hi)
        f_mid = mdot_of_p0(p_mid) - m_dot

        if abs(f_mid) <= rtol * m_dot:
            return float(p_mid)

        if f_mid > 0:
            p_hi = p_mid
        else:
            p_lo = p_mid

    return float(0.5*(p_lo + p_hi))

def build_inverse_interpolators(fu_p_map, ox_p_map, pc_axis, of_axis, fill_value=np.nan):
    """
    Builds inverse interpolators:
      (p_fu, p_ox) -> pc
      (p_fu, p_ox) -> OF
    """
    PC, OF = np.meshgrid(pc_axis, of_axis, indexing="ij")

    pf = fu_p_map.ravel()
    po = ox_p_map.ravel()
    pc_vals = PC.ravel()
    of_vals = OF.ravel()

    mask = np.isfinite(pf) & np.isfinite(po) & np.isfinite(pc_vals) & np.isfinite(of_vals)
    pts = np.column_stack([pf[mask], po[mask]])

    pc_interp = LinearNDInterpolator(pts, pc_vals[mask], fill_value=fill_value)
    of_interp = LinearNDInterpolator(pts, of_vals[mask], fill_value=fill_value)
    return pc_interp, of_interp


def pressure_map(minOF, maxOF, maxPc, resolution, cstar_eff,
                 throat_area, fuel_CdA, ox_CdA,
                 rho_fuel=800.0):
    
    engine = CEA_Obj(oxName="GOX", fuelName="H2")

    pc_scale = np.linspace(1.0, maxPc, resolution) # psi
    OF_scale = np.linspace(minOF, maxOF, resolution)  # n/a

    fu_p_map = np.full((resolution, resolution), np.nan, dtype=float)
    ox_p_map = np.full((resolution, resolution), np.nan, dtype=float)

    for i, pc in enumerate(pc_scale):
        for j, OF in enumerate(OF_scale):
            cstar = engine.get_Cstar(Pc=pc, MR=OF)  # RocketCEA returns ft/s
            p_c_pa = pc * psi_to_pa
            cstar_real = cstar * cstar_eff

            # mdot_total = Pc * At / c*  (convert c* to m/s by *ft_to_m)
            m_dot_total = p_c_pa * throat_area / (cstar_real * ft_to_m)
            m_dot_fu = m_dot_total / (1.0 + OF)
            m_dot_ox = m_dot_total - m_dot_fu

            # Fuel feed pressure: incompressible orifice model
            dp_fu = (m_dot_fu / fuel_CdA) ** 2 / (2.0 * rho_fuel)
            fu_p_map[i, j] = required_gox_feed_pressure_orifice(m_dot_fu, fuel_CdA, T=293.15, gamma=1.41, p_chamber=p_c_pa, R=R_H2, max_iter=80, rtol=1e-8) / psi_to_pa

            # Ox feed pressure (compressible flow)
            ox_p_map[i, j] = required_gox_feed_pressure_orifice(m_dot_ox, ox_CdA, T=293.15, gamma=1.4, p_chamber=p_c_pa, R=R_O2, max_iter=80, rtol=1e-8) / psi_to_pa

    print("m_dot_fu:", m_dot_fu)
    print("m_dot_ox:", m_dot_ox)
    print("Pc:", pc)
    print("cstar: ", cstar)
    print("R_H2:", R_H2)
    print("fuel_CdA: ", fuel_CdA)
    print("ox_CdA: ", ox_CdA)

    # Build inverse interpolators
    pc_interp, of_interp = build_inverse_interpolators(fu_p_map, ox_p_map, pc_scale, OF_scale)

    # Auto pressure grid for inverse maps (same resolution)
    pf_min, pf_max = finite_minmax(fu_p_map)
    po_min, po_max = finite_minmax(ox_p_map)
    pf_axis = np.linspace(pf_min, pf_max, resolution)
    po_axis = np.linspace(po_min, po_max, resolution)
    PF, PO = np.meshgrid(pf_axis, po_axis, indexing="ij")

    PC_pfpo = pc_interp(PF, PO)  # psi
    OF_pfpo = of_interp(PF, PO)  # -

    # Removes duplicates from the matrix inversion process
    bad = (PC_pfpo >= 0.999 * maxPc) | (OF_pfpo < 1.001 * minOF)
    PC_pfpo[bad] = np.nan
    OF_pfpo[bad] = np.nan

    return pc_scale, OF_scale, fu_p_map, ox_p_map, PC_pfpo, OF_pfpo