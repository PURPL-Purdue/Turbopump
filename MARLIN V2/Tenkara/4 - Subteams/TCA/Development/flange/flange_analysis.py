## Flange Analysis Code for Integral Flanges
# Companion to flange_sizing.py: that one searches for a valid design,
# this one checks an already-defined design and reports where it fails.

# Sources: ASME BPVC.VIII-1-2021 Mandatory Appendix 2

import numpy as np
from types import SimpleNamespace
from pathlib import Path

from flange_sizing import FlangeSizer, basic_seating_width, IN2M, PSI2PA, N2LBF
from refBolt import bolt_lookup
from margins import margins


def analyze_flange(B, N, G, C, A, bolt_name, bolt_count, t, g0, g1,
                    P, m, y, Sa, Sb, Sf, E, safetyFactor,
                    facing_sketch="1a", column="II", w=None, T=None,
                    thread_type="UNF"):
    """
    Checks an existing flange design against Appendix 2 2-7/2-8/2-14.
    Prints a pass/fail breakdown per check and makes a geometry plot.
    """

    # F, V, f hardcoded to g0/g1 = 1 values (Table 2-7.1); only valid for g1 == g0
    if not np.isclose(g0, g1):
        raise NotImplementedError(
            "g1 != g0: hardcoded F=0.908920, V=0.550103, f=1 are only "
            "valid for uniform-thickness hubs (g1 == g0). Implement the "
            "full Table 2-7.1 polynomial factors to support tapered hubs."
        )

    if not (B < G < C < A):
        print(f"WARNING: expected B < G < C < A, got B={B/IN2M:.3f} G={G/IN2M:.3f} C={C/IN2M:.3f} A={A/IN2M:.3f} in")

    # basic/effective gasket seating width, Table 2-5.2
    bo = basic_seating_width(N, facing_sketch, column, w, T)
    bo_mm = bo * 1000
    if bo_mm <= 6:
        b = bo
        gasket_ID = G - N
        gasket_OD = G + N
    else:
        Cb = 2.5
        b = Cb * np.sqrt(bo_mm) / 1000
        gasket_OD = G + 2 * b
        gasket_ID = gasket_OD - 2 * N

    # allowables actually used for sizing; raw yields kept for margin reporting
    Sa_used = Sa / safetyFactor
    Sb_used = Sb / safetyFactor
    Sf_used = Sf / safetyFactor

    # gasket loads
    H = 0.785 * G**2 * P
    Hp = 2 * b * np.pi * G * m * P
    Wm1 = H + Hp
    Wm2 = np.pi * b * G * y

    # required bolt area
    Am1 = Wm1 / Sb_used
    Am2 = Wm2 / Sa_used
    Am = max(Am1, Am2)
    seating_drives_area = Am2 > Am1

    # actual bolt properties
    boltArea = bolt_lookup("area", bolt_name, "name", thread_type) * IN2M**2
    boltDiam = bolt_lookup("diameter", bolt_name, "name", thread_type) * IN2M
    Ab = boltArea * bolt_count

    # bolt spacing
    Bs_actual = np.pi * C / bolt_count
    Bs_max = 2 * boltDiam + (6 * t) / (m + 0.5)

    # flange moments
    R = (C - B) / 2 - g1
    HD = 0.785 * (B**2) * P
    hD = R + 0.5 * g1
    MD = HD * hD

    HG = Wm1 - H
    hG = (C - G) / 2
    MG = HG * hG

    HT = H - HD
    hT = (R + g1 + hG) / 2
    MT = HT * hT

    W = (Am + Ab) * Sa_used / 2
    Mseating = W * (C - G) / 2

    Moperating = MD + MT + MG
    Mo = max(Mseating, Moperating)
    seating_drives_moment = Mseating > Moperating

    if Bs_actual > (2 * boltDiam + t):
        BSC = (Bs_actual / (2 * boltDiam + t)) ** 0.5
    else:
        BSC = 1
    Mo *= BSC

    # flange stresses
    K = A / B
    Tf = (K**2 * (1 + 8.55246 * np.log10(K)) - 1) / ((1.04720 + 1.9448 * K**2) * (K - 1))
    U = (K**2 * (1 + 8.55246 * np.log10(K)) - 1) / (1.36136 * (K**2 - 1) * (K - 1))
    Y = 1 / (K - 1) * (0.66845 + 5.71690 * (K**2 * np.log10(K)) / (K**2 - 1))
    Z = (K**2 + 1) / (K**2 - 1)

    # values for g0/g1 = 1 (validity enforced above)
    F = 0.908920
    V = 0.550103
    f = 1

    ho = (B * g0) ** 0.5
    e = F / ho
    d = (U / V) * ho * g0**2
    L = (t * e + 1) / Tf + (t**3) / d

    SH = (f * Mo) / (L * g1**2 * B)
    SR = (1.33 * t * e + 1) * Mo / (L * t**2 * B)
    ST = (Y * Mo) / (t**2 * B) - Z * SR

    # Table 2-14 rigidity index for integral-type flanges: KI = 0.3
    KI = 0.3
    J = (52.14 * V * Mo) / (L * E * g0**2 * KI * ho)

    maxStress = max(SH, SR, ST)
    SFmargin = margins(1, Sf, maxStress)

    # 2-8(a) stress limits + bolt area + bolt spacing + 2-14 rigidity
    checks = [
        ("Bolt area Am <= Ab",       Am / IN2M**2,       Ab / IN2M**2,           Am <= Ab,                 "in^2"),
        ("Bolt spacing Bs <= Bs_max", Bs_actual / IN2M,   Bs_max / IN2M,          Bs_actual <= Bs_max,      "in"),
        ("Radial stress SR <= Sf",    SR / PSI2PA,        Sf_used / PSI2PA,       SR <= Sf_used,            "psi"),
        ("Tangential stress ST <= Sf", ST / PSI2PA,       Sf_used / PSI2PA,       ST <= Sf_used,            "psi"),
        ("Hub stress SH <= 1.5*Sf",   SH / PSI2PA,        1.5 * Sf_used / PSI2PA, SH <= 1.5 * Sf_used,      "psi"),
        ("(SH+SR)/2 <= Sf",           (SH + SR) / 2 / PSI2PA, Sf_used / PSI2PA,   (SH + SR) / 2 <= Sf_used, "psi"),
        ("(SH+ST)/2 <= Sf",           (SH + ST) / 2 / PSI2PA, Sf_used / PSI2PA,   (SH + ST) / 2 <= Sf_used, "psi"),
        ("Rigidity index J <= 1.0",   J,                  1.0,                    J <= 1.0,                 ""),
    ]
    allPass = all(c[3] for c in checks)

    print(f"Operating Load Wm1: {Wm1*N2LBF:.3f} lbf")
    print(f"Seating Load Wm2: {Wm2*N2LBF:.3f} lbf")
    print(f"Required bolt area Am: {Am/IN2M**2:.3f} in^2 ({'seating' if seating_drives_area else 'operating'} governs)")
    print(f"Actual bolt area Ab ({bolt_count}x {bolt_name} {thread_type}): {Ab/IN2M**2:.3f} in^2")
    print(f"Total moment Mo: {Mo:.1f} N*m ({'seating' if seating_drives_moment else 'operating'} governs, BSC={BSC:.3f})")
    print(f"(practical min spacing guideline, 1.5x bolt diam: {1.5*boltDiam/IN2M:.3f} in; actual: {Bs_actual/IN2M:.3f} in)")
    print()
    print("Design checks:")
    for name, actual, allow, ok, unit in checks:
        status = "PASS" if ok else "FAIL"
        print(f"  [{status}] {name}: {actual:.3f} {unit} vs {allow:.3f} {unit}")
    print()

    if allPass:
        print("Design satisfies all 2-8(a) stress limits and the 2-14 rigidity check.")
    else:
        failing = [c[0] for c in checks if not c[3]]
        print(f"Design FAILS: {', '.join(failing)}")

    print(f"Stress margin vs yield (SF={safetyFactor}): {SFmargin:.2f}%")

    # plot geometry (reuses FlangeSizer.plot_geometry directly, no FlangeSizer instance needed)
    path = Path("Development/flange/analysis_figures")
    path.mkdir(parents=True, exist_ok=True)
    ns = SimpleNamespace(
        G=G, gasket_ID=gasket_ID, gasket_OD=gasket_OD, N=N,
        facing_sketch=facing_sketch, bo=bo, t=t, safetyFactor=safetyFactor,
        path=path,
    )
    FlangeSizer.plot_geometry(ns, A, B, C, bolt_count, boltDiam, bolt_name, SFmargin, J)

    return {
        "Wm1": Wm1, "Wm2": Wm2, "Am": Am, "Ab": Ab, "Mo": Mo,
        "SH": SH, "SR": SR, "ST": ST, "J": J,
        "checks": checks, "allPass": allPass, "SFmargin": SFmargin,
        "gasket_ID": gasket_ID, "gasket_OD": gasket_OD, "bo": bo, "b": b,
    }


def main():

    # Flange geometry (from drawing)
    B = 3.826 * IN2M   # Flange ID [m]
    N = 3/8 * IN2M      # Gasket width [m]
    G = 4.326 * IN2M    # Gasket reaction diameter [m]
    C = 5.5 * IN2M    # Bolt circle diameter [m]
    A = 6.25 * IN2M    # Flange OD [m]
    t = 3/8 * IN2M       # Flange thickness [m]

    # Hub thicknesses
    g0 = 0.337 * IN2M
    g1 = g0 # equal for straight integral flange

    # Bolts
    bolt_name = "1/4"
    bolt_count = 16
    thread_type = "UNF"

    # working chamber pressure
    P = 440 * PSI2PA

    # gasket properties: vermiculite with SS insert
    m = 2.0
    y = 2500 * PSI2PA

    # bolt yield stress
    Sa = 30000 * PSI2PA
    Sb = Sa

    # flange yield at temp
    Sf = 36000 * PSI2PA

    # modulus of elasticity of flange material at design temperature
    E = 29e6 * PSI2PA

    analyze_flange(B, N, G, C, A, bolt_name, bolt_count, t, g0, g1,
                    P, m, y, Sa, Sb, Sf, E, safetyFactor=1.25,
                    thread_type=thread_type)
    return


if __name__ == "__main__":
    main()
