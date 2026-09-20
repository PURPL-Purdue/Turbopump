## FEA Results Analysis for Flange Assembly
# Checks bolt preload stress, pressure-induced load, gasket bearing stress,
# and flange body stress against FEA-derived reaction/stress values.
#
# Reproduces the calc chain from the reference screenshots:
#   Image 1: bolt stress from total preload force -> S/Sy
#   Image 2: pressure force on chamber ID (info) + gasket bearing stress = R / A_gasket
# and adds a flange-body margin check using margins.py.
#
# Units: inches, psi, lbf (no SI conversion, consistent with refBolt.py)

import numpy as np
from refBolt import bolt_lookup
from margins import margins
from flange_sizing import basic_seating_width


def gasket_seating_check(d1, gasket_width, m, P, R_fea,
                          facing_sketch="1a", column="II", w=None, T=None,
                          G=None, d2=None):
    """
    Required vs. actual gasket contact stress under OPERATING conditions.

    Required side (per BPVC.VIII-1 App. 2, Wm1 = H + Hp):
        Hp = 2*b*pi*G*m*P              (joint-contact force, Appendix 2)
        required_stress = Hp / A_gasket
    Actual side (from FEA):
        actual_stress = R_fea / A_gasket

    PASS if actual_stress >= required_stress, i.e. the sim shows at least
    as much residual contact stress on the gasket as the code requires to
    resist blow-out at operating pressure.

    Parameters
    ----------
    d1            : gasket ID [in]
    gasket_width  : physical/radial gasket width, N [in]
    m             : gasket factor []
    P             : operating (design) pressure [psi]
    R_fea         : FEA reaction force at the gasket, under operating load [lbf]
    facing_sketch, column, w, T : passed to basic_seating_width (Table 2-5.2)
    G, d2         : optional overrides. If not given, this assumes the
                     bo <= 6 mm case (b = bo, G = d1 + gasket_width,
                     d2 = d1 + 2*gasket_width), same as flange_analysis.py.
                     For a wide gasket (bo > 6 mm) you must supply G and d2
                     from the actual flange design (see flange_analysis.py),
                     since b and G are no longer simply derived from N.
    """
    IN2MM = 25.4

    bo = basic_seating_width(gasket_width, facing_sketch, column, w, T)  # in
    bo_mm = bo * IN2MM

    if bo_mm <= 6:
        b = bo
        if G is None:
            G = d1 + gasket_width
        if d2 is None:
            d2 = d1 + 2 * gasket_width
    else:
        if G is None or d2 is None:
            raise ValueError(
                f"bo = {bo_mm:.2f} mm > 6 mm: b and G are not simply "
                "derived from d1/gasket_width in this regime. Supply G "
                "and d2 explicitly (see flange_analysis.py Table 2-5.2 logic)."
            )
        Cb = 2.5
        b = (Cb * np.sqrt(bo_mm)) / IN2MM

    A_gasket = (np.pi / 4) * (d2**2 - d1**2)

    Hp = 2 * b * np.pi * G * m * P
    required_stress = Hp / A_gasket
    actual_stress = R_fea / A_gasket

    margin_pct = (actual_stress / required_stress - 1) * 100 if required_stress > 0 else float("inf")
    passed = actual_stress >= required_stress

    print("=" * 60)
    print("GASKET SEATING CHECK (operating condition, Hp basis)")
    print("=" * 60)
    print(f"Gasket ID d1: {d1} in, width N: {gasket_width} in, OD d2: {d2:.4f} in")
    print(f"Basic seating width bo: {bo:.4f} in ({bo_mm:.2f} mm), effective b: {b:.4f} in")
    print(f"Gasket reaction diameter G: {G:.4f} in")
    print(f"Gasket area A = pi/4*(d2^2-d1^2): {A_gasket:.4f} in^2")
    print(f"Gasket factor m: {m}, operating pressure P: {P} psi")
    print(f"Required joint-contact force Hp = 2*b*pi*G*m*P: {Hp:.3f} lbf")
    print(f"Required stress Hp/A: {required_stress:.3f} psi")
    print(f"FEA reaction force R: {R_fea} lbf")
    print(f"Actual (sim) stress R/A: {actual_stress:.3f} psi")
    print(f"[{'PASS' if passed else 'FAIL'}] actual >= required, margin: {margin_pct:.2f}%")
    print("=" * 60)

    return {
        "bo": bo, "b": b, "G": G, "d2": d2, "A_gasket": A_gasket,
        "Hp": Hp, "required_stress": required_stress,
        "actual_stress": actual_stress, "margin_pct": margin_pct, "passed": passed,
    }


def analyze_fea_results(B, gasket_width, P,
                         bolt_count, bolt_name, Sy_bolt, F_preload,
                         R_fea, max_stress_flange, Sy_flange,
                         safetyFactor=1, thread_type="UNF",
                         preload_allowable_frac=0.7):
    """
    Parameters
    ----------
    B                    : chamber/flange ID [in]
    gasket_width         : physical gasket contact width, N [in]
    P                    : chamber design pressure [psi]
    bolt_count           : number of bolts, N (Image 1 notation)
    bolt_name            : bolt size string for refBolt lookup (e.g. "5/16")
    Sy_bolt              : bolt material yield stress [psi]
    F_preload            : total measured/target bolt preload force, F [lbf]
    R_fea                : FEA reaction force used for gasket bearing check [lbf]
    max_stress_flange    : max stress on flange from FEA [psi]
    Sy_flange            : flange material yield stress [psi]
    safetyFactor         : safety factor used in margins() calls
    thread_type          : "UNF" or "UNC", passed to bolt_lookup
    preload_allowable_frac : fraction of Sy_bolt used as an informal preload
                              stress limit (0.7 in Image 1's reference calc)

    Returns
    -------
    dict of all computed values
    """

    # ---------- Bolt preload stress (Image 1) ----------
    A_bolt = bolt_lookup("area", bolt_name, "name", thread_type)  # in^2, single bolt
    F_bolt = F_preload / bolt_count                                # lbf, per bolt
    S_bolt = F_bolt / A_bolt                                       # psi
    bolt_ratio = S_bolt / Sy_bolt                                  # S/Sy, as in screenshot
    bolt_margin = margins(safetyFactor, Sy_bolt, S_bolt)           # %, vs given SF

    preload_allow = preload_allowable_frac * Sy_bolt
    preload_ok = S_bolt <= preload_allow

    # ---------- Pressure-induced load (Image 2, informational) ----------
    A_inj = (np.pi / 4) * B**2
    F_inj = P * A_inj
    F_inj_per_bolt = F_inj / bolt_count

    # ---------- Gasket bearing stress (Image 2) ----------
    gasket_OD = B + 2 * gasket_width
    A_gasket = (np.pi / 4) * (gasket_OD**2 - B**2)
    S_gasket = R_fea / A_gasket

    # ---------- Flange body stress margin ----------
    flange_margin = margins(safetyFactor, Sy_flange, max_stress_flange)
    flange_ratio = max_stress_flange / Sy_flange
    flange_ok = max_stress_flange * safetyFactor <= Sy_flange

    # ---------- Print report ----------
    print("=" * 60)
    print("BOLT PRELOAD STRESS")
    print("=" * 60)
    print(f"Bolt: {bolt_count}x {bolt_name} {thread_type}, A_bolt = {A_bolt:.4f} in^2")
    print(f"Total preload force F: {F_preload:.3f} lbf")
    print(f"Force per bolt F_bolt = F/N: {F_bolt:.3f} lbf")
    print(f"Bolt stress S = F_bolt/A_bolt: {S_bolt:.3f} psi")
    print(f"S/Sy: {bolt_ratio:.4f}   (Sy = {Sy_bolt:.0f} psi)")
    print(f"Margin vs yield (SF={safetyFactor}): {bolt_margin:.2f}%")
    print(f"[{'PASS' if preload_ok else 'FAIL'}] S <= {preload_allowable_frac:.1f}*Sy "
          f"({preload_allow:.0f} psi): actual {S_bolt:.0f} psi")
    print()

    print("=" * 60)
    print("PRESSURE-INDUCED LOAD (info)")
    print("=" * 60)
    print(f"Chamber ID area A_inj = pi/4*B^2: {A_inj:.4f} in^2  (B = {B} in)")
    print(f"Pressure P: {P:.1f} psi")
    print(f"Total pressure force F_inj = P*A_inj: {F_inj:.3f} lbf")
    print(f"Pressure force per bolt F_inj/N: {F_inj_per_bolt:.3f} lbf")
    print(f"(compare to preload per bolt F_bolt = {F_bolt:.3f} lbf -- "
          f"preload should exceed this to keep joint from separating)")
    print()

    print("=" * 60)
    print("GASKET BEARING STRESS")
    print("=" * 60)
    print(f"Gasket width N: {gasket_width} in, ID: {B} in, OD: {gasket_OD:.4f} in")
    print(f"Gasket contact area A_gasket = pi/4*(OD^2-ID^2): {A_gasket:.4f} in^2")
    print(f"FEA reaction force R: {R_fea:.3f} lbf")
    print(f"Gasket bearing stress R/A_gasket: {S_gasket:.3f} psi")
    print()

    print("=" * 60)
    print("FLANGE BODY STRESS")
    print("=" * 60)
    print(f"Max FEA stress on flange: {max_stress_flange:.3f} psi")
    print(f"Flange yield Sy: {Sy_flange:.0f} psi")
    print(f"Stress/Sy: {flange_ratio:.4f}")
    print(f"[{'PASS' if flange_ok else 'FAIL'}] margin vs yield (SF={safetyFactor}): {flange_margin:.2f}%")
    print("=" * 60)

    return {
        "A_bolt": A_bolt, "F_bolt": F_bolt, "S_bolt": S_bolt,
        "bolt_ratio": bolt_ratio, "bolt_margin": bolt_margin, "preload_ok": preload_ok,
        "A_inj": A_inj, "F_inj": F_inj, "F_inj_per_bolt": F_inj_per_bolt,
        "gasket_OD": gasket_OD, "A_gasket": A_gasket, "S_gasket": S_gasket,
        "flange_margin": flange_margin, "flange_ratio": flange_ratio, "flange_ok": flange_ok,
    }


def main():
    # --- Fill these in with your values / FEA results ---

    # Geometry
    B = 3.826            # chamber/flange ID [in]
    gasket_width = 3/8   # gasket contact width, N [in]

    # Loads
    P = 660              # chamber design pressure [psi]
    F_preload = 16423.189  # total bolt preload force [lbf]
    R_fea = 9085.1        # FEA reaction force for gasket check [lbf]
    max_stress_flange = 26288  # <-- REPLACE: max FEA stress on flange body [psi]

    # Bolts
    bolt_count = 14
    bolt_name = "5/16"
    thread_type = "UNF"
    Sy_bolt = 30000       # bolt yield stress [psi]

    # Flange material
    Sy_flange = 36000     # flange yield stress [psi]

    # Safety factor for margin reporting
    safetyFactor = 1

    analyze_fea_results(B, gasket_width, P,
                         bolt_count, bolt_name, Sy_bolt, F_preload,
                         R_fea, max_stress_flange, Sy_flange,
                         safetyFactor=safetyFactor, thread_type=thread_type)

    # --- Gasket seating check (operating condition) ---
    # d1, gasket_width, R_fea below are from your latest screenshot.
    # m, P are NOT shown in that screenshot -- these are pulled from
    # flange_analysis.py / flange_sizing.py's main() (m=2, P=440*1.5=660 psi).
    # Replace with your actual operating m/P if different.
    d1 = 3.826
    R_fea_gasket = 8494.6
    m = 2
    P_operating = 660

    gasket_seating_check(d1, gasket_width, m, P_operating, R_fea_gasket)


if __name__ == "__main__":
    main()
