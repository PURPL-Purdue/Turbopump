"""Tests of inducer_sizing.py against the chapter 6 sample calculations in
Huzel & Huang, "Design of Liquid Propellant Rocket Engines" (NASA SP-125).

Page numbers are the PRINTED book pages. In the PDF, add 9.

Sample calculations used (the others in chapter 6 don't touch this script):
    6-2  p. 191        pump specific speed, eq. (6-7)
    6-3  p. 193        suction specific speed, eq. (6-10)
    6-7  pp. 215-218   LOX pump inducer, A-1 stage (centrifugal pump)
    6-10 pp. 233-236   LH2 pump inducer, A-2 stage (axial-flow pump)

Run any of:
    python -m unittest test_inducer_vs_huzel_huang -v
    python -m pytest test_inducer_vs_huzel_huang.py -v
    python test_inducer_vs_huzel_huang.py        # tests + end-to-end tables

Test groups
    1. Component tests. Each script function gets the book's own inputs and
       must reproduce the book's own answer. Tolerances cover the book's
       slide-rule rounding; each test says how big that rounding is.
    2. Known differences (unittest.expectedFailure). Places where the script
       and the book disagree. They are EXPECTED to fail. If one ever passes,
       unittest reports an "unexpected success" so you notice the change.
    3. Internal checks. Not from the book: they check the script's constants
       and solver against their own derivations.
    4. End-to-end. size_inducer() on the book's design data, printed next to
       the book's answers. Not pass/fail, because the script and the book
       size the inducer by different methods (see notes in the table).
"""

import math
import unittest
import warnings

import bladeParams as ind


def dm(deg: float, minutes: float = 0.0) -> float:
    """Degrees + arc-minutes -> decimal degrees (the book prints 5°42' etc.)."""
    return deg + minutes / 60.0


class BookTestCase(unittest.TestCase):
    def assertRel(self, actual, expected, rel, msg=""):
        err = abs(actual - expected) / abs(expected)
        self.assertLessEqual(
            err, rel,
            f"{msg} script={actual:.6g} book={expected:.6g} "
            f"diff={100 * (actual - expected) / expected:+.2f}% (tol {100 * rel:.2f}%)",
        )

    def assertDeg(self, actual_rad, expected_deg, tol_deg, msg=""):
        actual_deg = math.degrees(actual_rad)
        self.assertLessEqual(
            abs(actual_deg - expected_deg), tol_deg,
            f"{msg} script={actual_deg:.3f} deg book={expected_deg:.3f} deg "
            f"(tol {tol_deg} deg)",
        )


# =============================================================================
# 1. COMPONENT TESTS
# =============================================================================

class TestSpecificSpeed_SampleCalc_6_2(BookTestCase):
    """p. 191, eq. (6-7): N_s = N Q^0.5 / dH^0.75.

    Book arithmetic slips (oxidizer): it takes (12 420)^0.5 = 111.7 (true
    111.45) and (2930)^0.75 = 395 (true 398.2). Together these put its 1980
    about 1.1% high; the correct value is 1959. Tolerance 1.5% for the
    oxidizer, 1% for the fuel.
    """

    def test_oxidizer_pump(self):
        # N = 7000 rpm, Q = 12 420 gpm, dH = 2930 ft  ->  N_s = 1980 (book)
        self.assertRel(ind.specific_speed(7000, 12420, 2930), 1980, 0.015)

    def test_fuel_pump(self):
        # N = 7000 rpm, Q = 7960 gpm, dH = 4790 ft  ->  N_s = 1083
        self.assertRel(ind.specific_speed(7000, 7960, 4790), 1083, 0.01)


class TestSuctionSpecificSpeed_SampleCalcs_6_3_and_6_10(BookTestCase):
    """eq. (6-10), p. 192: N_ss = N Q^0.5 / (NPSH)_c^0.75.

    Book rounding: 111.7 for (12 420)^0.5 again, and (58)^0.75 = 21.
    Tolerance 1%.
    """

    def test_6_3_oxidizer_pump(self):
        # p. 193: N = 7000, Q = 12 420, (NPSH)c = 58 ft  ->  N_ss = 37 230
        self.assertRel(ind.suction_specific_speed(7000, 12420, 58), 37230, 0.01)

    def test_6_3_fuel_pump(self):
        # p. 193: N = 7000, Q = 7960, (NPSH)c = 70 ft  ->  N_ss = 25 790
        self.assertRel(ind.suction_specific_speed(7000, 7960, 70), 25790, 0.01)

    def test_6_10_speed_from_nss(self):
        # p. 234: book solves eq. (6-10) for N = 27 000 rpm from
        # N_ss = 53 400, (NPSH)c = 135 ft, Q = 6080 gpm. Run it forward.
        self.assertRel(ind.suction_specific_speed(27000, 6080, 135), 53400, 0.01)


class TestTipSpeed_Eq_6_65(BookTestCase):
    """p. 214, eqs. (6-65) and (6-65a): u = pi N d / 720 (d in inches).
    Script: u = pi D n / 60 (D in ft). Tolerance 0.5%.
    """

    def test_6_7_mean_tip(self):
        # p. 215: d_t = 11.62 in at 7000 rpm  <->  u_t = 355 ft/s
        self.assertRel(ind.tip_speed(11.62 / 12, 7000), 355, 0.005)

    def test_6_7_inlet_tip(self):
        # p. 217: d_0t = 12.19 in at 7000 rpm  ->  u_0t = 372.5 ft/s
        self.assertRel(ind.tip_speed(12.19 / 12, 7000), 372.5, 0.005)

    def test_6_10_tip(self):
        # p. 234: d_t = 7 in at 27 000 rpm  ->  u_t = 826 ft/s
        self.assertRel(ind.tip_speed(7 / 12, 27000), 826, 0.005)


class TestTipDiameter_K_DT(BookTestCase):
    """K_DT form, derived from eqs. (6-59) p. 213 and (6-64), (6-65a), (6-68)
    p. 214. Round trip: feed the book's Q_ind, N, phi and inlet hub/tip ratio;
    the script must give back the book's inlet tip diameter.

    The book's 3.12 constant is 448.83/144 = 3.117 rounded, so expect ~0.1%.
    Tolerance 0.5%.
    """

    def test_K_DT_value(self):
        exact = (4 * 60 / (448.83 * math.pi ** 2)) ** (1 / 3)
        self.assertAlmostEqual(ind.K_DT, exact, delta=1e-4)

    def test_6_7_inlet_tip_diameter(self):
        # pp. 215-217: Q_ind = 13 040 gpm, N = 7000, phi = 0.0998,
        # d_0h = 2.33 in, d_0t = 12.19 in
        d = ind.tip_diameter(13040, 7000, 2.33 / 12.19, 0.0998) * 12
        self.assertRel(d, 12.19, 0.005)

    def test_6_10_tip_diameter(self):
        # pp. 234-235: Q_ind = 6450 gpm, N = 27 000, phi = 0.0784,
        # d_0h = 2.9 in, d_t = 7 in
        d = ind.tip_diameter(6450, 27000, 2.9 / 7, 0.0784) * 12
        self.assertRel(d, 7.0, 0.005)


class TestVelocityTriangles(BookTestCase):
    """Velocity-triangle relations used in Sample Calcs 6-7 and 6-10."""

    def test_6_7_inlet_flow_angle(self):
        # p. 217 (unnumbered): tan(beta'_0t) = c_m0 / u_0t = 37.2 / 372.5
        # -> 5°42'. Script's gamma_1 = book's beta'_0t.
        self.assertDeg(ind.inlet_flow_angle(37.2, 372.5), dm(5, 42), 0.05)

    def test_6_7_outlet_flow_angle(self):
        # p. 216 (unnumbered): tan(beta'_1) = c_m1 / (u_1 - c'_u1)
        # = 53.1 / (258.5 - 29.2) -> 13°3'. Script's gamma_2 = book's beta'_1.
        self.assertDeg(ind.outlet_flow_angle(53.1, 258.5, 29.2), dm(13, 3), 0.05)

    def test_6_7_tangential_velocity(self):
        # p. 216, eq. (6-66): c'_u1 = dH_ind g / u_1 = 235 x 32.2 / 258.5
        # = 29.2 ft/s. Eq. (6-66) has no efficiency, so eta = 1.
        # Book g = 32.2, script g = 32.174 (0.08%). Tolerance 0.5%.
        self.assertRel(ind.tangential_velocity(235, 258.5, 1.0), 29.2, 0.005)

    def test_6_10_tangential_velocity(self):
        # p. 235, eq. (6-66): c'_u1 = 6500 x 32.2 / 774 = 270 ft/s
        self.assertRel(ind.tangential_velocity(6500, 774, 1.0), 270, 0.005)


class TestBladeGeometry_Eqs_6_50_6_51(BookTestCase):
    """p. 212, eqs. (6-50) and (6-51), used on p. 217."""

    def test_6_7_blade_spacing(self):
        # eq. (6-50): P_i = pi d_t / z = pi x 11.62 / 3 = 12.18 in
        # (book's 12.18 vs exact 12.17: rounding, 0.1%). Tolerance 0.5%.
        s = ind.blade_spacing(11.62 / 12, 3) * 12
        self.assertRel(s, 12.18, 0.005)

    def test_6_7_chord_from_solidity(self):
        # eq. (6-51): S_v = C_i / P_i = 26.57 / 12.18 = 2.18.
        # Script: chord = solidity x spacing. Tolerance 0.5%.
        chord = 2.18 * ind.blade_spacing(11.62 / 12, 3) * 12
        self.assertRel(chord, 26.57, 0.005)


# =============================================================================
# 2. KNOWN DIFFERENCES (expected to fail)
# =============================================================================

class TestKnownDifferences(BookTestCase):

    @unittest.expectedFailure
    def test_nss_relation_vs_eq_6_49_worked_value(self):
        """p. 217 evaluates eq. (6-49) with 8150 and gets (N_ss)_ind = 75 700
        at phi = 0.0998, r_d = 0.3. The script's constant is K_NSS = 3574,
        which gives ~33 650. See also test_script_relation_matches_fig_6_40
        and test_8150_reproduces_book_worked_value below.
        """
        phi, r_d = 0.0998, 0.3
        n_ss = ind.K_NSS * (1 - 2 * phi ** 2) ** 0.75 / phi * math.sqrt(1 - r_d ** 2)
        self.assertRel(n_ss, 75700, 0.02)

    @unittest.expectedFailure
    def test_flow_to_blade_angle_ratio(self):
        """p. 217: flow angle 5°42', chosen vane angle 9° (ratio 0.633; the
        angle of attack is 3°18', inside the 4° limit given on p. 215).
        Script: blade = flow / 0.575 -> 9.91°. The 0.575 is not in this book.
        """
        gamma = ind.inlet_flow_angle(37.2, 372.5)
        beta = gamma / ind.FLOW_TO_BLADE_ANGLE
        self.assertDeg(beta, 9.0, 0.25)

    @unittest.expectedFailure
    def test_tip_diameter_uses_Q_not_Q_ind(self):
        """size_inducer() passes the pump flow Q to tip_diameter(). The book
        uses Q_ind = Q + Q_ee + Q_e/2 (eq. 6-63, p. 214; 13 040 vs 12 420 gpm
        on p. 216). With Q, the 6-7 tip comes out ~1.6% small.
        """
        d = ind.tip_diameter(12420, 7000, 2.33 / 12.19, 0.0998) * 12
        self.assertRel(d, 12.19, 0.005)

    @unittest.expectedFailure
    def test_6_7_flow_coefficient_from_nss(self):
        """Script: phi from N_ss via the K_NSS relation. Book (p. 217): phi
        comes from the geometry (tip sized from the head coefficient, eq. 6-66
        on p. 215), giving 0.0998, then eq. (6-49) is only a check.
        """
        n_ss_corr = ind.corrected_suction_specific_speed(37230, 0.3)
        self.assertRel(ind.flow_coefficient(n_ss_corr), 0.0998, 0.02)

    @unittest.expectedFailure
    def test_6_10_flow_coefficient_from_nss(self):
        """Book (p. 235): tip set equal to the impeller tip (7 in), giving
        phi = 0.0784 (spec: 0.09 max). Script's N_ss-based phi is ~0.061.
        """
        n_ss_corr = ind.corrected_suction_specific_speed(53400, 2.9 / 7)
        self.assertRel(ind.flow_coefficient(n_ss_corr), 0.0784, 0.02)


# =============================================================================
# 3. INTERNAL CHECKS (not from the book's worked examples)
# =============================================================================

class TestInternalChecks(BookTestCase):

    def test_8150_reproduces_book_worked_value(self):
        """Same relation with the book's 8150 instead of K_NSS reproduces the
        p. 217 value. So the form matches; only the constant differs.
        The book's arithmetic has a slip: it writes (1 - 2 phi^2) = 0.9601,
        but 1 - 2(0.0998)^2 = 0.9801. That makes its 75 700 about 1.4% low
        against the correct 76 740. Tolerance 2%.
        """
        phi, r_d = 0.0998, 0.3
        n_ss = 8150 * (1 - 2 * phi ** 2) ** 0.75 / phi * math.sqrt(1 - r_d ** 2)
        self.assertRel(n_ss, 75700, 0.02)

    def test_K_NSS_is_unit_constant_over_3_to_the_3_4(self):
        """Brumfield-criterion derivation (mine, not the book's):
        K = (60/pi) sqrt(448.83 pi / 4) (2g)^(3/4) / 3^(3/4) = 3574.
        The numerator alone is ~8147, which matches the book's 8150.
        """
        unit = 60 / math.pi * math.sqrt(448.83 * math.pi / 4) * (2 * ind.G) ** 0.75
        self.assertRel(unit, 8150, 0.001)
        self.assertRel(ind.K_NSS, unit / 3 ** 0.75, 0.001)

    def test_script_relation_matches_fig_6_40(self):
        """p. 212, Fig. 6-40: dashed curve 'one-dimensional theory (equation
        6-49)'. Values read BY EYE off the scanned figure, so tolerance 5%.
        """
        for phi, fig_value in [(0.10, 36000), (0.13, 26500), (0.19, 17500)]:
            with self.subTest(phi=phi):
                n = ind.K_NSS * (1 - 2 * phi ** 2) ** 0.75 / phi
                self.assertRel(n, fig_value, 0.05)

    def test_flow_coefficient_round_trip(self):
        """flow_coefficient(exact=True) must invert the K_NSS relation."""
        for phi in (0.06, 0.10, 0.15, 0.20):
            with self.subTest(phi=phi):
                n = ind.K_NSS * (1 - 2 * phi ** 2) ** 0.75 / phi
                self.assertRel(ind.flow_coefficient(n, exact=True), phi, 1e-6)

    def test_linearized_phi_over_table_6_5_range(self):
        """Table 6-5 (p. 213) gives phi = 0.06-0.20. Over that range the
        closed-form (exact=False) phi should be within 0.1% of exact."""
        for phi in (0.06, 0.10, 0.15, 0.20):
            with self.subTest(phi=phi):
                n = ind.K_NSS * (1 - 2 * phi ** 2) ** 0.75 / phi
                self.assertRel(ind.flow_coefficient(n, exact=False), phi, 0.001)


# =============================================================================
# 4. END-TO-END: size_inducer() on the book's design data
# =============================================================================

# Inputs taken from the book. blade_efficiency = 1.0 because the book's
# eq. (6-66), used for c'_u1 in both examples, has no efficiency term.
CASE_6_7 = ind.InducerInputs(
    flow_gpm=12420,        # p. 215
    speed_rpm=7000,        # p. 215
    npsh_ft=58,            # p. 215, (NPSH)c
    inducer_head_ft=235,   # p. 215, dH_ind from eq. (6-47)
    hub_tip_ratio=0.3,     # p. 215, r_d
    blade_efficiency=1.0,
    num_blades=3,          # p. 217
    solidity=2.2,          # p. 215 (given); p. 217 computes 2.18
)

CASE_6_10 = ind.InducerInputs(
    flow_gpm=6080,         # p. 234
    speed_rpm=27000,       # p. 234
    npsh_ft=135,           # p. 234, (NPSH)c
    inducer_head_ft=6500,  # p. 235, dH_ind from eq. (6-66)
    hub_tip_ratio=2.9 / 7,  # p. 235: d_0h = 2.9 in, d_t = 7 in
    blade_efficiency=1.0,
    num_blades=3,          # rotor blade count not given in 6-10
    solidity=2.2,          # rotor solidity not given in 6-10
)

# (label, script attribute, unit scale, book value or None, note)
ROWS_6_7 = [
    ("N_ss", "suction_specific_speed", 1, 37230, "p. 193/215"),
    ("phi", "flow_coefficient", 1, 0.0998, "p. 217, from geometry"),
    ("D_t inlet tip [in]", "tip_diameter_ft", 12, 12.19, "p. 215 (mean tip 11.62)"),
    ("D_h [in]", "hub_diameter_ft", 12, 2.33, "p. 215 inlet hub (mean hub 3.49)"),
    ("u tip [ft/s]", "tip_speed_fps", 1, 372.5, "p. 217, inlet tip"),
    ("C_m [ft/s]", "meridional_velocity_fps", 1, 37.2, "p. 216, inlet"),
    ("C_u [ft/s]", "tangential_velocity_fps", 1, 29.2, "p. 216, at mean outlet d_1"),
    ("gamma_1 [deg]", "inlet_flow_angle_rad", "deg", dm(5, 42), "p. 217, inlet tip"),
    ("beta_1 blade [deg]", "inlet_blade_angle_rad", "deg", 9.0, "p. 217, inlet tip"),
    ("gamma_2 [deg]", "outlet_flow_angle_rad", "deg", dm(13, 3), "p. 216, at mean outlet d_1"),
    ("spacing [in]", "blade_spacing_ft", 12, 12.18, "p. 217, at mean tip 11.62 in"),
    ("chord [in]", "chord_ft", 12, 26.57, "p. 217"),
]

ROWS_6_10 = [
    ("N_ss", "suction_specific_speed", 1, 53400, "p. 234"),
    ("phi", "flow_coefficient", 1, 0.0784, "p. 235, from geometry"),
    ("D_t tip [in]", "tip_diameter_ft", 12, 7.0, "p. 234, set = impeller tip"),
    ("D_h inlet [in]", "hub_diameter_ft", 12, 2.9, "p. 235"),
    ("u tip [ft/s]", "tip_speed_fps", 1, 826, "p. 234"),
    ("C_m [ft/s]", "meridional_velocity_fps", 1, 64.8, "p. 235, inlet"),
    ("C_u [ft/s]", "tangential_velocity_fps", 1, 270, "p. 235, at mean outlet d_1"),
    ("gamma_1 [deg]", "inlet_flow_angle_rad", "deg", None, "not given"),
    ("beta_1 blade [deg]", "inlet_blade_angle_rad", "deg", None, "not given"),
    ("gamma_2 [deg]", "outlet_flow_angle_rad", "deg", None, "not given"),
]


def _value(res, attr, scale):
    v = getattr(res, attr)
    return math.degrees(v) if scale == "deg" else v * scale


def comparison_rows(case, rows):
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        res = ind.size_inducer(case)
    out = []
    for label, attr, scale, book, note in rows:
        s = _value(res, attr, scale)
        diff = None if book is None else 100 * (s - book) / book
        out.append((label, s, book, diff, note))
    return out


def print_comparison():
    for title, case, rows in [
        ("Sample Calc 6-7 (A-1 LOX pump inducer, pp. 215-218)", CASE_6_7, ROWS_6_7),
        ("Sample Calc 6-10 (A-2 LH2 pump inducer, pp. 233-236)", CASE_6_10, ROWS_6_10),
    ]:
        print(f"\n{title}")
        print(f"{'quantity':<20}{'script':>12}{'book':>12}{'diff':>9}  source / note")
        for label, s, book, diff, note in comparison_rows(case, rows):
            b = "-" if book is None else f"{book:.4g}"
            d = "" if diff is None else f"{diff:+.1f}%"
            print(f"{label:<20}{s:>12.4g}{b:>12}{d:>9}  {note}")


class TestEndToEnd(BookTestCase):
    """Only N_ss is pass/fail here: it uses the same eq. (6-10) inputs as the
    book. Everything else is method-dependent; see print_comparison()."""

    def test_6_7_runs_and_nss_matches(self):
        res = ind.size_inducer(CASE_6_7)
        self.assertRel(res.suction_specific_speed, 37230, 0.01)

    def test_6_10_runs_and_nss_matches(self):
        res = ind.size_inducer(CASE_6_10)
        self.assertRel(res.suction_specific_speed, 53400, 0.01)


if __name__ == "__main__":
    import sys
    print_comparison()
    print()
    unittest.main(argv=[sys.argv[0], "-v"])
