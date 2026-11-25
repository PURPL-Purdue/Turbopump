# Where every number comes from

Path: **CEA sizing → contour → wall temp solver → channel solver** (IPA / LOX).

**How to read the tables**

- **In code** is the variable or function name. Search for it in the file named in the section heading.
- **Source** is one of:
  - **Book**: Heister et al., *Rocket Propulsion*. I checked the page against the PDF in `reference/`.
  - **You**: in your own first commit, before Claude touched the repo.
  - **No source**: someone picked it, possibly Claude, and nothing backs it up. **This needs your decision.**
  - **CEA** or **thermo**: output of that library.
  - **Math**: follows directly from the rows above it.

---

## 1. Engine inputs: `regen.yaml` (`engine`, `chamber`, `contour`), read by `engine/engine_inputs.py`

Values updated 2026-11-21 to match `regen.yaml`. Where a value changed after Oct 11, the note says so. The reason for those changes wasn't written down.

| In code | Value | Source | Note |
|---|---|---|---|
| `reac_names` | IPA / LOX | You | |
| `T_reactant` | fuel 293 K, LOX 90.17 K | You | 90.17 K is LOX's boiling point. In reality the fuel reaches the injector hot from the jacket (~500 K), not at 293 K. |
| `fuel_dH_vap` | 45.5 kJ/mol = 7.57e5 J/kg | NIST WebBook | IPA heat of vaporization (ΔHf gas −272.6, liquid −318.1 kJ/mol). Subtracted in `_solve()` because CEA's library only has IPA as a gas. Without it Tc is 30–115 K too high. |
| `F_target_lbf` | 1500 lbf | You | Only printed for comparison. Nothing is sized to it. |
| `p0_psi` | 440 psi | You | What the heatsink engine runs at right now. |
| `pa_psi` | 18 psi | You (commit msg) | Pulled from heatsink engine parameters. |
| `OF` | 1.5 | You | Completly arbitrary for now. Have not decided on a specific value yet.  |
| `m_dot_pump_max` | 2.0 kg/s per propellant | You | |
| `c_star_eff` | 0.875 | **Prior test data** | Book p.328 says ~95% is typical for LOX/RP, so this is low. |
| `c_thrust_eff` | 0.97 | **Quick lookup** | Book p.328: 97–99% for large nozzles, "low 90s" for small engines. |
| `contraction_ratio` | 6 | **No source** | You started at 4. It later became 5, then 9, then 6. Book p.305: anything from 1.3 to 10, with small engines toward 10. |
| `lstar` | 1.3 m | You | You said it came from a paper. **Add which one.** Book Table 8.1 gives LOX/RP-1 1.0–1.3 m and doesn't list IPA. |
| `theta_conv_deg` | 30° | **No source** | |
| `bell_percent` | 80 | You | Chosen 2026-11-21 in place of the 15° cone. Heister §4.3 p.99 describes the Rao bell and what "80%" means, but doesn't say 80% is the standard choice. See section 4. |
| `R_up_over_Rt` | 1.65 | **No source** | Was 1.5 on Oct 11. The code comment says Huzel & Huang. I couldn't check: that book isn't in `reference/`. |
| `R_down_over_Rt` | 0.382 | **No source** | The code comment says Rao. I couldn't check it. The team's bell script uses the same 0.382 Rt. |

## 2. CEA settings: `engine/cea_interface.py`

The library is `cea` 3.3.4 (NASA CEA). Its chemistry and transport data are its own.

| What | Setting | Source | Note |
|---|---|---|---|
| Combustor | infinite-area | CEA default | Never set explicitly in our code |
| Nozzle expansion | shifting equilibrium | CEA default | So T, Mach, γ, c*, Cf and ε are all equilibrium values |
| Exit plane | where pressure = `pa_psi` | Math | |
| cp and Pr fed to Bartz | **frozen** values | **No source** | A code comment claims Bartz and Huzel & Huang use frozen properties. Not verified. |
| Unit conversions | bar→Pa, millipoise→Pa·s, kJ→J | Math | Checked |

## 3. Engine sizing: `size_engine()` in `engine/cea_interface.py`

| Step | Equation | Source |
|---|---|---|
| Mass flow | O/F > 1, so oxidizer is pump-limited: `ṁ_ox = 2`, `ṁ_f = ṁ_ox / OF` | You |
| Real c* | `c* = c*_CEA × c_star_eff` | Math |
| Real Cf | `Cf = Cf_CEA × c_thrust_eff` | Math |
| Throat area | `At = ṁ c* / p0` | Book Eq 8.44, p.328 |
| Exit area | `Ae = ε At` (ε from CEA) | CEA |
| Chamber area | `Ac = CR × At` | Book Eq 8.45, p.328 |
| Chamber volume | `Vc = L* × At` | Book Eq 8.24, p.304 |
| Thrust | `F = Cf p0 At` | Math. Printed only. Comes out at 1597 lbf. |

Current results (2026-11-21, O/F 1.4): Rt 23.12 mm, Rc 56.62 mm, Re 46.69 mm (ε 4.079), ṁ_fuel 1.429 kg/s, T0 2951 K.

## 4. Contour: `build_contour()` in `engine/contour.py`

| What | Source | Note |
|---|---|---|
| Chamber shape: cylinder, 30° cone, arc, throat | **No source** | Chosen because it's easy to sketch in CAD. The cylinder→cone joint is a sharp corner (no blend radius). |
| Cylinder length set so the injector→throat volume equals `Vc` | Math | Includes the cone and the upstream arc |
| Divergent: 0.382 Rt arc to θn, then a parabola to the exit at θe | Book §4.3 p.99 | Rao's approximate method: a parabola from the end of the downstream throat arc (at θn) to the exit (at θe) |
| Parabola drawn as a quadratic Bézier N→Q→E, Q where the two end tangents cross | Team code, `engine/rao_bell.py` `bell_nozzle()` | That script cites the Aspire Space note "The Thrust Optimised Parabolic Nozzle" (eqns 2–10). The note isn't in `reference/` and hasn't been checked. A quadratic Bézier is a parabola, so this matches the book's description. |
| θn, θe from a table of (ε, bell %) | Team code, `find_wall_angles()` | Said to be read off Rao (1958) charts. **Not checked against Rao.** At ε 4.079 and 80%, interpolating between ε 4 and 5 gives θn 21.62°, θe 13.92°. The table starts at ε 4. Below that it uses the ε 4 values without warning. Its 60% row is non-monotonic (θn drops at ε 50), which looks like a typo. |
| Bell exit at x = 0.8 (√ε − 1) Rt / tan 15° past the throat | Book p.99 | 80% of the length of a 15° cone with the same ε. Comes out to 70.38 mm, against about 88 mm for the old 15° cone. |
| 100 stations along the length | Book p.213 | "at least 100 points" |
| `output/contour.dxf` | Math | The wall as exact entities (lines, arcs, degree-2 spline with control points N, Q, E), mm. The channel floor and top are offset **perpendicular** to the wall by t_w and t_w + h(x), which is how the solver uses t_w and h. Floor lines and arcs are exact, the bell floor is a fitted spline, and the sharp cylinder→cone corner is filled with an arc of radius t_w. Checked: floor 0.7000 mm thick everywhere. The channel top is within 0.3 µm (median) of t_w + h; between stations, h is whatever the spline does. |
| `output/channel_heights.csv`, `channel_floor.dat`, `channel_top.dat` | Math | Still a **radial** shift (r + t_w, r + t_w + h). On sloped walls that gives a wall t_w cos θ thick (0.61 mm on the 30° cone), so it doesn't match the DXF there. The solver's heated width `b = 2π(r + t_w)/N` is also radial. |

## 5. Gas-side heat transfer: `heat_transfer/gas_side.py`

| In code | Equation | Source | Note |
|---|---|---|---|
| `gas_environment()` | Gas state along the nozzle from CEA at A/At = (r/Rt)² | CEA | |
| `recovery_temperature()` | `T_r = T (1 + r (γ−1)/2 M²)` | Book Eq 6.15, p.196 | |
| same | recovery factor `r = Pr^(1/3)` | Book Table 6.2, p.196 | Turbulent boundary layer |
| `station_fluxes()` | `q = h_g (T_r − T_wg)` | Book Eq 6.16, p.196 | |
| `bartz_hg()` | Bartz 1957 | Book Eq 6.13, p.195 | |
| `bartz_hg()` | `T_am = (T + T_wg)/2` | Book p.195 | |
| `bartz_hg()` | density ratio `ρ_am/ρ = T/T_am` | **No source** | The book doesn't say how to evaluate this. Assumes ideal gas. |
| `bartz_hg()` | viscosity ratio `μ_am/μ0 = (T_am/T0)^ω`, ω fitted to CEA data | **No source** | The book doesn't say how to evaluate this either |
| `gas["G"]` | `ρV = ṁ / (π r²)` | Math | |
| gas temperature | ideal CEA, not reduced by `c_star_eff` | **No source** | Conservative: the gas is hotter than it would be with the efficiency applied |

## 6. IPA properties: `coolant/ipa_props.py`

Comes from the `thermo` Python package (Caleb Bell). It is tabulated at 850 psi from 270 to 480 K. Above 480 K the 480 K values are reused.

`thermo` picks a published correlation for each property:

| Property | Correlation `thermo` uses | Valid up to | We use it up to |
|---|---|---|---|
| density | DIPPR (Perry's Handbook 8th ed.) | 508 K | 480 K |
| **viscosity** | **DIPPR (Perry's 8th ed.)** | **355 K** | **480 K: everything above 355 K is extrapolated** |
| conductivity | thermo's own 2023 fit | 470 K | 480 K |
| cp | Zabransky | 473 K | 480 K |

Also:

- The code claims the table was "checked against handbook values". It doesn't say which handbook.
- IPA's critical point is 508 K and 691 psi. The jacket runs at 850 → ~717 psi, which is above the critical pressure. So near 508 K, IPA isn't a normal liquid, and none of these correlations apply.

## 7. Coolant-side heat transfer and pressure drop: `heat_transfer/coolant_side.py`, `heat_transfer/friction.py`, `jacket/geometry.py`, `jacket/march.py`

The correlation is picked in `regen.yaml` (`jacket.coolant_htc`): `sieder_tate` or `gnielinski`. The Re check in the report uses each correlation's own lower limit.

"VDI G1" below is VDI Heat Atlas, 2nd ed. (Springer, 2010), chapter G1 "Heat Transfer in Pipe Flow" by V. Gnielinski, `reference/LB_G1 Heat Transfer in Pipe Flow.pdf`.

| In code | Equation | Source | Note |
|---|---|---|---|
| `sieder_tate()` | Sieder-Tate `Nu = a Re^m Pr^n (μ/μ_w)^b` | Book Eq 6.17, p.197 | |
| `sieder_tate()` | `a, m, n = 0.023, 0.8, 0.4` | **No source** | These are the book's **kerosene** constants (Table 6.3). IPA isn't in the table. |
| `sieder_tate()` | `b = 0.114` | Book p.197 | |
| `sieder_tate()` | properties at `(T_l + T_wl)/2` | Book p.198 | "the mean temperature, which is reflective of the average between the wall and the bulk fluid". Above ~480 K this mean is past the IPA property table, so Re and Pr use clamped values (see Problems 1-2). |
| `gnielinski_nu()` | `Nu = (ξ/8) Re Pr / (1 + 12.7 √(ξ/8) (Pr^(2/3) − 1))` | VDI G1 Eq. 26, p.696 | 2010 form (uses Re, not Re − 1000). Checked against VDI's worked example (Table 1, Eq. 37): 82.53 vs 82.57, difference is VDI rounding ξ to 0.0308. |
| `gnielinski_nu()` | local entrance factor `× [1 + (1/3)(d_i/x)^(2/3)]` | VDI G1 Eq. 28, p.696 | Only with `entrance_effects: true` (x = distance from the coolant inlet). Off: factor = 1, fully developed, same as before (Nov 24). +17% at 5 mm, +4% at 45 mm from the inlet (D_h 1.8 mm). |
| `laminar_nu()` | `Nu₁ = 4.354`, `Nu₂ = 1.302 (Re Pr d_i/x)^(1/3)`, `Nu₃ = 0.462 Pr^(1/3) (Re d_i/x)^(1/2)` | VDI G1 Eqs. 14, 15, 22, p.695 | Local, constant heat flux, velocity and temperature developing together. Pr > 0.7. |
| `laminar_nu()` | `Nu = (Nu₁³ + (Nu₂ − 1)³ + Nu₃³)^(1/3)` | VDI G1 Eq. 23, p.695 | Implemented **as printed**. The related Eq. 16 has a `+ 1` term (`Nu₁³ + 1 + (Nu₂ − 1)³`) that Eq. 23 lacks; with it, the fully developed limit is exactly 4.354; without it, 4.335. Effect ≤ 0.5%. With `entrance_effects: false` the code returns Eq. 14 (4.354) directly. |
| `laminar_nu()` | no wall-property correction | VDI G1 §3.2.3 | G1 gives none for constant heat flux (no experimental data). |
| `vdi_nu()` | `Nu = (1 − γ) Nu_lam,2300 + γ Nu_turb,10⁴`, `γ = (Re − 2300)/(10⁴ − 2300)` | VDI G1 Eqs. 29–30, p.696 | Re ≤ 2300: `laminar_nu` at Re. Re ≥ 10⁴: `gnielinski_nu`. Continuous at both ends (checked). |
| `vdi_nu()` | both ends of Eq. 29 evaluated **locally** at the station's d_i/x | **No source** (our step) | G1 writes Eq. 29 with mean (whole-tube) values (Eq. 25 / Eq. 26 with d_i/l). The march needs local h, so local Eqs. 23 and 28 are used instead. |
| `konakov_smooth_darcy()` | `ξ = (1.8 log10 Re − 1.5)^−2`, smooth pipe, Darcy | VDI G1 Eq. 27, p.696 (Konakov) | Always the smooth value, whatever `eps_rough` is; the pressure drop uses `friction_factor()` separately. |
| `gnielinski_nu()` | liquids: `× (Pr/Pr_w)^0.11` | VDI G1 Eq. 40, p.697 | Valid 0.1 < Pr/Pr_w < 10. Pr_w uses all properties at T_wl, so above 480 K it's from clamped values. |
| `gnielinski()` | Re, Pr, k at the bulk temperature `T_l` | VDI G1 §4.1, p.696 | G1 refers properties to the mean of inlet and outlet bulk temperature; in a station march that's the local bulk temperature. |
| `gnielinski` | turbulent part valid 10⁴ ≤ Re ≤ 10⁶, 0.1 ≤ Pr ≤ 1000; laminar part 0.7 < Pr < 1000 | VDI G1 §3.2.2, §4.1 | Laminar and transition implemented Nov 24 (`vdi_nu`), so the report no longer flags Re < 10⁴; it counts transition and laminar stations instead. Pr and Re upper limits are not checked. |
| `inlet_distance()` in `jacket/march.py` | x = arc length along the wall from the coolant inlet (nozzle exit for counterflow) | **No source** (our choice) | The inlet station uses half its own segment length (~5 mm), not 0, where the developing-flow terms would be infinite. Assumes heating and flow both start at the manifold. |
| `util/report.print_checks()` | warning when Re < 10,000 (both correlations) | Book p.199; VDI G1 §4.1 | |
| `channel_geometry()` | `D_h = 4A / perimeter` | Book Eq 6.18, p.198 | |
| `friction_factor()` | smooth-wall curve fits | Book Eq 6.37, p.209 | |
| `friction_factor()` | which fit to use where they overlap | **No source** | The book's ranges overlap and leave a gap |
| `advance_coolant()` | `Δp = 4 cf (ds/D_h) ρv²/2` | Book Eq 6.36, p.209 | |

## 8. The march (Book §6.5.2, p.211–213): `jacket/march.py`, `jacket/geometry.py`

| In code | What | Source | Note |
|---|---|---|---|
| `march()` | Heat balance at each station: `h_g(T_r − T_wg) = k(T_wg − T_wl)/t_w = h_l(T_wl − T_l)` | Book Eq 6.29, p.204 | |
| `solve_station()` | Solves for T_wg with a root-finder | Book step 4–9 | The book just says "iterate" |
| `run()` | `ṁ_chan = ṁ_fuel / N`, all fuel through the jacket | Book Eq 6.31, p.206 | |
| `channels_from_land()` | N = as many channels as fit around the throat with `land_min` between them | **No source** | The book's formula (Eq 6.30) is for brazed tubes. This rule for milled channels was made up. |
| `advance_coolant()` | Coolant heats up: `ΔT = q b ds / (ṁ cp)` | Book Eq 6.39, p.212 | |
| `channel_geometry()` | `b` (heated width per channel) measured at the channel floor radius | **No source** | Arguably it should be at the gas-side radius. The difference overstates coolant heating by ~4%. |
| `march()` | Outlet pressure = last station's static pressure | Book step 11, p.213 | |
| — | Not modelled: curvature, entrance effects, heat loss out the back | — | The book mentions all of these. Left out on purpose. |

### 8a. Fin conduction through the lands: `heat_transfer/fins.py`, `coolant_area_factor()` in `jacket/march.py`

On/off with `jacket.fin_conduction`. Off: φ = 1 and the march is the book's (checked: channel and wall temp solver output identical to before, Nov 23). "Blong" is Blong, Martinez, Heister (Purdue, 2024), "Investigation of Heat Transfer in a Water-Cooled GOx/RP-1 RDRE", `reference/investigation_of_heat_transfer_in_a_water-cooled_RDRE-DESENSITIZED-blong-2024.pdf`.

| In code | Equation | Source | Note |
|---|---|---|---|
| `fin_width()` | `w_fin = 2π(r + t_w + h_c/2)/N − w_c` | Blong Eq. 5, **interpreted** | Land width at mid-height (Blong: "width at the midpoint of the fin"). As printed, Eq. 5 reads `(r + 2π(t_w + 0.5h_c)·is_not_centerbody)/N − w_c`. Taken literally that gives −0.56 mm at our throat (r 23.45, t_w 1.0, h 2.51 mm, N 85), so the 2π is read as multiplying the whole radius: pitch circumference at mid-fin over N. For a chamber wall (not a centerbody) the flag is 1. The centerbody case isn't used here. |
| `fin_factor()` | `A_wl = dx w_c`, `A_wg = dx (w_c + w_fin)`, `A_fin = 2 dx L_fin` | Blong Eqs. 22–24 | |
| `fin_factor()` | `L_fin = h_c + w_c/2` | Blong Eq. 27 | Fin = land plus half the closeout, up to the channel centreline (symmetry line), so the tip is adiabatic. Wetted area `A_wl + A_fin` is then the whole channel perimeter. |
| `fin_efficiency()` | `η = tanh(mL)/(mL)` | Blong Eq. 25 | Adiabatic-tip straight fin (Incropera, cited by Blong) |
| `fin_efficiency()` | `m = √(2 h_l / (k_w w_fin))` | Incropera/Bergman straight fin, `m = √(hP/(kA_c))` with P = 2 dx (two wetted side faces), A_c = w_fin dx. **Deviates from Blong Eq. 26** (changed Nov 24, your approval) | Blong's `m = √(2 h_l (w_fin + dx)/(k_w dx w_fin))` counts the front and back faces of a dx-long slice of land as wetted; they're cuts through a continuous land. That term made results depend on the station count: 100 → 3000 stations took throat φ 1.81 → 1.15, max T_wg 850 → 996 K, Δp 58 → 257 psi, and as dx → 0, η → 0. Without it, results converge: φ 1.93, 850 K, Δp 46.4 → 45.7 psi over the same range. At 100 stations Blong's form gave m about 13% high. |
| `fin_factor()` | `q_l = h_l (T_wl − T_l)(η A_fin + A_wl)/A_wg`, returned as `φ = (η A_fin + A_wl)/A_wg` | Blong Eq. 21 | `station_fluxes()` uses `q_l = h_l φ (T_wl − T_l)` |
| `fins.py` | 1D fin, constant k_w, same h_l on every wetted face | Blong | Blong reports fins carried at least 60% of the liquid-side heat flux in his tests. |
| `advance_coolant()` | coolant heating still `q b ds / (ṁ cp)` | Book Eq 6.39 | Not Blong Eq. 29 (`q dx (w_c + 2 h_c η)`). See Problems 9. |
| `fins.py` | Hand check at the throat (N 85, w 1.0 mm, T_l 314 K): h = 2.45 mm, η = 0.684, φ = 2.66, h_l φ = 46.2 kW/m²K = required | Scratch calculation, Nov 23 | Without fins the same station needs h = 0.84 mm. |
| — | Structural check (book step 12) | — | Removed. Left to FEA. |

## 9. Jacket inputs: `regen.yaml` (`jacket`), read by `jacket/jacket_inputs.py`

| In code | Value | Source | Note |
|---|---|---|---|
| `counterflow` | True (coolant enters at the nozzle exit) | **No source** | Confirm you want this |
| `p_l_in` | 797 psi | You | Your note: "pump outlet before regen cooling". Was 850 on Oct 11. |
| `T_l_in` | 293 K | **No source** | |
| `t_w` | 1.0 mm | **No source** | Was 2.0 mm, then 0.7 mm; 1.0 mm since Nov 21. |
| `k_wall` | 284 W/m-K | Blong 2024, section on wall conduction (Fig. 4.6) | L-PBF GRCop-42, average over Blong's observed wall temperatures (curve is fairly flat). Copper only: an Al chamber needs its own value. Was 350, 300, 320. Also used as the fin conductivity (section 8a). |
| `T_wg_allow` | 850 K | **No source** | Was 800 K, then 750 K, 700 K, back to 750 K (all Nov 21), 850 K by Nov 23. Needs a GRCop-42 limit (strength/creep at temperature). |
| `land_min` | 0.8 mm | **No source** | Should come from whoever machines or prints it |
| `fin_conduction` | true | You | Switch for section 8a. |
| `entrance_effects` | true | You | true = developing flow from the coolant inlet (VDI G1 Eqs. 23, 28, section 7); false = fully developed everywhere (laminar Nu 4.354), the pessimistic bound. |
| `eps_rough` | 0.08 mm | Blong 2024 (vendor-reported, printed copper) | Blong quotes 0.0035" = 0.089 mm (Outerbody-1, AM) and 0.002" = 0.051 mm (Quadrus). 0.08 mm sits near the AM value. Vendor absolute roughness, possibly Ra rather than sand-grain. Was 0, 0.20 (Nov 21), 0.05, 0.06, 0.07 mm. Needs a printed-surface value; equivalent sand-grain roughness, not Ra. >0 switches friction to Colebrook (`fluids`, exact). |

## 10. Channel solver: `jacket/sizing.py`, run by `scripts/channel_solver.py`

| In code | What | Source | Note |
|---|---|---|---|
| `height_for_target()` | Holds T_wg at the target and solves for the channel height at each station | **No source** | Not a book method. The book only says height can vary along the length (p.211). The final heights are re-run through the normal march, so the temperatures it reports are from the book procedure. |
| `height_for_target()` | Picks the tallest height that meets the target | **No source** | Tallest = lowest pressure drop |
| `height_for_target()` | Condition `h_l(h) φ(h) = q/(T_wl − T_l)`; scans `N_SCAN = 60` log-spaced heights, takes the tallest that meets it, refines with brentq between that point and the next | **No source** | Needed with fins: `h_l φ` can peak inside the range, so a single bracket on [h_min, h_max] isn't safe. With fins off the result matches the old single brentq (station table identical). |
| `height_for_target()` | Status `best`: no height meets the target, the peak of `h_l φ` is used | **No source** | Only possible with fins on |
| `T_wg_target` | 750 K | **No source** | Same as `T_wg_allow` (`null` in regen.yaml) |
| `w_chan` | 1.0 mm | **No source** | Was 2.0 mm, then 0.5 mm; 1.0 mm since Nov 21. |
| `N_channels` | None = use `channels_from_land()` | **No source** | See section 8 |
| `h_min`, `h_max` | 1.0 mm, 8.0 mm | **No source** | Your note: "currently arbitrarily chosen". h_min was 0.25 mm. h_max was 9.5, 4.0, 6.0, 9.0 mm. |

---

## Problems found (not fixed)

1. **Viscosity is extrapolated almost everywhere.** In the current wall temp run, the coolant-side wall is above 480 K at every station. So wall viscosity is the 480 K value throughout, and that value is already extrapolated beyond the correlation's 355 K limit.
2. **The coolant hits IPA's critical temperature** (508 K) at the injector end. For the first 60 mm every property is frozen at its 480 K value. That's why `h_l` is flat there.
3. **The coolant is above IPA's critical pressure.** Near-critical behaviour isn't modelled.
4. **Kerosene constants are used for IPA** in Sieder-Tate.
5. **CEA sees fuel at 293 K**, but in reality it arrives ~500 K from the jacket. Small effect, conservative.
6. **The c* efficiency shrinks the throat but doesn't cool the gas.** Conservative.
7. **`b` is measured at the channel floor**, overstating coolant heating by ~4%.
8. **The recovery temperature goes above the stagnation temperature downstream of the throat.** (Numbers below are from the Oct 11 inputs and the old cone; not rechecked.) T0 is 3339 K, but T_r reaches 3518 K at the exit. The book says this can't happen (p.196). Cause: the ideal-gas T_r formula is fed shifting-equilibrium T, M and γ from CEA. About +6 K at the throat and about +13% driving ΔT at the exit. Options: `T_r = T + r(T0 − T)`, a frozen CEA expansion, or accept it as conservative.

9. **Fin area conventions differ (fins on).** Blong references the flux to `A_wg = dx (w_c + w_fin)` (pitch at the land's mid-height). The coolant heating uses `b` at the channel floor (Book Eq 6.39). About 2% apart at the throat. Blong's own temperature rise (Eq. 29) uses `q dx (w_c + 2 h_c η)`, which isn't the same area as his Eq. 21. The code keeps the book's energy balance (heat in from the gas = heat into the coolant).
10. **Low-Re stations** (Re < 10⁴, near the coolant inlet). Was: Gnielinski extrapolated there. Fixed Nov 24 with the VDI transition and laminar range (section 7). Remaining: G1 is for round tubes; laminar Nu in rectangular channels heated on one side differs (Shah & London; VDI G2). Rough walls likely trip turbulence below Re 2300, which no correlation here captures.

## What you need to decide or cite

- Contraction ratio
- Which paper L* came from
- Why pa = 18 psi
- Sources for the c* efficiency (0.875) and Cf efficiency (0.97)
- O/F: keep 1.65, or optimize
- Contour: 30° inlet, arc radii 1.65 Rt and 0.382 Rt, sharp cylinder→cone corner. The 80% bell's θn/θe table still needs checking against Rao (1958).
- Keep or drop the radial CSV / `.dat` outputs now that the DXF uses perpendicular offsets
- Wall material: conductivity and allowable temperature, from a datasheet
- Manufacturing limits: wall thickness, channel width, minimum land, height range, roughness
- Sieder-Tate constants for IPA
- An IPA property source that's valid up to ~500 K
- Counterflow and coolant inlet temperature
