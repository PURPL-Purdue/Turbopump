import numpy as np
import matplotlib.pyplot as plt
from matplotlib.patches import Circle
from matplotlib.lines import Line2D
from scipy.optimize import minimize_scalar
import yaml

# ─────────────────────────────────────────────
# UNIT CONVERSIONS
# ─────────────────────────────────────────────
psi_into_pa         = 6894.76
pa_to_psi           = 1 / psi_into_pa
meters_into_inches  = 39.37
inches_into_meters  = 0.0254
kg_into_lbm         = 2.20462


# Import TCA Parameters
with open(r'Inputs/TCA_params.yaml') as file:
    tca_params = yaml.safe_load(file)

# ─────────────────────────────────────────────
# DESIGN PARAMETERS
# ─────────────────────────────────────────────
mdot_total   = tca_params["tp_mdot"]           # Total propellant mass flow [kg/s]
OF_Ratio     = tca_params["of_ratio"]          # O/F ratio
mdot_f       = mdot_total / (1 + OF_Ratio)     # Fuel mass flow [kg/s]
rho_f        = tca_params['densities']['fuel'] # Fuel density [kg/m3]
Pc           = tca_params['chamber_pressure']  # Chamber pressure [psi]
P_in_psi = Pc * (1 + tca_params["stiffness"]["fuel"]/100) # fuel inlet pressure

# Injector element layout
n_type1      = 6            # large mini manifolds
holes_type1  = 4            # holes per large mini manifold

n_type2      = 6            # small mini manifolds
holes_type2  = 2            # holes per small mini manifold

n_holes_total = n_type1 * holes_type1 + n_type2 * holes_type2   # total fuel hole count

# ─────────────────────────────────────────────
# ASYMMETRIC INLET FLAG
# ─────────────────────────────────────────────
inlet_between_types = True   # <-- SET THIS: True = inlet sits between a type1 & type2 manifold

if n_holes_total != tca_params["hole_number"]["fuel"]:
    print("WARNING: HOLE COUNT DOES NOT MATCH TCA HARDWARE DEFINITION!")
    print(f"YAML: {tca_params['hole_number']['fuel']}")
    print(f"Inputs: {n_holes_total}")
    print("Review Inputs.")
else:
    mdot_per_hole = mdot_f / n_holes_total

    # Manifold geometry  ← FIXED OUTER, ECCENTRIC INNER
    R_outer      = 2.35 * inches_into_meters    # Outer circle radius [m]  ← fixed
    h_manifold   = 0.238 * 2.8 * inches_into_meters    # Manifold height (axial depth) [m]  ← constant

    # Dynamic pressure assumption (sets velocity)
    dyn_pressure_fraction = 0.005

    # ─────────────────────────────────────────────
    # DERIVED FLOW QUANTITIES
    # ─────────────────────────────────────────────
    dyn_pressure = dyn_pressure_fraction * P_in_psi * psi_into_pa
    v_manifold   = np.sqrt(2 * dyn_pressure / rho_f)   # target constant velocity, both branches

    mdot_branch  = mdot_f / 2.0
    A_inlet      = mdot_branch / (rho_f * v_manifold)
    width_inlet  = A_inlet / h_manifold

    print("=" * 55)
    print("DUALLY / AVERAGED ECCENTRIC MANIFOLD DESIGN")
    print("=" * 55)
    print(f"Fuel mass flow (total):   {mdot_f*kg_into_lbm:.4f} lbm/s  |  {mdot_f:.4f} kg/s")
    print(f"Mass flow per branch:     {mdot_branch*kg_into_lbm:.4f} lbm/s  |  {mdot_branch:.4f} kg/s")
    print(f"Target manifold velocity: {v_manifold:.3f} m/s  ({v_manifold*3.28084:.3f} ft/s)")
    print(f"Inlet area (per branch):  {A_inlet*1e6:.3f} mm2  ({A_inlet*meters_into_inches**2:.4f} in2)")
    print(f"Inlet width:              {width_inlet*1e3:.3f} mm  ({width_inlet*meters_into_inches:.4f} in)")
    print(f"Inlet condition:          {'BETWEEN type1 & type2 (asymmetric)' if inlet_between_types else 'ALIGNED with one type (symmetric)'}")

    # ─────────────────────────────────────────────
    # ASYMMETRIC BRANCH HOLE SEQUENCES
    # (each branch alternates type1/type2; asymmetric inlet means the
    #  two branches encounter opposite first types)
    # ─────────────────────────────────────────────
    n_type1_branch = n_type1 // 2
    n_type2_branch = n_type2 // 2
    assert n_type1_branch == n_type2_branch, \
        "This layout logic assumes equal per-branch type counts for 1:1 alternation."

    def build_branch(start_with, n1, n2, h1, h2):
        seq = []
        for i in range(max(n1, n2)):
            if start_with == "type2":
                if i < n2: seq.append(h2)
                if i < n1: seq.append(h1)
            else:
                if i < n1: seq.append(h1)
                if i < n2: seq.append(h2)
        return seq

    if inlet_between_types:
        branch_A_holes = build_branch("type2", n_type1_branch, n_type2_branch, holes_type1, holes_type2)
        branch_B_holes = build_branch("type1", n_type1_branch, n_type2_branch, holes_type1, holes_type2)
    else:
        branch_A_holes = build_branch("type2", n_type1_branch, n_type2_branch, holes_type1, holes_type2)
        branch_B_holes = build_branch("type2", n_type1_branch, n_type2_branch, holes_type1, holes_type2)

    n_steps = len(branch_A_holes)
    theta_elements = np.linspace(0, np.pi, n_steps + 2)[1:-1]

    def required_widths(holes_per_step, v_target):
        """Same as original single-manifold code, but also returns areas (mm2 info)."""
        mdot_remaining = mdot_branch
        areas, mdots = [], []
        for n_holes in holes_per_step:
            areas.append(mdot_remaining / (rho_f * v_target))
            mdots.append(mdot_remaining)
            mdot_remaining -= n_holes * mdot_per_hole
        areas = np.array(areas)
        mdots = np.array(mdots)
        widths = areas / h_manifold
        return widths, mdots, areas

    widths_A, mdot_at_A, areas_A = required_widths(branch_A_holes, v_manifold)
    widths_B, mdot_at_B, areas_B = required_widths(branch_B_holes, v_manifold)

    def print_station_table(label, holes, theta, mdots, areas, widths):
        """Matches original code's station table: Station, th, mdot, Area, Width."""
        print(f"\n{label}  (hole sequence: {holes})")
        print(f"{'Station':>8} {'th (deg)':>9} {'mdot (kg/s)':>12} {'Area (mm2)':>12} {'Width (mm)':>12}")
        print("-" * 58)
        for i in range(len(holes)):
            print(f"{i+1:>8} {np.degrees(theta[i]):>9.1f} {mdots[i]:>12.5f} "
                  f"{areas[i]*1e6:>12.3f} {widths[i]*1e3:>12.3f}")

    print_station_table("Branch A (inlet-adjacent: type2)", branch_A_holes, theta_elements, mdot_at_A, areas_A, widths_A)
    print_station_table("Branch B (inlet-adjacent: type1)", branch_B_holes, theta_elements, mdot_at_B, areas_B, widths_B)

    # ─────────────────────────────────────────────
    # ECCENTRIC FIT — same method as original code
    # (fixed outer, eccentric inner, single free parameter e,
    #  mean_gap tied to inlet width via mean_gap + e = width_inlet)
    # ─────────────────────────────────────────────
    def gap_model(theta, mean_gap, e):
        """Channel width at angle theta for eccentric inner circle."""
        return mean_gap + e * np.cos(theta)

    def fit_eccentric(widths, theta):
        def fit_error(e):
            mean_gap = width_inlet - e
            if mean_gap - e <= 0:       # gap would close at 180
                return 1e10
            if mean_gap - e > R_outer:  # inner radius would be negative
                return 1e10
            predicted = gap_model(theta, mean_gap, e)
            return np.sqrt(np.mean((predicted - widths) ** 2))
        result = minimize_scalar(fit_error, bounds=(0, width_inlet * 0.99), method='bounded')
        e_best = result.x
        mean_gap = width_inlet - e_best
        R_inner = R_outer - mean_gap
        assert R_inner > 0, "Inner radius went negative — reduce R_outer or h_manifold"
        assert mean_gap - e_best > 0, "Gap closes before 180 — increase R_outer"
        return e_best, mean_gap, R_inner, result.fun

    def print_geometry_header(R_inner, e_best, mean_gap):
        print(f"Outer circle radius R_outer:  {R_outer*1e3:.3f} mm  (fixed)")
        print(f"Inner circle radius R_inner:  {R_inner*1e3:.3f} mm")
        print(f"Eccentricity e:               {e_best*1e3:.3f} mm  (inner center shift)")
        print(f"Mean gap:                     {mean_gap*1e3:.3f} mm")
        print(f"Gap at inlet  (th=0 deg):     {gap_model(0, mean_gap, e_best)*1e3:.3f} mm  "
              f"[target: {width_inlet*1e3:.3f} mm]")
        print(f"Gap at outlet (th=180 deg):   {gap_model(np.pi, mean_gap, e_best)*1e3:.3f} mm")

    def print_station_error_table(label, widths_req, theta, mean_gap, e_best):
        print(f"\n{label}")
        print(f"{'Station':>8} {'th (deg)':>9} {'Required (mm)':>15} {'Model (mm)':>12} {'Error (mm)':>12}")
        print("-" * 60)
        sq_err = []
        for i in range(n_steps):
            w_model = gap_model(theta[i], mean_gap, e_best)
            err = (w_model - widths_req[i]) * 1e3
            sq_err.append(err**2)
            print(f"{i+1:>8} {np.degrees(theta[i]):>9.1f} "
                  f"{widths_req[i]*1e3:>15.3f} "
                  f"{w_model*1e3:>12.3f} "
                  f"{err:>12.4f}")
        rms = np.sqrt(np.mean(sq_err))
        print(f"RMS fit error: {rms:.4f} mm")
        return rms

    # ═══════════════════════════════════════════════════════
    # PART 1: DUALLY ECCENTRIC MANIFOLD (fit computed, not printed —
    # feeds into the Part 2 average below)
    # ═══════════════════════════════════════════════════════
    e_A, mean_gap_A, R_inner_A, err_A = fit_eccentric(widths_A, theta_elements)
    e_B, mean_gap_B, R_inner_B, err_B = fit_eccentric(widths_B, theta_elements)

    # ═══════════════════════════════════════════════════════
    # PART 2: AVERAGE (SINGLE) ECCENTRIC MANIFOLD
    # Average of the two independent fits, applied to both branches.
    # ═══════════════════════════════════════════════════════
    e_avg        = (e_A + e_B) / 2
    mean_gap_avg = (mean_gap_A + mean_gap_B) / 2
    R_inner_avg  = R_outer - mean_gap_avg

    print(f"\n{'='*55}")
    print("PART 2: AVERAGE (SINGLE) ECCENTRIC GEOMETRY")
    print(f"{'='*55}")
    print_geometry_header(R_inner_avg, e_avg, mean_gap_avg)

    # One shared wall — geometry printed once above; only the fit against
    # each branch's required widths differs, so only the station tables repeat.
    print_station_error_table("vs. Branch A required widths", widths_A, theta_elements, mean_gap_avg, e_avg)
    print_station_error_table("vs. Branch B required widths", widths_B, theta_elements, mean_gap_avg, e_avg)

    # ─────────────────────────────────────────────
    # PLOTS — same two-panel style as original code
    # ─────────────────────────────────────────────
    def plot_manifold_design(title, filename, geoms, widths_list, labels, colors):
        """
        geoms: list of (mean_gap, e, R_inner) tuples, one per wall drawn
        widths_list: list of required-width arrays, one per branch (for scatter)
        labels: list of branch labels (for scatter legend)
        colors: list of colors, one per geom / branch pair
        """
        fig, axes = plt.subplots(1, 2, figsize=(14, 6))
        fig.suptitle(title, fontsize=13, fontweight='bold')

        # ── Plot 1: Width profile ──
        ax1 = axes[0]
        theta_full = np.linspace(0, np.pi, 300)
        for (mean_gap, e, R_inner), color, label in zip(geoms, colors, labels):
            widths_model = gap_model(theta_full, mean_gap, e)
            ax1.plot(np.degrees(theta_full), widths_model * 1e3, color=color,
                     linewidth=2, label=f'{label} fit wall')
        for w_req, color, label in zip(widths_list, colors, labels):
            ax1.scatter(np.degrees(theta_elements), w_req * 1e3, color=color,
                        zorder=5, s=70, label=f'{label} required')
        ax1.set_xlabel("Angular position th [deg]", fontsize=12)
        ax1.set_ylabel("Manifold channel width [mm]", fontsize=12)
        ax1.set_title("Channel Width vs. Angular Position", fontsize=12)
        ax1.legend(fontsize=8)
        ax1.grid(True, alpha=0.3)
        ax1.set_xlim(0, 180)

        # ── Plot 2: Top-view cross section ──
        ax2 = axes[1]
        ax2.set_aspect('equal')
        theta_circle = np.linspace(0, 2 * np.pi, 500)

        ax2.plot(R_outer * np.cos(theta_circle) * 1e3, R_outer * np.sin(theta_circle) * 1e3,
                 'k-', linewidth=2.5, label='Outer wall (fixed)')
        ax2.plot(0, 0, 'k+', markersize=12, markeredgewidth=2)

        for (mean_gap, e, R_inner), color, label in zip(geoms, colors, labels):
            inner_cx = -e * 1e3
            ax2.plot(inner_cx + R_inner * np.cos(theta_circle) * 1e3,
                     R_inner * np.sin(theta_circle) * 1e3,
                     color=color, linewidth=2.5,
                     label=f'{label} inner (R={R_inner*1e3:.2f}mm, e={e*1e3:.2f}mm)')
            ax2.plot(inner_cx, 0, '+', color=color, markersize=12, markeredgewidth=2)

            # Element stations on inner wall (upper half)
            for th in theta_elements:
                x_pt = inner_cx + R_inner * np.cos(th) * 1e3
                y_pt = R_inner * np.sin(th) * 1e3
                ax2.plot(x_pt, y_pt, 'o', color=color, markersize=4)

        # Inlet marker
        ax2.annotate('Inlet (max width)', xy=(R_outer * 1e3, 0), xytext=(R_outer * 1e3 + 5, 8),
                    fontsize=8, color='darkred',
                    arrowprops=dict(arrowstyle='->', color='darkred'))

        margin = 0.1 * R_outer * 1e3
        ax2.set_xlim(-(R_outer * 1e3 + margin), R_outer * 1e3 + margin)
        ax2.set_ylim(-(R_outer * 1e3 + margin), R_outer * 1e3 + margin)
        ax2.set_xlabel("x [mm]", fontsize=12)
        ax2.set_ylabel("y [mm]", fontsize=12)
        ax2.set_title("Top View - Eccentric Circle Geometry", fontsize=12)
        ax2.legend(fontsize=7, loc='lower right')
        ax2.grid(True, alpha=0.25)

        plt.tight_layout()
        plt.savefig(filename, dpi=150)
        plt.close(fig)
        print(f"Plot saved: {filename}")

    # Figure 1: Dually eccentric (both branch walls shown together)
    plot_manifold_design(
        title="Dually Eccentric Manifold Design (independent per-branch wall)",
        filename=r"Development\Injector_codes\dual_eccentric_design.png",
        geoms=[(mean_gap_A, e_A, R_inner_A), (mean_gap_B, e_B, R_inner_B)],
        widths_list=[widths_A, widths_B],
        labels=["Branch A", "Branch B"],
        colors=["tab:red", "tab:blue"],
    )

    # Figure 2: Averaged single wall (both branches' requirements shown against one wall)
    plot_manifold_design(
        title="Average (Single) Eccentric Manifold Design",
        filename=r"Development\Injector_codes\average_eccentric_design.png",
        geoms=[(mean_gap_avg, e_avg, R_inner_avg)],
        widths_list=[widths_A, widths_B],
        labels=["Averaged wall"],
        colors=["tab:green"],
    )
    # add branch B scatter manually since plot_manifold_design zips 1:1 above
    # (re-plot with both scatter sets against the single averaged wall)
    fig, ax = plt.subplots(figsize=(7, 6))
    theta_full = np.linspace(0, np.pi, 300)
    ax.plot(np.degrees(theta_full), gap_model(theta_full, mean_gap_avg, e_avg) * 1e3,
            'g-', linewidth=2, label='Averaged fit wall')
    ax.scatter(np.degrees(theta_elements), widths_A * 1e3, color='tab:red', s=70, zorder=5, label='Branch A required')
    ax.scatter(np.degrees(theta_elements), widths_B * 1e3, color='tab:blue', s=70, zorder=5, label='Branch B required')
    ax.set_xlabel("Angular position th [deg]")
    ax.set_ylabel("Manifold channel width [mm]")
    ax.set_title("Average Wall vs. Both Branches' Required Widths")
    ax.legend(fontsize=9)
    ax.grid(True, alpha=0.3)
    plt.tight_layout()
    plt.savefig(r"Development\Injector_codes\average_eccentric_both_branches.png", dpi=150)
    plt.close(fig)
    print("Plot saved: average_eccentric_both_branches.png")

    print("\nDesign complete.")