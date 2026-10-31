# Pumps (MARLIN V2 turbopump)

Preliminary design code for the LOX and IPA pumps: 1D impeller sizing and H-Q curve, inducer sizing and radial equilibrium, the impeller meridional contour with an FEM potential-flow solve, impeller blade design (inlet angles and plan-view wrap) and a CadQuery model of the impeller, plus standalone inducer-sizing and volute scripts. Design-point inputs are in `Inputs/Inputs.yaml`: shaft speed, and mass flow and feed/discharge pressure for each propellant. Impeller geometry, meridional geometry, blade parameters and fluid properties are set in `Code/Main.py`.

The git repository root is `Turbopump/`, and this folder is `MARLIN V2/Tenkara/4 - Subteams/Pumps` in it.

## Running

- Use the local venv `env\Scripts\python.exe` (Python 3.13). It is not in git, and it has numpy, scipy, pint, plotly, matplotlib, CoolProp, PyYAML and cadquery. The system Python has numpy, scipy, pint and plotly but not CoolProp, PyYAML, matplotlib or cadquery, so `Main.py` fails on it. `requirements.txt` is UTF-16 encoded.
- The entry point is `Code/Main.py`: `MakeDesign(Propellant.LOX)` or `MakeDesign(Propellant.IPA)`, which returns the `Pump`. Running Main.py designs both, LOX first. `Code/README.md` says the code is designed to run from Main.py.
  - `Propellant` is the enum in `Pump.py`, imported into Main.py. `CreatePump` and `MakeDesign` take it, not a string, and Main.py's per-propellant dicts are keyed by it. `CreatePump` reads `Inputs[propellant.name]`, so the top-level keys of `Inputs.yaml` must match the enum names.
  - Imports need `Code/` as the script directory or on `sys.path`:
    - top-level files in `Code/` import each other by bare name (`import Pump`);
    - `Impeller/` and `Inducer/` are namespace packages (no `__init__.py`), imported as `Impeller.X` and `Inducer.X`;
    - files in `Impeller/` import each other relatively (`from .MeridionalContour import ...`), so they only work imported as `Impeller.X`, or run with `python -m` from `Code/`;
    - files in `Inducer/` import each other as `Inducer.X`. A script's own folder comes first on `sys.path`, and from there `import Inducer` would find `Inducer/Inducer.py`. So the runnable ones take their folder off `sys.path` and put `Code/` first: `RadialEquilibrium.py` and `MeanlineTriangle.py` at import, `Sizing.py` and `1DPerformance.py` in `__main__`.
  - `Impeller` and `Inducer` each name a package, a module in it, and a class. `Pump.py` and Main.py do `import Impeller.Impeller as Impeller`, so there `Impeller` is the module and `Impeller.Impeller` the class.
  - Main.py adds the Pumps folder to `sys.path` for `from Inputs import Inputs`.
- `MakeDesign` runs the pipeline in this order:
  1. `CreatePump` (overall efficiency 0.50), `pump.Impeller.PlotOutletVelocities` and `Pump1D.PlotHQ`. It prints the torque and the inducerless NPSH_r;
  2. `Pump1D.GetInducerHeadRequirement`, then `pump.AddInducer`. The impeller's `H_rise` becomes the total head rise minus the inducer head;
  3. `InducerRadialEquilibrium(pump.Inducer)` with its defaults;
  4. `MeridionalDesign`, which includes the FEM potential-flow solve, then `PotentialFlow.Plot` and `PlotMeridional`;
  5. `ImpellerBladeDesign` with Main.py's `blade_params`, then `PlotInletTriangles` and `PlotBladeLines`;
  6. `ImpellerCAD`, which builds the CadQuery model and opens it in the cadquery viewer.

  The inducer is sized on its own (Gülich tip diameter), so its tip and exit hub diameters don't match the impeller's `d_1` and `d_hub`. Only its outlet angular momentum reaches the impeller, as free-vortex preswirl. `Inducer/Sizing.py` matches the two, but the pipeline doesn't use it.

  Figures are shown by the `.show()` calls there, most of them commented out (currently all but `PlotBladeLines`). The cadquery viewer blocks until its window is closed, so the IPA design starts only after the LOX viewer is closed. A run takes a few seconds per propellant, about half of it in `ImpellerBladeDesign`; `MeridionalDesign`, FEM solve included, takes about 0.01 s.
- There is no test suite; `test.py` in the Pumps folder is a pint scratch file. To check changes without windows popping up:
  - set `plotly.graph_objects.Figure.show` to a no-op in a scratch script and use `fig.write_html(...)` instead. `MeridionalMesh.Plot()` calls `.show()` itself;
  - set `cadquery.vis.show` to a no-op before `MakeDesign` reaches the CAD step. Main.py imports `Impeller.ImpellerCAD` inside `MakeDesign`, and that module binds `show` at import (`from cadquery.vis import show`);
  - there is no kaleido, so plotly figures can't be exported as images;
  - set `MPLBACKEND=Agg` for the matplotlib scripts.
- Standalone scripts, not part of the `MakeDesign` pipeline:
  - `Inducer/Sizing.py`: sizes the inducer inlet from the impeller inlet. With `USE_PUMP = True`, its `__main__` takes the impeller inlet, flow, speed and NPSH_a from the pumps in Main.py, and lets the impeller inlet move by up to 20 %;
  - `Inducer/1DPerformance.py`: runs `OptimizeImpellerInlet` on hard-coded LOX and IPA cases and prints their 1D performance. Its name starts with a digit, so an `import` statement can't load it;
  - `Inducer/MeanlineTriangle.py`: prints and plots the inducer inlet meanline triangles;
  - `Inducer/BladeThickness.py` (formerly `Sizing/Blade_Param_Structural_Resilience.py`): asks for its inputs with `input()`, matplotlib;
  - `Inducer/RadialEquilibrium.py`: matplotlib `__main__`. It fails, because it still calls `Main.CreatePump('LOX')` with a string;
  - `Volute.py`: imports Main. Its `__main__` fails the same way (`Main.MakeDesign('LOX')`); fixed, it would run the whole pipeline, CAD viewer included. The `Volute` class is unfinished;
  - `VolutePlot.py`: matplotlib, hard-coded inputs. It saves `volute_geometry.png` to the working directory, and a copy is tracked in the Pumps folder;
  - `Impeller/JoinParabola.py`: plotly demo `__main__`.

  `Code/InducerSizing.py` is an older version of `SizeInletForImpeller` in `Inducer/Sizing.py`. It fails on import, because it imports from the old flat `Inducer.py`. `Code/Old/` is archived inducer code whose imports no longer resolve.

## Code layout (`Code/`)

Top level:

- `Main.py`: dicts keyed by `Propellant`: `basic_impeller_geometry` (`Impeller.InputGeometry`), `meridional_geometry` (`MeridionalDesign.MeridionalGeometry`), `blade_params` (`BladeParams`) and `fluid_props` (LOX vapor pressure from CoolProp). Also `CreatePump` and `MakeDesign`.
- `Pump.py`: `Propellant`, `PumpDesignPoint` and `Pump`, which holds `.Propellant`, `.Impeller` and, after `AddInducer`, `.Inducer`. It also has the design head rise and flow, `NPSH_a`, `n_s` and `n_q`.
- `Pump1D.py`: Euler head, power and torque, `PlotHQ` (empirical H-Q curve), `GetInducerlessNPSH_r` (CFturbo), `GetInducerHeadRequirement` (rounded up to a multiple of 5 m).
- `EmpiricalRelations.py`: flow-speed ratio, PUMPA slip and efficiency correlations, the volute K_v fit.
- `VelocityTriangle.py`: `VelocityTriangle(u, c_m, c_u)`, `VelocityTriangle.FromBeta`, `PlotTriangles`, `CONVENTION_NOTES`. An optional `incidence` (blade angle = β + i) makes `PlotTriangles` draw the blade direction and the arc i, as in Gülich Fig. 3.1; `PlotInletTriangles` uses it on the 1′ triangle.
- `Volute.py`, `VolutePlot.py`, `InducerSizing.py`: see Running.

`Impeller/`:

- `Impeller.py`: `InputGeometry` (`Z_blade`, `Beta2B`, `d_1`, `d_2`, `b_2`, `blade_thick`, `d_hub`), `DesignPoint` (`Q`, `N_shaft`, `H_rise`) and `Impeller`. `Impeller` covers Wiesner slip, exit blade blockage, the outlet velocities and triangles (`PlotOutletVelocities`), and the dimensionless numbers (`PrintDimensionlessParameters`).
- `MeridionalDesign.py`: `MeridionalGeometry` (contour fullness, outlet angles, hub taper ratio, shroud flat, `Z_e_factor`) and `MeridionalDesign(pump, geometry)`. `MeridionalDesign.__init__` builds:
  - `Z_e`;
  - `ShroudContour`, `HubContour` and `MeanContour`;
  - `PotentialFlow`;
  - the leading edge (`LeadingEdge`, a `BezierCurve` from hub to shroud);
  - the trailing edge (`TrailingEdge`, a straight `Polyline` from the hub outlet to the shroud outlet).
  
  Plotting: `MeridionalTraces`, `PlotMeridional`.
- `MeridionalContour.py`: meridional curves, the FEM domain and the mesh.
  - `MeridionalCurve` is the common base:
    - `curve(s)` returns (z, r) Quantities at normalized arclength;
    - `Points(s_m)` returns an (n, 2) array in m at true arclength;
    - `Length`;
    - `Tangent(s)` returns the unit (z, r) tangent at normalized arclength;
    - `StaticMoment(s_1, s_2)` returns ∫r dx;
    - `FindIntersect(other)` returns an `IntersectResults`: a list of `(z, r, s_1, s_2)` tuples ordered by `s_1` (normalized arclength on this curve, `s_2` on `other`), empty when the curves don't cross. A start or end of one curve lying on the other (within 1 µm) counts, with `s` exactly 0 or 1. Overlapping stretches and a curve with itself are not handled;
    - `Closest(p)` returns the distance from a point (Quantity or m) to the curve and the normalized arclength of the closest point;
    - `Crop(s_1, s_2)` returns the part between two normalized arclengths as a `SubCurve`; `Between(a, b)` the part between where the curve crosses `a` and `b` once each, e.g. `streamline.Between(LeadingEdge, TrailingEdge)` for the blade part of a streamline;
    - `Trace(l, **scatter)` returns the curve as a plotly line.
  - Subclasses (a new kind of curve only needs `length_m` and `Points`):
    - `Polyline`: straight lines through points; two points make a straight line. `Interpolate(values, s)` interpolates per-point data, such as a streamline's `v_m`, at normalized arclength the same way as the points;
    - `BezierCurve`: cubic Bezier by arclength, with `ControlPoints`. `BezierCurve.ThroughPoints(P_a, P_m, P_b)` is a quadratic through three points;
    - `SegmentChain`: any curves joined end to end, chains included. `Segments` holds the pieces and `s_m` the arclength (m) of each joint;
    - `SubCurve`: part of a `Parent` curve, from `Crop`/`Between`. `ToParent(s)` maps normalized arclength back to the parent, e.g. to interpolate a streamline's `v_m`;
    - `MeridionalContour(SegmentChain)`: a flat `Polyline`, then the bend `Bend` (a `BezierCurve`), then an outlet `Polyline`.
  - `MeanCurve(a, b)` returns the `Polyline` midway between two curves at equal normalized arclength.
  - `Cross(a, b)` is the 2D cross product, used by `FindIntersect` and the mesh.
  - `FlowDomain` extends the hub and shroud upstream and downstream for the FEM solve, with zones `inlet`/`impeller`/`outlet`. Its walls are `SegmentChain`s with the contour nested as one segment. `FlowDomain.Mesh()` builds a `MeridionalMesh`, by default 101 stations × 11 lines.
- `PotentialFlow.py`: the FEM stream-function solve (`PotentialFlow`, with `Plot`). `GetStreamline(phi)` returns a `Streamline(Polyline)` with:
  - `.x`: arclength from the impeller inlet;
  - `.v_m`: meridional velocity;
  - `.Impeller`: the impeller part only.
- `ImpellerBladeDesign.py`: `BladeParams` and `ImpellerBladeDesign(meridional_design, inducer_radial_equilibrium=None, blade_params=BladeParams())`.
  - Finds the inlet flow where the leading edge crosses streamline φ: r₁, c_m1 (the streamline's potential-flow v_m) and w_u1, using the inducer preswirl c_u = K/r (K = 0 without an inducer).
  - Adds leading-edge blade blockage and solves β₁′ explicitly. `_INCIDENCE_PADDING` (2°) is only drawn: `PlotInletTriangles` shows it as the incidence on the 1′ triangle, but the blade lines start at β₁′, at zero incidence.
  - Builds blade lines along 11 streamlines from φ = 0 to 1, each over its blade part (`Between(LeadingEdge, TrailingEdge)`) at `BladeParams.resolution` points. `BladeParams.use_cubic` picks the method:
    - `False`, Main.py's choice: `_CreateBladeLines`. β(x) is a `JoinParabola` from β₁′ to `Beta2B`, joined at x_j = sj·x_f, and θ = ∫ dx/(r tan β). The shroud is the reference, with `BladeParams.sj`. On the other streamlines sj is bisected so every trailing edge ends at the shroud's wrap (with a warning if no sj reaches it), and the leading edge is swept by `lead_sweep`;
    - `True`, the `BladeParams` default: `_CreateBladeLines2`. `ConformalMapping` puts every trailing edge at `target_wrap`, and the sweep is measured from the hub.
  - `GetBladeLines()` returns `(points, betas)`, one entry per streamline, hub first. `points` are (resolution, 3) arrays of (z, r, θ) in m, m and rad. Each entry of `betas` has `Get(x)` (β in rad at x in m from the leading edge), `Domain` and `TotalWrap`.
  - Plotting: `PlotInletTriangles`, and `PlotBladeLines`, which plots the plan view of all `Z_blade` blades and β against meridional arclength.
- `ConformalMapping.py`: `ConformalMapping(beta_1, beta_2, arclen, target_wrap, radii).SolveDistribution()` returns a `BetaCMM`. A cubic A(x) has A′ = 1/tan β at both ends, and θ = ∫ dA/r. Its end value A(x_f) is interpolated linearly for the target wrap within the range where β stays monotonic; outside that range it is clamped, with a warning.
- `JoinParabola.py`: `JoinParabola(y0, yf, xj, xf)` returns two parabolas as one function of x: flat at 0 and at xf, and tangent-continuous at xj. It is singular for xj = 0 or xj = xf, so the sj bisection stays between 0.01 and 0.99.
- `ImpellerCAD.py`: `ImpellerCAD(blade_design)` lofts the blade lines into one zero-thickness camber surface and copies it `Z_blade` times about the z axis. It revolves the hub and shroud contours into surfaces and shows the compound with `cadquery.vis.show`. It exports nothing; the STEP files in `Outputs/` came from an earlier version that has been deleted.
- `Bezier.py`: `bez(P, t)` and `dbez(P, t)` for cubic Beziers, by Bezier parameter t. Use `BezierCurve` when you need a curve evaluated by arclength.

`Inducer/` (NASA SP-8052, whose PDF is in this folder):

- `Inducer.py`: `InducerDesignPoint`, `InducerInletDesignConstraints` and `Inducer`, which does the inducer sizing (tip and hub diameter, N_ss, NPSH_3, hub taper, velocity triangles). `Inducer` defaults to the module-level `constraints` (hub-to-tip ratio 0.5, σ = 0.03), not to the dataclass defaults (0.3, 0.06).
- `RadialEquilibrium.py`: `InducerRadialEquilibrium`, which gives c_θ(r, z) and c_axial(r, z) behind the tapered-hub inducer. `GetOutletAngularMomentum()` gives the impeller inlet preswirl.
- `MeanlineTriangle.py`: `MeanlineTriangle`, the zero-blockage inlet triangle on the rms radius. It passes the incidence ratio to `InducerRadialEquilibrium`.
- `Sizing.py`: sizes the inducer inlet to match the impeller inlet.
  - `SizeInletForImpeller` takes the tip and exit hub diameters from the impeller's `d_1` and `d_hub`. It solves the inlet hub diameter from the hub taper, and the cavitation parameter from the Gülich tip diameter;
  - `OptimizeImpellerInlet` grid-searches for the largest impeller inlet, within a tolerance, that satisfies `InletSizingLimits`;
  - `SizeInducerFromPump` runs either one on a `Pump`;
  - `HeadRequirementModel` repeats `Pump1D.GetInducerHeadRequirement` for any impeller inlet;
  - `PlotContour` plots the inducer's meridional section.
- `1DPerformance.py`: `Calculate1DPerformance(sizing)` gives φ, N_ss, NPSH_3, the NPSH margin, ψ₂ and the annulus areas for a `Sizing.InletSizing`.
- `BladeThickness.py`: see Running.

`Old/`: archived inducer code, see Running.

Outside `Code/`: `Inputs/` holds `Inputs.yaml`, and `Serialize.py`, which turns every `[value, 'unit']` pair into a `Q_`. `Inputs/__init__.py` exports the result as `Inputs`. `Outputs/` holds `LOX_impeller.step` and `IPA_impeller.step` (see `ImpellerCAD.py`).

## Conventions

- File split: geometry goes in `Impeller/MeridionalContour.py`, flow theory goes in `Impeller/PotentialFlow.py`. `FlowDomain` and `MeridionalMesh` count as geometry because they contain no flow equations. The one deliberate exception is `QuasiOrthogonalFlow` (1D continuity), which stays with the mesh because `FlowDomain.Mesh(by_flow=True)` uses it; moving it would make the two files import each other.
- Try to reuse code
- Use pint `Q_` at interfaces. Curve, mesh and blade-line internals are plain floats in meters and radians, and so are the `ConformalMapping` and `JoinParabola` arguments; convert with `.to(...)` and `.m`.
- Meridional coordinates are always ordered (z, r): z is axial and is the plot x-axis, r is radial. Blade lines are (z, r, θ). The CAD model uses x = r cos θ, y = r sin θ and z, in mm.
- Normalized arclength s runs from 0 at the impeller inlet to 1 at the outlet on the contours and the meanline, and from the hub (0) to the shroud (1) on the leading and trailing edges. Streamlines span the whole FEM domain; use `Streamline.Impeller` for the impeller part. Along a blade line, x is the meridional arclength in m from the leading edge.
- `phi` in the code is the flow fraction (stream function / its shroud value: 0 at the hub, 1 at the shroud). It is not Milos's velocity potential φ. In the inducer files, `phi` is the inducer flow coefficient instead.
- `psi` is the stream function in `PotentialFlow`, but the swirl ratio c_u/(ω r) in `RadialEquilibrium.py`; `psi_1` and `psi_2` are inducer head coefficients.
- In `Inducer/Sizing.py`, `SizingConstants.lambda_w` and `lambda_c` (0.03 and 1.1) are σ and K from `Inducer.py`. `HeadRequirementModel.lambda_w` and `lambda_c` (1.0 and 1.1) are the impeller NPSH coefficients of `Pump1D.GetInducerlessNPSH_r`.
- Velocities:
  - `c` is the absolute-frame velocity in the sizing and velocity-triangle code (`c_m`, `c_u`). The inducer files also use SP-8052's `c_theta` and `c_axial`;
  - `PotentialFlow` deliberately keeps `v_z`, `v_r`, `v_m`.
- Blade and flow angles β are measured from the tangential direction, as in `VelocityTriangle`: tan β = c_m / (u − c_u). Some references measure them from the meridional direction instead, such as NASA TN D-1562.
- Style:
  - match the file's indentation. Most of `Code/` uses tabs. These files use spaces: `EmpiricalRelations.py`, `Volute.py`, `VolutePlot.py`, `Impeller/Bezier.py`, `Old/`, and `Inducer.py`, `RadialEquilibrium.py`, `MeanlineTriangle.py` and `BladeThickness.py` in `Inducer/`. The other two files in `Inducer/`, `Sizing.py` and `1DPerformance.py`, use tabs;
  - CamelCase method names (`GetX`, `PlotX`) and object attributes (`HubContour`). The inducer files and `Volute.py` also have snake_case and camelCase names (`c_theta`, `print_parameters`, `getKv`); follow the file;
  - cite the source equation or page in docstrings.
- Plots use plotly, mostly with `template='plotly_white'`; `PlotInletTriangles` uses `'simple_white'`. A length-unit argument `l` defaults to `'in'`; `PlotBladeLines` calls it `unit`, and the `Inducer/Sizing.py` plots default to `'mm'`. Meridional and plan views use equal axis scaling.

## Design method and references

The main reference is M. Milos, "Impeller Design Using CAD Techniques and Conformal Mapping Method", InTech, *Centrifugal Pumps* (2012): https://cdn.intechopen.com/pdfs/29651/InTech-Impeller_design_using_cad_techniques_and_conformal_mapping_method.pdf. Page numbers below are the printed ones.

To read the PDF: WebFetch can't parse it, no PDF library is installed, and the Read tool can't render PDFs here (no poppler). Download it and inflate its text streams with zlib in a small script.

- **Domain and mesh**, pp. 35-36: 101 × 11 nodes with straight quasi-orthogonals. `MeridionalDesign` currently calls `Mesh(by_flow=False)`, which divides each QO equally, as Milos does.
- **Stream-function FEM**, pp. 36-38: eq. 4, boundary conditions in eqs. 7-13, Galerkin method in eqs. 14-16, then streamlines by cubic-spline interpolation of ψ.
  - Milos computes velocities at element centres (eq. 17).
  - Our solve differs: it uses the symmetric divergence form, a direct sparse solve, ψ from 0 to 1, nodal finite-difference velocities, and PCHIP instead of a cubic spline.
- **Inlet edge**, pp. 40-41, Figs. 7-8: Milos places it by judgment from the velocity plots.
- **Static moment**, p. 45, eqs. 34-35: S = ∫r dx feeds Pfleiderer's slip coefficient p = ψ′ r₂²/(z S).
- **Leading edge (our choice, not Milos's):**
  - the hub, the geometric meanline and the shroud have equal S from the leading edge to the outlet, which gives equal slip from hub to shroud;
  - the shroud leading edge sits at the impeller inlet (`s_LE_shroud = 0`), and the other two are found by bisection;
  - the curve is a quadratic Bezier through the three points, stored as a cubic (`BezierCurve.ThroughPoints`);
  - the geometric meanline lies about 0.4-0.6 mm closer to the hub than the true φ = 0.5 streamline at the impeller inlet, and 0.8-1.1 mm from it where both cross the leading edge, for the current LOX and IPA geometry. Using `self.PotentialFlow.GetStreamline(0.5).Impeller` instead of `self.MeanContour` is a small swap in `MeridionalDesign.__init__`.
- **Leading-edge blockage** (`ImpellerBladeDesign._SolveInletFlow`): CFturbo; Hassan et al., Sci. Rep. 2024, eq. 17; Milos eq. 33.
- **Blade angles and wrap** (`ImpellerBladeDesign._CreateBladeLines`), p. 51, eqs. 51-52, as the code cites them: β runs from inlet to outlet as two parabolas joined at Milos's x_3, and the wrap angle is θ = ∫ dx/(r tan β).
- **Limits of the potential flow:** inviscid, irrotational, axisymmetric, no blades. It is fine for contour shaping and streamlines, but not for head, losses or off-design. The ~2 mm outlet width b₂ is dominated by viscosity.
- **NASA TN D-1562** (Stockman & Kramer, 1963), "Method for Design of Pump Impellers Using a High-Speed Digital Computer". It is an inverse method: from the hub, the hub relative velocity and the blade surface θ(r, z), it marches streamline by streamline to the *shroud*. Its β is measured from the meridional direction. We don't use the marching itself, because our shroud comes from Milos, and no code uses this report yet. The parts that apply to us are:
  - eq. 9: tan β = r dθ/dm along a streamline, which is the wrap-angle (plan-view) relation. It is the relation of Milos eqs. 51-52, with β from the meridional direction;
  - eq. 10: circumferential blade thickness t_θ. The blades have no thickness yet; `blade_thick` only enters the blockage;
  - App. C, eqs. C3-C7: approximate pressure- and suction-side velocities. C3 reduces to g Δh = (2πr − Z t_θ)/Z · w_m · (1/r) d(r c_u)/dm.
- **Other references:** Gülich, *Centrifugal Pumps* (axial extension Z_e p. 348, NPSH ch. 6, the inducer tip diameter, the volute pp. 403 and 420); Wiesner slip (CFturbo manual); PUMPA (slip correlation); Stepanoff (volute); NASA SP-8052, in `Code/Inducer/` (inducer radial equilibrium, meanline triangle and blade thickness).
