# Pumps (MARLIN V2 turbopump)

Preliminary design code for the LOX and IPA pumps: 1D impeller sizing and H-Q curve, inducer sizing and radial equilibrium, the impeller meridional contour with an FEM potential-flow solve, and impeller blade design, plus standalone volute scripts. Design-point inputs are in `Inputs/Inputs.yaml`: shaft speed, and mass flow and feed/discharge pressure for each propellant. Impeller and meridional geometry and the fluid properties are set in `Code/Main.py`.

## Running

- Use the local venv `env\Scripts\python.exe`. It is not in git, and it has numpy, scipy, pint, plotly, matplotlib, CoolProp, PyYAML and cadquery; the system Python lacks them. `requirements.txt` is UTF-16 encoded. cadquery is installed, but no current code uses it.
- The entry point is `Code/Main.py`: `MakeDesign('LOX')` or `MakeDesign('IPA')`. `Code/README.md` says the code is designed to run from Main.py.
  - Files in `Code/` import each other by bare name (`import Pump`), so `Code/` must be the script directory or on `sys.path`.
  - Main.py adds the Pumps folder to `sys.path` for `from Inputs import Inputs`.
- `MakeDesign` runs the pipeline in this order:
  1. `CreatePump`, then `Pump1D.PlotHQ`;
  2. `Pump1D.GetInducerHeadRequirement`, then `pump.AddInducer`. The impeller's `H_rise` becomes the total head rise minus the inducer head;
  3. `InducerRadialEquilibrium`;
  4. `MeridionalDesign`, which includes the FEM potential-flow solve and takes a few seconds;
  5. `ImpellerBladeDesign`.
  
  Figures are shown by the `.show()` lines there, most of which are commented out.
- There is no test suite. To check changes without windows popping up:
  - set `plotly.graph_objects.Figure.show` to a no-op in a scratch script and use `fig.write_html(...)` instead;
  - there is no kaleido, so plotly figures can't be exported as images;
  - set `MPLBACKEND=Agg` for the matplotlib scripts.
- Standalone scripts, not part of the `MakeDesign` pipeline:
  - `VolutePlot.py`: matplotlib, hard-coded inputs;
  - `Volute.py`: imports Main;
  - `Sizing/Blade_Param_Structural_Resilience.py`: asks for its inputs with `input()`;
  - `InducerRadialEquilibrium.py`: has a matplotlib `__main__`.

## Code layout (`Code/`)

- `Main.py`: per-propellant `Impeller.InputGeometry` and `MeridionalDesign.MeridionalGeometry`, fluid properties (LOX vapor pressure from CoolProp), `CreatePump`, `MakeDesign`.
- `Pump.py`: `PumpDesignPoint` and `Pump`, which holds `.Impeller` and, after `AddInducer`, `.Inducer`. It also has the design head rise and flow, `NPSH_a`, `n_s` and `n_q`.
- `Impeller.py`: `InputGeometry` (`Z_blade`, `Beta2B`, `d_1`, `d_2`, `b_2`, `blade_thick`, `d_hub`), `DesignPoint` (`Q`, `N_shaft`, `H_rise`) and `Impeller`. `Impeller` covers Wiesner slip, exit blade blockage, the outlet velocities and triangles, and the dimensionless numbers.
- `Pump1D.py`: Euler head, power and torque, `PlotHQ` (empirical H-Q curve), `GuessInducerlessNPSH_r` (CFturbo), `GetInducerHeadRequirement` (rounded to 5 m).
- `EmpiricalRelations.py`: flow-speed ratio, PUMPA slip and efficiency correlations, the volute K_v fit.
- `Inducer.py`: inducer sizing (tip and hub diameter, N_ss, NPSH_3, hub taper, velocity triangles).
- `InducerRadialEquilibrium.py`: c_θ(r, z) behind the tapered-hub inducer (NASA SP-8052). `GetOutletAngularMomentum()` gives the impeller inlet preswirl.
- `MeridionalDesign.py`: `MeridionalGeometry` (contour fullness, outlet angles, hub taper ratio, `Z_e_factor`) and `MeridionalDesign`. `MeridionalDesign.__init__` builds:
  - `Z_e`;
  - `ShroudContour`, `HubContour` and `MeanContour`;
  - `PotentialFlow`;
  - the leading edge (`LeadingEdge`, a `BezierCurve` from hub to shroud);
  - the trailing edge (`TrailingEdge`, a straight `Polyline` from the hub outlet to the shroud outlet).
  
  Plotting: `MeridionalTraces`, `PlotMeridional`.
- `ImpellerBladeDesign.py`: `ImpellerBladeDesign(meridional_design, inducer_radial_equilibrium)`.
  - Finds the inlet flow where the leading edge crosses streamline φ: r₁, c_m1, and w_u1, using the inducer preswirl c_u = K/r.
  - Adds leading-edge blade blockage and solves β₁′ explicitly. `_INCIDENCE_PADDING` adds 2° to the blade angle.
  - Plotting: `PlotInletTriangles`.
  - The plan view (blade wrap) is not implemented yet.
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
    - `MeridionalContour(SegmentChain)`: a flat `Polyline`, then the bend `Bend` (a `BezierCurve`), then an outlet `Polyline`;
    - `MeanCurve(a, b)` returns a `Polyline`.
  - `Cross(a, b)` is the 2D cross product, used by `FindIntersect` and the mesh.
  - `FlowDomain` extends the hub and shroud upstream and downstream for the FEM solve, with zones `inlet`/`impeller`/`outlet`. Its walls are `SegmentChain`s with the contour nested as one segment. `FlowDomain.Mesh()` builds a `MeridionalMesh`, by default 101 stations × 11 lines.
- `PotentialFlow.py`: the FEM stream-function solve (`PotentialFlow`). `GetStreamline(phi)` returns a `Streamline(Polyline)` with:
  - `.x`: arclength from the impeller inlet;
  - `.v_m`: meridional velocity;
  - `.Impeller`: the impeller part only.
- `VelocityTriangle.py`: `VelocityTriangle(u, c_m, c_u)`, `VelocityTriangle.FromBeta`, `PlotTriangles`, `CONVENTION_NOTES`. An optional `incidence` (blade angle = β + i) makes `PlotTriangles` draw the blade direction and the arc i, as in Gülich Fig. 3.1; `PlotInletTriangles` uses it on the 1′ triangle.
- `Bezier.py`: `bez(P, t)` and `dbez(P, t)` for cubic Beziers, by Bezier parameter t. Use `BezierCurve` when you need a curve evaluated by arclength.
- `Volute.py`, `VolutePlot.py`, `Sizing/`: see Running.

`Inputs/`: `Inputs.yaml`, and `Serialize.py`, which turns every `[value, 'unit']` pair into a `Q_`. `Inputs/__init__.py` exports the result as `Inputs`.

## Conventions

- File split: geometry goes in `MeridionalContour.py`, flow theory goes in `PotentialFlow.py`. `FlowDomain` and `MeridionalMesh` count as geometry because they contain no flow equations. The one deliberate exception is `QuasiOrthogonalFlow` (1D continuity), which stays with the mesh because `FlowDomain.Mesh(by_flow=True)` uses it; moving it would make the two files import each other.
- Try to reuse code
- Use pint `Q_` at interfaces. Curve and mesh internals are plain floats in meters; convert with `.to(...)` and `.m`.
- Meridional coordinates are always ordered (z, r): z is axial and is the plot x-axis, r is radial.
- Normalized arclength s runs from 0 at the impeller inlet to 1 at the outlet on the contours and the meanline. Streamlines span the whole FEM domain; use `Streamline.Impeller` for the impeller part.
- `phi` in the code is the flow fraction (stream function / its shroud value: 0 at the hub, 1 at the shroud). It is not Milos's velocity potential φ.
- Velocities:
  - `c` is the absolute-frame velocity in the sizing and velocity-triangle code (`c_m`, `c_u`);
  - `PotentialFlow` deliberately keeps `v_z`, `v_r`, `v_m`.
- Blade and flow angles β are measured from the tangential direction, as in `VelocityTriangle`: tan β = c_m / (u − c_u). Some references measure them from the meridional direction instead, such as NASA TN D-1562.
- Style:
  - match the file's indentation: tabs in most of `Code/`, spaces in `Bezier.py`, `EmpiricalRelations.py`, `Inducer.py`, `InducerRadialEquilibrium.py`, `Volute.py`, `VolutePlot.py` and `Sizing/`;
  - CamelCase method names (`GetX`, `PlotX`) and object attributes (`HubContour`);
  - cite the source equation or page in docstrings.
- Plots use plotly, mostly with `template='plotly_white'`; `PlotInletTriangles` uses `'simple_white'`. A length-unit argument `l` defaults to `'in'`. Meridional views use equal axis scaling.

## Design method and references

The main reference is M. Milos, "Impeller Design Using CAD Techniques and Conformal Mapping Method", InTech, *Centrifugal Pumps* (2012): https://cdn.intechopen.com/pdfs/29651/InTech-Impeller_design_using_cad_techniques_and_conformal_mapping_method.pdf. Page numbers below are the printed ones.

To read the PDF: WebFetch can't parse it and no PDF library is installed. Download it and inflate its text streams with zlib in a small script.

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
  - the geometric meanline lies about 0.5-0.8 mm closer to the hub than the true φ = 0.5 streamline at the inlet. Using `self.PotentialFlow.GetStreamline(0.5).Impeller` instead is a one-line swap.
- **Leading-edge blockage** (`ImpellerBladeDesign._SolveInletFlow`): CFturbo; Hassan et al., Sci. Rep. 2024, eq. 17; Milos eq. 33.
- **Limits of the potential flow:** inviscid, irrotational, axisymmetric, no blades. It is fine for contour shaping and streamlines, but not for head, losses or off-design. The ~2 mm outlet width b₂ is dominated by viscosity.
- **NASA TN D-1562** (Stockman & Kramer, 1963), "Method for Design of Pump Impellers Using a High-Speed Digital Computer". It is an inverse method: from the hub, the hub relative velocity and the blade surface θ(r, z), it marches streamline by streamline to the *shroud*. Its β is measured from the meridional direction. We don't use the marching itself, because our shroud comes from Milos. The parts we use are:
  - eq. 9: tan β = r dθ/dm along a streamline, which is the wrap-angle (plan-view) relation;
  - eq. 10: circumferential blade thickness t_θ;
  - App. C, eqs. C3-C7: approximate pressure- and suction-side velocities. C3 reduces to g Δh = (2πr − Z t_θ)/Z · w_m · (1/r) d(r c_u)/dm.
- **Other references:** Gülich, *Centrifugal Pumps* (axial extension Z_e p. 348, NPSH ch. 6); Wiesner slip (CFturbo manual); NASA SP-8052 (inducer radial equilibrium).
