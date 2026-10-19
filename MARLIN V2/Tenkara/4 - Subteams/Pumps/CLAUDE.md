# Pumps (MARLIN V2 turbopump)

Preliminary design code for the LOX and IPA centrifugal pump impellers, plus standalone inducer, volute and NPSH scripts. Design-point inputs are in `Inputs/Inputs.yaml`: shaft speed, and mass flow and feed/discharge pressure for each propellant.

## Running

- Use the local venv `env\Scripts\python.exe`. It is not in git, and it has cadquery, CoolProp, PyYAML and scipy; the system Python lacks them. `requirements.txt` doesn't list scipy or cadquery.
- Entry points are `Code/Impeller/SizingIPA.py` and `SizingLOX.py`, run directly; they put the Pumps folder on `sys.path`. They print the sizing and open plotly figures in the browser. `ShroudHubRevolve()` opens the cadquery viewer.
- `import Code.Impeller` builds both impellers at import time (`_IPA_impeller.py`, `_LOX_impeller.py`), including the FEM potential-flow solve, so it takes a few seconds. Any error in `Impeller.__init__` or `BuildMeridian` breaks every import of the package.
- There is no test suite. To check changes without windows popping up:
  - set `plotly.graph_objects.Figure.show` to a no-op in a scratch script and use `fig.write_html(...)` instead;
  - there is no kaleido, so figures can't be exported as images.
- `Code/Inducer*.py`, `NPSHrScaling.py` and `Volute.py` are separate matplotlib scripts with bare imports, not part of the `Code.Impeller` package.

## Code/Impeller layout

- `Impeller.py`
  - `InputGeometry`, `DesignPoint` and `Impeller`.
  - `Impeller` does the 1D sizing: Wiesner slip, velocity triangles, H-Q curve, NPSH_i.
  - `BuildMeridian` builds `ShroudContour`, `HubContour`, `MeanContour`, `PotentialFlow`, the leading edge (`LeadingEdge`, a `BezierCurve` from hub to shroud), and the trailing edge (`TrailingEdge`, a straight `Polyline` from the shroud outlet to the hub outlet).
  - Plotting: `PlotMeridional`, `PlotPerformanceHQ`, `ShroudHubRevolve`.
- `_IPA_impeller.py`, `_LOX_impeller.py`: geometry inputs for each pump, with the design point taken from `Inputs.yaml`.
- `MeridionalContour.py`: meridional curves, the FEM domain and the mesh.
  - `MeridionalCurve` is the common base:
    - `curve(s)` returns (z, r) Quantities at normalized arclength;
    - `Points(s_m)` returns an (n, 2) array in m at true arclength;
    - `Length`;
    - `StaticMoment(s_1, s_2)` returns ∫r dx;
    - `FindIntersect(other)` returns an `IntersectResults`: a list of `(z, r, s_1, s_2)` tuples ordered by `s_1` (normalized arclength on this curve, `s_2` on `other`), empty when the curves don't cross. A start or end of one curve lying on the other (within 1 µm) counts, with `s` exactly 0 or 1. Overlapping stretches and a curve with itself are not handled;
    - `Closest(p)` returns the distance from a point (Quantity or m) to the curve and the normalized arclength of the closest point;
    - `Trace(l, **scatter)` returns the curve as a plotly line.
  - Subclasses (a new kind of curve only needs `length_m` and `Points`):
    - `Polyline`: straight lines through points; two points make a straight line. `Interpolate(values, s)` interpolates per-point data, such as a streamline's `v_m`, at normalized arclength the same way as the points;
    - `BezierCurve`: cubic Bezier by arclength, with `ControlPoints`. `BezierCurve.ThroughPoints(P_a, P_m, P_b)` is a quadratic through three points;
    - `SegmentChain`: any curves joined end to end, chains included. `Segments` holds the pieces and `s_m` the arclength (m) of each joint;
    - `MeridionalContour(SegmentChain)`: a flat `Polyline`, then the bend `Bend` (a `BezierCurve`), then an outlet `Polyline`;
    - `MeanCurve(a, b)` returns a `Polyline`.
  - `Cross(a, b)` is the 2D cross product, used by `FindIntersect` and the mesh.
  - `FlowDomain` extends the hub and shroud upstream and downstream for the FEM solve, with zones `inlet`/`impeller`/`outlet`. Its walls are `SegmentChain`s with the contour nested as one segment. `FlowDomain.Mesh()` builds a `MeridionalMesh`, by default 101 stations × 11 lines.
- `PotentialFlow.py`: the FEM stream-function solve (`PotentialFlow`). `GetStreamline(phi)` returns a `Streamline(Polyline)` with:
  - `.x`: arclength from the impeller inlet;
  - `.v_m`: meridional velocity;
  - `.Impeller`: the impeller part only.
- `Bezier.py`: `bez(P, t)` and `dbez(P, t)` for cubic Beziers, by Bezier parameter t. Use `BezierCurve` when you need a curve evaluated by arclength.
- `VelocityTriangle.py`, `EmpiricalRelations.py`.

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
- Style:
  - match the file's indentation: tabs in most of `Code/Impeller`, spaces in `Bezier.py` and the inducer scripts;
  - CamelCase method names (`GetX`, `PlotX`) and object attributes (`HubContour`);
  - cite the source equation or page in docstrings.
- Plots use plotly with `template='plotly_white'`. A length-unit argument `l` defaults to `'in'`. Meridional views use equal axis scaling.

## Design method and references

The main reference is M. Milos, "Impeller Design Using CAD Techniques and Conformal Mapping Method", InTech, *Centrifugal Pumps* (2012): https://cdn.intechopen.com/pdfs/29651/InTech-Impeller_design_using_cad_techniques_and_conformal_mapping_method.pdf. Page numbers below are the printed ones.

To read the PDF: WebFetch can't parse it and no PDF library is installed. Download it and inflate its text streams with zlib in a small script.

- **Domain and mesh**, pp. 35-36: 101 × 11 nodes with straight quasi-orthogonals. `BuildMeridian` currently calls `Mesh(by_flow=False)`, which divides each QO equally, as Milos does.
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
- **Limits of the potential flow:** inviscid, irrotational, axisymmetric, no blades. It is fine for contour shaping and streamlines, but not for head, losses or off-design. The ~2 mm outlet width b₂ is dominated by viscosity.
- **Other references:** Gülich, *Centrifugal Pumps* (axial extension Z_e p. 348, NPSH ch. 6); Wiesner slip (CFturbo manual).