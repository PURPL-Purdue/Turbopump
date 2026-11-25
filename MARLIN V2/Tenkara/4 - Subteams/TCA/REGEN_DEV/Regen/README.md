# Regen

Regenerative cooling analysis for the IPA / LOX engine, following Heister et al.,
*Rocket Propulsion*, section 6.5.2. Where every input and equation comes from:
[TRACEABILITY.md](TRACEABILITY.md).

## Layout

Physics lives in packages that never read the config. The scripts read
`regen.yaml` and pass the values in. Printing and files are in `util/`.

| Path | What's in it |
|---|---|
| `regen.yaml` | Every input: engine, chamber, contour, jacket, models, both solvers' channel settings |
| **Run scripts** | |
| `scripts/channel_solver.py` | Design: give w and a target T_wg, get the channel height h(x). Writes `output/channel_geometry.csv`, and the DXF / CSV for CAD |
| `scripts/wall_temp_solver.py` | Analysis: give the channels (constant, or a geometry file with h(x)), get wall / coolant temperatures and pressure |
| `scripts/engine_sizing.py` | CEA sizing sanity check |
| `studies/studies.py` | Parameter studies (currently CEA chamber temperature vs O/F) |
| **Physics** | |
| `engine/engine_inputs.py` | Engine definition from `regen.yaml`, converted to SI |
| `engine/cea_interface.py` | Engine sizing (`size_engine`) and hot-gas state along the contour (`gas_side_states`) from CEA |
| `engine/contour.py` | Wall profile, injector face -> exit, with an 80% Rao bell |
| `engine/rao_bell.py` | Team's bell contour script (heatsink engine repo), vendored unmodified; only its Rao angle lookup is used |
| `heat_transfer/gas_side.py` | Recovery temperature, Bartz h_g, hot-gas conditions along the contour |
| `heat_transfer/coolant_side.py` | Coolant-side h_l correlations (Sieder-Tate), picked by name in `regen.yaml` |
| `heat_transfer/friction.py` | Channel friction factor (smooth-wall fits, Colebrook when rough) |
| `coolant/` | Coolant property modules (`ipa_props.py`), picked by name in `regen.yaml` |
| `jacket/jacket_inputs.py` | Jacket section of `regen.yaml` -> the `Jacket` the physics takes |
| `jacket/geometry.py` | Channel dimensions, channel-count rule, reading a channel geometry file |
| `jacket/march.py` | The 6.5.2 march: station balance, coolant update, forward march |
| `jacket/sizing.py` | Inverse: channel height for a target T_wg |
| **Support** | |
| `util/config.py` | Loads `regen.yaml` (or the file named in `REGEN_CONFIG`) |
| `util/report.py` | Printed summary, checks, station table, plots |
| `util/export.py` | `contour.dxf`, `channel_heights.csv` + NX points, `channel_geometry.csv` |
| `util/provenance.py` | The stamp on every output: DXF header, CSV comment lines, plot footer, `inputs_used.yaml` |
| `output/` | Generated files (gitignored). `inputs_used.yaml` records what made them: time, git commit, yaml hash, and an exact copy of the yaml |
| `../Regen Cad/` | NX parts (outside this repo) |
| `reference/` | Textbook, papers, RPA output, CEA example. Local only, not in git (copyrighted) |

Optimization studies, the Blong 2024 paper models, the old sizing scripts, the throat
solver and the ethanol inputs are on branch `archive/full-2026-10-01`.

## Setup

Python 3.13. From this folder:

```
python -m pip install -r requirements.txt
```

Versions are pinned exactly. CEA and the property packages ship their own data,
so a different version can change the numbers.

## Running

```
python scripts/channel_solver.py      # size h(x); writes output/channel_geometry.csv
python scripts/wall_temp_solver.py    # analyse a channel (set geometry_file in regen.yaml)
python scripts/engine_sizing.py       # sizing sanity check
python studies/studies.py
```

To analyse the channel the solver just sized, or an edited copy of it: set
`wall_temp_solver.geometry_file: output/channel_geometry.csv` (or your copy).
The file must match the current contour; if the engine inputs changed, re-run
the channel solver first.

## Adding a model

- Coolant-side correlation: a function in `heat_transfer/coolant_side.py`, an entry
  in `CORRELATIONS`, a TRACEABILITY row, then `jacket.coolant_htc` in `regen.yaml`.
- Coolant: a module with the `ipa_props` interface in `coolant/`, an entry in
  `FLUIDS`, then `jacket.coolant`.
