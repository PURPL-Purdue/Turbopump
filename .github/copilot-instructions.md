# Copilot / AI agent instructions — TurboPump

Purpose: Help AI coding agents quickly become productive in this mechanical-engineering codebase by describing architecture, workflows, conventions, and safe edit boundaries.

1) Big picture
- Repo is an engineering workspace combining analysis scripts, notebooks, CAD archives, and parameter YAMLs. Primary code lives under `TCA/` and `Systems design software/`; many calculations reference top-level YAMLs such as `params.yaml`, `pumps_params.yaml`, and `rotordynamics_params.yaml`.
- The repo is NOT a single packaged application — changes should be focused to scripts, notebooks, and small utilities rather than binary CAD assets in `Archive/`, `Pump Assembly/`, `Fall 2025/`.

2) Key locations (examples)
- Parameters: [params.yaml](params.yaml) and [pumps_params.yaml](pumps_params.yaml)
- Thrust chamber / injector analysis: [TCA/injector_calculator.py](TCA/injector_calculator.py) and [TCA/injecotr_interactive.ipynb](TCA/injecotr_interactive.ipynb)
- Sizing code: [TCA/TCA Sizing Code Files/Final_Contouring_Script.py](TCA/TCA%20Sizing%20Code%20Files/Final_Contouring_Script.py)
- Systems utilities: [Systems design software/InducerSuctionPerformance.py](Systems%20design%20software/InducerSuctionPerformance.py)

3) Conventions and patterns discovered
- Parameter-first: Many scripts read configuration YAML/CSV files instead of hard-coded constants. Prefer adding new parameters to the YAMLs and wiring them into scripts.
- Notebook-first exploration: Jupyter notebooks are used for interactive analysis; convert notebook logic into small, importable `.py` modules when you need repeatable CLI runs.
- Unit awareness: Engineering code assumes physical units; preserve units and clarify conversions when changing numeric routines.

4) Developer workflows (how to run common tasks)
- Run a Python analysis script (use the workspace Python; examples used Python 3.13):
  - `python3 TCA/TCA\ Sizing\ Code\ Files/1.py`
  - `python3 TCA/injector_calculator.py`
- Open or run notebooks interactively with Jupyter Lab / Notebook: `jupyter lab` in the repo root.
- MATLAB files and CAD binaries are present but out-of-scope for code edits; do not attempt to programmatically modify `.prt`, `.SLDASM`, or MATLAB `.m` files unless asked.

5) Integration points & external dependencies
- Python scripts rely on common scientific packages (numpy, matplotlib, scipy). Check local environment and prefer `python3` from a modern interpreter (3.11+ / 3.13 is used in terminal history).
- Data inputs often come from CSV or YAML at repo root (e.g., `injector_parameters.csv`, `params.yaml`). Modifying these will affect many scripts — add new keys carefully.

6) Safe edit rules for AI agents
- Do: Edit `.py` and `.ipynb` files that implement analysis, add small helper modules, and update top-level YAML parameters.
- Do not: Modify or compress large binary CAD/assembly files inside `Archive/`, `Pump Assembly/`, or `Fall 2025/`.
- Back up any changed parameter YAMLs and notify maintainers when altering canonical inputs.

7) Suggested PR contents when code changes are made
- Short description of intent and which scripts are affected.
- A small runnable example or brief instructions to reproduce the change (one or two commands).
- If a parameter is added, include the YAML diff and an explanation of units/expected range.

8) When in doubt
- Run the related script or notebook locally to confirm results before making large changes.
- Ask for clarification when edits impact CAD/archival directories or broad parameter changes.

If this covers the missing context, I can iterate with more specific examples (e.g., show how to add a parameter and wire it through `TCA/injector_calculator.py`). What would you like changed or expanded?
