#!/usr/bin/env python3
"""
run.py  —  Full pipeline entry point for TCA heat transfer analysis
--------------------------------------------------------------------
Reads TCA_params.yaml (engine design parameters) and heat_transfer.yaml
(solver settings), runs the Bartz gas-property calculator, then runs
the 2-D transient heat transfer solver.

Usage (from anywhere):
    python /path/to/run.py
Or from this directory:
    python run.py

All paths resolve relative to this script's directory unless they are
absolute.  The keys `tca_params_path` and `bartz_module_dir` in
heat_transfer.yaml control where TCA_params.yaml and the Bartz_Values
module are located; both default to legacy sibling-folder paths but you
can override them to anything (including absolute paths).

After the run, results are saved to results.npz.  Inspect with:
    python view_results.py results.npz
"""

import os
import sys
import yaml


def _resolve(path, base):
    """Return abs path: keep absolute paths, resolve relative ones against base."""
    if not path:
        return None
    return path if os.path.isabs(path) else os.path.normpath(os.path.join(base, path))


# ---------------------------------------------------------------------------
# Resolve paths
# ---------------------------------------------------------------------------
SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
HT_YAML    = os.path.join(SCRIPT_DIR, 'heat_transfer.yaml')
PROPS_CSV  = os.path.join(SCRIPT_DIR, 'turbopump_properties.csv')

print('=' * 60)
print('Loading heat_transfer.yaml ...')
if not os.path.exists(HT_YAML):
    raise FileNotFoundError(f'heat_transfer.yaml not found at:\n  {HT_YAML}')
with open(HT_YAML) as f:
    ht_params = yaml.safe_load(f) or {}

# TCA params — path comes from heat_transfer.yaml, default fallback if absent
TCA_REL    = ht_params.get('tca_params_path', '../../TCA_params.yaml')
TCA_YAML   = _resolve(TCA_REL, SCRIPT_DIR)
if not os.path.exists(TCA_YAML):
    _alt = os.path.join(os.getcwd(), 'TCA_params.yaml')
    if os.path.exists(_alt):
        TCA_YAML = _alt
if not os.path.exists(TCA_YAML):
    raise FileNotFoundError(
        f"TCA_params.yaml not found.\n"
        f"  Tried: {TCA_YAML}\n"
        f"  Fix: set 'tca_params_path' in heat_transfer.yaml to point to it "
        f"(relative or absolute).")

print(f'Loading {os.path.basename(TCA_YAML)} ...')
with open(TCA_YAML) as f:
    tca_params = yaml.safe_load(f) or {}

# Bartz module location
BARTZ_REL = ht_params.get('bartz_module_dir', '../../TCA Sizing Code Files')
BARTZ_DIR = _resolve(BARTZ_REL, SCRIPT_DIR)
if not BARTZ_DIR or not os.path.isdir(BARTZ_DIR):
    # Fall back to the local copy that ships with this repo
    if os.path.exists(os.path.join(SCRIPT_DIR, 'Bartz_Values.py')):
        BARTZ_DIR = SCRIPT_DIR
        print(f'Bartz: using local copy in {SCRIPT_DIR}')
    else:
        raise FileNotFoundError(
            f"Bartz_Values module not found.\n"
            f"  Tried: {BARTZ_DIR}\n"
            f"  Fix: set 'bartz_module_dir' in heat_transfer.yaml.")

# ---------------------------------------------------------------------------
# Validate critical inputs early
# ---------------------------------------------------------------------------
_problems = []
for _k in ('oxidizer_fuel_ratio', 'tca_chamber_pressure', 'turbopump_mdot'):
    if _k not in tca_params:
        _problems.append(f'TCA_params.yaml missing required key: {_k}')
if tca_params.get('tca_chamber_pressure', 0) <= 0:
    _problems.append(f'tca_chamber_pressure must be > 0 (got {tca_params.get("tca_chamber_pressure")})')
if tca_params.get('turbopump_mdot', 0) <= 0:
    _problems.append(f'turbopump_mdot must be > 0 (got {tca_params.get("turbopump_mdot")})')
if tca_params.get('oxidizer_fuel_ratio', 0) <= 0:
    _problems.append(f'oxidizer_fuel_ratio must be > 0 (got {tca_params.get("oxidizer_fuel_ratio")})')
if _problems:
    raise ValueError('Invalid TCA_params.yaml settings:\n  - ' + '\n  - '.join(_problems))

print(f"  O/F ratio         : {tca_params['oxidizer_fuel_ratio']}")
print(f"  Chamber pressure  : {tca_params['tca_chamber_pressure']} psia")
print(f"  Mass flow         : {tca_params['turbopump_mdot']} lbm/s")
print('=' * 60)

# ---------------------------------------------------------------------------
# Step 1: Bartz gas-property calculator
# ---------------------------------------------------------------------------
print('\n[Step 1/2]  Running Bartz gas-property calculator...')
sys.path.insert(0, BARTZ_DIR)
import Bartz_Values

Bartz_Values.run(tca_params, SCRIPT_DIR)

# ---------------------------------------------------------------------------
# Step 2: 2-D heat transfer solver
# ---------------------------------------------------------------------------
print('\n[Step 2/2]  Running 2-D heat transfer solver...')
sys.path.insert(0, SCRIPT_DIR)
import HeatTransferSolver_2D

HeatTransferSolver_2D.run(ht_params, tca_params, PROPS_CSV, SCRIPT_DIR)

print('\nPipeline complete.')
