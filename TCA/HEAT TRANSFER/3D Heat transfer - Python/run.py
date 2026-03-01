#!/usr/bin/env python3
"""
run.py  —  Full pipeline entry point for TCA heat transfer analysis
--------------------------------------------------------------------
Reads TCA_params.yaml (engine design parameters) and heat_transfer.yaml
(solver settings), runs the Bartz gas-property calculator, then runs
the 2-D transient heat transfer solver.

Usage (from repo root):
    python "TCA/HEAT TRANSFER/3D Heat transfer - Python/run.py"

Or from this directory:
    python run.py

Input files read:
    TCA/TCA_params.yaml                      — engine design parameters
    heat_transfer.yaml  (same dir as run.py) — solver / geometry settings

Intermediate file written:
    turbopump_properties.csv  (same dir as run.py)

Output files written (same dir as run.py):
    nozzle_heat_transfer.mp4  (video, if record_video: true)
    live matplotlib figures
"""

import os
import sys
import yaml

# ---------------------------------------------------------------------------
# Resolve paths
# ---------------------------------------------------------------------------
SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))

# TCA_params.yaml: two levels up from this script -> .../TCA/TCA_params.yaml
TCA_YAML = os.path.normpath(os.path.join(SCRIPT_DIR, '..', '..', 'TCA_params.yaml'))
if not os.path.exists(TCA_YAML):
    # Fallback: search from cwd (repo root -> TCA/TCA_params.yaml)
    TCA_YAML = os.path.join(os.getcwd(), 'TCA', 'TCA_params.yaml')

HT_YAML    = os.path.join(SCRIPT_DIR, 'heat_transfer.yaml')
PROPS_CSV  = os.path.join(SCRIPT_DIR, 'turbopump_properties.csv')
BARTZ_DIR  = os.path.normpath(os.path.join(SCRIPT_DIR, '..', '..', 'TCA Sizing Code Files'))

# ---------------------------------------------------------------------------
# Load YAML files
# ---------------------------------------------------------------------------
print('=' * 60)
print('Loading TCA_params.yaml ...')
if not os.path.exists(TCA_YAML):
    raise FileNotFoundError(f'TCA_params.yaml not found at:\n  {TCA_YAML}')
with open(TCA_YAML) as f:
    tca_params = yaml.safe_load(f) or {}

print('Loading heat_transfer.yaml ...')
if not os.path.exists(HT_YAML):
    raise FileNotFoundError(f'heat_transfer.yaml not found at:\n  {HT_YAML}')
with open(HT_YAML) as f:
    ht_params = yaml.safe_load(f) or {}

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
import HeatTranferSolver_2D

HeatTranferSolver_2D.run(ht_params, tca_params, PROPS_CSV, SCRIPT_DIR)

print('\nPipeline complete.')
