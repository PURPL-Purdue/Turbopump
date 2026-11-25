"""
Loads regen.yaml (repo root). Every input lives there; this file only reads it.
Unit conversions to SI happen in engine_inputs.py / jacket_inputs.py and the solvers.

Use a different file:  set the REGEN_CONFIG environment variable to its path.
"""
import hashlib
import os
import pathlib
import yaml

ROOT = pathlib.Path(__file__).resolve().parents[1]
PATH = pathlib.Path(os.environ.get("REGEN_CONFIG", ROOT / "regen.yaml"))

# Keep the exact bytes that were loaded, so outputs can record what made them
# (util/provenance.py). The hash is of the raw file, so it matches
# `sha256sum regen.yaml` whatever the line endings are.
_RAW = PATH.read_bytes()
SHA256 = hashlib.sha256(_RAW).hexdigest()
TEXT = _RAW.decode("utf-8")
CFG = yaml.safe_load(TEXT)

engine = CFG["engine"]
chamber = CFG["chamber"]
contour = CFG["contour"]
jacket = CFG["jacket"]
wall_temp_solver = CFG["wall_temp_solver"]
channel_solver = CFG["channel_solver"]
