"""Load the original standalone correlations without its CEA/pyfluids drivers."""
import ast
from pathlib import Path
from types import SimpleNamespace

import numpy as np


def load_correlations():
    path = Path(__file__).resolve().parent.parent / 'Development' / 'CFD BCs' / 'Constant_BCs.py'
    source = ast.parse(path.read_text(encoding='utf-8-sig'), filename=str(path))
    names = {'friction_factor', 'gnielinski', 'dittus_boelter'}
    functions = [node for node in source.body
                 if isinstance(node, ast.FunctionDef) and node.name in names]
    if {node.name for node in functions} != names:
        raise ValueError(f'Missing required correlation functions in {path}')
    # Execute the original function definitions; do not run the unrelated
    # rocketcea/pyfluids imports or hardcoded torch/chamber calculations.
    namespace = {'np': np}
    exec(compile(ast.Module(body=functions, type_ignores=[]), str(path), 'exec'), namespace)
    return SimpleNamespace(**{name: namespace[name] for name in names})
