"""CEA gas-property calculations for a rocket thrust-chamber contour.

The GUI can import this module without running CEA. Geometry may have an
unmeasured throat at x=0. Wall temperature is supplied to Bartz separately.
"""
from pathlib import Path
import numpy as np
import pandas as pd
import yaml

PROJECT_ROOT = Path(__file__).resolve().parent.parent


# Read x/y from a CSV and normalize their units to millimeters.
def read_geometry_csv(path):
    df = pd.read_csv(path)

    def axis(name):
        # Supported column spellings and conversion factors to millimeters.
        options = [
            (f'{name}_mm', 1.0), (f'{name}_cm', 10.0),
            (f'{name}_m', 1000.0), (f'{name} [mm]', 1.0),
            (f'{name} [cm]', 10.0), (f'{name} [m]', 1000.0),
        ]
        found = [(col, scale) for col, scale in options if col in df.columns]
        if len(found) != 1:
            raise ValueError(
                f'Expected exactly one unit-labeled {name} column, '
                f'for example {name}_mm or {name}_m.'
            )
        col, scale = found[0]
        return pd.to_numeric(df[col], errors='raise').to_numpy(dtype=float) * scale

    geometry = pd.DataFrame({'x_mm': axis('x'), 'y_mm': axis('y')})
    return validate_geometry(geometry)


# Check the contour is finite, spans x=0, and converges/diverges as expected.
def validate_geometry(geometry):
    if not {'x_mm', 'y_mm'}.issubset(geometry.columns):
        raise ValueError('Geometry requires x_mm and y_mm columns.')

    df = geometry[['x_mm', 'y_mm']].copy()
    for name in ('x_mm', 'y_mm'):
        df[name] = pd.to_numeric(df[name], errors='raise')
    df = df.sort_values('x_mm').reset_index(drop=True)

    if not np.isfinite(df.to_numpy()).all() or (df.y_mm <= 0).any():
        raise ValueError('Geometry must have finite coordinates and positive radii.')
    # Adjacent contour segments can repeat their shared endpoint. Merge only
    # matching radii, allowing negligible floating-point rounding in millimeters.
    duplicates = df.loc[df.x_mm.duplicated(keep=False)]
    for xpos, group in duplicates.groupby('x_mm'):
        radii = group.y_mm.to_numpy(dtype=float)
        if not np.allclose(radii, radii[0], rtol=1e-9, atol=1e-9):
            raise ValueError(
                f'Conflicting radii at x={xpos:.12g} mm; '
                'each axial station must have a single radius.'
            )
    df = df.drop_duplicates(subset='x_mm', keep='first').reset_index(drop=True)
    if len(df) < 4:
        raise ValueError('Provide at least four distinct axial stations.')
    if not (df.x_mm.min() < 0 < df.x_mm.max()):
        raise ValueError('Geometry must include points on both sides of x=0.')

    # Tolerance in millimeters; permits small measurement/CAD sampling noise.
    tolerance_mm = 0.01
    before = df.loc[df.x_mm < 0, 'y_mm'].to_numpy(dtype=float)
    after = df.loc[df.x_mm > 0, 'y_mm'].to_numpy(dtype=float)
    if np.any(np.diff(before) > tolerance_mm):
        raise ValueError('Radius must generally decrease or stay constant before x=0.')
    if np.any(np.diff(after) < -tolerance_mm):
        raise ValueError('Radius must generally increase or stay constant after x=0.')
    return df


# Convert either a YAML species string or a custom YAML record to a CEA reactant.
def make_reactant(spec, cea):
    if isinstance(spec, str):
        return spec
    if not isinstance(spec, dict):
        raise TypeError('Propellant must be a species string or a custom reactant mapping.')
    required = ('name', 'formula', 'molecular_weight', 'enthalpy', 'enthalpy_units', 'temperature')
    missing = [key for key in required if key not in spec]
    if missing:
        raise ValueError(f'Custom propellant missing YAML fields: {missing}')
    return cea.Reactant(
        name=str(spec['name']),
        formula={str(k): float(v) for k, v in spec['formula'].items()},
        molecular_weight=float(spec['molecular_weight']),
        enthalpy=float(spec['enthalpy']),
        enthalpy_units=str(spec['enthalpy_units']),
        temperature=float(spec['temperature']),
    )


# Run CEA on each geometry station and save the resulting gas properties.
def generate_properties(geometry, params_path=PROJECT_ROOT / 'Inputs' / 'TCA_params.yaml',
                        output_path=None, log=print):
    """Return gas-property DataFrame; also write output_path.
    Chamber pressure is interpreted as psi (absolute), CEA Cp as kJ/kg/K,
    viscosity as millipoise, and conductivity as W/cm/K (as in the original).
    These source-unit assumptions should be verified for your CEA installation.
    """
    import cea  # Delayed so existing-properties mode needs no CEA installation.

    g = validate_geometry(geometry)
    with open(params_path, encoding='utf-8') as f:
        p = yaml.safe_load(f)

    x = g.x_mm.to_numpy(dtype=float) / 1000.0
    y = g.y_mm.to_numpy(dtype=float) / 1000.0

    # Use minimum *measured* radius as a conservative estimate of r_t. Since
    # x=0 is not sampled, this is not necessarily the true physical throat.
    # It guarantees that no measured station has A/At below unity.
    rt = float(np.min(y))
    area_ratios = (y / rt) ** 2
    log(f'Using minimum measured throat radius {rt * 1000:.6f} mm; '
        'the actual throat radius at x=0 may differ.')

    fuel = make_reactant(p['propellants']['fuel'], cea)
    oxidizer = make_reactant(p['propellants']['oxidizer'], cea)
    species = [fuel, oxidizer]
    reac = cea.Mixture(species)
    prod = cea.Mixture(species, products_from_reactants=True)
    solver = cea.RocketSolver(prod, reactants=reac, transport=True)
    solution = cea.RocketSolution(solver)

    weights = reac.of_ratio_to_weights(
        np.array([0.0, 1.0]), np.array([1.0, 0.0]), float(p['of_ratio'])
    )
    temperatures = p['reactant_T']
    T_reactant = np.array([
        float(temperatures['fuel']), float(temperatures['oxidizer'])
    ])
    hc = reac.calc_property(cea.ENTHALPY, weights, T_reactant) / cea.R

    pc_bar = float(p['chamber_pressure']) * 0.0689475729
    ac_at = float(p['contraction_ratio'])
    if pc_bar <= 0 or ac_at <= 0:
        raise ValueError('Chamber pressure and contraction ratio must be positive.')

    rows = []
    throat_index = int(np.argmin(y))
    x_throat = float(x[throat_index])
    for i, xpos in enumerate(x):
        # The throat reference at x=0 isn't a station in the nozzle branch solve.
        if i == throat_index or np.isclose(xpos, 0.0, rtol=0, atol=1e-12):
            continue
        ar = float(area_ratios[i])
        if ar < 1.0 - 1e-12:
            raise ValueError(f'Invalid A/At={ar:.8g} at x={xpos:.8g} m.')
        kwargs = {'subar': [ar]} if xpos < 0 else {'supar': [ar]}
        try:
            solver.solve(solution, weights, pc_bar, ac_at=ac_at,
                         iac=False, hc=hc, **kwargs)
        except Exception as exc:
            raise RuntimeError(
                f'CEA failed at x={xpos:.9g} m, A/At={ar:.9g}. '
                'Check throat-radius estimate and contour geometry.'
            ) from exc
        rows.append({
            'x_m': x[i], 'y_m': y[i],
            'gamma': float(solution.gamma_s[-1]),
            'Cp [(KJ/kg-K)]': float(solution.cp[-1]),
            'k [W/m-K]': float(solution.conductivity_eq[-1]) * 100.0,
            'MW [kg/kmol]': float(solution.MW[-1]),
            'R [J/kg-K]': float(cea.R / solution.MW[-1]),
            'T_local [K]': float(solution.T[-1]),
            'P_local [bar]': float(solution.P[-1]),
            'ae_at': float(solution.ae_at[-1]),
            'Viscosity [Pa*s]': float(solution.viscosity[-1]) * 1e-4,
            'Prandtl Number': float(solution.Pr_eq[-1]),
            'Mach': float(solution.Mach[-1]),
        })
        if len(rows) % 10 == 0:
            log(f'CEA solved {len(rows)} stations')

    props = pd.DataFrame(rows)
    if props.empty:
        raise ValueError('No CEA stations were solved.')
    # No automatic file writes. GUI export handles timestamped outputs.
    if output_path is not None:
        raise ValueError('Direct output_path writes are disabled. Export through the GUI.')
    log('CEA properties calculated in memory; nothing saved.')
    return props


if __name__ == '__main__':
    result = generate_properties(read_geometry_csv(PROJECT_ROOT / 'Outputs' / 'contour.csv'))
    print(result.head())
    print('No files saved. Use thermal_gui.py to export.')
