"""Bartz gas-side thermal calculation callable from thermal_gui.py.

Preserves the original script's array-level h_g and recovery-temperature
expressions. Wall temperature comes from the function argument, not the CSV.
"""
from pathlib import Path
import numpy as np
import pandas as pd
import yaml
from scipy.interpolate import CubicSpline


def normalize_properties(frame):
    df = frame.copy()
    aliases = {
        'Temperature [K]': 'T_local [K]',
        'T_chamber [K]': 'T_local [K]',
        'P_chamber [bar]': 'P_local [bar]',
        'Pressure [bar]': 'P_local [bar]',
        'Viscosity [Pa*s]': 'Viscosity_local [Pa*s]',
    }
    for old, new in aliases.items():
        if old in df and new not in df:
            df = df.rename(columns={old: new})
    # Convert geometry columns to metres when necessary.
    for name in ('x', 'y'):
        target = f'{name}_m'
        if target not in df:
            if f'{name}_mm' in df:
                df[target] = pd.to_numeric(df[f'{name}_mm']) / 1000
            elif f'{name} [mm]' in df:
                df[target] = pd.to_numeric(df[f'{name} [mm]']) / 1000
            else:
                raise ValueError(f'Missing {target} or {name}_mm')
    # Allow alternate explicit Cp and viscosity units.
    if 'Cp [(KJ/kg-K)]' not in df and 'Cp [J/kg-K]' in df:
        df['Cp [(KJ/kg-K)]'] = pd.to_numeric(df['Cp [J/kg-K]']) / 1000
    if 'Viscosity_local [Pa*s]' not in df:
        for source, factor in [('Viscosity [cP]', 0.001), ('Viscosity [mPa*s]', 0.001)]:
            if source in df:
                df['Viscosity_local [Pa*s]'] = pd.to_numeric(df[source]) * factor
                break
    # Normalize the local specific gas constant for velocity calculations.
    if 'R [J/kg-K]' not in df and 'R [kJ/kg-K]' in df:
        df['R [J/kg-K]'] = pd.to_numeric(df['R [kJ/kg-K]'], errors='raise') * 1000
    # Local gas density [kg/m^3], using station pressure and temperature.
    if 'Density_free [kg/m3]' not in df and 'Density [kg/m3]' in df:
        df['Density_free [kg/m3]'] = df['Density [kg/m3]']
    if 'Density_free [kg/m3]' not in df:
        if 'P_local [bar]' not in df or 'T_local [K]' not in df:
            raise ValueError('Density needs local pressure and temperature or Density [kg/m3].')
        if 'R [J/kg-K]' in df:
            gas_constant = pd.to_numeric(df['R [J/kg-K]'], errors='raise')
        elif 'R [kJ/kg-K]' in df:
            gas_constant = pd.to_numeric(df['R [kJ/kg-K]'], errors='raise') * 1000
        else:
            raise ValueError('Density needs R [J/kg-K] or R [kJ/kg-K].')
        pressure = pd.to_numeric(df['P_local [bar]'], errors='raise') * 1e5
        temperature = pd.to_numeric(df['T_local [K]'], errors='raise')
        if (not np.isfinite(gas_constant).all() or (gas_constant <= 0).any()
                or not np.isfinite(pressure).all() or (pressure <= 0).any()):
            raise ValueError('Local pressure and gas constant must be finite and positive.')
        df['Density_free [kg/m3]'] = pressure / (gas_constant * temperature)
    required = ['x_m', 'y_m', 'gamma', 'Cp [(KJ/kg-K)]',
                'T_local [K]', 'Viscosity_local [Pa*s]', 'Prandtl Number', 'Mach', 'Density_free [kg/m3]', 'R [J/kg-K]']
    missing = [col for col in required if col not in df]
    if missing:
        raise ValueError('Missing properties CSV headers: ' + ', '.join(missing))
    for col in required:
        df[col] = pd.to_numeric(df[col], errors='raise')
    df = df.sort_values('x_m').reset_index(drop=True)
    if len(df) < 4 or not np.isfinite(df[required].to_numpy()).all():
        raise ValueError('At least 4 finite property stations are needed.')
    for col in ['gamma', 'R [J/kg-K]', 'Density_free [kg/m3]', 'y_m', 'Cp [(KJ/kg-K)]', 'T_local [K]', 'Viscosity_local [Pa*s]', 'Prandtl Number']:
        if (df[col] <= 0).any():
            raise ValueError(f'{col} must be positive.')
    if not (df.x_m.min() <= 0 <= df.x_m.max()):
        raise ValueError('Properties must span the throat at x=0.')
    # Repeated samples at the same axial station are common in exported CEA files.
    # Remove only duplicates whose required thermal inputs agree numerically.
    dup_count = int(df.x_m.duplicated().sum())
    if dup_count:
        for xpos, group in df.loc[df.x_m.duplicated(keep=False)].groupby('x_m'):
            for column in required:
                values = group[column].to_numpy(dtype=float)
                if not np.allclose(values, values[0], rtol=1e-9, atol=1e-12):
                    raise ValueError(
                        f'Conflicting values at x={xpos:.9g} m for {column}; '
                        'resolve these stations in the CSV before continuing.'
                    )
        df = df.drop_duplicates(subset='x_m', keep='first').reset_index(drop=True)
    if len(df) < 4:
        raise ValueError('At least 4 distinct axial stations are required.')
    return df


def find_converging_start(x, y):
    """Return the first station where radius falls below the chamber plateau.

    Assumes an initial constant-radius chamber and a decreasing-radius
    converging section upstream of the throat (x=0).
    """
    upstream = np.flatnonzero(x < 0)
    if len(upstream) < 2:
        raise ValueError('Need upstream chamber stations to detect convergence.')
    chamber_radius = float(y[upstream[0]])
    tolerance = max(1e-9, chamber_radius * 1e-5)
    falling = upstream[y[upstream] < chamber_radius - tolerance]
    if not len(falling):
        raise ValueError(
            'Cannot detect converging section from the geometry. '
            'Provide a contour with a constant-radius chamber followed by decreasing radius.'
        )
    # Include the preceding point so the thermal domain starts at the transition.
    return int(max(0, falling[0] - 1))


def bartz_correlation_rpe(*, x, y, tw, gamma, cp, T, mu_local, mu_am, pr, mach,
                             pc, cstar, rt, dt, dc, rc, ff, sigma, rho_free, rho_wall, rho_am, u):
    """Return the Bartz gas-side heat transfer coefficient from Rocket Propulsion Elements by Sutton [W/(m^2 K)].

    SI units: x/y/rt/dt/dc/rc [m], tw/T [K], cp [J/(kg K)], mu_local [Pa s],
    rho_free is free-stream gas density [kg/m^3], one value per station.
    rho_wall is gas density at the wall [kg/m^3], one value per station.
    rho_am is gas density at the arithmetic mean temperature [kg/m^3].
    u is local gas velocity [m/s], one value per station; rho * u is mass flux.
    pc [Pa], cstar [m/s]. gamma/pr/mach/ff/sigma are dimensionless.
    Station inputs are NumPy arrays; tw, pc, cstar, rt, dt, dc, rc are scalars.
    dc is the chamber diameter from the original upstream chamber radius.
    Use np for array operations. Return h_g [W/(m^2 K)] as one value per
    station or a constant. None hides the unfinished curve.
    """
    # Implement the proposed RPE comparison equation with D = throat diameter.
    # Stagnation viscosity uses the same approximate power law as mu_am.
    T0 = T * ff
    mu0 = mu_local * (T0 / T)**0.6
    return ((0.026 / dt**0.2) * (cp * mu_local**0.2 / pr**0.6)
            * (rho_free * u)**0.8 * (rho_am / rho_free)
            * (mu_am / mu0)**0.2)

def calculate_bartz(properties, wall_temperature, params_path='Inputs/TCA_params.yaml',
                    output_dir='Outputs', log=print, contour_geometry=None,
                    chamber_roughness_m=6.35e-6):
    """Return (heat_df, props_df, (x_ansys, hg_ansys, Tr_ansys))."""
    props = normalize_properties(pd.read_csv(properties) if isinstance(properties, (str, Path)) else properties)
    tw = float(wall_temperature)
    if not np.isfinite(tw) or tw <= 0:
        raise ValueError('Wall temperature must be a positive Kelvin value.')
    # Gas density at the wall, using the supplied wall temperature [K].
    props['T_average [K]'] = (props['T_local [K]'] + tw) / 2
    props['Density_wall [kg/m3]'] = (props['Density_free [kg/m3]'] * props['T_local [K]'] / tw).to_numpy(dtype=float)
    props['Density_am [kg/m3]'] = (props['Density_free [kg/m3]'] * props['T_local [K]'] / props['T_average [K]']).to_numpy(dtype=float)
    props['Viscosity_average [Pa*s]'] = (props['Viscosity_local [Pa*s]'] * (props['T_average [K]'] / props['T_local [K]'])**0.6).to_numpy(dtype=float)
    with open(params_path, encoding='utf-8') as f:
        params = yaml.safe_load(f)
    pc = float(params['chamber_pressure']) * 6894.757293  # psi -> Pa
    cstar = float(params['cstar'])  # m/s
    if pc <= 0 or cstar <= 0:
        raise ValueError('Chamber pressure and cstar must be positive.')
    x = props.x_m.to_numpy()
    y = props.y_m.to_numpy()
    # Use the original full geometry for throat curvature when available:
    # the CEA properties output excludes the throat station.
    if contour_geometry is not None:
        geom = contour_geometry.sort_values('x_mm')
        xs = geom.x_mm.to_numpy(dtype=float) / 1000
        ys = geom.y_mm.to_numpy(dtype=float) / 1000
    else:
        xs, ys = x, y
    xu, ind = np.unique(xs, return_index=True)
    if len(xu) < 4:
        raise ValueError('Need at least 4 distinct geometry stations.')
    spline = CubicSpline(xu, ys[ind])
    rt = float(spline(0.0))
    if rt <= 0 or not np.isfinite(rt):
        raise ValueError('Invalid throat radius.')
    curvature = abs(float(spline(0., 2)))
    if curvature < 1e-12:
        raise ValueError('Throat curvature is near zero; cannot calculate rc.')
    rc = (1 + float(spline(0., 1))**2)**1.5 / curvature
    dt = 2 * rt
    # Do not apply the Bartz calculation to the straight combustion chamber.
    # Detect the beginning of convergence from the original contour if supplied,
    # otherwise from the geometry embedded in the properties CSV.
    start_index = find_converging_start(xs, ys)
    # Use the upstream chamber radius before restricting the thermal domain.
    dc = 2 * float(ys[np.flatnonzero(xs < 0)[0]])  # chamber diameter [m]
    start_x = float(xs[start_index])
    active = x >= start_x
    if np.count_nonzero(active) < 4:
        raise ValueError('Fewer than four stations remain from convergence onward.')
    log(f'Gnielinski upstream of x={start_x:.6g} m; Bartz from convergence onward.')
    gamma = props.gamma.to_numpy()
    cp = props['Cp [(KJ/kg-K)]'].to_numpy() * 1000  # kJ -> J
    T = props['T_local [K]'].to_numpy()
    mu_local = props['Viscosity_local [Pa*s]'].to_numpy()
    # am means arithmetic mean, same nomenclature as seen in RPE
    mu_am = props['Viscosity_average [Pa*s]'].to_numpy()
    pr = props['Prandtl Number'].to_numpy()
    mach = props['Mach'].to_numpy()
    rho_wall = props['Density_wall [kg/m3]'].to_numpy(dtype=float)
    rho_am = props['Density_am [kg/m3]'].to_numpy(dtype=float)
    rho_free = props['Density_free [kg/m3]'].to_numpy(dtype=float)
    R = props['R [J/kg-K]'].to_numpy(dtype=float)
    u = mach * np.sqrt(gamma * R * T)  # local gas velocity [m/s]
    if not np.isfinite(u).all() or (u < 0).any():
        raise ValueError('Local gas velocity must be finite and nonnegative.')
    ff = 1 + ((gamma - 1) / 2) * mach**2
    sigma = 1 / ((0.5 * (tw / T) + 0.5)**0.68 * ff**0.12)
    # Preserved from your original array-level Bartz equation: (At/A)^0.9
    hg = (0.026 / dt**0.2) * (mu_local**0.2 * cp / pr**0.6) * \
         (pc / cstar)**0.8 * (dt / rc)**0.1 * (rt / y)**1.8 * sigma
    tr = T * (1 + (gamma - 1) / 2 * pr**(1 / 3) * mach**2)
    chamber = ~active
    reynolds = np.full(len(x), np.nan)
    friction = np.full(len(x), np.nan)
    nusselt = np.full(len(x), np.nan)
    if chamber.any():
        from constant_bc_loader import load_correlations
        bc = load_correlations()
        roughness = float(chamber_roughness_m)
        mdot = float(params['tp_mdot'])  # total combustion-gas flow [kg/s]
        if not np.isfinite(roughness) or roughness < 0 or not np.isfinite(mdot) or mdot <= 0:
            raise ValueError('Chamber roughness must be nonnegative and tp_mdot positive (SI units).')
        for i in np.flatnonzero(chamber):
            reynolds[i] = 4 * mdot / (np.pi * dc * mu_local[i])
            f = bc.friction_factor(roughness, dc, reynolds[i])
            nu = bc.gnielinski(f, reynolds[i], pr[i]) if f is not None else None
            if nu is None or not np.isfinite(nu) or nu <= 0:
                raise ValueError(f'Chamber Gnielinski outside validity at x={x[i]:.6g} m: '
                                 f'Re={reynolds[i]:.6g}, Pr={pr[i]:.6g}, roughness/D={roughness/dc:.6g}.')
            friction[i], nusselt[i] = f, nu
            # k = mu * cp / Pr, using local CEA transport properties.
            hg[i] = nu * mu_local[i] * cp[i] / (pr[i] * dc)
        log(f'Chamber Gnielinski uses tp_mdot={mdot:g} kg/s, D={dc:g} m, roughness={roughness:g} m.')
    if not np.isfinite(hg).all() or not np.isfinite(tr).all():
        raise ValueError('Non-finite gas-side results; check your input data.')
    heat = pd.DataFrame({'x_m': x, 'Tw_K': np.full(len(x), tw),
                         'hg_W_m2K': hg, 'Tr_K': tr,
                         'correlation': np.where(chamber, 'Gnielinski', 'Bartz'),
                         'Re_chamber': reynolds, 'Darcy_f_chamber': friction,
                         'Nu_chamber': nusselt})
    hg_custom = bartz_correlation_rpe(
        x=x, y=y, tw=tw, gamma=gamma, cp=cp, T=T, mu_local=mu_local, mu_am=mu_am, pr=pr,
        mach=mach, pc=pc, cstar=cstar, rt=rt, dt=dt, dc=dc, rc=rc,
        ff=ff, sigma=sigma, rho_free=rho_free, rho_wall=rho_wall, rho_am=rho_am, u=u)
    if hg_custom is not None:
        hg_custom = np.asarray(hg_custom, dtype=float)
        if hg_custom.ndim == 0:
            hg_custom = np.full(len(x), float(hg_custom))
        if hg_custom.shape != x.shape:
            raise ValueError('Custom Bartz must return one h_g per station or a scalar.')
        if not np.isfinite(hg_custom).all() or (hg_custom <= 0).any():
            raise ValueError('Custom Bartz h_g values must be finite and positive.')
        heat['hg_custom_W_m2K'] = np.where(active, hg_custom, np.nan)
        log('RPE comparison correlation added to the comparison graph.')
    else:
        log('Custom Bartz equation is empty; comparison curve skipped.')
    # Calculate in memory only. No files are written by this function.
    x_ansys = np.linspace(x.min(), x.max(), 500)
    # Interpolate each region separately so a correlation jump cannot cause
    # cubic overshoot or blend the two models across the chamber transition.
    hg_ansys = np.empty_like(x_ansys)
    for mask, query in [(chamber, x_ansys < start_x), (active, x_ansys >= start_x)]:
        if query.any() and mask.any():
            hg_ansys[query] = np.interp(x_ansys[query], x[mask], hg[mask])
    tr_ansys = np.interp(x_ansys, x, tr)
    log('Bartz and ANSYS arrays calculated in memory; nothing saved.')
    return heat, props, (x_ansys, hg_ansys, tr_ansys)

if __name__ == '__main__':
    heat, _, _ = calculate_bartz('Outputs/properties.csv', 800)
    print(heat.head())
    print('No files saved. Use thermal_gui.py to export.')
