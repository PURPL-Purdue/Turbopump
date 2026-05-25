import numpy as np
from rocketcea.cea_obj import CEA_Obj
from scipy import interpolate

# Unit conversion factors
psi_to_pa = 6894.76
foot_to_meters = 0.3048

def valid_range(values):
    values = np.asarray(values) # Convert input to Numpy array (so it works with lists or arrays)
    mask = np.isfinite(values) # Filter out non-finite values (NaN, inf)
    if not np.any(mask): # If no valid values are found, raise an error to avoid returning invalid ranges
        raise ValueError("No valid values found.")
    # Return the minimum and maximum of the valid values as floats (no NaN or inf)
    return float(np.min(values[mask])), float(np.max(values[mask]))

def build_inverse_interpolator(fuel_feed_pressure_map, ox_feed_pressure_map, chamber_pressure_axis, mixture_ratio_axis, mdot_total_forward, 
                               fill_value=np.nan):
    # Reconstruct the original grid that matches the foward map dimensions
    PC_grid, OF_grid = np.meshgrid(chamber_pressure_axis, mixture_ratio_axis, indexing='ij')
    # Flatten grid and both feed-pressure maps into 1D vectors. Each (i, j) cell becomes one sample in a scaterred dataset.
    p_fuel_feed = np.ravel(fuel_feed_pressure_map) # input coordinate 1 (fuel feed pressure)
    p_ox_feed = np.ravel(ox_feed_pressure_map) # input coordinate 2 (oxidizer feed pressure)
    Pc_values = np.ravel(PC_grid) # output value 1 (chamber pressure)
    OF_values = np.ravel(OF_grid) # output value 2 (mixture ratio) 
    mdot_values = np.ravel(mdot_total_forward) # output value 3 (mass flow rate)

    # Build a validity mask we can only interpolate using finite points
    valid = (np.isfinite(p_fuel_feed) & np.isfinite(p_ox_feed) & np.isfinite(Pc_values) & np.isfinite(OF_values) & np.isfinite(mdot_values))
    # Assemble the 2D input point cloud:
    feed_pressure_points = np.column_stack((p_fuel_feed[valid], p_ox_feed[valid]))
    # Fit an inverse interpolator for chamber pressure.
    chamber_pressure_interpolator = interpolate.LinearNDInterpolator(feed_pressure_points, Pc_values[valid], fill_value=fill_value)
    # Fit an inverse interpolator for mixture ratio.
    mixture_ratio_interpolator = interpolate.LinearNDInterpolator(feed_pressure_points, OF_values[valid], fill_value=fill_value)
    # Fit an inverse interpolator for mass flow rate 
    mdot_total_interpolator = interpolate.LinearNDInterpolator(feed_pressure_points, mdot_values[valid], fill_value=fill_value)
    return chamber_pressure_interpolator, mixture_ratio_interpolator, mdot_total_interpolator

def pressure_map(minOF, maxOF, maxPc, resolution, cstar_eff, throat_area, fuel_CdA, ox_CdA, rho_fuel, rho_ox):
    cea = CEA_Obj(oxName='LOX', fuelName='Isopropanol') # Create a RocketCEA object to compute c*(Pc, MR) for LOX / Kerosene
    pc_axis_psi = np.linspace(1.0, maxPc, resolution) # Define the forward grid in chamber conditions
    of_axis = np.linspace(minOF, maxOF, resolution)
    # For each (Pc, OF), store the required FEED pressure of fuel and oxidizer. 
    fuel_p_feed_map_psi = np.full((resolution, resolution), np.nan, dtype=np.float64) 
    ox_p_feed_map_psi = np.full((resolution, resolution), np.nan, dtype=np.float64)
    mdot_total_map = np.full((resolution, resolution), np.nan, dtype=np.float64) # Optional: store total mass flow rate for reference
    for i, pc in enumerate(pc_axis_psi):
        for j, of in enumerate(of_axis):
            # Combustion performance
            cstar = cea.get_Cstar(Pc=pc, MR=of)
            cstar_eff_real = cstar_eff * cstar
            p_c_pa = pc * psi_to_pa
            # Calculate mass flow rate
            mdot_total = (p_c_pa * throat_area) / (cstar_eff_real * foot_to_meters)
            mdot_fuel = mdot_total / (1 + of)
            mdot_ox = mdot_total - mdot_fuel
            # Incompressible flow assumption required for pressure drop calculations
            dp_fuel = (mdot_fuel / fuel_CdA) ** 2 / (2.0 * rho_fuel)
            dp_ox = (mdot_ox / ox_CdA) ** 2 / (2.0 * rho_ox)
            # required feed pressures
            fuel_p_feed_map_psi[i, j] = (p_c_pa + dp_fuel) / psi_to_pa
            ox_p_feed_map_psi[i, j] = (p_c_pa + dp_ox) / psi_to_pa
            mdot_total_map[i, j] = mdot_total
    # Building inverse interpolators
    Pc_from_feed, OF_from_feed, mdot_from_feed = build_inverse_interpolator(fuel_p_feed_map_psi, ox_p_feed_map_psi, pc_axis_psi, of_axis, mdot_total_map)
    # Autogenerate a feed pressure grid for testing
    pf_min, pf_max = valid_range(fuel_p_feed_map_psi)
    po_min, po_max = valid_range(ox_p_feed_map_psi)
    fuel_feed_axis_psi = np.linspace(pf_min, pf_max, resolution)
    ox_feed_axis_psi = np.linspace(po_min, po_max, resolution)
    PF, PO = np.meshgrid(fuel_feed_axis_psi, ox_feed_axis_psi, indexing='ij')
    # Evaluating inverse maps on the feed grid
    Pc_map_psi = Pc_from_feed(PF, PO)
    OF_map = OF_from_feed(PF, PO)
    mdot_map = mdot_from_feed(PF, PO)
    # Apply validity mask to filter out any points where interpolation failed or returned non-physical values. 
    invalid = ( ~np.isfinite(Pc_map_psi) | ~np.isfinite(OF_map) | ~np.isfinite(mdot_map) | (Pc_map_psi <= 0) | (Pc_map_psi > maxPc) | 
        (OF_map < minOF) | (OF_map > maxOF) )
    Pc_map_psi[invalid] = np.nan
    OF_map[invalid] = np.nan
    mdot_map[invalid] = np.nan
    return (pc_axis_psi, of_axis, fuel_p_feed_map_psi, ox_p_feed_map_psi, Pc_map_psi, OF_map, mdot_map, fuel_feed_axis_psi, ox_feed_axis_psi)
