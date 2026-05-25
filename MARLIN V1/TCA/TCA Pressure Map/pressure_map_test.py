import numpy as np
import matplotlib.pyplot as plt
from pressure_map import pressure_map
from CoolProp.CoolProp import PropsSI
from matplotlib.lines import Line2D

# Unit conversion constants 
ft_to_m = 0.3048
meters_to_inches = 39.3701 
psi_to_pa = 6894.76
kg_to_lbm = 2.20462

def finite_minmax(values): #Return (min, max) cpnsidering only finites entries. 
    values = np.asarray(values)
    # true only where values are real numbers
    finite_mask = np.isfinite(values)
    if not np.any(finite_mask):
        raise ValueError("No valid values found.")
    return float(np.min(values[finite_mask])), float(np.max(values[finite_mask]))

def plot_pc_heatmap_with_of_contours(fuel_pressure_map, ox_pressure_map, PC_pfpo, OF_pfpo, mdot_pfpo,
                                     of_levels = np.linspace(1.3, 2.0, 8), pc_vmin = None, pc_vmax = None):
    """
    Plot Pc heatmap in (p_ox, p_fuel) space with OF contour lines.
    Assuming that PC_pfpo / OF_pfpo are indexed as [fuel_index, ox_index]
    """
    # infer pressure axes range (Pa) from the forward maps, ignoring any NaN values
    #p_fuel_axis = np.asarray(fuel_pressure_map)
    #p_ox_axis = np.asarray(ox_pressure_map)
    pf_min, pf_max = finite_minmax(fuel_pressure_map)
    po_min, po_max = finite_minmax(ox_pressure_map)
    # building axes that matches the inverse map indexing
    p_fuel_axis = np.linspace(pf_min, pf_max, PC_pfpo.shape[0])
    p_ox_axis = np.linspace(po_min, po_max, PC_pfpo.shape[1])
    # Convert window bounds -> index ranges
    i0 = np.searchsorted(p_fuel_axis, 200.0, side="left")
    i1 = np.searchsorted(p_fuel_axis, 1200.0, side="right")
    j0 = np.searchsorted(p_ox_axis, 200.0, side="left")
    j1 = np.searchsorted(p_ox_axis, 1200.0, side="right")
    # doing zoom in the graph 
    PC_view = PC_pfpo[i0:i1, j0:j1]
    OF_view = OF_pfpo[i0:i1, j0:j1]
    mdot_view = mdot_pfpo[i0:i1, j0:j1]
    p_fuel_view = p_fuel_axis[i0:i1]
    p_ox_view   = p_ox_axis[j0:j1]
    PF, PO = np.meshgrid(p_fuel_view, p_ox_view, indexing="ij")
    # Autochoose OF contour levels safely
    if of_levels is None:
        of_finite = OF_pfpo[np.isfinite(OF_pfpo)]
        if of_finite.size > 0:
            lo, hi = np.percentile(of_finite, [5, 95])
            # if percentile collapses, fall back to min/max range with a small margin
            if not np.isfinite(lo) or not np.isfinite(hi) or (lo >= hi):
                lo, hi = float(np.min(of_finite)), float(np.max(of_finite))
            of_levels = np.round(np.linspace(1.3, 2.0, 8))
        else:
            of_levels = np.round(np.linspace(1.3, 2.0, 8))
    # Plotting
    fig, ax = plt.subplots(figsize=(10, 8), constrained_layout=True)    
    # heatmap where the x-axis is oxidizer pressure, y axis will be fuel pressure
    extent_view = [p_ox_view.min(), p_ox_view.max(), p_fuel_view.min(), p_fuel_view.max()]
    im = ax.imshow(PC_view, origin = "lower", aspect='auto', extent=extent_view, vmin=pc_vmin, vmax=pc_vmax)
    # Force axes to match the cropped window exactly
    ax.set_xlim(p_ox_view.min(), p_ox_view.max())
    ax.set_ylim(p_fuel_view.min(), p_fuel_view.max())
    cb = fig.colorbar(im, ax=ax)
    cb.set_label('Chamber Pressure (psi)')
    mdot_masked = np.ma.masked_invalid(mdot_view)
    OF_masked = np.ma.masked_invalid(OF_view)
    mdot_target = 23.630/kg_to_lbm # convert target mass flow rate to lbm/s for contouring
    cs_mdot = ax.contour(PO, PF, mdot_masked, levels= [mdot_target], colors='r', linewidths=1)
    mdot_handle = Line2D([0], [0], color='r', lw=2, label=f'ṁ = {mdot_target:.1f} kg/s')
    #mdot_labels = ax.clabel(cs_mdot, fmt=lambda v: f"ṁ = {v:.1f} kg/s", inline=True, fontsize=8)
    cs_mdot.set_clip_on(True) # Ensure contour lines are clipped to the axes limits
    cs_mdot.set_clip_path(ax.patch) # Clip to the axes patch (the plotting area)
    #for t in mdot_labels:
        #t.set_clip_on(True)
        #t.set_clip_path(ax.patch)
    ps = ax.contour(PO, PF, OF_masked, levels = of_levels, colors = "k", linewidths = 1)
    of_labels = ax.clabel(ps, fmt = lambda v: f"OF = {v:.2f}", inline=True, fontsize=8)
    ps.set_clip_on(True)
    ps.set_clip_path(ax.patch)
    for t in of_labels:
        t.set_clip_on(True)
        t.set_clip_path(ax.patch)
    pc_regulated = 500.0106 # psi
    pc_common = ax.contour(PO, PF, PC_view, levels=[pc_regulated], colors='b', linewidths=1.5)
    pc_handle = Line2D([0], [0], color='b', lw=2, label=f'Pc = {pc_regulated:.0f} psi')
    handles = [mdot_handle, pc_handle]
    ax.legend(handles=handles, loc='upper left', frameon=True)
    #ax.clabel(pc_common, fmt=lambda v: f"Pc = {v:.0f} psi", inline=True, fontsize=10)
    ax.set_xlabel('Oxidizer Feed Pressure (psi)')
    ax.set_ylabel('Fuel Feed Pressure (psi)')
    ax.set_title('Chamber Pressure Heatmap with OF Contours')
    plt.show()

def plot_maps(pc_axis_psi, of_axis, fuel_p_feed_map_psi, ox_p_feed_map_psi, Pc_map_psi, OF_map):
    # extents for imshow: [x_min, x_max, y_min, y_max]
    extent_forward = [of_axis.min(), of_axis.max(), pc_axis_psi.min(), pc_axis_psi.max()]
    # inverse maps are indexed as [fuel_index, ox_index]
    pf_min, pf_max = finite_minmax(fuel_p_feed_map_psi)
    po_min, po_max = finite_minmax(ox_p_feed_map_psi)
    extent_inv = [po_min, po_max, pf_min, pf_max] 

    fig, axs = plt.subplots(2, 2, figsize=(12, 10), constrained_layout=True)

    # (0,0) Forward fuel map
    im0 = axs[0,0].imshow(np.ma.masked_invalid(fuel_p_feed_map_psi), origin='lower', aspect='auto', extent=extent_forward)

    axs[0,0].set_title('Fuel Feed Pressure Map (psi)')
    axs[0,0].set_xlabel('OF Ratio')
    axs[0,0].set_ylabel('Chamber Pressure (psi)')
    cb0 = fig.colorbar(im0, ax=axs[0,0])
    cb0.set_label('Fuel Feed Pressure (psi)')

    # (0,1) Forward oxidizer map
    im1 = axs[0,1].imshow(np.ma.masked_invalid(ox_p_feed_map_psi), origin='lower', aspect='auto', extent=extent_forward)
    axs[0,1].set_title('Oxidizer Feed Pressure Map (psi)')
    axs[0,1].set_xlabel('OF Ratio')
    axs[0,1].set_ylabel('Chamber Pressure (psi)')
    cb1 = fig.colorbar(im1, ax=axs[0,1])
    cb1.set_label('Ox Feed Pressure (psi)')

    # (1,0) Inverse chamber pressure map
    im2 = axs[1,0].imshow(np.ma.masked_invalid(Pc_map_psi), origin='lower', aspect='auto', extent=extent_inv)
    axs[1,0].set_title('Chamber Pressure from Inverse Map (psi)')
    axs[1,0].set_xlabel('Oxidizer Feed Pressure (psi)')
    axs[1,0].set_ylabel('Fuel Feed Pressure (psi)')
    cb2 = fig.colorbar(im2, ax=axs[1,0])
    cb2.set_label('Pc (psi)')

    im3 = axs[1,1].imshow(np.ma.masked_invalid(OF_map), origin='lower', aspect='auto', extent=extent_inv)
    axs[1,1].set_title('OF Ratio from Inverse Map')
    axs[1,1].set_xlabel('Oxidizer Feed Pressure (psi)')
    axs[1,1].set_ylabel('Fuel Feed Pressure (psi)')
    cb3 = fig.colorbar(im3, ax=axs[1,1])
    cb3.set_label('O/F')

    plt.show()

# === Parameters for the pressure map ===
min_OF = 1.0 # Minimum O/F ratio to consider in the map
max_OF = 2.0 # Maximum O/F ratio to consider in the map
max_pc = 700 # Maximum chamber pressure (psi) to consider in the map
resolution = 400 # Number of points along each axis (OF and Pc) for the forward maps. Higher means smoother but slower.
cstar_eff = 0.875 # C* efficiency to apply to the ideal C* calculation, accounting for non-idealities in the combustion process. Should be between 0 and 1, where 1 means ideal performance. Adjusting this can help match the model to experimental data if available.

tca_throat_diameter_inches = 3.0746 # Throat diameter in inches, used to calculate throat area. Adjust this value based on your actual nozzle throat size.
tca_throat_diameter_meters = tca_throat_diameter_inches / meters_to_inches
throat_area = np.pi * (tca_throat_diameter_meters / 2) ** 2 # Throat area in m^2

# Effective injector / orifice CdA (m^2)
# These represent the *flow-limiting elements* on each side.
# Consider the injector holes, if we have a lot of holes. CdA_total = N_holes * CdA_per_hole.
number_of_fuel_holes = 80
Cd = 0.8 # Discharge coefficient for the injector holes, accounting for non-ideal flow. 
diameter_fuel_hole_meters = 0.0015113 # Diameter of each fuel injector hole in meters. [1.2 mm]
diameter_ox_hole_meters = 0.0017018 # Diameter of each oxidizer injector hole in meters. [1.7018 mm]
fuel_CdA = number_of_fuel_holes * Cd * np.pi * (diameter_fuel_hole_meters / 2) ** 2 # Effective flow area for fuel side (m^2)
ox_CdA = number_of_fuel_holes * Cd * np.pi * (diameter_ox_hole_meters / 2) ** 2 # Effective flow area for oxidizer side (m^2)

# Fluid properties (density). Using approximation inlet conditions. 
P_ref = 700*psi_to_pa # Reference pressure for density calculation (Pa)
rho_fuel = 789 # PropsSI('D', 'T', 298, 'P', P_ref, 'Isopropanol') fuel density at 298 K [kg/m^3]
rho_ox = PropsSI('D', 'T', 90, 'P', P_ref, 'Oxygen') # kg/m^3 

mdot_target = 23.630/kg_to_lbm #[kg/s]

pc_axis_psi, of_axis, fuel_p_feed_map_psi, ox_p_feed_map_psi, Pc_map_psi, OF_map, mdot_total_map, fuel_feed_axis_psi, ox_feed_axis_psi = pressure_map(min_OF, max_OF, max_pc, resolution, cstar_eff, throat_area, fuel_CdA, ox_CdA, rho_fuel, rho_ox)
plot_pc_heatmap_with_of_contours(fuel_p_feed_map_psi, ox_p_feed_map_psi, Pc_map_psi, OF_map, mdot_total_map, of_levels = np.arange(min_OF + 0.2, max_OF, 0.2), pc_vmin=0, pc_vmax=max_pc)
plot_maps(pc_axis_psi, of_axis, fuel_p_feed_map_psi, ox_p_feed_map_psi, Pc_map_psi, OF_map)

