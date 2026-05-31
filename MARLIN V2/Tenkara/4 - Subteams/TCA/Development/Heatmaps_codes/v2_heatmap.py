import numpy as np
import pint
import matplotlib.pyplot as plt

from v2_heatmap_functions import get_isp_heatmap, get_Tc_heatmap, get_isp_density_heatmap

ureg = pint.UnitRegistry()

### Engine Parameters/Inputs
of_range = (0.55, 6.0)
Pc_range = (15, 60)  * ureg.bar # bar
Pamb = 1.01325 *ureg.bar # bar
Pexit = 18 * ureg.psi # currently unused: see later note on pressure ratio
rho_ox = 1141 * ureg.kg / ureg.m**3 # LOX density
rho_f = 786 * ureg.kg / ureg.m**3 #IPA and ethanol density

### CEA Setup ### 

#reac_names = ["C3H8O,2propanol", "O2(L)"]
reac_names = ["C2H5OH(L)", "O2(L)"]
reac_temps = np.array([298.15, 90.17]) * ureg.K
fuel_weights = np.array([1.0, 0.0])
oxidant_weights = np.array([0.0, 1.0])
Prat = (Pc_range / Pamb).to_base_units().magnitude # currently unused:
# perfect expansion assumed for parity with initial code. will want to change in design consideration
# need to implement the heatmap function to calculate pressure ratio for every case, just have to pick a design exit pressure 
cont_ratio = 1.8


## create isp heatmap grid
of_set, Pc_set, isp_grid = get_isp_heatmap(ureg, of_range, Pc_range, reac_names, oxidant_weights, fuel_weights, reac_temps, Prat, cont_ratio)

## plot isp heatmap
fig, axes = plt.subplots(1, 3, figsize=(24, 6))
mesh1 = axes[0].pcolormesh(Pc_set, of_set, isp_grid, shading='auto', cmap='viridis')
#axes[0].axhline(mr_stoich, color='red', linestyle='-', linewidth=1, label=f'Stoich O/F = {mr_stoich:.2f}')
fig.colorbar(mesh1, ax=axes[0], label='Ambient Isp (s)')
axes[0].set_xlabel('Chamber Pressure (bar)')
axes[0].set_ylabel('O/F Ratio')
axes[0].set_title(f'Specific Impulse Heatmap ({reac_names[0]} / {reac_names[1]})')
axes[0].legend()


## create chamber temp heatmap grid
of_set, Pc_set, Tc_grid = get_Tc_heatmap(ureg, of_range, Pc_range, reac_names, oxidant_weights, fuel_weights, reac_temps, Prat, cont_ratio)

## plot chamber temp heatmap
mesh2 = axes[1].pcolormesh(Pc_set, of_set, Tc_grid, shading='auto', cmap='viridis')
#axes[1].axhline(mr_stoich, color='red', linestyle='-', linewidth=1, label=f'Stoich O/F = {mr_stoich:.2f}')
fig.colorbar(mesh2, ax=axes[1], label='Chamber Temperature (K)')
axes[1].set_xlabel('Chamber Pressure (bar)')
axes[1].set_ylabel('O/F Ratio')
axes[1].set_title(f'Chamber Temperature Heatmap ({reac_names[0]} / {reac_names[1]})')
axes[1].legend()


## create isp density heatmap grid
of_set, Pc_set, isp_density_grid = get_isp_density_heatmap(ureg, rho_ox, rho_f, of_range, Pc_range, reac_names, oxidant_weights, fuel_weights, reac_temps, Prat, cont_ratio)

## plot isp denisty heatmap
mesh3 = axes[2].pcolormesh(Pc_set, of_set, isp_density_grid, shading='auto', cmap='viridis')
#axes[2].axhline(mr_stoich, color='red', linestyle='-', linewidth=1, label=f'Stoich O/F = {mr_stoich:.2f}')
fig.colorbar(mesh2, ax=axes[2], label='Isp Density (kg*s/m^3)')
axes[2].set_xlabel('Chamber Pressure (bar)')
axes[2].set_ylabel('O/F Ratio')
axes[2].set_title(f'Isp Density Heatmap ({reac_names[0]} / {reac_names[1]})')
axes[2].legend()

plt.show()


