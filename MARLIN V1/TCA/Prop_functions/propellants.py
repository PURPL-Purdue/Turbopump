from rocketcea.cea_obj import add_new_fuel
import matplotlib.pyplot as plt
from funct import get_isp_heatmap
from funct import get_tempstar_heatmap


##### To add acetone - https://atct.anl.gov/Thermochemical%20Data/version%201.130/species/?species_number=54 for Enthalpy and https://pubchem.ncbi.nlm.nih.gov/compound/Acetone for density
acetone_card = """
fuel Acetone  C 3 H 6 O 1
h,kj/mol=-216.94
t(k)=298.15
rho=0.79
"""
add_new_fuel('Acetone', acetone_card)
#####

Fuel = 'JetA' 
rho_fuel = 790  # kg/m^3 (0.79 g/cm^3 acetone density)
rho_ox = 1141  # kg/m^3 (LOX density)
Ox = 'LOX'
mr_range = (0.55, 6.0)
p_range = (15, 60)  # bar
pamb = 1.01325  # bar

Mr, Pin, isp_grid, mr_stoich, density = get_isp_heatmap(
    Fuel=Fuel,
    Ox=Ox,
    mr_range=mr_range,
    p_range=p_range,
    pamb=pamb,
    rho_fuel=rho_fuel,
    rho_ox=rho_ox,
)
Mr, Pin, tempstar_grid, mr_stoich = get_tempstar_heatmap(
    Fuel=Fuel,
    Ox=Ox,
    mr_range=mr_range,
    p_range=p_range,
    pamb=pamb,
)
fig, axes = plt.subplots(1, 3, figsize=(24, 6))
mesh1 = axes[0].pcolormesh(Pin, Mr, isp_grid, shading='auto', cmap='viridis')
axes[0].axhline(mr_stoich, color='red', linestyle='-', linewidth=1, label=f'Stoich O/F = {mr_stoich:.2f}')
fig.colorbar(mesh1, ax=axes[0], label='Ambient Isp (s)')
axes[0].set_xlabel('Chamber Pressure (bar)')
axes[0].set_ylabel('O/F Ratio')
axes[0].set_title(f'Specific Impulse Heatmap ({Ox} / {Fuel})')
axes[0].legend()

mesh2 = axes[1].pcolormesh(Pin, Mr, tempstar_grid, shading='auto', cmap='viridis')
axes[1].axhline(mr_stoich, color='red', linestyle='-', linewidth=1, label=f'Stoich O/F = {mr_stoich:.2f}')
fig.colorbar(mesh2, ax=axes[1], label='Temperature Star (K)')
axes[1].set_xlabel('Chamber Pressure (bar)')
axes[1].set_ylabel('O/F Ratio')
axes[1].set_title(f'Throat Temperature Heatmap ({Ox} / {Fuel})')
axes[1].legend()

mesh3 = axes[2].pcolormesh(Pin, density, isp_grid, shading='auto', cmap='viridis')
axes[2].set_xlabel('Chamber Pressure (bar)')
axes[2].set_ylabel('Density (kg/m^3)')
axes[2].set_title(f'Isp Heatmap ({Ox} / {Fuel})')

plt.tight_layout()
plt.show()
