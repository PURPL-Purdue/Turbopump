"""
Calls contour_script.py to generate a CSV and DXF files in Output directory.
"""

import yaml
from pint import UnitRegistry
from matplotlib import contour

from contour_script import bell_nozzle
from contour_script import export_nozzle_csv
from contour_script import export_nozzle_dxf

ureg = UnitRegistry()

with open('Inputs/TCA_params.yaml') as f:
    p = yaml.safe_load(f)

eratio = p['expansion_ratio']
Rt = p['throat_diameter']/2 * ureg.inch
l_percent = p['l_percent'] 
cratio = p['contraction_ratio']
alpha = p['alpha_divergence'] * ureg.deg
Lc = p['chamber_length'] * ureg.m



angles, contour, R2 = bell_nozzle(eratio, Rt.to(ureg.mm).magnitude, l_percent, cratio, alpha.magnitude, Lc.to(ureg.mm).magnitude)
export_nozzle_csv(contour, 'Outputs/contour.csv')
export_nozzle_dxf(contour, 'Outputs/contour.dxf')
