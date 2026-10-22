import numpy as np
from Impeller import DesignPoint, InputGeometry, Impeller
from pint import Quantity as Q_
#import CoolProp.CoolProp as CP
from pathlib import Path
import yaml

script_dir = Path(__file__).resolve().parent
yaml_path = script_dir.parent / "Inputs/Inputs.yaml"	# Pumps/Inputs

with open(yaml_path, "r") as file:
	inputs = yaml.safe_load(file)
dp = inputs['design_point']

g = Q_(9.81, 'm/s^2')				# gravitational acceleration

opt_m_dot = Q_(
	dp['IPA']['mass_flow'], 'kg/s')	# mass flow rate at BEP
rho = Q_(786, 'kg/m^3')				# fluid density

feed_pressure = Q_(					# inlet feed pressure, from tank pressure
	dp['IPA']['feed_pressure'], 'bar')

discharge_pressure = Q_(
	dp['IPA']['discharge_pressure'], 'bar')	# ideal discharge pressure

opt_dP = discharge_pressure - feed_pressure	# total pressure rise at BEP

design_point: DesignPoint = DesignPoint(
	Q = opt_m_dot/rho,				# volumetric flow rate at BEP
	H = (opt_dP / (g * rho)),    	# developed head at BEP
	N_shaft = Q_(dp['N_shaft'], 'rpm'),	# shaft speed
	n_hyd_BEP = 0.5,
	rho=rho
)
ipa_geometry: InputGeometry = InputGeometry(
	Z_blade = 6,
	Beta2B = Q_(25,'deg').to('rad'),# blade angle at exit, relative to tangent
    d_1 = Q_(1.2, 'in'),			# impeller inlet diameter
	d_2 = Q_(2.8, 'in'),			# impeller outlet diameter
	b_2 = Q_(0.08, 'in'),			# impeller outlet height
	blade_thick = Q_(0.04, 'in'),			# blade thickness at exit
    d_hub=Q_(0.7, 'in'),			# hub diameter
    
	Z_e_factor = 1.5,				# axial length extension factor, modifying Gulich formula
    shroud_fullness=(1, 0.5),
    hub_fullness=(0.7, 0.7),
    shroud_outlet_angle=Q_(5, 'deg')	# shroud contour angle from radial at outlet
)

ipa_imp = Impeller(ipa_geometry, design_point)