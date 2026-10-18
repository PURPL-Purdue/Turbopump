import numpy as np
from Impeller import DesignPoint, InputGeometry, Impeller
from pint import Quantity as Q_
import CoolProp.CoolProp as CP
from pathlib import Path
import yaml

script_dir = Path(__file__).resolve().parent
yaml_path = script_dir.parents[1] / "Inputs/Inputs.yaml"	# Pumps/Inputs

with open(yaml_path, "r") as file:
	inputs = yaml.safe_load(file)
dp = inputs['design_point']

g = Q_(9.81, 'm/s^2')					# gravitational acceleration

opt_m_dot = Q_(
	dp['LOX']['mass_flow'], 'kg/s')		# mass flow rate at BEP
rho = Q_(1141, 'kg/m^3')				# density of LOX at saturation, atm pressure (lol?)

feed_pressure = Q_(150, 'psi')			# inlet feed pressure, from tank pressure
discharge_pressure = Q_(
	dp['LOX']['discharge_pressure'], 'bar')	# ideal discharge pressure
opt_dP = discharge_pressure - feed_pressure	# total pressure rise at BEP

design_point: DesignPoint = DesignPoint(
	Q = opt_m_dot/rho,					# volumetric flow rate at BEP
	H = (opt_dP / (g * rho)).to('m'),	# developed head at BEP
	N_shaft = Q_(dp['N_shaft'], 'rpm'),	# shaft speed
	n_hyd_BEP = 0.5,
	rho=rho
)

lox_geometry: InputGeometry = InputGeometry(
	Z_blade = 6,
	Beta2B = Q_(20,'deg').to('rad'),# blade angle at exit, relative to tangent
	d_1 = Q_(1, 'in'),				# impeller inlet diameter
	d_2 = Q_(2.1, 'in'),			# impeller outlet diameter
	b_2 = Q_(0.07, 'in'),			# impeller outlet height
	thk2 = Q_(0.04, 'in'),			# blade thickness at exit
	d_hub=Q_(0.6, 'in'),			# hub diameter
	
	Z_e_factor = 1.5,				# axial length extension factor, modifying Gulich formula
	shroud_fullness=(1, 0.5),
	hub_fullness=(0.7, 0.7),
	shroud_outlet_angle=Q_(4, 'deg')	# shroud contour angle from radial at outlet
)

lox_imp = Impeller(lox_geometry, design_point)