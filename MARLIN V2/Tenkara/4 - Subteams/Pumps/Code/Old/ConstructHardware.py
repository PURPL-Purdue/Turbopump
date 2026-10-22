import numpy as np
import plotly.graph_objects as go
from pint import Quantity as Q_
from Inducer import Inducer, InducerDesignPoint, InducerInletDesignConstraints, InducerGeometry
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

g = Q_(9.81, 'm/s^2')					# gravitational acceleration

# design mass flow
opt_m_dot_LOX = Q_(
	dp['LOX']['mass_flow'], 'kg/s'
)
opt_m_dot_IPA = Q_(
	dp['IPA']['mass_flow'], 'kg/s'
)

rho_IPA = Q_(786, 'kg/m^3')		# density of IPA at room temp, 10 bar ish
rho_LOX = Q_(1141, 'kg/m^3')	# density of LOX at saturation, atm pressure (lol? should really be the density at tank pressure)

feed_pressure = Q_(150, 'psi')			# inlet feed pressure, from tank pressure
discharge_pressure_LOX = Q_(
	dp['LOX']['discharge_pressure'], 'bar')	# ideal discharge pressure
discharge_pressure_IPA = Q_(
	dp['IPA']['discharge_pressure'], 'bar')

shaft_speed = Q_(dp['N_shaft'], 'rpm')	# common shaft

impeller_design_point_LOX: DesignPoint = DesignPoint(
	Q = opt_m_dot_LOX/rho_LOX,					# volumetric flow rate at BEP
	H = ((discharge_pressure_LOX - feed_pressure) / (g * rho_LOX)).to('m'),	# developed head at BEP
	N_shaft = shaft_speed,
	n_hyd_BEP = 0.5,
	rho=rho_LOX
)

impeller_design_point_IPA: DesignPoint = DesignPoint(
	Q = opt_m_dot_IPA/rho_IPA,
	H = ((discharge_pressure_IPA - feed_pressure) / (g * rho_IPA)).to('m'),
	N_shaft = shaft_speed,
	n_hyd_BEP = 0.5,
	rho = rho_IPA,
)

inducer_design_point_LOX = InducerDesignPoint(
	Q=Q_(1.69, 'L/s'), N=Q_(35000, 'rpm'),
	NPSH_a=Q_(105.46, 'm'), NPSH_r=Q_(150.0, 'm')
)

inducer_design_point_IPA = InducerDesignPoint(
	Q=Q_(2.739, 'L/s'), N=Q_(35000, 'rpm'),
	NPSH_a=Q_(105.46, 'm'), NPSH_r=Q_(150.0, 'm')
)

inducer_inlet_constraints = InducerInletDesignConstraints(
    hub_tip_ratio=0.5, sigma=0.03, K=1.1, psi_1=0.08, axial_length_ratio=1.0,
)

inducer_geometry = InducerGeometry(
	Z_blade=3, thk=Q_(1, 'mm'),
)

impeller_geometry_LOX: InputGeometry = InputGeometry(
	Z_blade = 6,
	Beta2B = Q_(20,'deg').to('rad'),# blade angle at exit, relative to tangent
	d_1 = Q_(1, 'in'),				# impeller inlet diameter
	d_2 = Q_(2.1, 'in'),			# impeller outlet diameter
	b_2 = Q_(0.07, 'in'),			# impeller outlet height
	blade_thick = Q_(0.04, 'in'),			# blade thickness at exit
	d_hub=Q_(0.6, 'in'),			# hub diameter
	
	Z_e_factor = 1.5,				# axial length extension factor, modifying Gulich formula
	shroud_fullness=(1, 0.5),
	hub_fullness=(0.7, 0.7),
	shroud_outlet_angle=Q_(4, 'deg')	# shroud contour angle from radial at outlet
)

impeller_geometry_IPA: InputGeometry = InputGeometry(
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

def ImpellerLOX() -> Impeller:
	return Impeller(impeller_geometry_LOX, impeller_design_point_LOX)

def ImpellerIPA() -> Impeller:
	return Impeller(impeller_geometry_IPA, impeller_design_point_IPA)

def InducerLOX() -> Inducer:
	return Inducer(inducer_design_point_LOX, inducer_inlet_constraints, inducer_geometry)

def InducerIPA() -> Inducer:
	return Inducer(inducer_design_point_IPA, inducer_inlet_constraints, inducer_geometry)

