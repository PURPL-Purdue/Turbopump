import numpy as np
from pint import Quantity as Q_
import Pump
import Impeller
import Inducer
from CoolProp import CoolProp as CP
import MeridionalDesign
import sys
from pathlib import Path

sys.path.append(str(Path(__file__).resolve().parents[1]))	# Pumps folder, so the Inputs package imports when run directly
from Inputs import Inputs

#print(Inputs)

atm_pressure = Q_(1, 'atm').to('Pa')

basic_impeller_geometry = {
	'LOX': Impeller.InputGeometry(
		eta_hyd = 0.50,					# hydraulic efficiency at BEP
		Z_blade = 6,					# integer number of blades
		Beta2B = Q_(20,'deg').to('rad'),# blade angle at exit, relative to tangent
		d_1 = Q_(1, 'in'),				# impeller inlet diameter
		d_2 = Q_(2.1, 'in'),			# impeller outlet diameter
		b_2 = Q_(0.07, 'in'),			# impeller outlet height
		blade_thick = Q_(0.04, 'in'),			# blade thickness at exit
		d_hub=Q_(0.6, 'in'),			# hub diameter
	),

	'IPA': Impeller.InputGeometry(
		eta_hyd = 0.50,					# hydraulic efficiency at BEP
		Z_blade = 6,					# integer number of blades
		Beta2B = Q_(25,'deg').to('rad'),# blade angle at exit, relative to tangent
		d_1 = Q_(1.2, 'in'),			# impeller inlet diameter
		d_2 = Q_(2.8, 'in'),			# impeller outlet diameter
		b_2 = Q_(0.08, 'in'),			# impeller outlet height
		blade_thick = Q_(0.04, 'in'),			# blade thickness at exit
		d_hub=Q_(0.7, 'in'),			# hub diameter
	)
}

meridional_geometry = {
	'LOX': MeridionalDesign.MeridionalGeometry(
		Z_e_factor = 1.2,				# axial length extension factor, modifying Gulich formula
		shroud_fullness=(1, 0.5),
		hub_fullness=(0.7, 0.7),
		shroud_outlet_angle=Q_(4, 'deg')	# shroud contour angle from radial at outlet
	),
	'IPA': MeridionalDesign.MeridionalGeometry(
		Z_e_factor = 1.5,
		shroud_fullness=(1, 0.5),
		hub_fullness=(0.7, 0.7),
		shroud_outlet_angle=Q_(5, 'deg')
	)
}

# LO2 is initially saturated at atmospheric pressure and then later pressurized
# It can be assumed that the temperature of LO2 after tank pressurization is approximately equal to the saturation temperature before pressurization
lox_temp = CP.PropsSI('T', 'P', atm_pressure.m, 'Q', 0, 'O2')
lox_vapor_pressure = Q_(CP.PropsSI('P', 'T', lox_temp, 'Q', 0, 'O2'), 'Pa')

fluid_props = {
	'LOX': dict(
		rho = Q_(1146, 'kg/m^3'),
		vapor_pressure = lox_vapor_pressure
	),

	'IPA': dict(
		rho = Q_(786, 'kg/m^3'),
		vapor_pressure = Q_(0.64, 'psi') # TODO make this procedural based on temperature
	)
}

color = dict(LOX='blue', IPA='green')

def CreatePump(propellant: str='LOX'):
	return Pump.Pump(
		design_point=Pump.PumpDesignPoint(
			feed_pressure=Inputs[propellant]['feed_pressure'],
			discharge_pressure=Inputs[propellant]['discharge_pressure'],
			mass_flow=Inputs[propellant]['mass_flow'],
			N_shaft=Inputs['N_shaft'],
			rho=fluid_props[propellant]['rho'],
			inlet_vapor_pressure=fluid_props[propellant]['vapor_pressure']
		),

		impeller_geometry=basic_impeller_geometry[propellant]
	)

def MakeDesign(propellant: str='LOX'):
	pump = CreatePump(propellant)

	pump.Impeller.PlotOutletVelocities().show()

	import Pump1D
	Pump1D.PlotHQ(pump)#.show()

	inducer_head = Pump1D.GetInducerHeadRequirement(pump)
	print(f'{propellant} inducer desired head:', inducer_head)

	pump.AddInducer(inducer_head)
	pump.Impeller.DP.H_rise = pump.DesignHeadRise - inducer_head
	from InducerRadialEquilibrium import InducerRadialEquilibrium
	inducer_eq = InducerRadialEquilibrium(pump.Inducer)

	#print(inducer_eq.GetOutletAngularMomentum())

	pump.Impeller.PrintDimensionlessParameters()
	pump_meridional = MeridionalDesign.MeridionalDesign(pump, meridional_geometry[propellant])
	#pump_meridional.PotentialFlow.Plot().show()
	pump_meridional.PlotMeridional(color=color[propellant])#.show()

	from ImpellerBladeDesign import ImpellerBladeDesign

	imp_blade_design = ImpellerBladeDesign(pump_meridional, inducer_eq)
	imp_blade_design.PlotInletTriangles().show()
	return pump

if __name__ == "__main__":
	MakeDesign('LOX')
	MakeDesign('IPA')