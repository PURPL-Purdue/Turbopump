import numpy as np
from pint import Quantity as Q_
import Pump
import Impeller
from CoolProp import CoolProp as CP
#import CoolProp.CoolProp as CP
import sys
from pathlib import Path

sys.path.append(str(Path(__file__).resolve().parents[2]))	# Pumps folder, so the Inputs package imports when run directly
from Inputs import Inputs

#print(Inputs)

atm_pressure = Q_(1, 'atm').to('Pa')

def LOXPump():

	# LO2 is initially saturated at atmospheric pressure and then later pressurized
	# It can be assumed that the temperature of LO2 after tank pressurization is approximately equal to the saturation temperature before pressurization
	lox_temp = CP.PropsSI('T', 'P', atm_pressure.m, 'Q', 0, 'O2')
	lox_vapor_pressure = Q_(CP.PropsSI('P', 'T', lox_temp, 'Q', 0, 'O2'), 'Pa')

	return Pump.Pump(
		design_point=Pump.PumpDesignPoint(
			feed_pressure=Inputs['LOX']['feed_pressure'],
			discharge_pressure=Inputs['LOX']['discharge_pressure'],
			mass_flow=Inputs['LOX']['mass_flow'],
			N_shaft=Inputs['N_shaft'],
			rho=Q_(1146, 'kg/m^3'),
			inlet_vapor_pressure=lox_vapor_pressure
		),

		impeller_geometry=Impeller.InputGeometry(
			eta_hyd = 0.50,				# hydraulic efficiency at BEP
			Z_blade = 6,					# integer number of blades
			Beta2B = Q_(20,'deg').to('rad'),# blade angle at exit, relative to tangent
			d_1 = Q_(1, 'in'),				# impeller inlet diameter
			d_2 = Q_(2.1, 'in'),			# impeller outlet diameter
			b_2 = Q_(0.07, 'in'),			# impeller outlet height
			thk2 = Q_(0.04, 'in'),			# blade thickness at exit
			d_hub=Q_(0.6, 'in'),			# hub diameter
			
			# Z_e_factor = 1.5,				# axial length extension factor, modifying Gulich formula
			# shroud_fullness=(1, 0.5),
			# hub_fullness=(0.7, 0.7),
			# shroud_outlet_angle=Q_(4, 'deg')	# shroud contour angle from radial at outlet
		)
	)

if __name__ == "__main__":
	pump = LOXPump()
	import Pump1D
	print(Pump1D.GuessInducerlessNPSH_r(pump))
	print(pump.NPSH_a)
	print(Pump1D.GetInducerHeadRequirement(pump))
	#pump.Impeller.PrintDimensionlessParameters()