# conda activate ct-env
import cantera as ct
import numpy as np

gas1 = ct.Solution('gri30.yaml')

gas1()

gas1.TP = 1200, 101325

gas1()

#Need 2 properties -> all equivalent initializations
gas1.TP = 1200, 101325          # temperature, pressure
gas1.TD = 1200, 0.020473        # temperature, density
gas1.HP = 1.3295e7, 101325      # specific enthalpy, pressure
gas1.UV = 8.3457e6, 1/0.020473  # specific internal energy, specific volume
gas1.SP = 85222, 101325         # specific entropy, pressure
gas1.SV = 85222, 1/0.020473     # specific entropy, specific volume

# read off property
print(gas1.T)

#mole fractions X, mass fractions Y
gas1.X = 'CH4:1, O2:2, N2:7.52'

# can use dictionary as well
phi = 0.8
gas1.X = {'CH4':1, 'O2':2/phi, 'N2':2*3.76/phi}


# control what properties are held constant by passing None to property setter
gas1.SV = None, 2.1
gas1.TPX = None, None, 'CH4:1.0, O2:0.5'
