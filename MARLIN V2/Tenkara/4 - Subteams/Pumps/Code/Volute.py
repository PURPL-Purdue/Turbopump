#   TODO: 
#   Integrate this with the impeller sizing (automatically wooooo)
#   Select cross section shapes
#   integrate auto cross section dimension calcs

from pint import Quantity as Q_
from matplotlib import pyplot as plt
import numpy as np
import math

# Constants
g = Q_(9.81, 'm/s**2')

# System inputs
Q_design = Q_(1.78, 'liters/second')
H_design = Q_(294.8, 'meters')

# Impeller inputs
C_u2 =  Q_(34.516, 'm/s')   # tangential velocity component
r2 = Q_(1.1, 'inches')      # impeller outlet radius
C_2 = Q_(68.2, 'm/s')       # Absolute velocity

# Volute inputs
e_sp = Q_(2 * math.pi, 'radians')          # volute wrap angle (2pi for single volutes)
K_v = 0.55                                 # experimental average velocity coefficient from h&h equation 6-70 pg 220
                                           # and from https://ijsea.com/archive/volume8/issue8/IJSEA08081016.pdf 

# Throat area calculations
avg_V = K_v * np.sqrt(2 * g.to('ft/s**2') * H_design.to('feet'))    # use empirical relation from h&h pg 220
A_throat = Q_design / avg_V
throat_diameter = np.sqrt((4 / math.pi) * A_throat)                 # if it was a circle
HH_throatV = Q_design.to('gallons per minute') / (A_throat.to('in**2') * 3.12) # flow velocity at throat outlet according to H&H pg 222

print(f"Throat area: {A_throat.to('cm**2'):.4f}")
print(f"Throat diameter: {throat_diameter.to('cm'):.4f}")
print(f"Average velocity: {avg_V.to('m/s'):.4f}")
print(f"Throat velocity: {HH_throatV.to('m/s'):.4f}")

# From gulich pg 420:
#X_sp = (Q_design / (math.pi * C_u2 * r2)) * ((e_sp / Q_(2 * math.pi, 'radians')).to('dimensionless'))
