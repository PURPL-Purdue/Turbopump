import numpy as np
from pint import Quantity as Q_

torque = Q_(8.0, 'N*m')

area = np.pi * Q_(49.43583207, 'mm^2')
r_o = Q_(5.08099784, 'mm')
r_i = Q_(3.17500000, 'mm')
r_mean = (r_o + r_i) / 2

friction = torque / r_mean
coeff = 0.5
normal_force = friction / coeff

axial_stress = (normal_force / area).to('MPa')

print('Axial stress: ', axial_stress)