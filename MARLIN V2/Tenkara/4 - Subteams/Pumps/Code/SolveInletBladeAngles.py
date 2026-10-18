import numpy as np
from Impeller import Impeller
from pint import Quantity as Q_

# inlet blockage factor
def _tau_1(impeller: Impeller, Beta1B: Q_[np.ndarray]) -> np.ndarray:
	# Gulich Table 0.1
	"""
	Returns the inlet blockage factor tau_1, where A_1' = A_1 * tau_1
	Beta1B is the metallic blade angle at the inlet
	lambda_LE is the angle measured in the meridional view between the shroud contour and the leading edge projection
	"""
	# TODO the diameter or radius at which the circumference is evaluated should vary along with Beta1B
	tau_1: Q_ = 1 / (
		1 - impeller.Z_blade * impeller.thk2 / (np.pi * impeller.d_hub * np.sin(Beta1B.to('rad').m) * np.sin(impeller.lambdaLE.to('rad').m))
		)
	return tau_1.to('dimensionless').m

from _LOX_impeller import lox_imp

print(_tau_1(lox_imp, Q_(30, 'deg')))