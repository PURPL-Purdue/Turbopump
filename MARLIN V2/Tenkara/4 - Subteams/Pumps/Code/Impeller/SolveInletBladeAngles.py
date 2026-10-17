import numpy as np
from Impeller import Impeller
from pint import Quantity as Q_

# inlet blockage factor
def _tau_1(impeller: Impeller, Beta1B: Q_[np.ndarray], lambda_LE: Q_[float]) -> np.ndarray:
	# Gulich Table 0.1
	"""
	Returns the inlet blockage factor tau_1, where A_1' = A_1 * tau_1
	Beta1B is the metallic blade angle at the inlet
	lambda_LE is the angle measured in the meridional view between the shroud contour and the leading edge projection
	"""
	tau_1: Q_ = 1 / (
		1 - impeller.Z_blade * impeller.thk2 / (np.pi * impeller.d_1 * np.sin(Beta1B.to('rad').m) * np.sin(lambda_LE.to('rad').m))
		)
	return tau_1.to('dimensionless').m

