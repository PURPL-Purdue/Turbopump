import numpy as np
from Impeller import Impeller
from pint import Quantity as Q_

# inlet blockage factor
def _tau_1(impeller: Impeller, Beta1B: Q_[np.ndarray]) -> np.ndarray:
	# Gulich Table 0.1 and Milo's inlet sizing method
	# TODO the diameter or radius at which the circumference is evaluated should vary along with Beta1B

	rho_1 = .8 # initial guess
	Beta1B_prime = np.atan(np.tan(Beta1B) / rho_1)
	for _ in range(50):
		rho_1: Q_ = 1 - impeller.Z_blade * impeller.blade_thick / (np.pi * impeller.d_hub * np.sin(Beta1B_prime.to('rad').m) * np.sin(impeller.lambdaLE.to('rad').m))
		Beta1B_prime = np.atan(np.tan(Beta1B_prime) / rho_1)
		print(Beta1B_prime.to('deg'))
		print(rho_1)
		
		c_m1 = (impeller.DP.Q / (impeller.Area1 * rho_1)).to('m/s')
		#print(c_m1)
		c_u = Q_(0, 'm/s')
		u_1 = impeller.DP.N_shaft * impeller.d_hub/2
		beta = Q_(np.atan((c_m1 / (u_1 - c_u)).to('dimensionless').m), 'rad')
		print(beta.to('deg'))
	
	return rho_1, Beta1B_prime

from _LOX_impeller import lox_imp

hub_triangle = lox_imp.GetVelocities(1, 0)
beta1 = hub_triangle.beta
_tau_1(lox_imp, beta1)
#lox_imp.SweepVelocities([0]).show()