import numpy as np

# Inputs
d_1 = 0.03048	# [m] tip diameter
npsh_max = 80	# [m] maximum allowed NPSH required
phi_min = 0.06	# [-] minimum flow coefficient (cm_1 / u1)
lambda_w = 0.03
lambda_c = 1.1
Q = 0.001771	# [m^3/s]
n = 35000		# [rpm]
g = 9.81

u1 = np.pi * n * d_1 / 60

# phi = cm_1 / u1 rises with the hub size, so the smallest hub is the one that gives phi_min
rho_sq = 1 - (4 * Q) / (np.pi * d_1**2 * u1 * phi_min)
if rho_sq < 0:
	print("phi_min is already met with no hub, using d_h = 0")
	rho_sq = 0
rho = np.sqrt(rho_sq)	# hub to tip ratio
d_h = rho * d_1			# hub diameter at the inlet

# Check, same equations as the original code
cm_1 = (4 * Q) / (np.pi * (d_1**2 - d_h**2))
w_1 = np.sqrt(cm_1**2 + u1**2)
Npsh = lambda_c * (cm_1**2 / (2 * g)) + lambda_w * (w_1**2 / (2 * g))
flow_coeff = cm_1 / u1
if Npsh > npsh_max:
	raise ValueError(f"NPSH of {Npsh:.1f} m exceeds npsh_max at this d_1, even with the smallest hub")
print(d_h, rho)
print(Npsh, flow_coeff)