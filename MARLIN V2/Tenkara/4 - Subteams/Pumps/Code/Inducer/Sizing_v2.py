import numpy as np
lambda_w = 0.3
lambda_c = 1.1
Q = 0.001771
n = 35000
g = 9.81
d_1 = 0.03
d_h = d_1 / 2
cm_1 = (4*Q)/(np.pi*(d_1**2))
u1 = (np.pi * n * d_1)/60
w_1 = np.sqrt(cm_1**2 + u1**2)
Npsh = lambda_c*(cm_1**2/(2*g)) + lambda_w*(w_1**2/(2*g))
print(Npsh)