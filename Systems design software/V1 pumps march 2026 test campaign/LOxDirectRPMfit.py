import CoolProp
import numpy as np
import matplotlib.pyplot as plt
from CoolProp.CoolProp import PropsSI
import pandas as pd
import numpy as np
import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
from mpl_toolkits.mplot3d import Axes3D


#Pull pump curves
vol_flow, head = np.genfromtxt('Systems design software\\V1 pumps march 2026 test campaign\\CFTurbo performance curves\\LOx pump N2\\LOx_LN2_25000_55000_combo.csv', delimiter='	', unpack=True, skip_header=1, dtype=float)

#Define the 3D model function
def func(data, a, b, c, d, e, f):
    x, y = data
    return a + b*(x) + c * (y) + d*x**2 + e*y**2 + f*(x)*(y)

# Generate data based on the non-linear model
rpm_25000 = np.ones(99, int) * 25000
rpm_30000 = np.ones(100, int) * 30000
rpm_35000 = np.ones(100, int) * 35000
rpm_40000 = np.ones(100, int) * 40000
rpm_45000 = np.ones(100, int) * 45000
rpm_50000 = np.ones(100, int) * 50000
rpm_55000 = np.ones(100, int) * 55000

#Combinging all rpms into one vector
hstack_1 = np.hstack((rpm_25000, rpm_30000))
hstack_2 = np.hstack((hstack_1, rpm_35000))
hstack_3 = np.hstack((hstack_2, rpm_40000))
hstack_4 = np.hstack((hstack_3, rpm_45000))
hstack_5 = np.hstack((hstack_4, rpm_50000))
hstack_6 = np.hstack((hstack_5, rpm_55000))

rpm_vector = hstack_6

#Fit the model 
popt, _ = curve_fit(func, (vol_flow, head), rpm_vector, p0=[1, 1, 1, 1, 1, 1])

########Plotting slop ###########
x_flat = vol_flow
y_flat = head
z_data = rpm_vector
x_range = np.linspace(min(x_flat), max(x_flat), 30)
y_range = np.linspace(min(y_flat), max(y_flat), 30)
X_grid, Y_grid = np.meshgrid(x_range, y_range)
Z_fit_grid = func((X_grid, Y_grid), *popt)
X, Y = np.meshgrid(vol_flow, head)
fig = plt.figure(figsize=(10, 7))
ax = fig.add_subplot(111, projection='3d')
x_data = vol_flow
y_data = head
ax.scatter(x_data, y_data, z_data, label='Original Data', color='blue', alpha=0.5)
ax.set_xlabel('X axis')
ax.set_ylabel('Y axis')
ax.set_zlabel('Z axis')
ax.set_title('3D Regression Plot of LOx pump curves')
ax.plot_surface(X_grid, Y_grid, Z_fit_grid, color='red', alpha=0.3) # Surface
plt.show()

#Output
a_fit, b_fit, c_fit, d_fit, e_fit, f_fit = popt
print(f"Fitted parameters: a = {a_fit:.2f}, b = {b_fit:.2f}, c = {c_fit:.2f}, d = {d_fit:.2f}, e_fit = {e_fit:.2f}, f_fit = {f_fit:.2f}")

Q_TCA, H_TCA = 86.92, 835.21

prediction = func((Q_TCA, H_TCA), a_fit, b_fit, c_fit, d_fit, e_fit, f_fit)
print(f'We need {prediction:.2f} RPM to meet the LOx setpoint.\n')

A = e_fit
B = c_fit + f_fit*Q_TCA
C = a_fit + b_fit*Q_TCA + d_fit*Q_TCA**2 - 40372

disc = B**2 - 4*A*C

y1 = (-B + np.sqrt(disc)) / (2*A)
y2 = (-B - np.sqrt(disc)) / (2*A)

adjusted_H = np.minimum(y1, y2)

print(f'Now, adjusting for the kero prediction of 40372 RPM, H will be {adjusted_H:.2f} ft to meet the LOx setpoint')