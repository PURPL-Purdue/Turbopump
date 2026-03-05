import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

data_folder = "./"

# load data
mdot_range = np.loadtxt(data_folder + "turbine_off_design_mdot_range.csv", delimiter=",")
rpm_range = np.loadtxt(data_folder + "turbine_off_design_rpm_range.csv", delimiter=",")
HP = np.loadtxt(data_folder + "turbine_hp_surface.csv", delimiter=",")

MDOT, RPM = np.meshgrid(mdot_range, rpm_range)

# hp surface
fig = plt.figure()
ax = fig.add_subplot(111, projection="3d")

ax.plot_surface(RPM, MDOT, HP)

ax.set_xlabel("RPM")
ax.set_ylabel("Mass Flow Rate (lbm/s)")
ax.set_zlabel("Horsepower (HP)")
ax.set_title("Turbine Performance Surface")

plt.show()

# hp contour
plt.figure()

contour = plt.contourf(RPM, MDOT, HP, levels=30)
plt.colorbar(contour)

plt.xlabel("RPM")
plt.ylabel("Mass Flow Rate (lbm/s)")
plt.title("Turbine Horsepower Map")

plt.show()

