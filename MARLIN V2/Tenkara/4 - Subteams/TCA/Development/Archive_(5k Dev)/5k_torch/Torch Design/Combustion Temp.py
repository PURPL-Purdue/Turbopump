# Plot combustion temp vs. OF ratio

from pyfluids import Fluid, FluidsList, Input
import matplotlib.pyplot as plt
from rocketcea.cea_obj import CEA_Obj
import numpy as np

R2K = 0.555556 # rankine to Kelvin

p_chamber = 300 # chamber pressure [psi]
of_array = []
temp_array = []
ox = 'GOX'
fuel = 'GH2'
ispObj = CEA_Obj(oxName=ox, fuelName=fuel)

for i in range(0, 4000):
    of_array.append(i*0.01)
    temp_array.append(ispObj.get_Tcomb(Pc= p_chamber, MR= i*0.01) * R2K)

a = of_array.index(2)
T = temp_array[a]
print(of_array[a], T)

Cp_chamber, Cp_Throat, x = ispObj.get_HeatCapacities(Pc= p_chamber, MR=2, eps=1)
print(Cp_chamber, Cp_Throat)

print(ispObj.get_full_cea_output(Pc= p_chamber, MR=2.2, eps=1))

fig, ax = plt.subplots()
ax.plot(of_array, temp_array)
plt.xlabel('OF Ratio')
plt.ylabel('Combustion Temperature [K]')
plt.title('Combustion Temperature vs. OF Ratio: ' + ox + '/' + fuel)
plt.grid(True)
plt.show()

