import matplotlib.pyplot as plt
import numpy as np

thermal_expansion_coefficients = {
    'PTFE': 67e-6,
    'Inconel 718': 7.2e-6,
    'Monel': 7.8e-6,
    'Kel-F': 33e-6,
    '316 Stainless': 8.8e-6,
    '304 Stainless': 9.6e-6,
    'Brass': 10.4e-6,
    'Copper': 9.4e-6,
    'Bronze': 9.5e-6,
    '7075 Aluminum': 13.0e-6
}   # all units in number / degree F


#inside
mat1 = "316 Stainless"
rad1 = .7 #radius at room temp (inches)S
thermCo1 = thermal_expansion_coefficients.get(mat1)

#outside 
mat2 = 'Kel-F'
rad2 = .701 #radius at room temp (inches)
thermCo2 = thermal_expansion_coefficients.get(mat2)

roomtemp = 70
temp1 = -300
temp2 = 200

#length calcs
lengths1 = []
lengths2 = []
clearances = []
temps = list(range(temp1, temp2))
for i in temps:
    lgth1 = rad1*thermCo1*(i-roomtemp) + rad1
    lengths1.append(lgth1)
    lgth2 = rad2*thermCo2*(i-roomtemp) + rad2
    lengths2.append(lgth2)
    clearances.append(lgth2-lgth1)

#plots
fig, (ax1, ax2) = plt.subplots(2,1)

ax1.plot(temps, lengths1, color='b', label = "Inner: "+ mat1)
ax1.plot(temps, lengths2, color='r', label = "Outer: " + mat2)
ax1.set_xlabel('Temp (F)')
ax1.set_ylabel('Diameter (in)')
ax1.set_title('Diameters at through temperature range')
ax1.legend() 
ax1.grid(True)


ax2.plot(temps, clearances, color='k')
ax2.set_title('overlap at different times (negative is overlap)')
ax2.set_xlabel('Temp (F)')
ax2.set_ylabel('Clearances (in)')
ax2.grid(True)
plt.tight_layout()
plt.show()

