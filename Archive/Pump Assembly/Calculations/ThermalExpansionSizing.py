#this code is for finding sizing parts on the shaft
# and their wear rings at differnent temps
# PTFE might not work

# %% setup
# You should get the jupyter notebooks extension for VS code, this works a little better with that
import matplotlib.pyplot as plt
import numpy as np

thermal_expansion_coefficients = {
    'PTFE': 67e-6,
    'Inconel 718': 7.2e-6,
    'Inconel 625': 7.1e-6,
    'Monel': 7.8e-6,
    'Kel-F': 33e-6,
    '316 Stainless': 8.8e-6,
    '304 Stainless': 9.6e-6,
    'Brass': 10.4e-6,
    'Copper': 9.4e-6,
    'Bronze': 9.5e-6,
    '7075 Aluminum': 13.0e-6,
    'Mild Steel' : 6.5e-6
}   # all units in number / degree F

def thermalSize(deltaT, length, CFT):
    return length + length*deltaT*CFT
# %% calculate clearance at different temps with set dimensions

#inside
mat1 = "316 Stainless"
rad1 = .4 #radius at room temp (inches)S
thermCo1 = thermal_expansion_coefficients.get(mat1)

#outside 
mat2 = 'Copper'
rad2 = .402 #radius at room temp (inches)
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

# %% find dimensions for disired clearance

# set either the od or ID to 0 if trying to find the proper size
shaftOD = .7 #set this to zero if trying to find dia
shaft_mat = 'Inconel 718'
shaft_CFT = thermal_expansion_coefficients.get(shaft_mat)

featureID = 0 #set this to zero if trying to find
feature_mat = 'Brass'
feature_CFT = thermal_expansion_coefficients.get(feature_mat)

clearance = .000 #negative for interferacne fit
operatingTemp = -300 #degrees F
roomtemp = 70

if shaftOD != 0 and featureID != 0:
    print('YOU MUST set either featureID or shaftOD to 0')
    
elif shaftOD != 0:
    atTempOD = shaftOD + shaft_CFT*shaftOD*(operatingTemp-roomtemp)
    atTempNewDia = atTempOD + clearance
    finalDia = atTempNewDia*-(operatingTemp-roomtemp)*feature_CFT+atTempNewDia
else:
    atTempID = featureID + featureID*(operatingTemp-roomtemp)*feature_CFT
    atTempNewDia =  atTempID - clearance
    finalDia = atTempNewDia*-(operatingTemp-roomtemp)*shaft_CFT+atTempNewDia

print(finalDia)
# %% one dimension calc
length = 2
mat = 'PTFE'
cft = thermal_expansion_coefficients.get(mat)
roomtemp = 70
operatingTemp = -300
deltaT = operatingTemp-700
newLen = thermalSize(deltaT, length, cft)
print(newLen)
# %%
