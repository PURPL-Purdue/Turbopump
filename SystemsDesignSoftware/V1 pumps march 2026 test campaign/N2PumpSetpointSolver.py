import numpy as np
import math
import matplotlib.pyplot as plt

minDiff = 1000
in3TOgal = 0.004329
minToSec = 60
g = 32.174 # standard gravity in ft/sec^2
Cd = 0.9 # orifice coefficient

fig, axes = plt.subplots(1, 2, figsize=(14, 6))
fig.suptitle("Pump Performance: Orifice Equation vs CFTurbo", fontsize=14, fontweight='bold')

## WATER IN KERO PUMP ##

density = 62.4 # density of water in lbm/ft3
A = 3.14159 * ((0.301/12)/2) ** 2 # area of outlet orifice (ft^2)

kero_vol_flow, kero_head = np.genfromtxt('SystemsDesignSoftware\\V1 pumps march 2026 test campaign\\CFTurbo performance curves\\Kero pump water\\Kero_water_25000.txt', delimiter='	', unpack=True, skip_header=1, dtype=float)

kero_calc_flows = []

for H in kero_head:
    deltaP = (density / 12**3) * (g * 12) * (H * 12)
    mdot = Cd * (A * 12**2) * math.sqrt(2 * (density / 12**3) * deltaP)
    calc_vol_flow = minToSec * in3TOgal * mdot / (density / 12**3)
    kero_calc_flows.append(calc_vol_flow)

    table_vol_flow = kero_vol_flow[np.where(kero_head == H)]
    diff = abs(calc_vol_flow - table_vol_flow)
    if diff < minDiff:
        minDiff = diff
        reqHead = H
        reqQ = table_vol_flow
        outletPressure = (reqHead * 12) * (density / (12**3)) * (g / 12)

print("----For Water in Kero Pump @ 25000 rpm----")
print("Optimizer Accuracy =", round(minDiff[0], 2), "gal/min")
print("H =", round(reqHead, 2), "ft")
print("Q =", round(reqQ[0], 2), "gal/min")
print("Outlet Pressure:", round(outletPressure, 2), "psi\n")

ax1 = axes[0]
ax1.plot(kero_head, kero_calc_flows, 'b-o', markersize=4, label='Orifice Equation')
ax1.plot(kero_head, kero_vol_flow, 'r-s', markersize=4, label='CFTurbo')
ax1.axvline(reqHead, color='gray', linestyle='--', linewidth=1, label=f'Best match H = {round(reqHead, 1)} ft')
ax1.set_xlabel("Head (ft)", fontsize=11)
ax1.set_ylabel("Volumetric Flow Rate (gal/min)", fontsize=11)
ax1.set_title("Water in Kero Pump @ 25,000 rpm", fontsize=12)
ax1.legend()
ax1.grid(True, alpha=0.3)

## WATER IN LOX PUMP ##

minDiff = 1000

density = 62.4 # density of water in lbm/ft3
A = 3.14159 * ((0.289/12)/2) ** 2 # area of outlet orifice (ft^2)

lox_vol_flow, lox_head = np.genfromtxt('SystemsDesignSoftware\\V1 pumps march 2026 test campaign\\CFTurbo performance curves\\LOx pump water\\LOx_water_25000.txt', delimiter='	', unpack=True, skip_header=1, dtype=float)

coeffs = np.polyfit(lox_head, lox_vol_flow, 2)

print(coeffs)

i = lox_head[0]
print(i)
while(i < 700):
    i += 1
    lox_head = np.append(i,lox_head)
    lox_vol_flow = np.append((coeffs[0]*i**2 + coeffs[1]*i + coeffs[2]),lox_vol_flow)


lox_calc_flows = []

for H in lox_head:
    deltaP = (density / 12**3) * (g * 12) * (H * 12)
    mdot = Cd * (A * 12**2) * math.sqrt(2 * (density / 12**3) * deltaP)
    calc_vol_flow = minToSec * in3TOgal * mdot / (density / 12**3)
    lox_calc_flows.append(calc_vol_flow)

    table_vol_flow = lox_vol_flow[np.where(lox_head == H)]
    diff = abs(calc_vol_flow - table_vol_flow)
    if diff < minDiff:
        minDiff = diff
        reqHead = H
        reqQ = table_vol_flow
        outletPressure = (reqHead * 12) * (density / (12**3)) * (g / 12)

print("----For Water in LOx Pump @ 25000 rpm----")
print("Optimizer Accuracy =", round(minDiff[0], 2), "gal/min")
print("H =", round(reqHead, 2), "ft")
print("Q =", round(reqQ[0], 2), "gal/min")
print("Outlet Pressure:", round(outletPressure, 2), "psi\n")

ax2 = axes[1]
ax2.plot(lox_head, lox_calc_flows, 'b-o', markersize=4, label='Orifice Equation')
ax2.plot(lox_head, lox_vol_flow, 'r-s', markersize=4, label='CFTurbo')
ax2.axvline(reqHead, color='gray', linestyle='--', linewidth=1, label=f'Best match H = {round(reqHead, 1)} ft')
ax2.set_xlabel("Head (ft)", fontsize=11)
ax2.set_ylabel("Volumetric Flow Rate (gal/min)", fontsize=11)
ax2.set_title("Water in LOx Pump @ 25,000 rpm", fontsize=12)
ax2.legend()
ax2.grid(True, alpha=0.3)

plt.tight_layout()
plt.show()