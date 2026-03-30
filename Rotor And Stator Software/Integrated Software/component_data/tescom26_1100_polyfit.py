import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import json

D_outlet = 0.25 # inches
A_outlet = np.pi * (D_outlet / 2) ** 2 # in^2

data_file = "Tescom26_1100_n2.csv"
data = pd.read_csv(data_file)

print(data.head())
headers = data.columns
print(headers)

flow_rate_scfm = data[headers[0]].values
static_pressure_gauge = data[headers[1]].values

poly_deg = 6
coefficients = np.polyfit(flow_rate_scfm, static_pressure_gauge, poly_deg)
poly_fit = np.poly1d(coefficients)

fig, ax = plt.subplots(figsize=(10, 6))
ax.scatter(flow_rate_scfm, static_pressure_gauge, label="Data Points")
x_fit = np.linspace(min(flow_rate_scfm), max(flow_rate_scfm), 100)
y_fit = poly_fit(x_fit)
ax.plot(x_fit, y_fit, color="red", label=f"Polynomial Fit (degree {poly_deg})")
ax.set_xlabel("Flow Rate (SCFM)")
ax.set_ylabel("Static Pressure (Gauge)")
ax.set_title("Polynomial Fit to Tescom 26-1100 Data")
ax.legend()
ax.grid()
plt.show()

save_fit = input("Save polynomial coefficients to file? (y/n): ")
if save_fit.lower() == 'y':
    coeff_file = "tescom26_1100_polyfit_coefficients"
    np.savetxt(coeff_file + ".txt", coefficients)
    data_dict = {
        "Coefficients": coefficients.tolist(), 
        "Degree": poly_deg,
        "Flow Rate Range (SCFM)": (min(flow_rate_scfm), max(flow_rate_scfm)),
        "Static Pressure Range (Gauge)": (min(static_pressure_gauge), max(static_pressure_gauge)),
        "Outlet Area (in^2)": A_outlet,
        "Outlet Diameter (in)": D_outlet,
        "R-squared": np.corrcoef(static_pressure_gauge, poly_fit(flow_rate_scfm))[0, 1] ** 2
    }
    with open(coeff_file + ".json", 'w') as f:
        json.dump(data_dict, f, indent=4)
    print(f"Coefficients saved to {coeff_file}(.txt and .json)")