import pandas as pd
import numpy as np
from pint import UnitRegistry
import yaml
from scipy.interpolate import PchipInterpolator
import matplotlib.pyplot as plt

df = pd.read_csv("Outputs/properties_torch.csv")

with open('Inputs/TCA_params.yaml') as f:
    p = yaml.safe_load(f)

ureg = UnitRegistry()

#Interpolated data for Nusselt number as a function of x/D, based on data from "MPINGEMENT OF A CIRCULAR JET WITH AND WITHOUT CROSS FLOW R. J. GOLDSTEIN and A. I. BEHBAHANI
X_D = np.array([-20, -15, -10, -7.5, -5, -3, -1.5, 0,
                1.5, 3, 5, 7.5, 10, 15, 20, 30, 40, 48], dtype=float)
Y = np.array([35, 50, 80, 105, 140, 175, 200, 215,
              195, 165, 135, 105, 80, 52, 40, 20, 13, 7], dtype=float)

_interp = PchipInterpolator(X_D, Y, extrapolate=False)
 

def re108100(x):
    """Interpolated value at x/D (scalar or array). NaN outside [-20, 48]."""
    return _interp(x)

De = p['torch_dimensions']['throat_diameter'] * ureg.inch
ke = df['k [W/m-K]'].iloc[-1] *ureg.watt / (ureg.meter * ureg.kelvin)
Me = 1
rho = df['Density [kg/m^3]'].iloc[-1] * ureg.kilogram / ureg.meter**3
mu = df['Viscosity [Pa*s]'].iloc[-1] * ureg.pascal * ureg.second
Te = df['Temperature [K]'].iloc[-1] * ureg.kelvin
R = df['R [J/kg-K]'].iloc[-1] * ureg.joule / (ureg.kilogram * ureg.kelvin)
gamma = df['gamma'].iloc[-1] 


print(f"torch exit diameter: {De.to(ureg.m)} ")
print(f"torch thermal conductivity: {ke}")
print(f"torch exit Mach number: {Me}")
print(f"torch exit density: {rho} ")
print(f"torch exit viscosity: {mu}")
print(f"torch exit temperature: {Te}")
print(f"torch gas constant: {R}")
print(f"torch specific heat ratio: {gamma}")

ue = Me * np.sqrt(gamma * R * Te) .to(ureg.m/ureg.s)
print(f"torch exit velocity: {ue}")

Re = (rho * ue * De.to(ureg.m) / mu).to(ureg.dimensionless)
print(f"torch exit Reynolds number: {Re}")

#Array of x/D values for plotting
D = np.linspace(0, 20, 20)
Nu = re108100(D)
Tr = 1.1*Te

print(f"torch Nusselt number: {Nu}")

h = (Nu * ke / De.to(ureg.m)).to(ureg.watt / (ureg.meter**2 * ureg.kelvin))
print(f"torch convective heat transfer coefficient: {h}")

r = D*De.to(ureg.m).magnitude

fig, ax = plt.subplots(figsize=(7, 4.5))
ax.plot(r, h.magnitude, "o-", color="tab:red", label="h")
ax.set_xlabel("x [in]")
ax.set_ylabel("h  [W/(m²·K)]")
ax.set_title(f"Torch jet convective coefficient (Re = {Re.magnitude:,.0f})")
ax.grid(True, alpha=0.3)
ax.set_ylim(bottom=0)

fig.savefig("Outputs/torch_h_vs_xD.png", dpi=200)
plt.show()

# ============================================================
# 2D circular h(x,y) field for ANSYS External Data
# ============================================================
R_DISK = 0.0762         # [m] disk radius -- set to your real value
DX     = 0.001          # [m] grid spacing (1 mm resolves the peak at the jet centre)
R_CUT  = 1.05 * R_DISK  # small margin beyond the disk edge

De_m = De.to(ureg.m).magnitude
ke_over_D = (ke / De.to(ureg.m)).to(ureg.watt / (ureg.meter**2 * ureg.kelvin)).magnitude

# Cartesian grid centred on the jet axis (includes the exact point (0, 0))
n = int(np.ceil(R_CUT / DX))
x = np.arange(-n, n + 1) * DX
X_grid, Y_grid = np.meshgrid(x, x)

# Radial distance in jet diameters
R_D_grid = np.sqrt(X_grid**2 + Y_grid**2) / De_m

# Keep only points inside the cut radius
mask = R_D_grid * De_m <= R_CUT
R_D_pts = R_D_grid[mask]

Nu_pts = re108100(R_D_pts)
if np.isnan(Nu_pts).any():
    raise ValueError(
        f"Nu is NaN for r/D up to {R_D_pts.max():.1f}; "
        "correlation only covers r/D in [0, 48]."
    )

h_pts = Nu_pts * ke_over_D   # [W/m2-K]

ansys_data = pd.DataFrame({
    "X [m]": X_grid[mask],
    "Y [m]": Y_grid[mask],
    "h [W/m2-K]": h_pts,
    "Tr [K]": Tr.magnitude,
})

ansys_data.to_csv("Outputs/ANSYS_jet_h.csv", index=False)
print(ansys_data)
print(f"h at centre: {h_pts[np.argmin(R_D_pts)]:.1f} W/m2-K, max: {h_pts.max():.1f}")