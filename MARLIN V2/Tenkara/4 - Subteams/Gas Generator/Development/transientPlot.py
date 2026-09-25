import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.widgets import Slider, Button
from scipy.integrate import solve_ivp
from scipy.interpolate import interp1d, PchipInterpolator

# ---------------------------------------------------------
# Establish Constants
# ---------------------------------------------------------
hptowatts = 745.7        # unit conversion
lbmtokg = 0.453592       # unit conversion
I_total = 0.000373116    # kg-m^2
dragTorque = 0           # Nm - ASSUMED
N0 = 0                   # RPM
N_switch = 27175         # RPM
t_grid = np.linspace(0, 3, 200)  # s
t_span = (0, 3)
N_span = (0, 40000)

mdotRead = np.linspace(0.1678, 1.2838, 100) * lbmtokg  # kg/s
RPMRead = np.linspace(25000, 75000, 50)               # RPM
RPMExpansion = np.linspace(0, 25000, 50)             # RPM

# ---------------------------------------------------------
# Power Surface Processing
# ---------------------------------------------------------
script_dir = os.path.dirname(os.path.abspath(__file__))
csv_path = os.path.join(script_dir, "hp_surface_n2.csv")

pmatrix_raw = pd.read_csv(csv_path, header=None).values * hptowatts
border = pmatrix_raw[0, :]
slopes = border / np.max(RPMExpansion)
pExt = np.outer(RPMExpansion, slopes)

pmatrix = np.vstack([pExt, pmatrix_raw])
RPMRead_full = np.concatenate([RPMExpansion, RPMRead])

# Create 2D Meshgrids and Transpose to match MATLAB layout
RPM_mesh, mdot_mesh = np.meshgrid(RPMRead_full, mdotRead)
RPM = RPM_mesh.T
mdot = mdot_mesh.T

# ---------------------------------------------------------
# Power Surface Polynomial Fit
# ---------------------------------------------------------
x = RPM.flatten() / 1000.0  # kRPM
y = mdot.flatten()
z = pmatrix.flatten()

# Construct Design Matrix A (3rd-degree 2D polynomial surface fit)
A = np.column_stack([
    np.ones_like(x),
    x, y,
    x**2, x*y, y**2,
    x**3, (x**2)*y, x*(y**2), y**3
])

# Least-squares fit
c, _, _, _ = np.linalg.lstsq(A, z, rcond=None)

def p_fit(rpm_val, m_val):
    x_val = np.asarray(rpm_val) / 1000.0
    y_val = np.asarray(m_val)
    return (c[0] + c[1]*x_val + c[2]*y_val +
            c[3]*(x_val**2) + c[4]*x_val*y_val + c[5]*(y_val**2) +
            c[6]*(x_val**3) + c[7]*(x_val**2)*y_val + c[8]*x_val*(y_val**2) + c[9]*(y_val**3))

Z_fit = p_fit(RPM, mdot)

# ---------------------------------------------------------
# Solve Integral Equations for Transients (ODE45 equivalent)
# ---------------------------------------------------------
num_sims = len(mdotRead)
M_grid, T_grid = np.meshgrid(mdotRead, t_grid)
RPM_grid = np.zeros_like(T_grid)
Nswitch_grid = np.full_like(RPM_grid, N_switch)

time_results = []
RPM_results = []

def make_dNdt(m_val):
    def dNdt(t, N):
        N_curr = N[0] if isinstance(N, np.ndarray) else N
        denom = max(N_curr * (np.pi / 30.0), 1.0)
        power = p_fit(max(N_curr, 1.0), m_val)
        return [(30.0 / np.pi) * (((power / denom) - dragTorque) / I_total)]
    return dNdt

for i in range(num_sims):
    mdot_i = mdotRead[i]
    sol = solve_ivp(make_dNdt(mdot_i), t_span, [N0], method='RK45', rtol=1e-3, atol=1e-3)
    time_results.append(sol.t)
    RPM_results.append(sol.y[0])

# Project ODE results onto structured grid
for i in range(num_sims):
    t_vec = time_results[i]
    N_vec = RPM_results[i]
    
    # Ensure monotonic sequence for interpolation
    _, unique_idx = np.unique(t_vec, return_index=True)
    interp_func = PchipInterpolator(t_vec[unique_idx], N_vec[unique_idx], extrapolate=True)
    RPM_grid[:, i] = interp_func(t_grid)

# ---------------------------------------------------------
# Display all Plots
# ---------------------------------------------------------
fig = plt.figure(figsize=(11, 7))

# Create tab navigation buttons at the top
ax_tab1 = plt.axes([0.15, 0.92, 0.22, 0.05])
ax_tab2 = plt.axes([0.39, 0.92, 0.22, 0.05])
ax_tab3 = plt.axes([0.63, 0.92, 0.22, 0.05])

btn_tab1 = Button(ax_tab1, 'Power Surface')
btn_tab2 = Button(ax_tab2, 'Transient Surface')
btn_tab3 = Button(ax_tab3, 'Transient Plot')

# Plot Axes setups
ax1 = fig.add_subplot(111, projection='3d')
ax2 = fig.add_subplot(111, projection='3d')
ax3 = fig.add_subplot(111)

# Adjust axes position so 3D titles and labels don't clip top/bottom
for ax in [ax1, ax2]:
    ax.set_position([0.10, 0.12, 0.76, 0.72])
    ax.view_init(elev=30, azim=-120)  # Reorient 90 degrees CCW from default (-30 -> -120)

ax3.set_position([0.12, 0.25, 0.78, 0.62])

# --- TAB 1: Power Surface ---
ax1.plot_surface(RPM, mdot, pmatrix, cmap='viridis', alpha=0.8)
ax1.plot_wireframe(RPM, mdot, Z_fit, color='k', rstride=5, cstride=5, linewidth=0.5, label='Fit Mesh')
ax1.set_xlabel("RPM", labelpad=10)
ax1.set_ylabel("mdot (kg/s)", labelpad=10)
ax1.set_zlabel("Power (W)", labelpad=10)
ax1.set_title("Power Surface Fit Reconstruction", pad=15)

# --- TAB 2: Transient Surface ---
surf2 = ax2.plot_surface(M_grid, T_grid, RPM_grid, cmap='jet', alpha=0.85, vmin=0, vmax=40000)
cbar = fig.colorbar(surf2, ax=ax2, label='RPM', shrink=0.55, pad=0.1)
cbar.mappable.set_clim(0, 40000)

ax2.plot_surface(M_grid, T_grid, Nswitch_grid, color='red', alpha=0.3)
ax2.set_xlabel('mdot (kg/s)', labelpad=10)
ax2.set_ylabel('Time (s)', labelpad=10)
ax2.set_zlabel('RPM', labelpad=10)
ax2.set_zlim(N_span)
ax2.set_title("Transient Spinup Profile at Various Mass Flows", pad=15)

# --- TAB 3: Interactive Transient Plot ---
ax3.grid(True)
ax3.set_xlabel("Time (s)")
ax3.set_ylabel("RPM")
ax3.set_ylim(N_span)
ax3.set_xlim(t_span)

mdotinitial = 0.1244
line_handle, = ax3.plot(time_results[0], RPM_results[0], 'b-', linewidth=2, label="Transient Plot")
ax3.plot(t_span, [N_switch, N_switch], 'r--', linewidth=1.5, label="Switchoff RPM")
ax3.legend(loc="lower right")

time_text = ax3.text(0.04, 0.88, '', transform=ax3.transAxes, fontsize=11, 
                     bbox=dict(boxstyle='square', facecolor='white', edgecolor='black'))

# Slider control placement
ax_slider = plt.axes([0.25, 0.08, 0.50, 0.04])
sld = Slider(ax_slider, 'mdot (kg/s): ', min(mdotRead), max(mdotRead), valinit=mdotinitial, valfmt='%.4f')

# Function to handle tab switches
def show_tab(tab_num):
    ax1.set_visible(tab_num == 1)
    ax2.set_visible(tab_num == 2)
    cbar.ax.set_visible(tab_num == 2)
    ax3.set_visible(tab_num == 3)
    ax_slider.set_visible(tab_num == 3)
    fig.canvas.draw_idle()

btn_tab1.on_clicked(lambda event: show_tab(1))
btn_tab2.on_clicked(lambda event: show_tab(2))
btn_tab3.on_clicked(lambda event: show_tab(3))

# Slider update callback
def expand_update(val):
    sol = solve_ivp(make_dNdt(val), t_span, [N0], method='RK45', rtol=1e-3, atol=1e-3)
    t_out = sol.t
    N_out = sol.y[0]

    if np.max(N_out) >= N_switch:
        N_unique, unique_idx = np.unique(N_out, return_index=True)
        t_unique = t_out[unique_idx]
        interp_switch = interp1d(N_unique, t_unique, kind='linear', fill_value='extrapolate')
        t_switch = float(interp_switch(N_switch))
        time_str = f'Switchoff Time: {t_switch:.3f} s'
    else:
        time_str = 'Switchoff Time: N/A (Does not reach RPM)'

    line_handle.set_xdata(t_out)
    line_handle.set_ydata(N_out)
    time_text.set_text(time_str)
    fig.canvas.draw_idle()

sld.on_changed(expand_update)

# Display initial tab
show_tab(1)
expand_update(mdotinitial)

plt.show()