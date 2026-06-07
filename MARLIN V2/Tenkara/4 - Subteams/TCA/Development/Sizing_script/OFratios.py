import numpy as np
import matplotlib.pyplot as plt
import cea
from pint import UnitRegistry

ureg = UnitRegistry()

# Inputs
pc = 800 * ureg.psi
pe = 12 * ureg.psi
pi_p = [pc / pe]

of_ratio = np.arange(0.5, 5.0, 0.01)

# CEA setup
reac_names = ["C2H5OH(L)", "O2(L)"]
T_reactant = np.array([298.15, 90.17]) * ureg.K

fuel_weights   = np.array([1.0, 0.0])
oxidant_weights = np.array([0.0, 1.0])

reac   = cea.Mixture(reac_names)
prod   = cea.Mixture(reac_names, products_from_reactants=True)
solver = cea.RocketSolver(prod, reactants=reac, transport=True)

Tc_list    = []
cstar_list = []

for of in of_ratio:
    solution = cea.RocketSolution(solver)
    weights  = reac.of_ratio_to_weights(oxidant_weights, fuel_weights, of)
    hc = (
        reac.calc_property(
            cea.ENTHALPY,
            weights,
            T_reactant.to(ureg.K).magnitude
        ) / cea.R
    )
    solver.solve(
        solution,
        weights,
        pc.to(ureg.bar).magnitude,
        pi_p,
        iac=True,
        hc=hc
    )

    Tc_list.append(float(np.asarray(solution.T[0]).flatten()[0]))
    cstar_list.append(float(np.asarray(solution.c_star).flatten()[0]))

# Convert to clean 1D arrays
Tc_array    = np.array(Tc_list)
cstar_array = np.array(cstar_list)

# Print max values
max_cstar_idx = int(np.argmax(cstar_array))
max_Tc_idx    = int(np.argmax(Tc_array))

print(f"Max c*: {cstar_array[max_cstar_idx]:.2f} m/s at O/F = {of_ratio[max_cstar_idx]:.2f}")
print(f"Max Tc: {Tc_array[max_Tc_idx]:.2f} K    at O/F = {of_ratio[max_Tc_idx]:.2f}")

# Plot
fig, ax1 = plt.subplots(figsize=(8, 5))

line1 = ax1.plot(of_ratio, Tc_array,    color="tab:blue",   label="Chamber Temperature (K)")[0]
ax1.set_xlabel("O/F Ratio")
ax1.set_ylabel("Chamber Temperature [K]")
ax1.grid(True)

# Mark max Tc
ax1.axvline(of_ratio[max_Tc_idx], color="tab:blue", linestyle="--", linewidth=0.8, alpha=0.6)

ax2 = ax1.twinx()
line2 = ax2.plot(of_ratio, cstar_array, color="tab:orange", label="c* (m/s)")[0]
ax2.set_ylabel("c* [m/s]")

# Mark max c*
ax2.axvline(of_ratio[max_cstar_idx], color="tab:orange", linestyle="--", linewidth=0.8, alpha=0.6)

# Combined legend
ax1.legend([line1, line2], [line1.get_label(), line2.get_label()], loc="upper right")

plt.title("Chamber Temperature and c* vs O/F Ratio")
plt.tight_layout()
plt.show()