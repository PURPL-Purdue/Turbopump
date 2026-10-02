"""
thing
"""
import numpy as np
import matplotlib.pyplot as plt
from pint import Quantity as Q_

class VelocityTriangle:
	u: Q_[float]
	c: Q_[float]
	c_m: Q_[float]
	c_u: Q_[float]
	w: Q_[float]

	def __init__(self, u: Q_[float], c_m: Q_[float], c_u: Q_[float], w: Q_[float]):
		self.u = u
		self.c_m = c_m
		self.c_u = c_u
		self.w = w
		self.c = np.sqrt(c_m**2 + c_u**2)

	def Plot(self, title='Velocity Triangle', unit='m/s', station=0) -> None:

		plt.figure()
		c_u, c_m, c = self.c_u.to(unit).magnitude, self.c_m.to(unit).magnitude, self.c.to(unit).magnitude
		u, w = self.u.to(unit).magnitude, self.w.to(unit).magnitude
		# C: absolute velocity, origin -> (C_u, C_m)
		plt.plot([0, c_u], [0, c_m], linewidth=2,
			label=rf'$C_{station}$ = {c:.1f} {unit}')

		# U: blade speed, origin -> (U, 0)
		plt.plot([0, u], [0, 0], linewidth=2,
			label=rf'$U_{station}$ = {u:.1f} {unit}')

		# W: relative velocity, (C_u, C_m) -> (U, 0)
		plt.plot([c_u, u], [c_m, 0], linewidth=2,
			label=rf'$W_{station}$ = {w:.1f} {unit}')

		# C_m: meridional leg, (C_u, 0) -> (C_u, C_m)
		plt.plot([c_u, c_u], [0, c_m], linewidth=2,
			label=rf'$C_{{m{station}}}$ = {c_m:.1f} {unit}')

		# Labels
		plt.text(c_u / 2, c_m / 2, rf'$C_{station}$', fontsize=11)
		plt.text(u / 2, -0.04 * c_m, rf'$U_{station}$', fontsize=11, ha='center')
		plt.text((u + c_u) / 2, c_m / 2, rf'$W_{station}$', fontsize=11)
		plt.text(c_u, c_m / 2, rf'$C_{{m{station}}}$', fontsize=11, ha='left', va='center')

		plt.axhline(0, linewidth=0.8)
		plt.axvline(0, linewidth=0.8)

		plt.xlabel(f'Tangential velocity [{unit}]')
		plt.ylabel(f'Meridional velocity [{unit}]')
		plt.title(title)
		plt.axis('equal')
		plt.grid(True, alpha=0.2)
		plt.legend()

		plt.tight_layout()
