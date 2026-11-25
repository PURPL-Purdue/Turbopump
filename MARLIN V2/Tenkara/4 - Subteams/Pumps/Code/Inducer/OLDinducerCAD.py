import numpy as np
from pint import Quantity as Q_
import cadquery as cq
import sys
from pathlib import Path

sys.path.append(str(Path(__file__).resolve().parents[1]))	# Pumps folder, so the Inputs package imports when run directly

# ---------------------------------------------------------------------------------------------------------------------
# Inputs
# ---------------------------------------------------------------------------------------------------------------------
d_tip = Q_(1.0, 'in')			# inducer tip diameter, constant from inlet to outlet (straight shroud)
d_hub_1 = Q_(0.4, 'in')		# hub diameter at the inlet
d_hub_2 = Q_(0.6, 'in')			# hub diameter at the outlet, straight taper from the inlet
beta_1 = Q_(12.86, 'deg')		# blade angle at the inlet, rms radius, from tangential (eq.beta1b of InducerRadialEquilibrium)
beta_2 = Q_(25.5, 'deg')		# blade angle at the outlet, rms radius, from tangential (eq.beta2b of InducerRadialEquilibrium)
blade_thickness = Q_(1, 'mm')	# normal blade thickness, reached once the leading edge wedge has grown to it
LE_sweep = Q_(60, 'deg')		# leading edge sweepback: how far the tip leading edge is wrapped behind the hub leading edge
Z_blade = 2						# number of blades
L_ax = 0.8 * d_tip				# axial length of the hub taper

LE_thickness = Q_(0.1, 'mm')	# blade thickness at the very leading edge, the suction face grows from here up to blade_thickness
wedge_angle = Q_(6, 'deg')	# angle of the leading edge wedge (suction face to the flat pressure face), None takes SP-8052 eqs. 69-70 from flow_coefficient (see WedgeAngle)
flow_coefficient = 0.108		# tip flow coefficient c_m1 / u_tip at the inlet, only used when wedge_angle is None
round_TE = True					# rounds the trailing edge, False for a blunt flat end
flip_z = True					# same orientation as InducerRadialEquilibrium.CreateCAD, inlet at z = 0 and flow towards -z
flip_theta = True				# False for both flips gives the raw (z, r, theta) orientation
name = 'Inducer'				# name of the exported STEP part
# ---------------------------------------------------------------------------------------------------------------------


def WedgeAngle(beta: Q_, flow_coefficient: float) -> Q_:
	"""SP-8052's blade wedge angle (eqs. 69-70): the blade angle less arctan(1.10 flow coefficient). beta is the tip blade angle at the inlet."""
	return beta - Q_(np.arctan(1.10 * flow_coefficient), 'rad')


class InducerCAD:
	"""
	Solid inducer from a straight tapered hub and constant lead blades, following InducerRadialEquilibrium.GetBladeLines:
	- the blade angle follows the parabolic camber beta(z) = beta_1 + (beta_2 - beta_1) (z / L)^2 at the rms radius
	- the lead lambda(z) = r_ms tan(beta) is the same at every radius, so tan(beta(r, z)) = lambda(z) / r and theta is the same at every r
	- blade lines are moved axially so every span line leaves the hub perpendicular to it
	- the leading edge is an Archimedean spiral, theta = 0 at the inlet hub radius to theta = LE_sweep at the tip
	- the leading edge is a one sided wedge: the pressure face stays flat (a constant offset from the blade line), only the suction face opens up,
	  from LE_thickness at the edge with the angle wedge_angle to the pressure face until the thickness reaches blade_thickness
	- the suction face is the one on the side of increasing theta (the blade is a helix that advances with theta, rotation is towards -theta), swap_sides flips it
	The pressure and suction surfaces are fitted through the span lines, joined by the edge faces into a solid and fused to the hub. Lengths are in mm.
	"""

	def __init__(self, d_tip: Q_, d_hub_1: Q_, d_hub_2: Q_, beta_1: Q_, beta_2: Q_, blade_thickness: Q_, LE_sweep: Q_,
				 Z_blade: int = 2, L_ax: Q_ = None, round_TE: bool = True, flip_z: bool = True, flip_theta: bool = True,
				 LE_thickness: Q_ = Q_(0.1, 'mm'), wedge_angle: Q_ = None, flow_coefficient: float = None,
				 num_sections: int = 9, resolution: int = 160, hub_embed: Q_ = Q_(0.5, 'mm'), swap_sides: bool = True):
		if L_ax is None:
			L_ax = 0.6 * d_tip
		self.R_t, self.R_h1, self.R_h2, self.L, self.t, self.embed = (v.to('m').magnitude for v in (d_tip / 2, d_hub_1 / 2, d_hub_2 / 2, L_ax, blade_thickness, hub_embed))
		if not (self.R_h1 < self.R_t and self.R_h2 < self.R_t):
			raise ValueError('Hub diameters must be smaller than the tip diameter.')
		self.beta_1, self.beta_2 = beta_1.to('rad').magnitude, beta_2.to('rad').magnitude
		self.sweep = LE_sweep.to('rad').magnitude
		self.Z_blade, self.round_TE = Z_blade, round_TE
		self.flip_z, self.flip_theta = flip_z, flip_theta
		self.num_sections, self.resolution = num_sections, resolution
		self.swap_sides = swap_sides
		self.slope = (self.R_h2 - self.R_h1) / self.L

		if wedge_angle is None:
			if flow_coefficient is None:
				raise ValueError('Give wedge_angle or flow_coefficient.')
			wedge_angle = WedgeAngle(Q_(np.degrees(np.arctan(self._Lead(0.0) / self.R_t)), 'deg'), flow_coefficient)
		self.wedge = wedge_angle.to('rad').magnitude
		if not 0 < self.wedge < np.pi / 2:
			raise ValueError(f'The wedge angle must be positive, it is {np.degrees(self.wedge):.2f} deg. Check the flow coefficient and blade angle.')
		self.t_LE = min(LE_thickness.to('m').magnitude, self.t)

		# wrap angle of the blade reference line, theta(z) = integral of dz / lambda
		self._z_ref = np.linspace(0, self.L, 4001)
		self._theta_ref = np.concatenate(([0], np.cumsum(0.5 * (1 / self._Lead(self._z_ref[1:]) + 1 / self._Lead(self._z_ref[:-1])) * np.diff(self._z_ref))))
		if self.sweep > self._theta_ref[-1]:
			raise ValueError('The leading edge wraps further than the blade, reduce LE_sweep.')

		self.Hub = self._MakeHub()
		self.Blade = self._MakeBlade()
		blades = [self.Blade.rotate(cq.Vector(0, 0, 0), cq.Vector(0, 0, 1), 360 * n / Z_blade) for n in range(Z_blade)]
		self.Solid = self.Hub.fuse(*blades).clean()
		self.Assembly = (cq.Assembly()
			.add(self.Hub, name='hub', color=cq.Color(0.6, 0.6, 0.6))
			.add(cq.Compound.makeCompound(blades), name='blades', color=cq.Color(0.2, 0.3, 0.7)))

	# ---- geometry in SI floats ----------------------------------------------------------------------------------------
	def _HubRadius(self, z):
		return self.R_h1 + self.slope * z

	def _Beta(self, z):
		return self.beta_1 + (self.beta_2 - self.beta_1) * (z / self.L)**2

	def _Lead(self, z):
		return np.sqrt((self.R_t**2 + self._HubRadius(z)**2) / 2) * np.tan(self._Beta(z))		# r_ms tan(beta)

	def _Theta(self, z):
		return np.interp(z, self._z_ref, self._theta_ref)

	def _LeadingEdgeTheta(self, r):
		return self.sweep * (r - self.R_h1) / (self.R_t - self.R_h1)

	def _Point(self, alpha, z_ref):
		"""Surface point (x, y, z) in m and its unit tangent along the blade line and unit normal, at span fraction alpha."""
		r_hub = self._HubRadius(z_ref)
		span = self.R_t - r_hub
		r = r_hub + alpha * span
		z = z_ref - alpha * span * self.slope		# span line perpendicular to the hub
		theta = self._Theta(z_ref)
		dtheta = 1 / self._Lead(z_ref)
		dr, dz = self.slope * (1 - alpha), 1 + alpha * self.slope**2
		T1 = np.stack([dr * np.cos(theta) - r * np.sin(theta) * dtheta, dr * np.sin(theta) + r * np.cos(theta) * dtheta, dz * np.ones_like(r)], axis=-1)
		T2 = np.stack([span * np.cos(theta), span * np.sin(theta), -span * self.slope], axis=-1)
		n = np.cross(T1, T2)
		n /= np.linalg.norm(n, axis=-1, keepdims=True)
		T1 /= np.linalg.norm(T1, axis=-1, keepdims=True)
		return np.stack([r * np.cos(theta), r * np.sin(theta), z], axis=-1), T1, n

	def _Vector(self, p) -> cq.Vector:
		x, y, z = np.asarray(p) * 1e3		# m to mm
		return cq.Vector(x, -y if self.flip_theta else y, -z if self.flip_z else z)

	def _Line(self, alpha: float, z_start: float, z_end: float):
		# Blade line at span fraction alpha, the unit normal d pointing to the suction side, and the offsets of the two faces along d.
		# The pressure face is flat at -t/2, the suction face starts at LE_thickness above it and opens with the wedge angle up to full thickness.
		z = np.linspace(z_start, z_end, self.resolution)
		p, T, d = self._Point(alpha, z)
		theta_hat = np.stack([-np.sin(self._Theta(z)), np.cos(self._Theta(z)), np.zeros_like(z)], axis=-1)
		d *= np.sign(np.sum(d * theta_hat, axis=-1, keepdims=True)) * (-1 if self.swap_sides else 1)
		s = np.concatenate(([0], np.cumsum(np.linalg.norm(np.diff(p, axis=0), axis=1))))		# distance from the leading edge
		off_p = np.full_like(s, -0.5 * self.t)
		u = self.t_LE + s * np.tan(self.wedge) - self.t		# wedge thickness less full thickness
		b = 0.15 * self.t		# half width of the blend where the wedge runs into the constant thickness, removes the kink the surface fit would ripple on
		g = np.where(u < -b, 0, np.where(u > b, u, (u + b)**2 / (4 * b)))
		off_s = off_p + (self.t + u - g)
		return p, d, off_p, off_s

	def _Limits(self, alpha: float) -> tuple[float, float]:
		# The line starts at the swept leading edge and stops short of the hub ends, so that with the thickness added (mostly axial, the blade is a
		# shallow helix) the blade ends at the hub inlet and outlet planes and does not poke through them.
		margin = 0.01e-3		# m, keeps the blade corners just inside the hub ends
		z_ref = np.linspace(0, self.L, 4001)
		r = self._HubRadius(z_ref) + alpha * (self.R_t - self._HubRadius(z_ref))
		z_start = np.interp(0, self._Theta(z_ref) - self._LeadingEdgeTheta(r), z_ref)		# first point behind the leading edge
		p, T, n = self._Point(alpha, z_ref)
		reach = p[:, 2] + 0.5 * self.t * np.abs(n[:, 2])		# axial reach of the full thickness (the faces never go beyond +-t/2)
		z_end = z_ref[len(z_ref) - 1 - np.argmax((reach <= self.L - margin)[::-1])]
		for _ in range(30):		# the wedge thickens along the line, so move the start in until nothing pokes out of the inlet plane
			p, d, off_p, off_s = self._Line(alpha, z_start, z_end)
			deficit = margin - min(np.min(p[:, 2] + off_p * d[:, 2]), np.min(p[:, 2] + off_s * d[:, 2]))
			if deficit <= 0:
				break
			z_start += 1.2 * deficit
		return z_start, z_end

	@staticmethod
	def _Boundary(face: cq.Face, grid: list) -> dict:
		# The four edges of a fitted surface, found from the corners of the point grid it was fitted to
		corners = {'root': (grid[0][0], grid[0][-1]), 'tip': (grid[-1][0], grid[-1][-1]), 'LE': (grid[0][0], grid[-1][0]), 'TE': (grid[0][-1], grid[-1][-1])}
		near = lambda v, w: (v - w).Length < 0.05
		out = {}
		for name, (a, b) in corners.items():
			for e in face.Edges():
				s, f = e.startPoint(), e.endPoint()
				if (near(s, a) and near(f, b)) or (near(s, b) and near(f, a)):
					out[name] = e
		return out

	def _MakeBlade(self) -> cq.Solid:
		# Pressure and suction surfaces fitted through the blade lines of all span lines, closed by ruled faces along the root, tip and the two edges
		# The first span line sits inside the hub so the fuse leaves no gap at the root
		alpha_root = -self.embed / (self.R_t - 0.5 * (self.R_h1 + self.R_h2))
		alphas = np.concatenate(([alpha_root], np.linspace(0, 1, self.num_sections)))
		lines = [self._Line(a, *self._Limits(a)) for a in alphas]
		grids = [[[self._Vector(x) for x in p + off[:, None] * d] for p, d, *offs in lines for off in [offs[k]]] for k in (1, 0)]		# suction, pressure
		suction, pressure = (cq.Face.makeSplineApprox(g, tol=2e-4) for g in grids)
		edges_s, edges_p = self._Boundary(suction, grids[0]), self._Boundary(pressure, grids[1])
		sides = []
		for name in ('root', 'tip', 'LE', 'TE'):
			e_s, e_p = edges_s[name], edges_p[name]
			if (e_s.startPoint() - e_p.startPoint()).Length > (e_s.startPoint() - e_p.endPoint()).Length:
				e_p = cq.Edge(e_p.wrapped.Reversed())
			sides.append(cq.Face.makeRuledSurface(e_s, e_p))
		blade = cq.Solid.makeSolid(cq.Shell.makeShell([suction, pressure, *sides]))
		return self._RoundTE(blade, [edges_s['TE'], edges_p['TE']]) if self.round_TE else blade

	def _RoundTE(self, blade: cq.Solid, edges: list) -> cq.Solid:
		# Fillet the two edges of the trailing edge face
		at = lambda e: e.Center()
		selected = [b for b in blade.Edges() if any((at(b) - at(e)).Length < 1e-3 for e in edges)]
		for fraction in (0.45, 0.4, 0.3, 0.2):
			try:
				rounded = blade.fillet(fraction * self.t * 1e3, selected)
				if rounded.isValid():
					return rounded
			except Exception:
				continue
		return blade		# blunt trailing edge if no fillet works

	def _MakeHub(self) -> cq.Solid:
		s = -1 if self.flip_z else 1
		pts = [(0, 0), (self.R_h1 * 1e3, 0), (self.R_h2 * 1e3, s * self.L * 1e3), (0, s * self.L * 1e3)]
		return cq.Workplane('XZ').polyline(pts).close().revolve(360, (0, 0, 0), (0, 1, 0)).val()

	# ---- public API ---------------------------------------------------------------------------------------------------
	def Export(self, path: str = 'Inducer.step') -> None:
		cq.exporters.export(self.Solid, path)

	def Show(self) -> None:
		from cadquery.vis import show
		show(self.Assembly)


def InputsFromRadialEquilibrium(eq, blade_thickness: Q_, LE_sweep: Q_ = Q_(60, 'deg'), **kwargs) -> dict:
	"""Inputs of InducerCAD from an InducerRadialEquilibrium: geometry from its inducer, blade angles from beta1b and beta2b,
	wedge angle from the inducer's FlowCoefficient unless wedge_angle or flow_coefficient is given."""
	inducer = eq.inducer
	if 'wedge_angle' not in kwargs and 'flow_coefficient' not in kwargs:
		phi = inducer.FlowCoefficient
		kwargs['flow_coefficient'] = phi.to('dimensionless').magnitude if hasattr(phi, 'to') else phi
	return dict(d_tip=2 * inducer.R_tip, d_hub_1=2 * inducer.R_hub_1, d_hub_2=2 * inducer.R_hub_2, L_ax=inducer.L_ax,
				beta_1=Q_(np.degrees(eq.beta1b), 'deg'), beta_2=Q_(np.degrees(eq.beta2b), 'deg'), Z_blade=inducer.Z_blade,
				blade_thickness=blade_thickness, LE_sweep=LE_sweep, **kwargs)


if __name__ == "__main__":
	inducer_cad = InducerCAD(d_tip, d_hub_1, d_hub_2, beta_1, beta_2, blade_thickness, LE_sweep, Z_blade, L_ax, round_TE, flip_z, flip_theta,
								 LE_thickness, wedge_angle, flow_coefficient)
	print(f"Solid valid: {inducer_cad.Solid.isValid()}, volume {inducer_cad.Solid.Volume():.1f} mm^3")
	# Save STEP file
	from Helpers import CAD
	do_save = True if input(f'Save {name} STEP? (y):\n') == 'y' else False
	if do_save:
		CAD.ExportPart(inducer_cad.Assembly, name)
	inducer_cad.Show()