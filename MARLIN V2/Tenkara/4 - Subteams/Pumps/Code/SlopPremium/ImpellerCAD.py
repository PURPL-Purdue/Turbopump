# 3D solid of the closed impeller, from the meridional contours and the mean blade lines, with CadQuery
import numpy as np
import cadquery as cq
from cadquery.vis import show
from cadquery.func import loft, fill, shell, solid
from pint import Quantity as Q_
from typing import Optional
from pathlib import Path
from ImpellerBladeDesign import ImpellerBladeDesign, BladeLine

class ImpellerCAD:
	"""
		Closed impeller as CadQuery solids, in mm, with the axis along +z and z = 0 at the impeller inlet. Built from an
		ImpellerBladeDesign:
			Hub:		the region behind the hub contour, from the bore (or the axis) out to the contour, closed by a back
						disk disk_thickness thick behind the outlet
			Shroud:		a cover of constant shroud_thickness outside the shroud contour
			Blade:		one blade, lofted from hub to shroud through sections on the streamlines phis (GetBladeLine)
			Blades:		the Z_blade copies of the blade around the axis
		The parts are separate, overlapping solids, not fused yet (Assembly). The blade sections are extended past the
		hub and shroud by overlap (a fraction of the span), so the blades will cut into both walls when fused. Near the
		shroud they also poke a little upstream of the inlet plane z = 0, which is not trimmed yet.
	"""

	def __init__(self, blade_design: ImpellerBladeDesign, shroud_thickness: Q_=Q_(0.04, 'in'),
			disk_thickness: Q_=Q_(0.06, 'in'), bore_diameter: Optional[Q_]=None,
			phis: tuple[float, ...]=tuple(np.linspace(0, 1, 9)), le_ratio: float=2, overlap: float=0.1, n: int=80):

		self._design = blade_design
		self._mer = blade_design.Meridional
		impeller = blade_design.Pump.Impeller

		self.Z_blade = impeller.Z_blade
		self._t = impeller.blade_thick.to('mm').m
		self._t_shroud = shroud_thickness.to('mm').m
		self._t_disk = disk_thickness.to('mm').m
		self._r_bore = 0 if bore_diameter is None else bore_diameter.to('mm').m / 2
		self._r_2 = impeller.d_2.to('mm').m / 2
		assert self._r_bore < impeller.d_hub.to('mm').m / 2, 'The bore must be smaller than the hub'

		# a revolved wall has a seam in the plane theta = 0, where the trailing edge of the first blade lies, and booleans
		# through the seam fail silently (empty result): turn the walls half a pitch to put the seam between two blades
		axis = (cq.Vector(0, 0, 0), cq.Vector(0, 0, 1))
		self.Hub = self._MakeHub().rotate(*axis, 180 / self.Z_blade)
		self.Shroud = self._MakeShroud().rotate(*axis, 180 / self.Z_blade)
		self.Blade = self._MakeBlade(phis, le_ratio, overlap, n)
		self.Blades = [self.Blade.rotate(*axis, 360 * k / self.Z_blade) for k in range(self.Z_blade)]

	@staticmethod
	def _ContourPoints(contour, n: int=201) -> np.ndarray:
		"""(n, 2) array of (r, z) in mm along a meridional contour, inlet to outlet, equally spaced in arclength"""
		z, r = contour(np.linspace(0, 1, n))
		return np.column_stack([r.to('mm').m, z.to('mm').m])

	@staticmethod
	def _Revolve(profile: cq.Workplane) -> cq.Solid:
		"""Solid of revolution of a closed profile sketched on the XZ workplane (sketch x = r, y = z) about the z axis"""
		return profile.revolve(360, (0, 0, 0), (0, 1, 0)).val()

	def _MakeHub(self) -> cq.Solid:
		hub = self._ContourPoints(self._mer.HubContour)
		z_back = hub[-1, 1] + self._t_disk
		return self._Revolve(cq.Workplane('XZ')
			.moveTo(self._r_bore, hub[0, 1]).lineTo(*hub[0])
			.spline([tuple(p) for p in hub[1:]], includeCurrent=True)
			.lineTo(hub[-1, 0], z_back).lineTo(self._r_bore, z_back)
			.close()
		)

	def _MakeShroud(self) -> cq.Solid:
		"""
			Cover between the shroud contour and its offset by shroud_thickness along the normal, away from the passage.
			That is the concave side of the bend, so the offset only stays smooth where the bend radius exceeds the
			thickness.
		"""
		inner = self._ContourPoints(self._mer.ShroudContour)
		tangent = np.gradient(inner, axis=0)
		ds = np.linalg.norm(tangent, axis=1, keepdims=True)
		tangent /= ds
		curvature = np.linalg.norm(np.gradient(tangent, axis=0), axis=1) / ds[:, 0]
		assert 1 / curvature.max() > self._t_shroud, 'The shroud bend is tighter than the shroud thickness'

		outer = inner + self._t_shroud * np.column_stack([tangent[:, 1], -tangent[:, 0]])
		return self._Revolve(cq.Workplane('XZ')
			.moveTo(*inner[0]).spline([tuple(p) for p in inner[1:]], includeCurrent=True)
			.lineTo(*outer[-1]).spline([tuple(p) for p in outer[-2::-1]], includeCurrent=True)
			.close()
		)

	@staticmethod
	def _Resample(line: BladeLine, m_bar: np.ndarray) -> np.ndarray:
		"""(n, 5) array of z, r (mm), theta, beta_B (rad) and m (mm) at fractions m_bar of the blade line's length"""
		m = line.m.to('mm').m
		x = m_bar * m[-1]
		columns = (line.z.to('mm').m, line.r.to('mm').m, line.Theta.to('rad').m, line.Beta.to('rad').m)
		return np.column_stack([np.interp(x, m, q) for q in columns] + [x])

	def _Section(self, line: np.ndarray, le_ratio: float) -> cq.Wire:
		"""
			Closed blade section on a stream surface, from a resampled blade line (_Resample). The half thickness t/2
			grows from zero at the leading edge as an ellipse of semi-axes le_ratio * t/2 along m and t/2 across
			(Milos eq. 61), and is laid off in the circumferential direction as t / (2 sin beta_B) (Milos eq. 28, no
			lean): pressure side ahead in the direction of rotation, suction side behind. The pressure and suction sides
			are separate splines that meet at the leading edge, so the loft joins leading edge to leading edge, and a
			straight line closes the blunt trailing edge.
		"""
		z, r, theta, beta, m = line.T
		a = le_ratio * self._t / 2
		half = self._t / 2 * np.sqrt(1 - (1 - np.clip(m / a, 0, 1))**2)
		dtheta = half / (r * np.sin(beta))

		def Cartesian(th: np.ndarray) -> list[cq.Vector]:
			return [cq.Vector(*p) for p in np.column_stack([r * np.cos(th), r * np.sin(th), z])]

		ps, ss = Cartesian(theta + dtheta), Cartesian(theta - dtheta)
		return cq.Wire.assembleEdges([cq.Edge.makeSpline(ps[::-1]), cq.Edge.makeSpline(ss), cq.Edge.makeLine(ss[-1], ps[-1])])

	def _MakeBlade(self, phis: tuple[float, ...], le_ratio: float, overlap: float, n: int) -> cq.Solid:
		"""
			Blade through the sections (_Section) on streamlines phis, plus one extended past each wall. The loft is
			ruled, straight between sections: a smooth loft approximates the pressure and suction sides separately, and
			errors of a fraction of a millimetre change the 1 mm thickness (about 7% thinner on average with 5 or 9
			sections). The end sections lie on curved stream surfaces, so they are capped with filled surfaces, not
			planes. The caps end up inside the walls.
		"""
		# stations clustered toward the leading edge to resolve its rounding, the same fractions of length on every line
		m_bar = 1 - np.cos(np.linspace(0, np.pi / 2, n))
		lines = [self._Resample(self._design.GetBladeLine(phi), m_bar) for phi in phis]

		def Extend(wall: np.ndarray, inner: np.ndarray, dphi: float) -> np.ndarray:
			"""Section overlap past the wall, extrapolated linearly in phi in z, r and theta"""
			extended = wall.copy()
			extended[:, :3] += overlap / dphi * (wall[:, :3] - inner[:, :3])
			return extended

		sections = [Extend(lines[0], lines[1], phis[1] - phis[0]), *lines, Extend(lines[-1], lines[-2], phis[-1] - phis[-2])]
		wires = [self._Section(line, le_ratio) for line in sections]
		blade = solid(shell(*loft(wires, ruled=True).Faces(), fill(wires[0]), fill(wires[-1])))
		assert blade.isValid(), 'The blade loft is not a valid solid'
		return blade

	def Assembly(self, shroud: bool=True) -> cq.Assembly:
		"""Hub, blades and, unless shroud=False, a see-through shroud, as named and colored parts"""
		assy = cq.Assembly(name='Impeller')
		assy.add(self.Hub, name='Hub', color=cq.Color(0.6, 0.6, 0.6))
		for k, blade in enumerate(self.Blades):
			assy.add(blade, name=f'Blade {k + 1}', color=cq.Color(0.85, 0.55, 0.1))
		if shroud:
			assy.add(self.Shroud, name='Shroud', color=cq.Color(0.3, 0.5, 0.9, 0.3))
		return assy

	def Show(self, shroud: bool=True):
		"""Opens the CadQuery viewer on the impeller (Assembly) and waits until it is closed"""
		show(self.Assembly(shroud), title='Impeller')

	def ExportSTEP(self, path: str):
		"""
			Writes the impeller parts to the STEP file path, in mm, as separate solids (not fused). CadQuery writes
			nothing and raises no error when path is a folder, hence the check.
		"""
		self.Assembly().export(path, exportType='STEP')
		assert Path(path).is_file(), f'No STEP file was written to {path}: give a file path, not a folder'
		print(f'Saved impeller CAD file to {path}')