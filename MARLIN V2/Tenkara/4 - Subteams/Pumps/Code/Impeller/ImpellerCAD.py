from pathlib import Path
import numpy as np
import cadquery.func
from cadquery.vis import show
from Pump import Propellant
from .MeridionalContour import MeridionalContour
from . import ImpellerBladeDesign

"""
@Charlie Zhu
Minimal AI
Impeller CAD generation
"""

def GetSTEPPath(propellant: Propellant) -> Path:
	"""Default STEP file of the impeller: <propellant>_impeller.step in the Impeller folder."""
	return Path(__file__).parent / f'{propellant.name}_impeller.step'

def _ToCartesian(points_zrt: np.ndarray, flip_z: bool=False) -> list[cadquery.Vector]:
	z, r, theta = points_zrt.T
	if flip_z: z = -z
	x = r * np.cos(theta)
	y = r * np.sin(theta)
	points_xyz = np.column_stack((x, y, z)) * 1e3 # meter to mm
	return [cadquery.Vector(*p) for p in points_xyz]

def _MeridionalContourToCartesian(contour: MeridionalContour, resolution: int=200, flip_z: bool=False) -> list[cadquery.Vector]:
	
	points_zr = contour.Points(np.linspace(0, contour.length_m, resolution))
	z, r = points_zr.T
	if flip_z: z = -z
	y=np.zeros(len(r))
	points_xyz = np.column_stack((r, y, z)) * 1e3 # meter to mm
	return [cadquery.Vector(*p) for p in points_xyz]

class ImpellerCAD:

	def __init__(self, blade_design: ImpellerBladeDesign.ImpellerBladeDesign, flip_z: bool):
		self._blade_design = blade_design
		self.Impeller = self._blade_design.Pump.Impeller
		RESOLUTION: int = 200

		# blade_points is a list of points (ndarrays)
		blade_points, _ = self._blade_design.GetBladeLines()

		t = np.linspace(0, 1, RESOLUTION)
		blade_edges = [cadquery.Edge.makeSpline(_ToCartesian(p, flip_z), parameters=list(t)) for p in blade_points]
		camber = cadquery.func.loft(blade_edges)

		Z_blade: int = self.Impeller.Z_blade
		# rotate the blades about the z-axis
		blades = camber.moved(*[cadquery.Location((0, 0, 0), (0, 0, 1),  360 * n / Z_blade) for n in range(Z_blade)])

		shroud_contour = self._blade_design.MeridionalDesign.ShroudContour
		shroud_edge = cadquery.Edge.makeSpline(
			_MeridionalContourToCartesian(shroud_contour, resolution=RESOLUTION, flip_z=flip_z)
			)
		shroud_surface = cadquery.func.revolve(shroud_edge, (0, 0, 0), (0, 0, 1))

		hub_contour = self._blade_design.MeridionalDesign.HubContour
		hub_edge = cadquery.Edge.makeSpline(
			_MeridionalContourToCartesian(hub_contour, resolution=RESOLUTION, flip_z=flip_z)
			)
		hub_surface = cadquery.func.revolve(hub_edge, (0, 0, 0), (0, 0, 1))

		# L = 1000  # mm, just bigger than the impeller
		# quadrant = cadquery.Solid.makeBox(L, L, 2 * L, pnt=cadquery.Vector(0, 0, -L))

		self._assembly = (
			cadquery.Assembly()
			.add(hub_surface, name='hub', color=cadquery.Color(0.6, 0.6, 0.6, 0.5))
			.add(shroud_surface, name='shroud', color=cadquery.Color(0.3, 0.5, 0.8, 0.3))
			.add(blades, name='blades', color=cadquery.Color(0.85, 0.35, 0.2, 1.0))
		)

		show(self._assembly)

	def SaveSTEP(self, path: str | Path | None = None):
		"""Saves the assembly as a STEP file in mm, keeping the part names and colors.
		By default it goes in the Impeller folder as <propellant>_impeller.step."""
		if path is None:
			path = GetSTEPPath(self._blade_design.Pump.Propellant)

		print(f'Saved STEP to: {path.name}')
		self._assembly.export(str(path), exportType='STEP', unit='MM')
