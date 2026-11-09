from pathlib import Path
import numpy as np
import cadquery.func
from cadquery.vis import show
from Helpers import CAD
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
	return Path(__file__).resolve().parents[2] / 'Outputs' / f'{propellant.name}_impeller.step'

def _MeridionalContourToCartesian(contour: MeridionalContour, resolution: int=200, flip_z: bool=False) -> list[cadquery.Vector]:
	
	points_zr = contour.Points(np.linspace(0, contour.length_m, resolution))
	z, r = points_zr.T
	if flip_z: z = -z
	y=np.zeros(len(r))
	points_xyz = np.column_stack((r, y, z)) * 1e3 # meter to mm
	return [cadquery.Vector(*p) for p in points_xyz]

class ImpellerCAD:

	def __init__(self, blade_design: ImpellerBladeDesign.ImpellerBladeDesign, flip_z: bool=True, flip_theta=True):
		self._blade_design = blade_design
		self.Impeller = self._blade_design.Pump.Impeller
		Z_blade_1: int = self.Impeller.Z_blade_1

		blades = []
		for is_splitter in (False, True) if self._blade_design._splitter else (False,):
			
			blade_surfaces, _ = self._blade_design._splitter_bladelines if is_splitter else self._blade_design._bladelines

			suction_blade_points, pressure_blade_points, _ = blade_surfaces

			# blade_points is a list of point arrays (ndarrays)
			for blade_points in (suction_blade_points, pressure_blade_points):

				t = np.linspace(0, 1, len(blade_points[0]))
				blade_edges = [cadquery.Edge.makeSpline(CAD.CylindricalToCartesian(p, flip_z, flip_theta), parameters=list(t)) for p in blade_points]
				blade_surface = cadquery.func.loft(blade_edges)
				# rotate the blades about the z-axis
				blades.append(blade_surface.moved(
					*[cadquery.Location((0, 0, 0), (0, 0, 1),  360 * n / Z_blade_1) for n in range(Z_blade_1)]
				))
				

		RESOLUTION = 200
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

		self.Assembly = (
			cadquery.Assembly()
			.add(hub_surface, name='hub', color=cadquery.Color(0.6, 0.6, 0.6, 0.5))
			.add(shroud_surface, name='shroud', color=cadquery.Color(0.3, 0.5, 0.8, 0.3))
		)
		
		self.Assembly.add(cadquery.func.compound(*blades), name='blades', color=cadquery.Color(0.85, 0.35, 0.2, 1.0))

	def Show(self):
		show(self.Assembly)