import ImpellerBladeDesign
import numpy as np
import cadquery
import cadquery.func
from cadquery.vis import show
from MeridionalContour import MeridionalContour

def _ToCartesian(points_zrt: np.ndarray) -> list[cadquery.Vector]:
	z, r, theta = points_zrt.T
	x = r * np.cos(theta)
	y = r * np.sin(theta)
	points_xyz = np.column_stack((x, y, z)) * 1e3 # meter to mm
	return [cadquery.Vector(*p) for p in points_xyz]

def _MeridionalContourToCartesian(contour: MeridionalContour, resolution: int=200) -> list[cadquery.Vector]:
	
	points_zr = contour.Points(np.linspace(0, contour.length_m, resolution))
	z, r = points_zr.T
	y=np.zeros(len(r))
	points_xyz = np.column_stack((r, y, z)) * 1e3 # meter to mm
	return [cadquery.Vector(*p) for p in points_xyz]

class ImpellerCAD:

	def __init__(self, blade_design: ImpellerBladeDesign.ImpellerBladeDesign):
		self._blade_design = blade_design
		self.Impeller = self._blade_design.Pump.Impeller
		RESOLUTION: int = 200

		# blade_points is a list of points (ndarrays)
		blade_points, _ = self._blade_design.GetBladeLines()

		t = np.linspace(0, 1, RESOLUTION)
		blade_edges = [cadquery.Edge.makeSpline(_ToCartesian(p), parameters=list(t)) for p in blade_points]
		camber = cadquery.func.loft(blade_edges)

		Z_blade: int = self.Impeller.Z_blade
		# rotate the blades about the z-axis
		blades = camber.moved(*[cadquery.Location((0, 0, 0), (0, 0, 1),  360 * n / Z_blade) for n in range(Z_blade)])

		shroud_contour = self._blade_design.MeridionalDesign.ShroudContour
		shroud_edge = cadquery.Edge.makeSpline(
			_MeridionalContourToCartesian(shroud_contour, resolution=RESOLUTION)
			)
		shroud_surface = cadquery.func.revolve(shroud_edge, (0, 0, 0), (0, 0, 1))

		hub_contour = self._blade_design.MeridionalDesign.HubContour
		hub_edge = cadquery.Edge.makeSpline(
			_MeridionalContourToCartesian(hub_contour, resolution=RESOLUTION)
			)
		hub_surface = cadquery.func.revolve(hub_edge, (0, 0, 0), (0, 0, 1))

		# L = 1000  # mm, just bigger than the impeller
		# quadrant = cadquery.Solid.makeBox(L, L, 2 * L, pnt=cadquery.Vector(0, 0, -L))

		impeller = cadquery.func.compound([hub_surface, shroud_surface, blades])

		show(impeller, alpha=0.5)
