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
	print(y.shape)
	print(r.shape)
	points_xyz = np.column_stack((r, y, z)) * 1e3 # meter to mm
	return [cadquery.Vector(*p) for p in points_xyz]

class ImpellerCAD:

	def __init__(self, blade_design: ImpellerBladeDesign.ImpellerBladeDesign):
		self._blade_design = blade_design
		RESOLUTION: int = 200

		# blade_points is a list of points (ndarrays)
		blade_points, _ = self._blade_design.GetBladeLines()

		t = np.linspace(0, 1, RESOLUTION)
		blade_edges = [cadquery.Edge.makeSpline(_ToCartesian(p), parameters=list(t)) for p in blade_points]
		camber = cadquery.func.loft(blade_edges)

		shroud_contour = self._blade_design.MeridionalDesign.ShroudContour
		shroud_edge = [cadquery.Edge.makeSpline(
			_MeridionalContourToCartesian(shroud_contour, resolution=RESOLUTION)
			)]

		hub_contour = self._blade_design.MeridionalDesign.HubContour
		hub_edge = [cadquery.Edge.makeSpline(
			_MeridionalContourToCartesian(hub_contour, resolution=RESOLUTION)
			)]
		
		show(blade_edges, shroud_edge + hub_edge, camber)
