import ImpellerBladeDesign
import numpy as np
import cadquery as cq
from cadquery.vis import show

def _ToCartesian(points_zrt: np.ndarray) -> np.ndarray:
	z, r, theta = points_zrt.T
	x = r * np.cos(theta)
	y = r * np.sin(theta)
	points_xyz = np.column_stack((x, y, z))
	return [cq.Vector(*p) for p in points_xyz]

class ImpellerCAD:

	def __init__(self, blade_design: ImpellerBladeDesign.ImpellerBladeDesign):
		self._blade_design = blade_design

		blade_points, _ = self._blade_design.CreateBladeLines()

		edges = [cq.Edge.makeSpline(_ToCartesian(p)) for p in blade_points]
		show(edges)
