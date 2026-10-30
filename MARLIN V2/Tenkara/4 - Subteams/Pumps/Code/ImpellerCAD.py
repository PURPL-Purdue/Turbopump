import ImpellerBladeDesign
import cadquery as cq

class ImpellerCAD:

	def __init__(self, blade_design: ImpellerBladeDesign):
		self._blade_design = blade_design

		#cq.Sketch().
