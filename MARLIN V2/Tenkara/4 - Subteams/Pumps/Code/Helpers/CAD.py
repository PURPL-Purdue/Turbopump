from pathlib import Path
import numpy as np
import cadquery.func
from cadquery.vis import show
from OCP.Interface import Interface_Static
from OCP.STEPControl import STEPControl_Writer, STEPControl_AsIs

def CylindricalToCartesian(points_zrt: np.ndarray, flip_z: bool=False, flip_theta: bool=False) -> list[cadquery.Vector]:
	z, r, theta = points_zrt.T
	if flip_z: z = -z
	if flip_theta: theta = -theta
	x = r * np.cos(theta)
	y = r * np.sin(theta)
	points_xyz = np.column_stack((x, y, z)) * 1e3 # meter to mm
	return [cadquery.Vector(*p) for p in points_xyz]

def GetStepExportPath(file_name: str) -> Path:
	return Path(__file__).resolve().parents[2] / 'Outputs' / f'{file_name}.step'

def ExportStep(*shapes: cadquery.Shape, name: str):
	"""
		Write the shapes to Outputs/<name>.step as the bodies of one part (a single STEP product), with no product structure.
		Siemens NX imports this as one .prt with multiple bodies, where an assembly STEP file becomes one part file per component.
		Nested compounds and moved copies (which share geometry) are fine. Solids, shells, faces, wires and edges are all written.
	"""
	writer = STEPControl_Writer()	# creating the first writer loads the STEP defaults, so the parameters are set after it
	# The default, 2 (auto), writes every compound as an assembly, which NX opens as one part file per body
	Interface_Static.SetIVal_s('write.step.assembly', 0)
	Interface_Static.SetCVal_s('xstep.cascade.unit', 'MM')
	Interface_Static.SetCVal_s('write.step.unit', 'MM')
	writer.Transfer(cadquery.func.compound(*shapes).wrapped, STEPControl_AsIs)
	path = str(GetStepExportPath(name))
	writer.Write(path)
	print(f'Saved to {path}')

def ExportPart(assembly: cadquery.Assembly, name: str):
	"""Every shape of the assembly, nested assemblies and locations resolved, as one part with ExportStep. Names and colors are dropped."""
	ExportStep(assembly.toCompound(), name=name)
