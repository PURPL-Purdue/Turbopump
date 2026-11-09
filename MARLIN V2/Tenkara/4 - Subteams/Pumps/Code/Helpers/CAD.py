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

def ExportPart(assembly: cadquery.Assembly, name: str):
	"""
		Export every body of the assembly as a single part: a flat compound of the solids and the loose faces, with no product structure.
		Siemens NX imports this as one .prt with multiple bodies, where an assembly STEP file becomes one part file per component.
		Names and colors are dropped.
	"""
	shape = assembly.toCompound()	# nested assemblies and locations are resolved
	solids = shape.Solids()
	solid_faces = {f.hashCode() for solid in solids for f in solid.Faces()}
	loose_faces = [f for f in shape.Faces() if f.hashCode() not in solid_faces]
	writer = STEPControl_Writer()	# the writer has to exist before its parameters can be set
	# Shapes that share geometry (the copies of one blade) are otherwise written as product instances, which NX opens as separate parts
	Interface_Static.SetIVal_s('write.step.assembly', 0)
	Interface_Static.SetCVal_s('xstep.cascade.unit', 'MM')
	Interface_Static.SetCVal_s('write.step.unit', 'MM')
	writer.Transfer(cadquery.Compound.makeCompound([*solids, *loose_faces]).wrapped, STEPControl_AsIs)
	path = str(GetStepExportPath(name))
	writer.Write(path)
	print(f'Saved to {path}')
