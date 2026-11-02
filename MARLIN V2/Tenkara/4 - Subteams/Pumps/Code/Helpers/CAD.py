from pathlib import Path
import numpy as np
import cadquery.func
from cadquery.vis import show

def CylindricalToCartesian(points_zrt: np.ndarray, flip_z: bool=False, flip_theta: bool=False) -> list[cadquery.Vector]:
	z, r, theta = points_zrt.T
	if flip_z: z = -z
	if flip_theta: theta = -theta
	x = r * np.cos(theta)
	y = r * np.sin(theta)
	points_xyz = np.column_stack((x, y, z)) * 1e3 # meter to mm
	return [cadquery.Vector(*p) for p in points_xyz]