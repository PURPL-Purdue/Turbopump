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

def _TrimToRadius(points: np.ndarray, max_radius: float) -> np.ndarray:
	"""
		Cuts a (z, r, θ) blade line where it first crosses r = max_radius, interpolating the crossing point linearly.
		Lines that never go past max_radius are returned unchanged.
	"""
	beyond = points[:, 1] > max_radius
	if not beyond.any():
		return points
	i = np.argmax(beyond)
	assert i > 0, 'Blade line starts beyond max_radius'
	alpha = (max_radius - points[i-1, 1]) / (points[i, 1] - points[i-1, 1])
	interpolated = points[i-1] + alpha * (points[i] - points[i-1])
	return np.vstack((points[:i], np.array([interpolated[0], max_radius, interpolated[2]])))

_SPLINE_TOL = 1e-3		# mm, approximation tolerance of the blade-line splines

def _ExtendBladeLine(outer: np.ndarray, inner: np.ndarray) -> np.ndarray:
	"""
		(z, r, θ) blade line one streamline step beyond outer, away from inner, by linear extrapolation point by point:
		2 outer - inner. A whole step keeps the spanwise spacing of the loft sections even, so the loft does not overshoot
		between the wall section and its neighbour. Both lines must have the same number of points.
	"""
	return 2 * outer - inner

def _WallTool(edge: cadquery.Edge, length: float) -> cadquery.Shape:
	"""
		Cutting surface for a meridional wall: edge (in the xz plane) with straight lines of length (mm) tangent to it
		at both ends, revolved about the z axis. The middle face revolves edge itself, so it is the same surface as
		cadquery.func.revolve(edge, ...); the end extensions keep the cuts away from the surface boundaries.
	"""
	start, end = edge.startPoint(), edge.endPoint()
	wire = cadquery.Wire.assembleEdges([
		cadquery.Edge.makeLine(start - edge.tangentAt(0) * length, start),
		edge,
		cadquery.Edge.makeLine(end, end + edge.tangentAt(1) * length),
	])
	return cadquery.func.revolve(wire, (0, 0, 0), (0, 0, 1))

class ImpellerCAD:

	def __init__(self, blade_design: ImpellerBladeDesign.ImpellerBladeDesign, flip_z: bool=True, flip_theta=True):
		self._blade_design = blade_design
		self.Impeller = self._blade_design.Pump.Impeller

		RESOLUTION = 200
		shroud_contour = self._blade_design.MeridionalDesign.ShroudContour
		shroud_edge = cadquery.Edge.makeSpline(
			_MeridionalContourToCartesian(shroud_contour, resolution=RESOLUTION, flip_z=flip_z)
			)
		shroud_wire = cadquery.Wire.assembleEdges([
			cadquery.Edge.makeLine(shroud_edge.startPoint(), shroud_edge.endPoint()),
			shroud_edge
		])
		
		shroud_solid = cadquery.func.revolve(cadquery.func.face(shroud_wire), (0, 0, 0), (0, 0, 1))
		if not (shroud_solid.Volume() > 0): print("INVALID SHROUD SOLID!!!")

		hub_contour = self._blade_design.MeridionalDesign.HubContour
		hub_edge = cadquery.Edge.makeSpline(
			_MeridionalContourToCartesian(hub_contour, resolution=RESOLUTION, flip_z=flip_z)
			)

		start_hub_edge = cadquery.Edge.makeLine((0, 0, hub_edge.startPoint().z), hub_edge.startPoint())
		back_shroud_extend = cadquery.Vector(0, 0, 2.54 * (-1 if flip_z else 1))
		back_hub_extend_edge = cadquery.Edge.makeLine(hub_edge.endPoint(), hub_edge.endPoint() + back_shroud_extend)
		back_hub_comeback = cadquery.Edge.makeLine(back_hub_extend_edge.endPoint(), (0, 0, back_hub_extend_edge.endPoint().z))
		hub_close = cadquery.Edge.makeLine(back_hub_comeback.endPoint(), start_hub_edge.startPoint())
		hub_wire = cadquery.Wire.assembleEdges([
			start_hub_edge, hub_edge, back_hub_extend_edge, back_hub_comeback, hub_close
		])
		
		hub_solid = cadquery.func.revolve(cadquery.func.face(hub_wire), (0, 0, 0), (0, 0, 1))
		if not (hub_solid.Volume() > 0): print("INVALID HUB SOLID!!!")

		# the blades are cut by these, so their ends lie on the same surfaces as hub_surface and shroud_surface
		self._wall_tools = cadquery.func.compound(_WallTool(hub_edge, 2), _WallTool(shroud_edge, 2))

		Z_blade_1: int = self.Impeller.Z_blade_1

		blades = []

		blade_sets = [self._blade_design.GetBladeSurfaceLines(('suction', 'pressure'))]
		if self._blade_design.HasSplitters:
			blade_sets.append(self._blade_design.GetBladeSurfaceLines(('suction', 'pressure'), splitter=True))
		
		self._debug_cad = []

		for suction_lines, pressure_lines in blade_sets:
			blade_surface = self._BladeSolid(suction_lines, pressure_lines, flip_z=flip_z, flip_theta=flip_theta)
			# rotate the blades about the z-axis
			blades.append(blade_surface.moved(
				*[cadquery.Location((0, 0, 0), (0, 0, 1),  360 * n / Z_blade_1) for n in range(Z_blade_1)]
			))

		self.Assembly = (
			cadquery.Assembly()
			.add(hub_solid, name='hub', color=cadquery.Color(0.6, 0.6, 0.6, 0.5))
			.add(shroud_solid, name='shroud', color=cadquery.Color(0.3, 0.5, 0.8, 0.3))
		)
		
		self.Assembly.add(cadquery.func.compound(*blades), name='blades', color=cadquery.Color(0.85, 0.35, 0.2, 1.0))

	def _BladeSolid(self, suction_lines: list[np.ndarray], pressure_lines: list[np.ndarray], flip_z: bool=False, flip_theta: bool=False) -> cadquery.Shape:
		"""
			One blade as a solid, ending exactly on the hub and shroud revolve surfaces.
			Each streamline gives a closed wire: suction line from the leading edge to the trailing edge, an arc on r_2 across the trailing edge to the pressure line, and the pressure line back to the leading edge.
			The suction face, the trailing-edge face and the pressure face come from one loft, so they share their edges exactly (no sewing tolerance).
			The blade lines are approximated to _SPLINE_TOL (they pass exactly through their end points). Interpolating every point gives each section its own knots, and the loft merges
			them into surfaces of ~10^4 poles, which make the split below take minutes.
			The blade is extended one streamline step past the hub and shroud (_ExtendBladeLine), the oversized blade is closed with filled surfaces on its first and last wire, and the
			hub and shroud revolve surfaces (self._wall_tools) then cut off both ends. The blade between them is kept, so its end faces are pieces of the hub and shroud surfaces.
			clean() merges the end faces that the seam of the revolve surfaces (θ = 0) splits in two.
			The leading-edge points must be identical on both sides.
		"""
		assert len(suction_lines) == len(pressure_lines), 'Suction and pressure blade lines should have equal number of streamlines'
		
		max_radius: float =self.Impeller.d_2.to('m').m/2
		s = self.Impeller.blade_thick.to('m').m
		Beta_2 = self.Impeller.Beta2B.to('rad').m
		theta_offset = s/np.sin(Beta_2)/max_radius

		# extend the blade past the hub and shroud, the walls cut it back below
		suction_lines, pressure_lines = [
			[_ExtendBladeLine(lines[0], lines[1]), *lines, _ExtendBladeLine(lines[-1], lines[-2])]
			for lines in (suction_lines, pressure_lines)
		]

		wires = []

		for points_s, points_p in zip(suction_lines, pressure_lines):

			points_p = _TrimToRadius(points_p, max_radius)
			# points_s does not exceed radius
			pressure_last_point: np.ndarray = points_p[-1, :] # at the max radius
			
			suction_last_point = np.array([pressure_last_point[0], max_radius, pressure_last_point[2] + theta_offset])
			points_s = np.vstack((points_s, suction_last_point))

			cq_points_s: list[cadquery.Vector] = CAD.CylindricalToCartesian(points_s, flip_z, flip_theta)
			cq_points_p: list[cadquery.Vector] = CAD.CylindricalToCartesian(points_p, flip_z, flip_theta)

			self._debug_cad.append(cq_points_s[-1])
			self._debug_cad.append(cq_points_p[-1])
			assert (cq_points_s[0] - cq_points_p[0]).Length < 1e-9, 'Leading edge points of suction and pressure side differ'
			cq_points_p[0] = cq_points_s[0]	# same vertex

			# trailing edge: arc on the r_2 circle, through the midpoint in theta
			te_mid_point = np.array([[pressure_last_point[0], max_radius, pressure_last_point[2] + theta_offset/2]])
			cq_te_mid_point = CAD.CylindricalToCartesian(te_mid_point, flip_z, flip_theta)[0]

			wires.append(cadquery.Wire.assembleEdges([
				cadquery.Edge.makeSplineApprox(cq_points_s, tol=_SPLINE_TOL, minDeg=3, maxDeg=5),
				cadquery.Edge.makeThreePointArc(cq_points_s[-1], cq_te_mid_point, cq_points_p[-1]),
				cadquery.Edge.makeSplineApprox(cq_points_p[::-1], tol=_SPLINE_TOL, minDeg=3, maxDeg=5),
			]))

		self._debug_cad.extend(wires)
		sides = cadquery.func.loft(wires)
		oversized = cadquery.func.solid(cadquery.func.shell(*sides.Faces(), cadquery.func.fill(wires[0]), cadquery.func.fill(wires[-1])))

		# hub stub, blade, shroud stub: fewer pieces means a wall did not cut all the way across the blade
		pieces = cadquery.func.split(oversized, self._wall_tools).Solids()
		assert len(pieces) == 3, f'The hub and shroud should cut the blade into 3 pieces, got {len(pieces)}'
		return max(pieces, key=lambda piece: piece.Volume()).clean()

	def Show(self, debug: bool=True):
		if debug:
			show(self.Assembly, *self._debug_cad)
		else:
			show(self.Assembly)
