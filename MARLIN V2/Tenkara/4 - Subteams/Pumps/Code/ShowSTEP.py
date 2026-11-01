if __name__ == '__main__':
	from Impeller import ImpellerCAD
	from Pump import Propellant
	import numpy as np

	step_path = ImpellerCAD.GetSTEPPath(Propellant.LOX)
	if step_path.exists():
		import trimesh
		import plotly.graph_objects as go

		scene = trimesh.load_scene(step_path).convert_units('mm')	# cascadio loads STEP in meters

		# one mesh per part (hub, shroud, blades), in the part's STEP color
		data = []
		for name, geometry in scene.geometry.items():
			if isinstance(geometry, trimesh.Trimesh):
				v = geometry.vertices
				f = geometry.faces
				rgba = geometry.visual.material.main_color / 255
				rgb = np.where(rgba[:3] <= 0.0031308, 12.92 * rgba[:3], 1.055 * rgba[:3]**(1 / 2.4) - 0.055)	# glTF linear to sRGB
				data.append(go.Mesh3d(x=v[:,0], y=v[:,1], z=v[:,2],
									i=f[:,0], j=f[:,1], k=f[:,2],
									name=name, showlegend=True,
									color=f'rgb{tuple(int(round(255 * c)) for c in rgb)}',
									opacity=rgba[3]))

		fig = go.Figure(data=data)
		fig.update_layout(
			title=step_path.name,
			template='plotly_white',
			scene=dict(aspectmode='data', xaxis_title='x (mm)', yaxis_title='y (mm)', zaxis_title='z (mm)'),
		)
		fig.show()