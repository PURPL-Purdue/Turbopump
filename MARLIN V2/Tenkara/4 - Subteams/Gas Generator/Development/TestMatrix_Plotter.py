import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.colors import Normalize
from matplotlib.text import Text


def separate_overlapping_labels(fig, ax, labels):
	fig.canvas.draw()
	renderer = fig.canvas.get_renderer()
	placed_boxes = []
	axes_box = ax.get_window_extent(renderer)
	directions = np.linspace(0, 2 * np.pi, 16, endpoint=False)

	def label_box(label):
		label.update_positions(renderer)
		return Text.get_window_extent(label, renderer)

	for label in labels:
		for radius in range(12, 301, 12):
			for angle in directions:
				offset = (radius * np.cos(angle), radius * np.sin(angle))
				label.set_position(offset)
				box = label_box(label)
				inside_axes = (
					box.x0 >= axes_box.x0
					and box.x1 <= axes_box.x1
					and box.y0 >= axes_box.y0
					and box.y1 <= axes_box.y1
				)
				if inside_axes and all(
					box.x1 + 2 < placed.x0
					or box.x0 - 2 > placed.x1
					or box.y1 + 2 < placed.y0
					or box.y0 - 2 > placed.y1
					for placed in placed_boxes
				):
					placed_boxes.append(box)
					break
			else:
				continue
			break
		else:
			placed_boxes.append(label_box(label))


def main():
	csv_path = 'TestMatrix_Outputs.csv'
	data = pd.read_csv(csv_path)

	fuel_column = 'stiff_fuel'
	ox_column = 'stiff_ox'
	for column in (ox_column, fuel_column):
		if column not in data.columns:
			raise ValueError(f"Column {column!r} is not in the CSV.")
	data['min stiff'] = data[[ox_column, fuel_column]].min(axis=1)

	print("Columns:", ", ".join(str(column) for column in data.columns))
	# x_column = input("X-axis column: ").strip()
	# y_column = input("Y-axis column: ").strip()
	x_column = 'P_fuel_psia'
	y_column = 'P_ox_psia'
	#color_column = input("Color-by column: ").strip()
	color_column = 'mdot_total_kg_s'
	label_column = input("Label-by column: ").strip()


	for column in (x_column, y_column, color_column, label_column):
		if column not in data.columns:
			raise ValueError(f"Column {column!r} is not in the CSV.")

	color_values = data[color_column].dropna()
	value_min = color_values.min()
	value_max = color_values.max()
	norm = Normalize(vmin=value_min, vmax=value_max)
	cmap = plt.get_cmap("viridis")

	fig, ax = plt.subplots()
	scatter = ax.scatter(
		data[x_column],
		data[y_column],
		c=data[color_column],
		cmap=cmap,
		norm=norm,
	)
	labels = [
		ax.annotate(
			str(label),
			(x_value, y_value),
			xytext=(5, 5),
			textcoords="offset points",
			arrowprops={
				"arrowstyle": "-",
				"color": "gray",
				"linewidth": 0.5,
				"connectionstyle": "arc3,rad=0",
			},
		)
		for x_value, y_value, label in zip(
			data[x_column], data[y_column], data[label_column]
		)
	]
	ticks = [value_min] if value_min == value_max else np.linspace(value_min, value_max, 10)
	fig.colorbar(scatter, ax=ax, ticks=ticks, format="%.2f", label=color_column)
	ax.set_title(y_column)
	ax.set_xlabel(x_column)
	ax.set_ylabel(y_column)
	ax.grid(True)
	fig.tight_layout()
	separate_overlapping_labels(fig, ax, labels)
	fig.canvas.draw()
	plt.show()


if __name__ == "__main__":
	main()
