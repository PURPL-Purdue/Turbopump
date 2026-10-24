from pathlib import Path
import yaml
from pint import Quantity as Q_

script_dir = Path(__file__).resolve().parent
yaml_path = script_dir / "Inputs.yaml"


def ToQuantities(node):
	"""Recursively replace every [value, 'unit'] pair in a parsed YAML tree with Q_(value, 'unit').
	Dicts keep their nesting; anything else is returned unchanged."""
	if isinstance(node, dict):
		return {key: ToQuantities(value) for key, value in node.items()}
	if isinstance(node, list) and len(node) == 2 and isinstance(node[0], (int, float)) and isinstance(node[1], str):
		return Q_(node[0], node[1])
	return node


def LoadInputs(path: Path = yaml_path) -> dict:
	"""Load a YAML inputs file as a nested dict of pint Quantities."""
	with open(path, "r") as file:
		return ToQuantities(yaml.safe_load(file))


inputs = LoadInputs()

if __name__ == "__main__":
	print(inputs)
