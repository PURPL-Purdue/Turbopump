"""Export configured inducer design parameters for LOX and IPA as JSON."""
import json
import sys
from pathlib import Path
from typing import Iterable

if not __package__:
    _here = Path(__file__).resolve().parent
    _code = _here.parent
    sys.path[:] = [p for p in sys.path if Path(p or ".").resolve() != _here]
    sys.path.insert(0, str(_code))

import Main
import Pump1D

OUTPUT_PATH = Path(__file__).with_name("design_parameters.json")

# Hub Parameters
LAMBDA_W = 0.03 # [0.03, 0.06] <= Gülich recommended range
LAMBDA_C = 1.1 # [1, 1.1] <= Gülich recommended range
PSI_MAX = 0.08 # Max head coefficient for straight hub inducer (from ???)
AXIAL_LENGTH_RATIO = 0.6 # NASA SP-125

# Blade parameters
Z_BLADE = 2 
INCIDENCE_RATIO = 0.425 # NASA SP-8052
ETA_BLADE = 0.85 # NASA SP-8052

def _quantity(value, unit: str) -> dict[str, float | str]:
    converted = value.to(unit)
    return {"value": float(converted.m), "unit": unit}


def BuildDesignParameters(propellants: Iterable[Main.Propellant] = (Main.Propellant.LOX, Main.Propellant.IPA)) -> dict:
    """Collect configured pump inputs and model constants, plus Pump1D's calculated required head."""
    records = {}

    for propellant in propellants:
        pump = Main.CreatePump(propellant)
        required_head = Pump1D.GetInducerHeadRequirement(pump)
        name = propellant.name

        records[name] = {
            "operating_point": {
                "mass_flow": _quantity(pump.DP.mass_flow, "kg/s"),
                "volumetric_flow_rate": _quantity(pump.DesignFlowRate, "m^3/s"),
                "shaft_speed": _quantity(pump.DP.N_shaft, "rpm"),
            },
            "impeller_inlet_geometry": {
                "tip_diameter": _quantity(pump.Impeller.d_1, "m"),
                "hub_diameter": _quantity(pump.Impeller.d_hub, "m"),
            },
            "inducer_design": {
                "required_head_rise": _quantity(required_head, "m"),
                "sizing_constants": {
                    "lambda_w": float(LAMBDA_W),
                    "lambda_c": float(LAMBDA_C),
                    "psi_max": float(PSI_MAX),
                    "axial_length_ratio": float(AXIAL_LENGTH_RATIO),
                },
                "blade_parameters": {
                    "blade_count": int(Z_BLADE),
                    "incidence_ratio": float(INCIDENCE_RATIO),
                    "blade_efficiency": float(ETA_BLADE),
                },
            },
        }

    return {"propellants": records}


def WriteDesignParameters(output_path: Path = OUTPUT_PATH) -> Path:
    """Write all configured design parameters and return the JSON path."""
    parameters = BuildDesignParameters()
    output_path = Path(output_path)
    with output_path.open("w", encoding="utf-8") as file:
        json.dump(parameters, file, indent=2)
        file.write("\n")
    return output_path


if __name__ == "__main__":
    print(f"Wrote inducer design parameters to {WriteDesignParameters()}")