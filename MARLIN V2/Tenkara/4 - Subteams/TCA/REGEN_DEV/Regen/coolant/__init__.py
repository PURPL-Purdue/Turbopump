"""
Coolant property modules. Each one provides init(p_Pa), rho(T), mu(T), k(T),
cp(T), and NAME, T_MAX, LIMIT_NOTE for the report's validity checks.

To add one (e.g. ethanol, PLAN item 1): write the module with that interface,
add it to FLUIDS, then pick it in regen.yaml (jacket.coolant).
"""
from coolant import ipa_props

FLUIDS = {
    "ipa": ipa_props,
}


def get(name):
    """The property module registered under `name`."""
    try:
        return FLUIDS[name]
    except KeyError:
        raise ValueError(f"Unknown coolant '{name}'. Available: {', '.join(FLUIDS)}")
