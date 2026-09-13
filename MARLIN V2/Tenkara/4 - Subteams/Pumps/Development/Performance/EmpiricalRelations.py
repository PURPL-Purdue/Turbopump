def FlowSpeedRatio(Q_now: float, N_now: float, Q_design: float, N_design: float) -> float:
    return (Q_now/N_now) / (Q_design/N_design)

# PUMPA pg 10, eq. 35
def SlipFactor(flowspeed_ratio: float, slipfactor_opt: float = 1) -> float:
    factor: float = 1.534988 - 0.6681668 * flowspeed_ratio + 0.077472 * flowspeed_ratio**2 + 0.0571508 * flowspeed_ratio**3
    return slipfactor_opt * factor

def HydraulicEfficiency(f: float) -> float:
    return 0.86387 + 0.3096 * f - 0.14086 * f**2 - 0.029265 * f**3