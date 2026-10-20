import numpy as np

def FlowSpeedRatio(Q_now, N_now, Q_design, N_design) -> float:
    return (Q_now/N_now) / (Q_design/N_design)

# PUMPA pg 10, eq. 35
def SlipFactor(flowspeed_ratio: float, slipfactor_opt: float = 1) -> float:
    factor: float = 1.534988 - 0.6681668 * flowspeed_ratio + 0.077472 * flowspeed_ratio**2 + 0.0571508 * flowspeed_ratio**3
    return slipfactor_opt * factor

def HydraulicEfficiency(f: float) -> float:
    return 0.86387 + 0.3096 * f - 0.14086 * f**2 - 0.029265 * f**3

def Kv():
    Ns = np.array([600, 800, 1000, 1500, 2000, 3000, 4000, 5000, 6000, 7000, 8000, 10000, 15000])
    K = np.array([0.515, 0.47, 0.438, 0.382, 0.347, 0.304, 0.277, 0.259, 0.243, 0.233, 0.221, 0.206, 0.179])
    b, log_a = np.polyfit(np.log10(Ns), np.log10(K), 1)
    a = 10**log_a
    return a, b