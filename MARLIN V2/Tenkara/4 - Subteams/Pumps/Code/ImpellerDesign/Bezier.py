import numpy as np


# Cubic Bezier
def bez(P: np.ndarray, t: np.ndarray | list[float]) -> np.ndarray:
    t = np.asarray(t, float)[:, None]
    return ((1 - t) ** 3) * P[0] + 3 * ((1 - t) ** 2) * t * P[1] \
        + 3 * (1 - t) * t ** 2 * P[2] + t ** 3 * P[3]


# Derivative of cubic Bezier
def dbez(P: np.ndarray, t: np.ndarray | list[float]) -> np.ndarray:
    t = np.asarray(t, float)[:, None]
    return 3 * ((1 - t) ** 2 * (P[1] - P[0]) + 2 * (1 - t) * t * (P[2] - P[1])
                + t ** 2 * (P[3] - P[2]))