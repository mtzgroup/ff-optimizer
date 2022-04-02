import numpy as np


def checkFloat(f1, f2, thresh=0.000001):
    f1 = float(f1)
    f2 = float(f2)
    return abs(f1 - f2) < thresh


def checkArray(a1, a2, thresh=0.000001):
    a1 = np.asarray(a1, dtype=np.float32)
    a2 = np.asarray(a2, dtype=np.float32)
    if a1.shape != a2.shape:
        return False
    return np.max(np.abs(a1 - a2)) < thresh
