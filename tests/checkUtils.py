import numpy as np


def checkFloat(f1, f2, thresh=0.000001):
    f1 = float(f1)
    f2 = float(f2)
    return abs(f1 - f2) < thresh


def checkArray(a1, a2, thresh=0.000001):
    a1 = np.asarray(a1, dtype=np.float32)
    a2 = np.asarray(a2, dtype=np.float32)
    if a1.shape != a2.shape:
        print(f"a1.shape: {a1.shape}")
        print(f"a2.shape: {a2.shape}")
        return False
    return np.max(np.abs(a1 - a2)) < thresh

def checkLists(list1, list2):
    if len(list1) != len(list2):
        print(len(list1), len(list2))
        return False
    for i in range(len(list1)):
        if list1[i] != list2[i]:
            print(f"List 1: {list1[i]}")
            print(f"List 2: {list2[i]}")
            return False
    return True
