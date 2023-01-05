import numpy as np


def checkFloat(f1, f2, thresh=0.00001):
    f1 = float(f1)
    f2 = float(f2)
    return abs(f1 - f2) < thresh


def checkArray(a1, a2, thresh=0.00001):
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


def checkListsFloats(list1, list2, thresh=0.000001):
    if len(list1) != len(list2):
        print(len(list1), len(list2))
        return False
    for i in range(len(list1)):
        try:
            isSame = checkFloat(float(list1[i]), float(list2[i]), thresh=thresh)
        except:
            isSame = list1[i] == list2[i]
        if not isSame:
            print(f"List 1: {list1[i]}")
            print(f"List 2: {list2[i]}")
            return False
    return True


def checkFiles(f1, f2):
    with open(f1, "r") as f:
        lines1 = f.readlines()
    with open(f2, "r") as f:
        lines2 = f.readlines()
    return checkLists(lines1, lines2)


def checkFileFloatsNoWhitespace(f1, f2, thresh=0.000001):
    with open(f1, "r") as f:
        lines1 = f.readlines()
    with open(f2, "r") as f:
        lines2 = f.readlines()
    if len(lines1) != len(lines2):
        return False
    for i in range(len(lines1)):
        isSame = checkListsFloats(lines1[i].split(), lines2[i].split(), thresh=thresh)
        if not isSame:
            return False
    return True
