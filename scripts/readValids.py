import os
import numpy as np
from pathlib import Path
from shutil import copyfile

def readValid(path, reweight = False):
    targets = 0
    with open(path, 'r') as f:
        for line in f.readlines():
            if "Setup for target" in line:
                targets += 1
            splitLine = line.split()
            if len(splitLine) > 0:
                if splitLine[0] == "valid_1":
                    v1 = float(splitLine[1])
            if "Objective Function Single Point" in line:
                obj = float(line.split()[6])
                if reweight:
                    obj *= targets
                    obj -= v1
                    obj /= targets - 1
                return obj
    raise RuntimeError(f"{path} failed")

def readValids(suffix, optdir="1_opt", nvalids=10): 
    valids = []
    cycle = 1
    while True:
        validTemp = []
        for i in range(nvalids):
            validFile = f"valid_{cycle}_{i}{suffix}.out"
            if i == 0:
                validFile = f"valid_{cycle}{suffix}.out"
            validPath = os.path.join(optdir, validFile)
            try:
                if i == 0:
                    validTemp.append(readValid(validPath))
                else:
                    validTemp.append(readValid(validPath, reweight=True))
            except:
                print(f"read of {validPath} failed")
        if len(validTemp) < 10:
            break
        else:   
            valids.append(validTemp)
        cycle += 1
    return np.asarray(valids)

def getStoppingPoint(vj):
    patience = 5
    inPatience = False
    cutoff = -1
    for j in range(vj.shape[0]):
        if not inPatience and vj[j] > cutoff:
            inPatience = True
            patienceCycle = j
        if inPatience and vj[j] < cutoff:
            inPatience = False
        if inPatience and j - patienceCycle > patience:
            lastCycle = j
            break
    return lastCycle

def runValid(optdir, inp, out, parmFolder, frcmod, mol2):
    home = os.getcwd()
    optdir = Path(optdir)
    parmFolder = Path(parmFolder)
    copyfile(parmFolder / frcmod, optdir / "forcefield" / frcmod)
    copyfile(parmFolder / mol2, optdir / "forcefield" / mol2)
    os.chdir(optdir)
    os.system(f"ForceBalance.py {inp} > {out}")
    v = readValid(out)
    os.chdir(home)
    return v

def getFinalValidations(i, lastCycle, optdir="1_opt", patience=5, frcmod="dasa.frcmod", mol2="dasa.mol2"):
    optdir = Path(optdir)
    vs = []
    # only check the parameters in the patience region
    for j in range(lastCycle-patience, lastCycle + 1):
        if i == 0:
            inp = f"valid_{lastCycle}.in"
        else:
            inp = f"valid_{lastCycle}_{i}.in"
        out = f"valid_{j}_{i}_final.out"
        try:
            vs.append(readValid(optdir / out))
        except:
            runValid(optdir, inp, out, f"{optdir}/result/opt_{j}", frcmod, mol2)
            vs.append(readValid(optdir / out))
    vs = np.asarray(vs)
    print(vs)
    return np.argmin(vs) + 1 + lastCycle - patience # account for one-indexing of validations
        
def readValidDiff(optdir="1_opt"):
    valid = readValids("", optdir)
    previous = readValids("_previous", optdir)
    ncycles = min(valid.shape[0], previous.shape[0])
    validPrev = valid[:ncycles,:] - previous[:ncycles,:]
    vj = validPrev / valid[:ncycles,:] * 100
    return vj

vj = readValidDiff()
for i in range(vj.shape[1]):
    finalCycle = getStoppingPoint(vj[:,i])
    best = getFinalValidations(i, finalCycle)
    print(best)
#np.savetxt("validPrev.txt", validPrev, fmt="%10.3e")
