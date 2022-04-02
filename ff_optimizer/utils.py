import numpy as np
from qcelemental.models import Molecule
from . import units


def readPDB(pdb: str):
    coords = []
    with open(pdb, "r") as f:
        for line in f.readlines():
            if line.startswith("ATOM") or line.startswith("HETATM"):
                splitLine = line.split()
                coords.append(splitLine[5])
                coords.append(splitLine[6])
                coords.append(splitLine[7])
    return coords


def convertPDBtoMolecule(pdb: str):
    coords = []
    symbols = []
    with open(pdb, "r") as f:
        for line in f.readlines():
            if line.startswith("ATOM") or line.startswith("HETATM"):
                splitLine = line.split()
                coords.append([splitLine[5], splitLine[6], splitLine[7]])
                symbols.append(splitLine[-1])
    coords = np.asarray(coords, dtype=np.float32)
    # Molecule class by default has coordinates in bohr
    coords *= units.ANGSTROM_TO_AU
    mol = Molecule(**{"symbols": symbols, "geometry": coords})
    return mol


def convertPDBtoXYZ(pdb):
    name = pdb.replace(".pdb", ".xyz")
    xyzLines = []
    with open(pdb, "r") as f:
        for line in f.readlines():
            splitLine = line.split()
            if len(splitLine) == 0:
                continue
            if splitLine[0] == "ATOM" or splitLine[0] == "HETATM":
                xyzLines.append(
                    f"{splitLine[10]}\t{splitLine[5]}\t{splitLine[6]}\t{splitLine[7]}\n"
                )
    with open(name, "w") as f:
        f.write(f"{str(len(xyzLines))}\n")
        f.write("Converted from {pdb}\n")
        for line in xyzLines:
            f.write(line)
    return name


def readGradFromTCout(outFile: str):
    isTCCloud = False
    with open(outFile, "r") as f:
        line = f.readline()
    if "TCCloud gradient output file" in line:
        isTCCloud = True
        energy = line.split()[-1]
    if isTCCloud:
        try:
            grads = np.loadtxt(outFile).flatten()
        except:
            os.remove(outFile)
            raise RuntimeError(f"Reading TCCloud format gradient file {outFile} failed")
    else:
        inGradient = False
        gradCounter = 0
        molSize = 0
        grads = []
        energy = 0
        finished = False
        with open(outFile, "r") as f:
            for line in f.readlines():
                if "Total atoms" in line:
                    molSize = int(line.split()[2])
                if "FINAL ENERGY" in line:
                    energy = float(line.split()[2])
                if "Gradient units" in line:
                    inGradient = True
                if gradCounter < molSize + 3 and gradCounter > 2:
                    for token in line.split():
                        grads.append(token)
                if inGradient:
                    gradCounter = gradCounter + 1
                if "Job finished" in line:
                    finished = True
        if not finished:
            print("File %s did not complete calculation" % outFile)
            grads = None
    return energy, grads


def readOpt(filename):
    inInitialParams = False
    inFinalParams = False
    inFinalObj = False
    params = []
    initialParams = []
    labels = []
    results = {}
    status = -1
    with open(filename, "r") as f:
        for line in f.readlines():
            if "-------" in line:
                inFinalParams = False
                inInitialParams = False
            if inFinalParams:
                if "=======" not in line:
                    splitLine = line.split()
                    params.append(splitLine[2])
                    labels.append(splitLine[5])
            if inInitialParams:
                if "=======" not in line:
                    initialParams.append(line.split()[2])
            if inFinalObj:
                results["obj"] = float(line.split()[5])
                inFinalObj = False
            if "Final physical parameters" in line:
                inFinalParams = True
            if "Starting parameter indices" in line:
                inInitialParams = True
            if "Final objective function" in line:
                inFinalObj = True
            if "Optimization Converged" in line:
                status = 0
            if "Maximum number of optimization steps reached" in line:
                status = 2

    # if status == -1:
    #    raise RuntimeError("ForceBalance optimization of " + filename + " failed")

    params = np.asarray(params, dtype=np.float32)
    initialParams = np.asarray(initialParams, dtype=np.float32)
    results["params"] = params
    results["labels"] = labels
    results["initialParams"] = initialParams
    return status, results


def readEsp(filename):
    lineCounter = 1
    espXYZ = []
    esp = []
    with open(filename, "r") as f:
        for line in f.readlines():
            splitLine = line.split()
            if lineCounter > 2 and len(splitLine) > 0:
                for i in range(1, 4):
                    espXYZ.append(splitLine[i])
                esp.append(splitLine[4])
            lineCounter += 1
    return espXYZ, esp
