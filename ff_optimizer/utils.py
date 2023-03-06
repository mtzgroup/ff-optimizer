import os

import numpy as np
import yaml
from qcelemental.models import Molecule
from scipy.io import netcdf_file

from . import units


def convertTCtoFB(
    tcout, coors, stride, start=None, end=None, qdata="qdata.txt", mdcrd="all.mdcrd"
):

    molSize = 0
    coords = []
    frame = []
    energies = []
    grads = []
    coordIndex = []
    gradIndex = []
    lineCounter = 0
    frameStart = -1
    coorEnergies = []

    with open(coors, "r") as f:
        maxFrame = -1
        for line in f.readlines():
            splitLine = line.split()
            if lineCounter == 0:
                molSize = int(splitLine[0])
            if lineCounter == 1:
                index = int(splitLine[2]) + 1  # coor files are 0-indexed
                energy = splitLine[0]
                if frameStart == -1:
                    frameStart = index
            elif lineCounter > 1 and lineCounter < molSize + 2:
                for token in splitLine[1:]:
                    frame.append(token)
            if lineCounter == molSize + 1:
                if index > maxFrame:
                    maxFrame = index
                    coords.append(frame)
                    coordIndex.append(index)
                    coorEnergies.append(energy)
                else:
                    j = len(coords) - 1
                    while coordIndex[j] > index:
                        j -= 1
                    coords[j] = frame
                    coordIndex[j] = index
                    coorEnergies[j] = energy
                frame = []
                lineCounter = -1
            lineCounter = lineCounter + 1

    gradCounter = molSize + 37
    index = 0
    frameStart = -1

    with open(tcout, "r") as f:
        maxFrame = -1
        for line in f.readlines():
            if "MD STEP" in line:
                index = int(line.split()[4])
                if frameStart == -1:
                    if len(grads) > 0:
                        frameStart = index - 1  # first step is unnumbered "step 0"
                        gradIndex[0] = index - 1
                    else:
                        frameStart = index
            if "FINAL ENERGY" in line:
                energy = float(line.split()[2])
            if "Gradient units" in line:
                gradCounter = 0
                frame = []
            if gradCounter < molSize + 2 and gradCounter > 2:
                for token in line.split():
                    frame.append(token)
            if gradCounter == molSize + 2:
                for token in line.split():
                    frame.append(token)
                if index > maxFrame:
                    maxFrame = index
                    grads.append(frame)
                    gradIndex.append(index)
                    energies.append(energy)
                else:
                    j = len(grads) - 1
                    while gradIndex[j] > index:
                        j -= 1
                    grads[j] = frame
                    gradIndex[j] = index
                    energies[j] = energy
            gradCounter = gradCounter + 1

    if start == None:
        start = gradIndex[0]
    if end == None:
        end = gradIndex[-1]

    jobCounter = 0
    indices = []
    precision = len(coorEnergies[0].split(".")[1])
    for i in range(len(gradIndex)):
        for j in range(len(coordIndex)):
            if gradIndex[i] == coordIndex[j]:
                eFormat = "%." + str(precision) + "f"
                gradEnergy = eFormat % energies[i]
                if coorEnergies[j] != gradEnergy:
                    # print("Mismatched energies in step " + str(gradIndex[i]))
                    # raise RuntimeError("Mismatched energies from " + tcout + " and " + coors)
                    break
                indices.append([i, j])

    usedIndices = []
    lastFrame = -stride - 37
    for i in range(len(indices)):
    # Get first frame from timestep stride, not timestep 0
        if gradIndex[indices[i][1]] >= start + stride - 1 and gradIndex[indices[i][1]] <= end:
            if gradIndex[indices[i][1]] - lastFrame >= stride:
                lastFrame = gradIndex[indices[i][1]]
                usedIndices.append(indices[i])

    with open(qdata, "w") as f:
        for j in usedIndices:
            f.write("\n")
            f.write("JOB " + str(jobCounter) + "\n")
            coordLine = "COORDS "
            for coord in coords[j[1]]:
                coordLine = coordLine + str(coord) + " "
            gradLine = "FORCES "
            for grad in grads[j[0]]:
                gradLine = gradLine + str(grad) + " "
            f.write(coordLine + "\n")
            f.write("ENERGY " + str(energies[j[0]]) + "\n")
            f.write(gradLine + "\n")
            jobCounter = jobCounter + 1

    with open(mdcrd, "w") as f:
        f.write("converted from TC with convertTCMDtoFB.py\n")
        for j in usedIndices:
            tokenCounter = 1
            for coord in coords[j[1]]:
                f.write("%8.3f" % float(coord))
                if tokenCounter == 10:
                    f.write("\n")
                    tokenCounter = 1
                else:
                    tokenCounter = tokenCounter + 1
            if tokenCounter != 1:
                f.write("\n")
    return jobCounter


# unused
def readPDB(pdb: str):
    coords = []
    with open(pdb, "r") as f:
        for line in f.readlines():
            if line.startswith("ATOM") or line.startswith("HETATM"):
                coords.append(line[30:38].replace(" ", ""))
                coords.append(line[38:46].replace(" ", ""))
                coords.append(line[46:54].replace(" ", ""))
    return np.asarray(coords, dtype=np.float32)


def readXYZ(xyz: str, readSymbols=False):
    if readSymbols:
        symbols = np.loadtxt(xyz, skiprows=2, usecols=(0), dtype=str)
        geometry = np.loadtxt(
            xyz, skiprows=2, usecols=(1, 2, 3), dtype=np.float32
        ).flatten()
        return geometry, symbols
    else:
        return np.loadtxt(
            xyz, skiprows=2, usecols=(1, 2, 3), dtype=np.float32
        ).flatten()


def writeXYZ(geometry: np.array, symbols: list, dest: str):
    natoms = int(geometry.shape[0] / 3)
    with open(dest, "w") as f:
        f.write(f"{natoms}\n")
        f.write("Written by ff_opt active learning\n")
        for i in range(natoms):
            f.write(
                "%3s %14.7f %14.7f %14.7f\n"
                % (
                    symbols[i],
                    geometry[3 * i],
                    geometry[3 * i + 1],
                    geometry[3 * i + 2],
                )
            )


# unused
def convertPDBtoMolecule(pdb: str):
    coords = []
    symbols = []
    with open(pdb, "r") as f:
        for line in f.readlines():
            if line.startswith("ATOM") or line.startswith("HETATM"):
                coords.append(
                    [
                        line[30:38].replace(" ", ""),
                        line[38:46].replace(" ", ""),
                        line[46:54].replace(" ", ""),
                    ]
                )
                symbols.append(line.split()[-1])
    coords = np.asarray(coords, dtype=np.float32)
    # Molecule class by default has coordinates in bohr
    coords *= units.ANGSTROM_TO_AU
    mol = Molecule(**{"symbols": symbols, "geometry": coords})
    return mol


# unused
def convertPDBtoXYZ(pdb: str):
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
            grads = -1
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


# frame is a 2D list
def writeRst(frame, natoms, dest):
    with open(dest, "w") as f:
        f.write(f"Written by ff_optimizer\n")
        f.write(f"{str(natoms)}\n")
        for i in range(len(frame)):
            f.write(
                "%12.7f%12.7f%12.7f"
                % (float(frame[i][0]), float(frame[i][1]), float(frame[i][2]))
            )
            if int(i / 2) * 2 != i:
                f.write("\n")


# unused
def writePDB(geometry, dest, template):
    with open(template, "r") as f:
        templateLines = f.readlines()

    i = 0
    with open(dest, "w") as f:
        for line in templateLines:
            if i > len(geometry):
                break
            if line.startswith("ATOM") or line.startswith("HETATM"):
                f.write(
                    "%s%8.3f%8.3f%8.3f%s"
                    % (
                        line[:30],
                        geometry[i],
                        geometry[i + 1],
                        geometry[i + 2],
                        line[54:],
                    )
                )
                i += 3

            else:
                f.write(line)


def convertNCtoXYZs(nc, symbols, offset=0):
    f = netcdf_file(nc, "r", mmap=False)
    try:
        coords = f.variables["coordinates"]
    except:
        # if no coords, must be vels file
        return 0
    natoms = coords.shape[1]
    for i in range(coords.shape[0]):
        with open(f"{i+1+offset}.xyz", "w") as f:
            f.write(f"{natoms}\n")
            f.write("Converted from {nc}, frame {i}\n")
            for j in range(natoms):
                f.write(
                    "%3s %14.7f %14.7f %14.7f\n"
                    % (
                        symbols[j],
                        coords[i, j, 0],
                        coords[i, j, 1],
                        coords[i, j, 2],
                    )
                )
    return coords.shape[0]


def getSymbolsFromPrmtop(prmtop):
    with open(os.path.join(os.path.dirname(__file__), "elements.yaml"), "r") as f:
        elements = yaml.safe_load(f)
    elementsByNumber = {}
    for element in elements.keys():
        elementsByNumber[elements[element]["atomic_number"]] = element

    atomicNumbers = []
    inAtomicNumbers = False
    almostInAtomicNumbers = False
    with open(prmtop, "r") as f:
        for line in f.readlines():
            if inAtomicNumbers:
                if "FLAG" in line:
                    break
                atomicNumbers += line.split()
            if almostInAtomicNumbers:
                inAtomicNumbers = True
                almostInAtomicNumbers = False
            if "FLAG ATOMIC_NUMBER" in line:
                almostInAtomicNumbers = True

    symbols = []
    for number in atomicNumbers:
        symbols.append(elementsByNumber[int(number)])
    return symbols
