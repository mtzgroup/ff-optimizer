import os
from pathlib import Path
from shutil import rmtree
from time import sleep

import numpy as np
import yaml
from scipy.io import netcdf_file


def rmrf(f: str | Path):
    """
    Recursively remove a file or directory.

    Args:
        f (str or Path): Path to file or directory to remove.

    Raises:
        OSError: If the directory cannot be removed after multiple attempts.
    """
    maxTries = 100
    f = Path(f)
    i = 0
    while f.exists() and i < maxTries:
        try:
            rmtree(f)
        except:
            sleep(0.1)
            i += 1
    if f.exists():
        raise OSError(f"Could not remove directory {f}")


def checkForAmber(raiseException: bool = True) -> bool:
    """
    Check if Amber is available in the environment path.

    Args:
        raiseException (bool, optional): Whether to raise an exception if Amber is not found. Defaults to True.

    Returns:
        bool: True if Amber is available, False otherwise.

    Raises:
        RuntimeError: If Amber is not available and raiseException is True.
    """
    value = os.getenv("AMBERHOME")
    if value is None:
        if raiseException:
            raise RuntimeError("Amber is not available in the environment path!")
        print("No Amber available!")
        return False
    return True


def convertTCtoFB(
    tcout: str,
    coors: str,
    stride: int,
    start: int | None = None,
    end: int | None = None,
    qdata: str = "qdata.txt",
    mdcrd: str = "all.mdcrd",
) -> int:
    """
    Convert TeraChem output to ForceBalance format.

    Args:
        tcout (str): Path to TeraChem output file.
        coors (str): Path to coordinates file.
        stride (int): Stride for frame selection.
        start (int, optional): Start frame. Defaults to None.
        end (int, optional): End frame. Defaults to None.
        qdata (str, optional): Output qdata filename. Defaults to "qdata.txt".
        mdcrd (str, optional): Output mdcrd filename. Defaults to "all.mdcrd".

    Returns:
        int: Number of jobs processed.
    """

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
        if (
            gradIndex[indices[i][1]] >= start + stride - 1
            and gradIndex[indices[i][1]] <= end
        ):
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
    """
    Read coordinates from PDB file.

    Args:
        pdb (str): Path to PDB file.

    Returns:
        numpy.ndarray: Array of coordinates.
    """
    coords = []
    with open(pdb, "r") as f:
        for line in f.readlines():
            if line.startswith("ATOM") or line.startswith("HETATM"):
                coords.append(line[30:38].replace(" ", ""))
                coords.append(line[38:46].replace(" ", ""))
                coords.append(line[46:54].replace(" ", ""))
    return np.asarray(coords, dtype=np.float32)


def readXYZ(
    xyz: str, readSymbols: bool = False
) -> tuple[np.ndarray, np.ndarray] | np.ndarray:
    """
    Read XYZ file and return geometry and optionally symbols.

    Args:
        xyz (str): Path to XYZ file.
        readSymbols (bool, optional): Whether to read atomic symbols. Defaults to False.

    Returns:
        tuple or np.ndarray: Geometry (and symbols if readSymbols is True).

    Raises:
        RuntimeError: If the XYZ file is not formatted correctly.
    """
    try:
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
    except:
        raise RuntimeError(f"XYZ file {xyz} is formatted incorrectly!")


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
def convertPDBtoXYZ(pdb: str) -> str:
    """
    Convert PDB file to XYZ format.

    Args:
        pdb (str): Path to PDB file.

    Returns:
        str: Path to the created XYZ file.
    """
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


def readGradFromTCout(outFile: str) -> tuple[float, np.ndarray]:
    """
    Read gradient from TeraChem output file.

    Args:
        outFile (str): Path to TeraChem output file.

    Returns:
        tuple: Energy and gradients.

    Raises:
        RuntimeError: If reading TCCloud format gradient file fails.
    """
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


def readOpt(filename: str) -> tuple[int, dict]:
    """
    Read optimization results from ForceBalance output file.

    Args:
        filename (str): Path to ForceBalance output file.

    Returns:
        tuple: Status and results dictionary.
    """
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


def readEsp(filename: str) -> tuple[list, list]:
    """
    Read ESP data from file.

    Args:
        filename (str): Path to ESP file.

    Returns:
        tuple: ESP XYZ coordinates and ESP values.
    """
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
def writeRst(frame: list, natoms: int, dest: str):
    """
    Write restart file.

    Args:
        frame (list): 2D list containing frame data.
        natoms (int): Number of atoms.
        dest (str): Destination file path.
    """
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


def writePDB(
    geometry: np.ndarray,
    dest: str,
    atoms: list | None = None,
    resname: str | None = None,
    template: str | None = None,
):
    """
    Write PDB file.

    Args:
        geometry (np.ndarray): Geometry data.
        dest (str): Destination file path.
        atoms (list, optional): Atom names. Defaults to None.
        resname (str, optional): Residue name. Defaults to None.
        template (str, optional): Template PDB file. Defaults to None.
    """
    if template is None:
        assert atoms is not None
        assert resname is not None
        atomTypes = {}
        atomNames = []
        for atom in atoms:
            setAtom = False
            for key in atomTypes.keys():
                if atom == key:
                    atomTypes[atom] += 1
                    token = atom + "%0" + str(3 - len(atom)) + "d"
                    atomNames.append(token % atomTypes[atom])
                    setAtom = True
            if setAtom == False:
                atomTypes[atom] = 1
                token = atom + "%0" + str(3 - len(atom)) + "d"
                atomNames.append(token % atomTypes[atom])
        with open(dest, "w") as f:
            f.write("COMPND   " + resname + "\n")
            for i in range(len(atoms)):
                f.write(
                    "ATOM  %5d  %3s %3s     1    %8.3f%8.3f%8.3f  1.00  0.00          %2s\n"
                    % (
                        i + 1,
                        atomNames[i],
                        resname,
                        float(geometry[3 * i]),
                        float(geometry[3 * i + 1]),
                        float(geometry[3 * i + 2]),
                        atoms[i],
                    )
                )

    else:
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


def convertNCtoXYZs(nc: str, symbols: list, offset: int = 0) -> int:
    """
    Convert NetCDF file to XYZ files.

    Args:
        nc (str): Path to NetCDF file.
        symbols (list): List of atomic symbols.
        offset (int, optional): Frame offset. Defaults to 0.

    Returns:
        int: Number of frames converted.
    """
    f = netcdf_file(nc, "r", mmap=False)
    try:
        coords = f.variables["coordinates"]
    except:
        # if no coords, must be vels file
        return 0

    # Count number of real atoms
    natoms = 0
    for symbol in symbols:
        if symbol != "X":
            natoms += 1
    for i in range(coords.shape[0]):
        with open(f"{i+1+offset}.xyz", "w") as f:
            f.write(f"{natoms}\n")
            f.write(f"Converted from {nc}, frame {i}\n")
            for j in range(coords.shape[1]):
                # skip not atoms
                if symbols[j] != "X":
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


def getSymbolsFromPrmtop(prmtop: str) -> list:
    """
    Get atomic symbols from AMBER prmtop file.

    Args:
        prmtop (str): Path to AMBER prmtop file.

    Returns:
        list: List of atomic symbols.
    """
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
        # Some MM simulations use fake particles (e.g. OPC water model)
        # Have to label those as not atoms
        if int(number) > 0:
            symbols.append(elementsByNumber[int(number)])
        else:
            symbols.append("X")
    return symbols


def getXYZs(folder: str | Path = ".") -> list:
    """
    Get list of XYZ files in a folder.

    Args:
        folder (str or Path, optional): Path to folder. Defaults to current directory.

    Returns:
        list: Sorted list of XYZ file paths.
    """
    if type(folder) == str:
        folder = Path(folder)
    xyzs = []
    for f in folder.iterdir():
        if f.name.endswith(".xyz") and not f.name.startswith("esp"):
            xyzs.append(f)
    return sorted(xyzs)


# Separate file name from extension
def getName(f: str | Path) -> str:
    """
    Separate file name from extension.

    Args:
        f (str or Path): File path.

    Returns:
        str: File name without extension.
    """
    if type(f) == str:
        return f.split(".")[0]
    else:
        return f.name.split(".")[0]
