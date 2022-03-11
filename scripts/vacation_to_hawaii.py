#!/usr/bin/env python

import os
import argparse
from random import randint
import matplotlib as mpl
import numpy as np
from ff_optimizer import qmengine
from shutil import rmtree

mpl.use("Agg")
import matplotlib.pyplot as plt
from time import sleep
from textwrap import dedent

# Some helper functions
def die():
    raise RuntimeError("die here")


def addTargetLines(inputFile, targetLines, initialTarget, newTarget):
    addedLines = False
    with open(inputFile, "r") as f:
        for line in f.readlines():
            splitLine = line.split()
            if len(splitLine) > 1:
                if splitLine[0] == "name" and splitLine[1] == newTarget:
                    addedLines = True
    if not addedLines:
        with open(inputFile, "a") as f:
            f.write("\n")
            for line in targetLines:
                f.write(line.replace(initialTarget, newTarget))


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
    lastFrame = -args.stride - 37
    for i in range(len(indices)):
        if gradIndex[indices[i][1]] >= start and gradIndex[indices[i][1]] <= end:
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


def changeParameter(inputFile, prmName, prmValue):
    changed = False
    lines = []
    with open(inputFile, "r") as f:
        for line in f.readlines():
            if prmName in line:
                line = line.replace(line.split()[1], prmValue)
                changed = True
            lines.append(line)
        if not changed:
            lines.insert(1, f"{prmName} {prmValue}")
            pass

    with open("temp.txt", "w") as f:
        for line in lines:
            f.write(line)
    os.system(f"mv temp.txt {inputFile}")


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


def determineAdaptiveDamping(
    testFile, upperThreshold=0.3, lowerThreshold=0.01, adaptiveDamping=0.5
):
    changeParameter(testFile, "adaptive_damping", str(adaptiveDamping))
    maxCycles = 100
    testOut = testFile.split(".")[0] + ".out"
    for j in range(maxCycles):
        os.system(f"ForceBalance.py {testFile} > {testOut}")
        status, results = readOpt(testOut)
        diff = np.abs(results["params"] - results["initialParams"]) / np.maximum(
            results["params"], results["initialParams"]
        )
        print(diff)
        if np.argwhere(diff > upperThreshold).shape[0] > 0:
            adaptiveDamping *= 2
            changeParameter(testFile, "adaptive_damping", str(adaptiveDamping))
        elif np.argwhere(diff > lowerThreshold).shape[0] == 0:
            adaptiveDamping *= 0.75
            changeParameter(testFile, "adaptive_damping", str(adaptiveDamping))
        else:
            return adaptiveDamping




def readRst(filename):
    coords = []
    lineCounter = 0
    with open(filename, "r") as f:
        for line in f.readlines():
            if lineCounter > 1:
                for coord in line.split():
                    coords.append(coord)
            lineCounter = lineCounter + 1
    return coords


def readValid(filename):
    with open(filename, "r") as f:
        for line in f.readlines():
            if "Objective Function Single Point" in line:
                return float(line.split()[6])
        raise RuntimeError(
            "ForceBalance single-point evaluation of " + filename + " did not converge"
        )


def writeRst(frame, index, molSize, dest):
    fileName = os.path.join(dest, str(index) + ".rst7")
    with open(fileName, "w") as f:
        f.write("sampled from frame " + str(index) + "\n")
        f.write(str(molSize) + "\n")
        for i in range(len(frame)):
            f.write(
                "%12.7f%12.7f%12.7f"
                % (float(frame[i][0]), float(frame[i][1]), float(frame[i][2]))
            )
            if int(i / 2) * 2 != i:
                f.write("\n")


def getFrame(index, coordFile, dest):
    frame = []
    atoms = []
    inFrame = False
    lineCounter = 1
    frameCounter = 1
    with open(coordFile, "r") as f:
        molSize = int(f.readline())
        for line in f.readlines():
            if inFrame:
                atoms.append(line.split()[0])
                frame.append(line.split()[1:])
                lineCounter = lineCounter + 1
            if "frame " in line:
                if frameCounter == index:
                    inFrame = True
                frameCounter += 1
            if lineCounter > molSize:
                inFrame = False
    writeRst(frame, index, molSize, dest)


# Summary stuff
summary = dedent(
    """\
Summary of necessary directories and files
Dynamics directory (--dynamicsdir, 0_dynamics) contains:
    -- QM dynamics coordinates (--coors, coors.xyz)
    -- TeraChem output file (--tcout, tc.out)
FB Optimization directory (--optdir, 1_opt) contains:
    -- FB optimization input file (--opt0, opt_0.in)
    -- FB validation input file (--valid0, valid_0.in)
    -- PDB conformation file (conf.pdb)
    -- Tleap setup file (setup.leap)
    -- Parameter frcmod file (*.frcmod)
    -- Parameter mol2 file (*.mol2)
MM sampling directory (--sampledir, 2_sampling) contains:
    -- Amber equilibration inputs (heat*in)
    -- Amber MM sampling input (md.in)
    -- Cpptraj input file (cpptraj.in)
    -- TeraChem input file for fire (tc_template.in)
    -- Backup TeraChem input file if job fails (tc_template_long.in)
    -- sbatch template file for fire (sbatch_template.sh)
"""
)
parser = argparse.ArgumentParser(
    epilog=summary, formatter_class=argparse.RawDescriptionHelpFormatter
)
parser.add_argument(
    "--dynamicsdir",
    help="Folder containing TC output file from MD trajectory and XYZ coord file from MD trajectory, default is 0_dynamics",
    type=str,
    default="0_dynamics",
)
parser.add_argument(
    "--tcout",
    help="TeraChem output file from MD trajectory, default is tc.out",
    type=str,
    default="tc.out",
)
parser.add_argument(
    "--coors",
    help="XYZ file containing coords from TeraChem MD trajectory, default is coors.xyz",
    type=str,
    default="coors.xyz",
)
parser.add_argument(
    "--start",
    help="Lower bound for sampling of frames from MD trajectory, default is the first frame",
    type=int,
    default=None,
)
parser.add_argument(
    "--end",
    help="Upper bound for sampling of frames from MD trajectory, default is the last frame",
    type=int,
    default=None,
)
parser.add_argument(
    "--stride",
    help="Stride between successive samples of frames from MD trajectory, default is 50",
    type=int,
    default=50,
)
parser.add_argument(
    "--optdir",
    help="Directory where ForceBalance optimization is performed, default is 1_opt",
    type=str,
    default="1_opt",
)
parser.add_argument(
    "--maxcycles",
    help="Maximum number of optimization/sampling cycles to be performed, default is 10",
    type=int,
    default=10,
)
parser.add_argument(
    "--sampledir",
    help="Directory where MM sampling is performed, default is 2_sampling",
    type=str,
    default="2_sampling",
)
parser.add_argument(
    "--split",
    help="If set, sample one geometry from before this frame and one from after",
    type=int,
    default=None,
)
parser.add_argument(
    "--opt0",
    help="ForceBalance input file for first optimization, with a single target named dynamics, default is opt_0.in",
    type=str,
    default="opt_0.in",
)
parser.add_argument(
    "--valid0",
    help="Template ($options section only) for ForceBalance input files for validation single points, default is valid_0.in",
    type=str,
    default="valid_0.in",
)
parser.add_argument(
    "--engine",
    help="Engine for performing QM calculations, either queue, debug, or tccloud",
    type=str,
    default="queue",
)
parser.add_argument(
    "--sbatch",
    help="Sbatch template file for submitting jobs to fire, default is sbatch_template.sh",
    type=str,
    default="sbatch_template.sh",
)
parser.add_argument(
    "--tctemplate",
    help="TC input template file for submitting jobs to fire, default is tc_template.in",
    type=str,
    default="tc_template.in",
)
parser.add_argument(
    "--tctemplate_long",
    help="TC input template file for resubmitting failed jobs to fire, default is tc_template_long.in",
    type=str,
    default="tc_template_long.in",
)
parser.add_argument(
    "--restart",
    help="Restart partially complete forcefield optimization",
    action="store_true",
)
parser.add_argument(
    "--resp",
    help="Restart partially complete forcefield optimization",
    action="store_true",
    default=False,
)
args = parser.parse_args()

# Set up the calculation
train = []
valid = []
validInitial = []
validPrevious = []

# Determine cycle for restart, set restart variables
restartCycle = -1
if args.restart:
    types = ["BONDSK", "BONDSB", "ANGLESK", "ANGLESB", "DIHS"]
    aliases = [
        "Bond strength",
        "Bond length",
        "Angle strength",
        "Equilibrium angle",
        "Dihedral strength",
    ]
    colors = [
        "blue",
        "green",
        "firebrick",
        "goldenrod",
        "orange",
        "purple",
        "lightskyblue",
        "olive",
    ]
    for i in range(args.maxcycles + 2):
        optOutput = os.path.join(args.optdir, "opt_" + str(i) + ".out")
        if os.path.isfile(optOutput):
            status, results = readOpt(optOutput)
            if status == 0:
                train.append(results["obj"])
                if i == 0:
                    params = np.zeros((args.maxcycles + 2, len(results["labels"])))
                    labels = results["labels"]
                    params[0, :] = results["initialParams"]
                else:
                    try:
                        valid.append(
                            readValid(
                                os.path.join(args.optdir, "valid_" + str(i) + ".out")
                            )
                        )
                        validPrevious.append(
                            readValid(
                                os.path.join(
                                    args.optdir, "valid_" + str(i) + "_previous.out"
                                )
                            )
                        )
                        validInitial.append(
                            readValid(
                                os.path.join(
                                    args.optdir, "valid_" + str(i) + "_initial.out"
                                )
                            )
                        )
                    except:
                        print("valids didn't complete")
                        # Chop valid lists to the size before trying to append
                        valid = valid[:i]
                        validPrevious = validPrevious[:i]
                        validInitial = validInitial[:i]
                        break
                for j in range(len(results["labels"])):
                    if labels[j] == results["labels"][j]:
                        params[i + 1, j] = results["params"][j]
                    else:
                        for k in range(j + 1, len(labels)):
                            if labels[k] == results["labels"][j]:
                                params[i + 1, k] = results["params"][j]
                                break
            else:
                break
        else:
            break
    restartCycle = i - 1
    print("Restarting optimization at cycle " + str(restartCycle + 1))

# Check for necessary folders and files
if restartCycle < 0:
    if not os.path.isdir(args.dynamicsdir):
        raise RuntimeError("Dynamics directory " + args.dynamicsdir + " does not exist")
    if not os.path.isfile(os.path.join(args.dynamicsdir, args.coors)):
        raise RuntimeError(
            "XYZ coordinates (from QM dynamics) "
            + args.coors
            + " does not exist in "
            + args.dynamicsdir
        )
    if not os.path.isfile(os.path.join(args.dynamicsdir, args.tcout)):
        raise RuntimeError(
            "TC output file (from QM dynamics "
            + args.tcout
            + " does not exist in "
            + args.dynamicsdir
        )
    if not os.path.isdir(args.optdir):
        raise RuntimeError(
            "ForceBalance optimization directory " + args.optdir + " does not exist"
        )
    if not os.path.isfile(os.path.join(args.optdir, "conf.pdb")):
        raise RuntimeError(
            "Prototype PDB coordinates conf.pdb does not exist in " + args.optdir
        )
    if not os.path.isfile(os.path.join(args.optdir, "setup.leap")):
        raise RuntimeError(
            "Tleap input file setup.leap does not exist in " + args.optdir
        )
    if not os.path.isfile(os.path.join(args.optdir, args.opt0)):
        raise RuntimeError(
            "Initial ForceBalance optimization input file "
            + args.opt0
            + " does not exist in "
            + args.optdir
        )
    if not os.path.isfile(os.path.join(args.optdir, args.valid0)):
        raise RuntimeError(
            "Initial ForceBalance validation input file "
            + args.valid0
            + " does not exist in "
            + args.optdir
        )
    if not os.path.isdir(args.sampledir):
        raise RuntimeError(
            "MM sampling directory " + args.sampledir + " does not exist"
        )
    if not os.path.isfile(os.path.join(args.sampledir, "heat1.in")):
        raise RuntimeError(
            "No sander input file for equilibration named heat1.in provided in "
            + args.sampledir
        )
    if not os.path.isfile(os.path.join(args.sampledir, "md.in")):
        raise RuntimeError(
            "No sander input file for sampling named md.in provided in "
            + args.sampledir
        )
    if not os.path.isfile(os.path.join(args.sampledir, "cpptraj.in")):
        raise RuntimeError("No cpptraj input file provided in " + args.sampledir)
    if args.engine != "queue" and args.engine != "debug" and args.engine != "tccloud":
        raise RuntimeError("Engine " + args.engine + " is not implemented")
    if args.engine == "queue":
        if not os.path.isfile(os.path.join(args.sampledir, args.sbatch)):
            raise RuntimeError(
                "Sbatch template "
                + args.sbatch
                + " does not exist in "
                + args.sampledir
            )
        if not os.path.isfile(os.path.join(args.sampledir, args.tctemplate)):
            raise RuntimeError(
                "TC input template "
                + args.tctemplate
                + " does not exist in "
                + args.sampledir
            )
        if not os.path.isfile(os.path.join(args.sampledir, args.tctemplate_long)):
            raise RuntimeError(
                "TC input template (for resubmitting jobs) "
                + args.tctemplate_long
                + " does not exist in "
                + args.sampledir
            )


# Set some miscellaneous variables
home = os.getcwd()
targetLines = []
validTargetLines = []
validInitialTargetLines = []
inTarget = False
with open(os.path.join(args.optdir, args.opt0), "r") as f:
    for line in f.readlines():
        if len(line.split()) == 0:
            continue
        if "$target" in line:
            inTarget = True
        if inTarget:
            targetLines.append(line)
            if line.split()[0] == "resp":
                continue
            validTargetLines.append(line)
            if line.split()[0] == "amber_leapcmd":
                line = line.replace(line.split()[1], "setup_valid_initial.leap")
            validInitialTargetLines.append(line)
if args.resp:
    respAdded = False
    respWeightAdded = False
    for line in targetLines:
        if line.split()[0] == "resp":
            respAdded = True
        if line.split()[0] == "w_resp":
            respWeightAdded = True
    if not respWeightAdded:
        targetLines.insert(1, "w_resp 1\n")
    if not respAdded:
        targetLines.insert(1, "resp 1\n")
with open(os.path.join(args.optdir, "setup.leap"), "r") as leapRead:
    with open(os.path.join(args.optdir, "setup_valid_initial.leap"), "w") as leapWrite:
        for line in leapRead.readlines():
            if "loadamberparams" in line:
                oldName = line.split()[1]
                newName = "initial_" + oldName
                line = line.replace(oldName, newName)
            if "loadmol2" in line:
                oldName = line.split()[3]
                newName = "initial_" + oldName
                line = line.replace(oldName, newName)
            leapWrite.write(line)
if not os.path.isdir(os.path.join(args.optdir, "forcefield")):
    os.mkdir(os.path.join(args.optdir, "forcefield"))
mol2 = None
frcmod = None
with open(os.path.join(args.optdir, args.opt0), "r") as f:
    for line in f.readlines():
        splitLine = line.split()
        if len(splitLine) > 1:
            if splitLine[0] == "forcefield":
                for i in range(1, 3):
                    if ".mol2" in splitLine[i]:
                        mol2 = splitLine[i]
                    elif ".frcmod" in splitLine[i]:
                        frcmod = splitLine[i]
            if splitLine[0] == "name":
                initialTarget = splitLine[1]
if mol2 == None:
    raise RuntimeError("No mol2 file specified for optimization in " + args.opt0)
if frcmod == None:
    raise RuntimeError("No frcmod file specified for optimization in " + args.opt0)
if restartCycle < 0:
    if not os.path.isfile(os.path.join(args.optdir, mol2)):
        raise RuntimeError(
            "Mol2 " + mol2 + " specified in " + args.opt0 + " is not in " + args.optdir
        )
    if not os.path.isfile(os.path.join(args.optdir, frcmod)):
        raise RuntimeError(
            "Frcmod "
            + frcmod
            + " specified in "
            + args.opt0
            + " is not in "
            + args.optdir
        )
mdFiles = []
for f in os.listdir(args.sampledir):
    if os.path.isfile(os.path.join(args.sampledir,f)):
        mdFiles.append(f)
# Set start, split, end indices for sampling from QM trajectory
coordIndex = []
with open(os.path.join(args.dynamicsdir, args.coors), "r") as f:
    for line in f.readlines():
        if "frame" in line:
            coordIndex.append(int(line.split()[2])+1)
startIndex = 0
if args.start != None:
    while coordIndex[startIndex] < args.start:
        startIndex += 1
endIndex = len(coordIndex) - 1
if args.end != None:
    while coordIndex[endIndex] > args.end:
        endIndex -= 1
splitIndex = 0
if args.split != None:
    while coordIndex[splitIndex] < args.split:
        splitIndex += 1
    if splitIndex == 0:
        raise RuntimeError("There must be frames before " + args.split)
    elif splitIndex >= len(coordIndex):
        raise RuntimeError("There must be frames after " + args.split)

# Make validation input for initial MM parameters
if not os.path.isfile(os.path.join(args.optdir, "valid_0_initial.in")):
    with open(os.path.join(args.optdir, args.valid0), "r") as srcValid:
        with open(os.path.join(args.optdir, "valid_0_initial.in"), "w") as destValid:
            for line in srcValid.readlines():
                destValid.write(
                    line.replace(frcmod, "initial_" + frcmod).replace(
                        mol2, "initial_" + mol2
                    )
                )
                if "$target" in line:
                    break

# Initialize QMEngine
if args.engine == 'debug':
    qmEngine = qmengine.DebugEngine(os.path.join(args.sampledir,args.tctemplate),os.path.join(args.sampledir,args.tctemplate_long),doResp=args.resp)
elif args.engine == 'queue':
    qmEngine = qmengine.SbatchEngine(os.path.join(args.sampledir,args.tctemplate),os.path.join(args.sampledir,args.tctemplate_long),os.path.join(args.sampledir,args.sbatch),os.getenv('USER'),doResp=args.resp)
elif args.engine == 'tccloud':
    qmEngine = qmengine.TCCloudEngine(os.path.join(args.sampledir,args.tctemplate),os.path.join(args.sampledir,args.tctemplate_long))
    if args.resp == True:
        print("WARNING! RESP not implemented for this QM engine")

# First optimization cycle is not necessary if restarting from somewhere later
if restartCycle < 0:

    # Create initial target data from dynamics
    with open(os.path.join(args.optdir, args.opt0)) as f:
        for line in f.readlines():
            splitLine = line.split()
            if len(splitLine) > 1:
                if splitLine[0] == "name":
                    initialTarget = splitLine[1]
    if not os.path.isdir(os.path.join(args.optdir, "targets")):
        os.mkdir(os.path.join(args.optdir, "targets"))
    path = os.path.join(args.optdir, "targets", initialTarget)
    if not os.path.isdir(path):
        os.mkdir(path)
    l = convertTCtoFB(
        os.path.join(args.dynamicsdir, args.tcout),
        os.path.join(args.dynamicsdir, args.coors),
        args.stride,
        args.start,
        args.end,
        os.path.join(path, "qdata.txt"),
        os.path.join(path, "all.mdcrd"),
    )

    src = os.path.join(args.optdir, "setup.leap")
    dest = os.path.join(path, ".")
    os.system(f"cp {src} {dest}")
    src = os.path.join(args.optdir, "conf.pdb")
    os.system(f"cp {src} {dest}")
    src = os.path.join(args.optdir, "setup_valid_initial.leap")
    os.system(f"mv {src} {dest}")
    src = os.path.join(args.optdir, frcmod)
    dest = os.path.join(args.optdir, "forcefield", ".")
    os.system(f"cp {src} {dest}")
    dest = os.path.join(args.optdir, "forcefield", "initial_" + frcmod)
    os.system(f"cp {src} {dest}")
    src = os.path.join(args.optdir, mol2)
    dest = os.path.join(args.optdir, "forcefield", ".")
    os.system(f"cp {src} {dest}")
    dest = os.path.join(args.optdir, "forcefield", "initial_" + mol2)
    os.system(f"cp {src} {dest}")
    if args.opt0 != "opt_0.in":
        src = os.path.join(args.optdir, args.opt0)
        dest = os.path.join(args.optdir, "opt_0.in")
        os.system(f"cp {src} {dest}")

    os.chdir(args.optdir)
    os.system("ForceBalance.py opt_0.in > opt_0.out")
    src = os.path.join("result", "opt_0", "*")
    dest = os.path.join("forcefield", ".")
    os.system(f"cp {src} {dest}")
    os.chdir(home)
    status, results = readOpt(os.path.join(args.optdir, "opt_0.out"))
    if status != 0:
        raise RuntimeError("ForceBalance optimization of opt_0.in failed")
    train.append(results["obj"])
    labels = results["labels"]
    params = np.zeros((args.maxcycles + 2, len(labels)))
    params[0, :] = np.asarray(results["initialParams"])
    for j in range(len(results["labels"])):
        if labels[j] == results["labels"][j]:
            params[1, j] = results["params"][j]
        else:
            for k in range(j + 1, len(labels)):
                if labels[k] == results["labels"][j]:
                    params[1, k] = results["params"][j]
                    break

# Begin sampling/optimization cycling
os.rename(
    os.path.join(args.optdir, args.valid0), os.path.join(args.optdir, "valid_0.in")
)
os.rename(os.path.join(args.optdir, args.opt0), os.path.join(args.optdir, "opt_0.in"))

print("%7s%15s%15s%20s%23s" % ("Epoch","Validation","Valid ratio","Current-Previous","Current-last Current"))
for i in range(1, args.maxcycles + 1):

    if i <= restartCycle:
        continue

    # Make sampling directory and copy files into it
    sampleName = str(i) + "_cycle_" + str(i)
    samplePath = os.path.join(args.sampledir, sampleName)
    if not os.path.isdir(samplePath):
        os.mkdir(samplePath)
    elif restartCycle == -1:
        rmtree(samplePath)
        os.mkdir(samplePath)
    src = os.path.join(args.optdir, "result", "opt_" + str(i - 1), "*")
    dest = os.path.join(samplePath, ".")
    os.system(f"cp {src} {dest}")
    src = os.path.join(args.optdir, "conf.pdb")
    os.system(f"cp {src} {dest}")
    src = os.path.join(args.optdir, "setup.leap")
    os.system(f"cp {src} {dest}")
    for f in mdFiles:
        src = os.path.join(args.sampledir, f)
        os.system(f"cp {src} {dest}")

    # Get new .rst7 files for MM sampling initial conditions
    rstCount = 0
    for f in os.listdir(samplePath):
        if ".rst7" in f:
            rstCount += 1
    if rstCount < 2:
        if args.split != None:
            coor1 = coordIndex[randint(startIndex, splitIndex-1)]
            coor2 = coordIndex[randint(splitIndex, endIndex)]
        else:
            coor1 = coordIndex[randint(startIndex, endIndex)]
            coor2 = coordIndex[randint(startIndex, endIndex)]
        getFrame(coor1, os.path.join(args.dynamicsdir, args.coors), samplePath)
        getFrame(coor2, os.path.join(args.dynamicsdir, args.coors), samplePath)

    # Setup MM dynamics, get some important file names
    os.chdir(samplePath)
    os.system("tleap -f setup.leap > leap.out")
    os.chdir(home)
    files = os.listdir(samplePath)
    rsts = []
    heatCounter = 0
    prmtop = None
    for f in files:
        if ".rst7" in f and "heat" not in f:
            rsts.append(f)
        if "heat" in f and ".in" in f:
            heatCounter += 1
        if ".prmtop" in f:
            prmtop = f
        if ".frcmod" in f:
            frcmod = f
    if prmtop == None:
        raise RuntimeError(
            "Tleap failed to create a new .prmtop file, check "
            + os.path.join(samplePath, "leap.out")
            + " for more information"
        )

    # Run sander
    os.chdir(samplePath)
    for rst in rsts:
        name = rst.split(".")[0]
        if not os.path.isfile(name + "_sample.nc"):
            src = rst
            dest = name + "_heat0.rst7"
            os.system(f"cp {src} {dest}")
            for j in range(1, heatCounter + 1):
                os.system(
                    "sander -O -p "
                    + prmtop
                    + " -i heat"
                    + str(j)
                    + ".in -o "
                    + name
                    + "_heat"
                    + str(j)
                    + ".out -c "
                    + name
                    + "_heat"
                    + str(j - 1)
                    + ".rst7 -x "
                    + name
                    + "_heat"
                    + str(j)
                    + ".nc -r "
                    + name
                    + "_heat"
                    + str(j)
                    + ".rst7"
                )
            os.system(
                "sander -O -p "
                + prmtop
                + " -c "
                + name
                + "_heat"
                + str(heatCounter)
                + ".rst7 -x "
                + name
                + "_sample.nc -v "
                + name
                + "_vel.nc -i md.in -o "
                + name
                + ".out"
            )

        # Set up QM calculations
        calcPath = name
        if not os.path.isdir(calcPath):
            os.mkdir(calcPath)
        src = name + "_sample.nc"
        dest = os.path.join(calcPath, ".")
        os.system(f"cp {src} {dest}")

        with open("../cpptraj.in", "r") as srcIn:
            with open(os.path.join(calcPath, "cpptraj.in"), "w") as destIn:
                for line in srcIn.readlines():
                    destIn.write(line.replace("XXX", name + "_sample"))
        os.chdir(calcPath)
        os.system(f"cpptraj -p ../{prmtop} -i cpptraj.in > cpptraj.out")

        # Run QM calculations
        files = os.listdir()
        pdbs = []
        for f in files:
            if '.pdb' in f:
                if len(f.split(".")) > 2:
                    label = f.split(".")[2]
                    fName = label + '.pdb'
                    os.system(f"mv {f} {fName}")
                    pdbs.append(fName)
                else:
                    pdbs.append(f)
        if i == restartCycle + 1:
            qmEngine.restart('.')
        else:
            qmEngine.getQMRefData(pdbs,'.')
        os.chdir(home)
        os.chdir(samplePath)
    os.chdir(home)

    # Set up new ForceBalance optimization

    # Copy new QM data into appropriate folders
    src = os.path.join(args.optdir, "targets", initialTarget, "*")
    trainFolder = os.path.join(args.optdir, "targets", "train_" + str(i))
    validFolder = os.path.join(args.optdir, "targets", "valid_" + str(i))
    if not os.path.isdir(trainFolder):
        os.mkdir(trainFolder)
    if not os.path.isdir(validFolder):
        os.mkdir(validFolder)
    dest = os.path.join(args.optdir, "targets", "train_" + str(i), ".")
    os.system(f"cp {src} {dest}")
    dest = os.path.join(args.optdir, "targets", "valid_" + str(i), ".")
    os.system(f"cp {src} {dest}")

    rstNumber = int(rsts[0].split(".")[0])
    rstNumber2 = int(rsts[1].split(".")[0])
    if args.split is not None:
        if rstNumber < args.split:
            src = os.path.join(
                args.sampledir, str(i) + "_cycle_" + str(i), str(rstNumber), "all.mdcrd"
            )
            dest = os.path.join(args.optdir, "targets", "train_" + str(i), ".")
            os.system(f"cp {src} {dest}")
            src = os.path.join(
                args.sampledir, str(i) + "_cycle_" + str(i), str(rstNumber), "qdata.txt"
            )
            os.system(f"cp {src} {dest}")
            src = os.path.join(
                args.sampledir,
                str(i) + "_cycle_" + str(i),
                str(rstNumber2),
                "all.mdcrd",
            )
            dest = os.path.join(args.optdir, "targets", "valid_" + str(i), ".")
            os.system(f"cp {src} {dest}")
            src = os.path.join(
                args.sampledir,
                str(i) + "_cycle_" + str(i),
                str(rstNumber2),
                "qdata.txt",
            )
            os.system(f"cp {src} {dest}")
        else:
            src = os.path.join(
                args.sampledir, str(i) + "_cycle_" + str(i), str(rstNumber), "all.mdcrd"
            )
            dest = os.path.join(args.optdir, "targets", "valid_" + str(i), ".")
            os.system(f"cp {src} {dest}")
            src = os.path.join(
                args.sampledir, str(i) + "_cycle_" + str(i), str(rstNumber), "qdata.txt"
            )
            os.system(f"cp {src} {dest}")
            src = os.path.join(
                args.sampledir,
                str(i) + "_cycle_" + str(i),
                str(rstNumber2),
                "all.mdcrd",
            )
            dest = os.path.join(args.optdir, "targets", "train_" + str(i), ".")
            os.system(f"cp {src} {dest}")
            src = os.path.join(
                args.sampledir,
                str(i) + "_cycle_" + str(i),
                str(rstNumber2),
                "qdata.txt",
            )
            os.system(f"cp {src} {dest}")
    else:
        src = os.path.join(
            args.sampledir, str(i) + "_cycle_" + str(i), str(rstNumber), "all.mdcrd"
        )
        dest = os.path.join(args.optdir, "targets", "train_" + str(i), ".")
        os.system(f"cp {src} {dest}")
        src = os.path.join(
            args.sampledir, str(i) + "_cycle_" + str(i), str(rstNumber), "qdata.txt"
        )
        os.system(f"cp {src} {dest}")
        src = os.path.join(
            args.sampledir, str(i) + "_cycle_" + str(i), str(rstNumber2), "all.mdcrd"
        )
        dest = os.path.join(args.optdir, "targets", "valid_" + str(i), ".")
        os.system(f"cp {src} {dest}")
        src = os.path.join(
            args.sampledir, str(i) + "_cycle_" + str(i), str(rstNumber2), "qdata.txt"
        )
        os.system(f"cp {src} {dest}")

    # Copy previous validation and optimization FB input files to current ones
    src = os.path.join(args.optdir, "valid_" + str(i - 1) + ".in")
    dest = os.path.join(args.optdir, "valid_" + str(i) + ".in")
    os.system(f"cp {src} {dest}")
    src = os.path.join(args.optdir, "valid_" + str(i - 1) + "_initial.in")
    dest = os.path.join(args.optdir, "valid_" + str(i) + "_initial.in")
    os.system(f"cp {src} {dest}")
    src = os.path.join(args.optdir, "opt_" + str(i - 1) + ".in")
    dest = os.path.join(args.optdir, "opt_" + str(i) + ".in")
    os.system(f"cp {src} {dest}")

    # Add new targets section to each FB input file
    addTargetLines(
        os.path.join(args.optdir, "opt_" + str(i) + ".in"),
        targetLines,
        initialTarget,
        "train_" + str(i),
    )
    addTargetLines(
        os.path.join(args.optdir, "valid_" + str(i) + ".in"),
        validTargetLines,
        initialTarget,
        "valid_" + str(i),
    )
    addTargetLines(
        os.path.join(args.optdir, "valid_" + str(i) + "_initial.in"),
        validInitialTargetLines,
        initialTarget,
        "valid_" + str(i),
    )

    # Run ForceBalance on each input
    os.chdir(args.optdir)
    # TODO: valid_previous calculation broken on restart

    #What is this code doing?
    if len(validPrevious) <= i:
        if i > 1:
            src = os.path.join("result", "opt_" + str(i - 1), "*")
            dest = os.path.join("forcefield",".")
            os.system(f"cp {src} {dest}")
        os.system(
            "ForceBalance.py valid_"
            + str(i)
            + ".in > valid_"
            + str(i)
            + "_previous.out"
        )
        validPrevious.append(readValid("valid_" + str(i) + "_previous.out"))
    if (len(train)) <= i:
        os.system("ForceBalance.py opt_" + str(i) + ".in > opt_" + str(i) + ".out")
        status, results = readOpt("opt_" + str(i) + ".out")
        if status == -1:
            raise RuntimeError(
                "ForceBalance optimization of "
                + os.path.join(args.optdir, "opt_" + str(i) + ".in")
                + " failed"
            )
        if status == 1:
            print("WARNING: large change in one of the parameters")
            print("Ethan should implement adaptive changing of adaptive_damping")
        src = os.path.join("result", "opt_" + str(i), "*")
        dest = os.path.join("forcefield", ".")
        os.system(f"cp {src} {dest}")
        train.append(results["obj"])
    if len(valid) <= i:
        os.system("ForceBalance.py valid_" + str(i) + ".in > valid_" + str(i) + ".out")
        valid.append(readValid("valid_" + str(i) + ".out"))
    if len(validInitial) <= i:
        os.system(
            "ForceBalance.py valid_"
            + str(i)
            + "_initial.in > valid_"
            + str(i)
            + "_initial.out"
        )
        validInitial.append(readValid("valid_" + str(i) + "_initial.out"))
    os.chdir(home)

    if i == 1:
        print("%7d%15.8f%15.8f%20.8f" % (i,valid[-1],valid[-1]/validInitial[-1],valid[-1]-validPrevious[-1]))
    else:
        print("%7d%15.8f%15.8f%20.8f%23.8f" % (i,valid[-1],valid[-1]/validInitial[-1],valid[-1]-validPrevious[-1],valid[-1]-valid[-2]))

    # Graph results so far
    x = range(1, i + 1)
    x0 = range(i + 1)
    plt.close()
    fig, ax = plt.subplots(figsize=(9, 6))
    ax.plot(x, valid, label="Validation, current parameters", marker="o")
    ax.plot(x, validPrevious, label="Validation, previous parameters", marker="o")
    ax.plot(x, validInitial, label="Validation, initial parameters", marker="o")
    ax.plot(x0, train, label="Training", marker="o")
    ax.set_xlabel("Optimization cycle", size=17)
    ax.set_ylabel("Objective function", size=17)
    ax.set_xticks(x)
    fig.set_dpi(200)
    plt.legend(fontsize=14)
    plt.savefig("ObjectiveFunction.png", bbox_inches="tight")
    plt.close()

    types = ["BONDSK", "BONDSB", "ANGLESK", "ANGLESB", "DIHS", "VDWS", "VDWT", "COUL"]
    aliases = [
        "Bond strength",
        "Bond length",
        "Angle strength",
        "Equilibrium angle",
        "Dihedral strength",
        "LJ sigma",
        "LJ epsilon",
        "Atomic charge",
    ]
    colors = [
        "blue",
        "green",
        "firebrick",
        "goldenrod",
        "orange",
        "purple",
        "lightskyblue",
        "olive",
    ]
    name = "opt_" + str(i) + ".out"
    for j in range(len(results["labels"])):
        if labels[j] == results["labels"][j]:
            params[i + 1, j] = results["params"][j]
        else:
            for k in range(j + 1, len(labels)):
                if labels[k] == results["labels"][j]:
                    params[i + 1, k] = results["params"][j]
                    break
    adds = []
    scs = []
    legendLabels = []
    sortedParams = []
    for k in range(len(types)):
        adds.append(False)
        sortedParams.append([])
    for k in range(len(labels)):
        label = labels[k]
        for j in range(len(types)):
            if types[j] in labels[k]:
                sortedParams[j].append(params[:, k])
                # sc = ax.scatter(range(1,i + 1),relError[:,i],label=types[j],c=colors[j])
                if not adds[j]:
                    # scs.append(sc)
                    # legendLabels.append(types[j])
                    adds[j] = True
                break

    # ax.set_ylim([-100,100])
    # plt.legend(scs,legendLabels,fontsize=14)
    # plt.savefig('params.png',bbox_inches='tight')

    fig, ax = plt.subplots(figsize=(9, 6))
    ax.set_xticks(range(0, i + 1))
    for j in range(len(types)):
        if len(sortedParams[j]) == 0:
            continue
        sortedParams[j] = np.asarray(sortedParams[j], dtype=np.float32)
        diff = (sortedParams[j] - np.roll(sortedParams[j], 1, axis=1))[:, 1:]
        weights = np.maximum(
            np.abs(sortedParams[j]), np.roll(np.abs(sortedParams[j]), 1, axis=1)
        )[:, 1:]
        # normalizedDiff = np.zeros(i+1)
        mrc = np.zeros(i + 1)
        for k in range(i + 1):
            # normalizedDiff[j] = np.sqrt(np.dot(diff[:,j],diff[:,j]) / np.dot(sortedParams[i][:,j],sortedParams[i][:,j]))
            mrc[k] = np.mean(np.abs(diff[:, k]) / weights[:, k]) * 100
        # plt.plot(range(1,i+1),normalizedDiff,label=aliases[i],marker='o')
        # TODO: code crashes here after restarting at cycle 1
        plt.plot(range(0, i + 1), mrc, label=aliases[j], marker="o")
    plt.legend(fontsize=14)
    ax.tick_params(labelsize=14)
    ax.set_xlabel("Optimization Cycle", size=17)
    # ax.set_ylabel('Normalized RMS parameter change',size=17)
    ax.set_ylabel("Mean relative parameter change / %", size=17)
    plt.savefig("ParameterChange.png", bbox_inches="tight")
    plt.close()
