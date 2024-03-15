#!/usr/bin/env python

import argparse
import os
from shutil import rmtree

import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np

mpl.use("Agg")

# Some helper functions
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


def readValid(filename):
    with open(filename, "r") as f:
        for line in f.readlines():
            if "Objective Function Single Point" in line:
                return float(line.split()[6])
        raise RuntimeError(
            "ForceBalance single-point evaluation of " + filename + " did not converge"
        )


def makeValidSplits(validFile):
    name = validFile.split(".")[0]
    with open(validFile, "r") as f:
        with open(f"{name}_split1.in", "w") as v1:
            with open(f"{name}_split2.in", "w") as v2:
                for line in f.readlines():
                    splitLine = line.split()
                    if len(splitLine) > 0:
                        if splitLine[0] == "name":
                            v1.write(
                                line.replace(splitLine[1], splitLine[1] + "_split1")
                            )
                            v2.write(
                                line.replace(splitLine[1], splitLine[1] + "_split2")
                            )
                        else:
                            v1.write(line)
                            v2.write(line)
                    else:
                        v1.write("\n")
                        v2.write("\n")


def makeValidWeaves(validFile):
    name = validFile.split(".")[0]
    with open(validFile, "r") as f:
        with open(f"{name}_weave1.in", "w") as v1:
            with open(f"{name}_weave2.in", "w") as v2:
                for line in f.readlines():
                    splitLine = line.split()
                    if len(splitLine) > 0:
                        if splitLine[0] == "name":
                            v1.write(
                                line.replace(splitLine[1], splitLine[1] + "_weave1")
                            )
                            v2.write(
                                line.replace(splitLine[1], splitLine[1] + "_weave2")
                            )
                        else:
                            v1.write(line)
                            v2.write(line)
                    else:
                        v1.write("\n")
                        v2.write("\n")


def splitValids(maxcycles):
    with open("valid_1/qdata.txt", "r") as f:
        for line in f.readlines():
            if "JOB" in line:
                njobs = int(line.split()[1])

        halfjobs = int(njobs / 2)
    with open("valid_1/all.mdcrd", "r") as f:
        nlines = len(f.readlines())

    for i in range(1, maxcycles + 1):
        folder = f"valid_{str(i)}"
        folder1 = f"valid_{str(i)}_split1"
        folder2 = f"valid_{str(i)}_split2"
        if os.path.isdir(folder1):
            rmtree(folder1)
        if os.path.isdir(folder2):
            rmtree(folder2)
        os.mkdir(folder1)
        os.mkdir(folder2)
        os.system(f"cp {os.path.join(folder,'setup.leap')} {os.path.join(folder1,'.')}")
        os.system(f"cp {os.path.join(folder,'conf.pdb')} {os.path.join(folder1,'.')}")
        os.system(f"cp {os.path.join(folder,'setup.leap')} {os.path.join(folder2,'.')}")
        os.system(f"cp {os.path.join(folder,'conf.pdb')} {os.path.join(folder2,'.')}")
        switched = False
        with open(os.path.join(folder, "qdata.txt"), "r") as f:
            with open(os.path.join(folder1, "qdata.txt"), "w") as f1:
                with open(os.path.join(folder2, "qdata.txt"), "w") as f2:
                    for line in f.readlines():
                        if f"JOB {str(halfjobs+1)}" in line:
                            switched = True
                        if switched:
                            if "JOB" in line:
                                number = int(line.split()[1]) - halfjobs
                                f2.write(f"JOB {str(number)}\n")
                            else:
                                f2.write(line)
                        else:
                            f1.write(line)
        oldMdcrd = os.path.join(folder, "all.mdcrd")
        newMdcrd = os.path.join(folder1, "all.mdcrd")
        os.system(f"head -n {str(int(nlines/2)+1)} {oldMdcrd} > {newMdcrd}")
        newMdcrd = os.path.join(folder2, "all.mdcrd")
        os.system(f"echo 'comment line' > {newMdcrd}")
        os.system(f"tail -n {str(int(nlines/2))} {oldMdcrd} >> {newMdcrd}")


def weaveValids(maxcycles):
    for i in range(1, maxcycles + 1):
        folder = f"valid_{str(i)}"
        folder1 = f"valid_{str(i)}_weave1"
        folder2 = f"valid_{str(i)}_weave2"
        if os.path.isdir(folder1):
            rmtree(folder1)
        if os.path.isdir(folder2):
            rmtree(folder2)
        os.mkdir(folder1)
        os.mkdir(folder2)
        os.system(f"cp {os.path.join(folder,'setup.leap')} {os.path.join(folder1,'.')}")
        os.system(f"cp {os.path.join(folder,'conf.pdb')} {os.path.join(folder1,'.')}")
        os.system(f"cp {os.path.join(folder,'setup.leap')} {os.path.join(folder2,'.')}")
        os.system(f"cp {os.path.join(folder,'conf.pdb')} {os.path.join(folder2,'.')}")
        write1 = False
        with open(os.path.join(folder, "qdata.txt"), "r") as f:
            with open(os.path.join(folder1, "qdata.txt"), "w") as f1:
                with open(os.path.join(folder2, "qdata.txt"), "w") as f2:
                    for line in f.readlines():
                        if "JOB" in line:
                            splitLine = line.split()
                            line = line.replace(
                                splitLine[1], str(int((int(splitLine[1]) + 1) / 2))
                            )
                            if write1:
                                write1 = False
                            else:
                                write1 = True
                        if write1:
                            f1.write(line)
                        else:
                            f2.write(line)
        oldMdcrd = os.path.join(folder, "all.mdcrd")
        mdcrd1 = os.path.join(folder1, "all.mdcrd")
        mdcrd2 = os.path.join(folder2, "all.mdcrd")
        firstLine = True
        write1 = True
        lineCounter = 1
        with open(oldMdcrd, "r") as f:
            with open(mdcrd1, "w") as f1:
                with open(mdcrd2, "w") as f2:
                    f1.write("Comment line\n")
                    f2.write("Comment line\n")
                    for line in f.readlines():
                        if firstLine:
                            firstLine = False
                            continue
                        if write1:
                            f1.write(line)
                        else:
                            f2.write(line)
                        lineCounter += 1
                        if lineCounter == 27:
                            lineCounter = 1
                            if write1:
                                write1 = False
                            else:
                                write1 = True


parser = argparse.ArgumentParser()
parser.add_argument(
    "--optdir",
    help="Directory where ForceBalance optimization is performed, default is 1_opt",
    type=str,
    default="1_opt",
)
parser.add_argument(
    "--vsplit", help="Compute split validation sets", action="store_true", default=False
)
parser.add_argument(
    "--vweave", help="Compute woven validation sets", action="store_true", default=False
)
args = parser.parse_args()

# Set up the calculation
train = []
valid = []
validInitial = []
validPrevious = []
validFinalTarget = []
valid1 = []
valid2 = []

# Determine cycle for restart, set restart variables
totalCycles = -1
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
optCounter = 0
for f in os.listdir(args.optdir):
    if "opt_" in f and ".out" in f:
        optCounter += 1

for i in range(optCounter + 1):
    optOutput = os.path.join(args.optdir, "opt_" + str(i) + ".out")
    if os.path.isfile(optOutput):
        status, results = readOpt(optOutput)
        if status == 0:
            if i == 0:
                params = np.zeros((optCounter + 2, len(results["labels"])))
                labels = results["labels"]
                params[0, :] = results["initialParams"]
            else:
                try:
                    v = readValid(os.path.join(args.optdir, "valid_" + str(i) + ".out"))
                    vPrev = readValid(
                        os.path.join(args.optdir, "valid_" + str(i) + "_previous.out")
                    )
                    vInit = readValid(
                        os.path.join(args.optdir, "valid_" + str(i) + "_initial.out")
                    )
                except:
                    print("valids didn't complete")
                    break
                valid.append(v)
                validPrevious.append(vPrev)
                validInitial.append(vInit)

            train.append(results["obj"])
            for j in range(len(results["labels"])):
                if labels[j] == results["labels"][j]:
                    params[i + 1, j] = results["params"][j]
                else:
                    for k in range(j + 1, len(labels)):
                        if labels[k] == results["labels"][j]:
                            params[i + 1, k] = results["params"][j]
                            break
        else:
            print("opt didn't complete " + str(status))
            break
    else:
        break
    totalCycles = i
i = totalCycles
if totalCycles <= 0:
    raise RuntimeError("No opts completed!")

print(
    "%7s%15s%15s%20s%20s%23s"
    % ("Epoch", "Validation", "Valid ratio", "Current-Previous", "Diff percentage", "Current-last Current")
)
print(
    "%7d%15.8f%15.8f%20.8f"
    % (1, valid[0], valid[0] / validInitial[0], valid[0] - validPrevious[0])
)
for i in range(1, totalCycles):
    print(
        "%7d%15.8f%15.8f%20.8f%20.8f%23.8f"
        % (
            i + 1,
            valid[i],
            valid[i] / validInitial[i],
            valid[i] - validPrevious[i],
            (valid[i] - validPrevious[i]) / valid[i] * 100,
            valid[i] - valid[i - 1],
        )
    )


if args.vsplit:
    os.chdir(os.path.join(args.optdir, "targets"))
    splitValids(totalCycles)
    os.chdir("../..")
    for i in range(1, totalCycles + 1):
        if os.path.isfile(os.path.join(args.optdir, "valid_" + str(i) + "_split1.out")):
            try:
                v1 = readValid(
                    os.path.join(args.optdir, "valid_" + str(i) + "_split1.out")
                )
                v2 = readValid(
                    os.path.join(args.optdir, "valid_" + str(i) + "_split2.out")
                )
                valid1.append(v1)
                valid2.append(v2)
            except:
                src = os.path.join(args.optdir, "result", "opt_" + str(i), "*")
                dest = os.path.join(args.optdir, "forcefield")
                os.system(f"cp {src} {dest}")
                os.chdir(args.optdir)
                makeValidSplits(f"valid_{str(i)}.in")
                os.system(
                    f"ForceBalance.py valid_{str(i)}_split1.in > valid_{str(i)}_split1.out"
                )
                os.system(
                    f"ForceBalance.py valid_{str(i)}_split2.in > valid_{str(i)}_split2.out"
                )
                os.chdir("..")
                v1 = readValid(
                    os.path.join(args.optdir, "valid_" + str(i) + "_split1.out")
                )
                v2 = readValid(
                    os.path.join(args.optdir, "valid_" + str(i) + "_split2.out")
                )
                valid1.append(v1)
                valid2.append(v2)
        else:
            src = os.path.join(args.optdir, "result", "opt_" + str(i), "*")
            dest = os.path.join(args.optdir, "forcefield")
            os.system(f"cp {src} {dest}")
            os.chdir(args.optdir)
            makeValidSplits(f"valid_{str(i)}.in")
            os.system(
                f"ForceBalance.py valid_{str(i)}_split1.in > valid_{str(i)}_split1.out"
            )
            os.system(
                f"ForceBalance.py valid_{str(i)}_split2.in > valid_{str(i)}_split2.out"
            )
            os.chdir("..")
            v1 = readValid(os.path.join(args.optdir, "valid_" + str(i) + "_split1.out"))
            v2 = readValid(os.path.join(args.optdir, "valid_" + str(i) + "_split2.out"))
            valid1.append(v1)
            valid2.append(v2)

validInitial1 = []
validInitial2 = []
# split initial
if args.vsplit:
    src = (
        os.path.join(args.optdir, "*mol2") + " " + os.path.join(args.optdir, "*frcmod")
    )
    os.system(f"cp {src} {os.path.join(args.optdir,'forcefield','.')}")
    for i in range(1, totalCycles + 1):
        if os.path.isfile(
            os.path.join(args.optdir, "valid_" + str(i) + "_initial_split1.out")
        ):
            try:
                v1 = readValid(
                    os.path.join(args.optdir, "valid_" + str(i) + "_initial_split1.out")
                )
                v2 = readValid(
                    os.path.join(args.optdir, "valid_" + str(i) + "_initial_split2.out")
                )
            except:
                os.chdir(args.optdir)
                os.system(
                    f"ForceBalance.py valid_{str(i)}_split1.in > valid_{str(i)}_initial_split1.out"
                )
                os.system(
                    f"ForceBalance.py valid_{str(i)}_split2.in > valid_{str(i)}_initial_split2.out"
                )
                os.chdir("..")
                v1 = readValid(
                    os.path.join(args.optdir, "valid_" + str(i) + "_initial_split1.out")
                )
                v2 = readValid(
                    os.path.join(args.optdir, "valid_" + str(i) + "_initial_split2.out")
                )
        else:
            os.chdir(args.optdir)
            os.system(
                f"ForceBalance.py valid_{str(i)}_split1.in > valid_{str(i)}_initial_split1.out"
            )
            os.system(
                f"ForceBalance.py valid_{str(i)}_split2.in > valid_{str(i)}_initial_split2.out"
            )
            os.chdir("..")
            v1 = readValid(
                os.path.join(args.optdir, "valid_" + str(i) + "_initial_split1.out")
            )
            v2 = readValid(
                os.path.join(args.optdir, "valid_" + str(i) + "_initial_split2.out")
            )
        validInitial1.append(v1)
        validInitial2.append(v2)

validPrevious1 = []
validPrevious2 = []
if args.vsplit:
    for i in range(1, totalCycles + 1):
        if os.path.isfile(
            os.path.join(args.optdir, f"valid_{str(i)}_previous_split1.out")
        ):
            try:
                v1 = readValid(
                    os.path.join(args.optdir, f"valid_{str(i)}_previous_split1.out")
                )
                v2 = readValid(
                    os.path.join(args.optdir, f"valid_{str(i)}_previous_split2.out")
                )
            except:
                os.chdir(args.optdir)
                src = os.path.join("result", f"opt_" + str(i - 1), "*")
                os.system(f"cp {src} forcefield/.")
                os.system(
                    f"ForceBalance.py valid_{str(i)}_split1.in > valid_{str(i)}_previous_split1.out"
                )
                os.system(
                    f"ForceBalance.py valid_{str(i)}_split2.in > valid_{str(i)}_previous_split2.out"
                )
                v1 = readValid(f"valid_{str(i)}_previous_split1.out")
                v2 = readValid(f"valid_{str(i)}_previous_split2.out")
                os.chdir("..")
        else:
            os.chdir(args.optdir)
            src = os.path.join("result", f"opt_" + str(i - 1), "*")
            os.system(f"cp {src} forcefield/.")
            os.system(
                f"ForceBalance.py valid_{str(i)}_split1.in > valid_{str(i)}_previous_split1.out"
            )
            os.system(
                f"ForceBalance.py valid_{str(i)}_split2.in > valid_{str(i)}_previous_split2.out"
            )
            v1 = readValid(f"valid_{str(i)}_previous_split1.out")
            v2 = readValid(f"valid_{str(i)}_previous_split2.out")
            os.chdir("..")
        validPrevious1.append(v1)
        validPrevious2.append(v2)

valid1 = np.asarray(valid1)
valid2 = np.asarray(valid2)
validInitial1 = np.asarray(validInitial1)
validInitial2 = np.asarray(validInitial2)
validPrevious1 = np.asarray(validPrevious1)
validPrevious2 = np.asarray(validPrevious2)

if args.vsplit:
    print("Valid split 1")
    print(
        "%7s%15s%15s%20s%23s"
        % (
            "Epoch",
            "Validation",
            "Valid ratio",
            "Current-Previous",
            "Current-last Current",
        )
    )
    print(
        "%7d%15.8f%15.8f%20.8f"
        % (1, valid1[0], valid1[0] / validInitial1[0], valid1[0] - validPrevious1[0])
    )
    for i in range(1, totalCycles):
        print(
            "%7d%15.8f%15.8f%20.8f%23.8f"
            % (
                i + 1,
                valid1[i],
                valid1[i] / validInitial1[i],
                valid1[i] - validPrevious1[i],
                valid1[i] - valid1[i - 1],
            )
        )

    print("Valid split 2")
    print(
        "%7s%15s%15s%20s%23s"
        % (
            "Epoch",
            "Validation",
            "Valid ratio",
            "Current-Previous",
            "Current-last Current",
        )
    )
    print(
        "%7d%15.8f%15.8f%20.8f"
        % (1, valid2[0], valid2[0] / validInitial2[0], valid2[0] - validPrevious2[0])
    )
    for i in range(1, totalCycles):
        print(
            "%7d%15.8f%15.8f%20.8f%23.8f"
            % (
                i + 1,
                valid2[i],
                valid2[i] / validInitial2[i],
                valid2[i] - validPrevious2[i],
                valid2[i] - valid2[i - 1],
            )
        )

validw1 = []
validw2 = []
if args.vweave:
    os.chdir(os.path.join(args.optdir, "targets"))
    weaveValids(totalCycles)
    os.chdir("../..")
    for i in range(1, totalCycles + 1):
        if os.path.isfile(os.path.join(args.optdir, "valid_" + str(i) + "_weave1.out")):
            try:
                v1 = readValid(
                    os.path.join(args.optdir, "valid_" + str(i) + "_weave1.out")
                )
                v2 = readValid(
                    os.path.join(args.optdir, "valid_" + str(i) + "_weave2.out")
                )
            except:
                src = os.path.join(args.optdir, "result", "opt_" + str(i), "*")
                dest = os.path.join(args.optdir, "forcefield")
                os.system(f"cp {src} {dest}")
                os.chdir(args.optdir)
                makeValidWeaves(f"valid_{str(i)}.in")
                os.system(
                    f"ForceBalance.py valid_{str(i)}_weave1.in > valid_{str(i)}_weave1.out"
                )
                os.system(
                    f"ForceBalance.py valid_{str(i)}_weave2.in > valid_{str(i)}_weave2.out"
                )
                os.chdir("..")
                v1 = readValid(
                    os.path.join(args.optdir, "valid_" + str(i) + "_weave1.out")
                )
                v2 = readValid(
                    os.path.join(args.optdir, "valid_" + str(i) + "_weave2.out")
                )
        else:
            src = os.path.join(args.optdir, "result", "opt_" + str(i), "*")
            dest = os.path.join(args.optdir, "forcefield")
            os.system(f"cp {src} {dest}")
            os.chdir(args.optdir)
            makeValidWeaves(f"valid_{str(i)}.in")
            os.system(
                f"ForceBalance.py valid_{str(i)}_weave1.in > valid_{str(i)}_weave1.out"
            )
            os.system(
                f"ForceBalance.py valid_{str(i)}_weave2.in > valid_{str(i)}_weave2.out"
            )
            os.chdir("..")
            v1 = readValid(os.path.join(args.optdir, "valid_" + str(i) + "_weave1.out"))
            v2 = readValid(os.path.join(args.optdir, "valid_" + str(i) + "_weave2.out"))
        validw1.append(v1)
        validw2.append(v2)

validInitialw1 = []
validInitialw2 = []
# split initial
if args.vweave:
    src = (
        os.path.join(args.optdir, "*mol2") + " " + os.path.join(args.optdir, "*frcmod")
    )
    os.system(f"cp {src} {os.path.join(args.optdir,'forcefield','.')}")
    for i in range(1, totalCycles + 1):
        if os.path.isfile(
            os.path.join(args.optdir, "valid_" + str(i) + "_initial_weave1.out")
        ):
            try:
                v1 = readValid(
                    os.path.join(args.optdir, "valid_" + str(i) + "_initial_weave1.out")
                )
                v2 = readValid(
                    os.path.join(args.optdir, "valid_" + str(i) + "_initial_weave2.out")
                )
            except:
                os.chdir(args.optdir)
                os.system(
                    f"ForceBalance.py valid_{str(i)}_weave1.in > valid_{str(i)}_initial_weave1.out"
                )
                os.system(
                    f"ForceBalance.py valid_{str(i)}_weave2.in > valid_{str(i)}_initial_weave2.out"
                )
                os.chdir("..")
                v1 = readValid(
                    os.path.join(args.optdir, "valid_" + str(i) + "_initial_weave1.out")
                )
                v2 = readValid(
                    os.path.join(args.optdir, "valid_" + str(i) + "_initial_weave2.out")
                )
        else:
            os.chdir(args.optdir)
            os.system(
                f"ForceBalance.py valid_{str(i)}_weave1.in > valid_{str(i)}_initial_weave1.out"
            )
            os.system(
                f"ForceBalance.py valid_{str(i)}_weave2.in > valid_{str(i)}_initial_weave2.out"
            )
            os.chdir("..")
            v1 = readValid(
                os.path.join(args.optdir, "valid_" + str(i) + "_initial_weave1.out")
            )
            v2 = readValid(
                os.path.join(args.optdir, "valid_" + str(i) + "_initial_weave2.out")
            )
        validInitialw1.append(v1)
        validInitialw2.append(v2)

validPreviousw1 = []
validPreviousw2 = []
if args.vsplit:
    for i in range(1, totalCycles + 1):
        if os.path.isfile(
            os.path.join(args.optdir, f"valid_{str(i)}_previous_weave1.out")
        ):
            try:
                v1 = readValid(
                    os.path.join(args.optdir, f"valid_{str(i)}_previous_weave1.out")
                )
                v2 = readValid(
                    os.path.join(args.optdir, f"valid_{str(i)}_previous_weave2.out")
                )
            except:
                os.chdir(args.optdir)
                src = os.path.join("result", f"opt_" + str(i - 1), "*")
                os.system(f"cp {src} forcefield/.")
                os.system(
                    f"ForceBalance.py valid_{str(i)}_weave1.in > valid_{str(i)}_previous_weave1.out"
                )
                os.system(
                    f"ForceBalance.py valid_{str(i)}_weave2.in > valid_{str(i)}_previous_weave2.out"
                )
                v1 = readValid(f"valid_{str(i)}_previous_weave1.out")
                v2 = readValid(f"valid_{str(i)}_previous_weave2.out")
                os.chdir("..")
        else:
            os.chdir(args.optdir)
            src = os.path.join("result", f"opt_" + str(i - 1), "*")
            os.system(f"cp {src} forcefield/.")
            os.system(
                f"ForceBalance.py valid_{str(i)}_weave1.in > valid_{str(i)}_previous_weave1.out"
            )
            os.system(
                f"ForceBalance.py valid_{str(i)}_weave2.in > valid_{str(i)}_previous_weave2.out"
            )
            v1 = readValid(f"valid_{str(i)}_previous_weave1.out")
            v2 = readValid(f"valid_{str(i)}_previous_weave2.out")
            os.chdir("..")
        validPreviousw1.append(v1)
        validPreviousw2.append(v2)

validw1 = np.asarray(validw1)
validw2 = np.asarray(validw2)
validInitialw1 = np.asarray(validInitialw1)
validInitialw2 = np.asarray(validInitialw2)
validPreviousw1 = np.asarray(validPreviousw1)
validPreviousw2 = np.asarray(validPreviousw2)

if args.vsplit:
    print("Valid weave 1")
    print(
        "%7s%15s%15s%20s%23s"
        % (
            "Epoch",
            "Validation",
            "Valid ratio",
            "Current-Previous",
            "Current-last Current",
        )
    )
    print(
        "%7d%15.8f%15.8f%20.8f"
        % (
            1,
            validw1[0],
            validw1[0] / validInitialw1[0],
            validw1[0] - validPreviousw1[0],
        )
    )
    for i in range(1, totalCycles):
        print(
            "%7d%15.8f%15.8f%20.8f%23.8f"
            % (
                i + 1,
                validw1[i],
                validw1[i] / validInitialw1[i],
                validw1[i] - validPreviousw1[i],
                validw1[i] - validw1[i - 1],
            )
        )

    print("Valid weave 2")
    print(
        "%7s%15s%15s%20s%23s"
        % (
            "Epoch",
            "Validation",
            "Valid ratio",
            "Current-Previous",
            "Current-last Current",
        )
    )
    print(
        "%7d%15.8f%15.8f%20.8f"
        % (
            1,
            validw2[0],
            validw2[0] / validInitialw2[0],
            validw2[0] - validPreviousw2[0],
        )
    )
    for i in range(1, totalCycles):
        print(
            "%7d%15.8f%15.8f%20.8f%23.8f"
            % (
                i + 1,
                validw2[i],
                validw2[i] / validInitialw2[i],
                validw2[i] - validPreviousw2[i],
                validw2[i] - validw2[i - 1],
            )
        )

# Graph results so far
x = range(1, len(valid) + 1)
x0 = np.arange(len(valid) + 1)
if i < 25:
    xticks = x0
else:
    xticks = np.arange(0, i + 1, 2)
fig, ax = plt.subplots(figsize=(9, 6))
ax.plot(x, valid, label="Validation, current parameters", marker="o")
# ax.plot(x, validPrevious, label="Validation, previous parameters", marker="o")
ax.plot(x, validInitial, label="Validation, initial parameters", marker="o")
# ax.plot(x0, train, label="Training", marker="o")
if args.vsplit:
    ax.plot(x, valid1, label="Validation, split1", marker="o")
    ax.plot(x, valid2, label="Validation, split2", marker="o")
if args.vweave:
    ax.plot(x, validw1, label="Validation, weave1", marker="o")
    ax.plot(x, validw2, label="Validation, weave2", marker="o")
ax.set_xlabel("Optimization cycle", size=17)
ax.set_ylabel("Objective function", size=17)
ax.set_xticks(xticks)
fig.set_dpi(200)
plt.legend(fontsize=14)
plt.savefig("Validations.png", bbox_inches="tight")
plt.close()

fig, ax = plt.subplots(figsize=(9, 6))
ax.plot(
    x, np.asarray(valid) / np.asarray(validInitial), label="Valid ratio", marker="o"
)
if args.vsplit:
    ax.plot(x, valid1 / validInitial1, label="Valid ratio split 1", marker="o")
    ax.plot(x, valid2 / validInitial2, label="Valid ratio split 2", marker="o")
if args.vweave:
    ax.plot(x, validw1 / validInitialw1, label="Valid ratio weave 1", marker="o")
    ax.plot(x, validw2 / validInitialw2, label="Valid ratio weave 2", marker="o")
ax.set_xlabel("Optimization cycle", size=17)
ax.set_ylabel("Objective function", size=17)
ax.set_xticks(xticks)
fig.set_dpi(200)
plt.legend(fontsize=14)
plt.savefig("ValidationRatios.png", bbox_inches="tight")
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
fig, ax = plt.subplots(figsize=(9, 6))
for k in range(len(labels)):
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
plt.close()
fig, ax = plt.subplots(figsize=(9, 6))
ax.set_xticks(xticks)
for j in range(len(types)):
    if len(sortedParams[j]) == 0:
        continue
    sortedParams[j] = np.asarray(sortedParams[j], dtype=np.float32)
    diff = (sortedParams[j] - np.roll(sortedParams[j], 1, axis=1))[:, 1:]
    weights = np.maximum(
        np.abs(sortedParams[j]), np.roll(np.abs(sortedParams[j]), 1, axis=1)
    )[:, 1:]
    normalizedDiff = np.zeros(i)
    mrc = np.zeros(i + 1)
    for k in range(i + 1):
        # normalizedDiff[j] = np.sqrt(np.dot(diff[:,j],diff[:,j]) / np.dot(sortedParams[i][:,j],sortedParams[i][:,j]))
        mrc[k] = np.mean(np.abs(diff[:, k]) / weights[:, k]) * 100
    # plt.plot(range(1,i+1),normalizedDiff,label=aliases[i],marker='o')
    plt.plot(range(i + 1), mrc, label=aliases[j], marker="o")
plt.legend(fontsize=14)
ax.tick_params(labelsize=14)
ax.set_xlabel("Optimization Cycle", size=17)
# ax.set_ylabel('Normalized RMS parameter change',size=17)
ax.set_ylabel("Mean relative parameter change / %", size=17)
plt.savefig("ParameterChange.png", bbox_inches="tight")
plt.close()
