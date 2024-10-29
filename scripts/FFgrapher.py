#!/usr/bin/env python

import argparse
import os

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


parser = argparse.ArgumentParser()
parser.add_argument(
    "--optdir",
    help="Directory where ForceBalance optimization is performed, default is 1_opt",
    type=str,
    default="1_opt",
)
parser.add_argument(
    "--final",
    help="Compute performance of parameters from each epoch on the final validation set",
    action="store_true",
    default=False,
)
args = parser.parse_args()

# Set up the calculation
train = []
valid = []
validInitial = []
validPrevious = []
validFinalTarget = []

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
                params = np.zeros((optCounter + 1, len(results["labels"])))
                labels = results["labels"]
                params[0, :] = results["initialParams"]
            else:
                try:
                    v = readValid(os.path.join(args.optdir, f"valid_{str(i)}.out"))
                    vPrev = readValid(
                        os.path.join(args.optdir, f"valid_{str(i)}_previous.out")
                    )
                    vInitial = readValid(
                        os.path.join(args.optdir, f"valid_{str(i)}_initial.out")
                    )
                except:
                    break
                valid.append(v)
                validPrevious.append(vPrev)
                validInitial.append(vInitial)
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
    "%7s%15s%15s%20s%23s"
    % ("Epoch", "Validation", "Valid ratio", "Current-Previous", "Current-last Current")
)
print(
    "%7d%15.8f%15.8f%20.8f"
    % (1, valid[0], valid[0] / validInitial[0], valid[0] - validPrevious[0])
)
for i in range(1, totalCycles):
    print(
        "%7d%15.8f%15.8f%20.8f%23.8f"
        % (
            i + 1,
            valid[i],
            valid[i] / validInitial[i],
            valid[i] - validPrevious[i],
            valid[i] - valid[i - 1],
        )
    )

if args.final:
    for i in range(totalCycles + 1):
        if os.path.isfile(
            os.path.join(args.optdir, "valid_" + str(i) + "_final_target.out")
        ):
            try:
                validFinalTarget.append(
                    readValid(
                        os.path.join(
                            args.optdir, "valid_" + str(i) + "_final_target.out"
                        )
                    )
                )
            except:
                src = os.path.join(args.optdir, "result", "opt_" + str(i), "*")
                dest = os.path.join(args.optdir, "forcefield")
                os.system(f"cp {src} {dest}")
                os.chdir(args.optdir)
                os.system(
                    "ForceBalance.py valid_"
                    + str(totalCycles)
                    + ".in > valid_"
                    + str(i)
                    + "_final_target.out"
                )
                os.chdir("..")
                validFinalTarget.append(
                    readValid(
                        os.path.join(
                            args.optdir, "valid_" + str(i) + "_final_target.out"
                        )
                    )
                )
        else:
            src = os.path.join(args.optdir, "result", "opt_" + str(i), "*")
            dest = os.path.join(args.optdir, "forcefield")
            os.system(f"cp {src} {dest}")
            os.chdir(args.optdir)
            os.system(
                "ForceBalance.py valid_"
                + str(totalCycles)
                + ".in > valid_"
                + str(i)
                + "_final_target.out"
            )
            os.chdir("..")
            validFinalTarget.append(
                readValid(
                    os.path.join(args.optdir, "valid_" + str(i) + "_final_target.out")
                )
            )

# Graph results so far
x = range(1, i + 2)
x0 = np.arange(i + 2)
xticks = np.arange(0, i + 1, int(i / 12))
fig, ax = plt.subplots(figsize=(9, 6))
ax.plot(x, valid, label="Validation, current parameters", marker="o")
ax.plot(x, validPrevious, label="Validation, previous parameters", marker="o")
ax.plot(x, validInitial, label="Validation, initial parameters", marker="o")
ax.plot(x0, train, label="Training", marker="o")
ax.set_xlabel("Optimization cycle", size=17)
ax.set_ylabel("Objective function", size=17)
ax.set_xticks(xticks)
fig.set_dpi(200)
plt.legend(fontsize=14)
plt.savefig("ObjectiveFunction.png", bbox_inches="tight")
plt.close()

if args.final:
    fig, ax = plt.subplots(figsize=(9, 6))
    ax.plot(
        x0, validFinalTarget, label="Performance on final validation set", marker="o"
    )
    ax.set_xlabel("Optimization cycle", size=17)
    ax.set_ylabel("Objective function", size=17)
    ax.set_xticks(xticks)
    fig.set_dpi(200)
    plt.legend(fontsize=14)
    plt.savefig("ObjectiveFunctionFinalTarget.png", bbox_inches="tight")
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
    ax.semilogy(range(i + 1), mrc, label=aliases[j], marker="o")
    # plt.plot(range(i + 1), mrc, label=aliases[j], marker="o")
plt.legend(fontsize=14)
ax.tick_params(labelsize=14)
ax.set_xlabel("Optimization Cycle", size=17)
# ax.set_ylabel('Normalized RMS parameter change',size=17)
ax.set_ylabel("Mean relative parameter change / %", size=17)
plt.savefig("ParameterChange.png", bbox_inches="tight")
plt.close()
