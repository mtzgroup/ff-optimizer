import os
from shutil import copyfile

import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np

from . import resp_prior

mpl.use("Agg")


class OptEngine:

    # We assume __init__ and all the functions it calls run from the top directory in the optimization
    def __init__(self, options):
        self.home = os.getcwd()
        self.optdir = options["optdir"]
        self.resp = options["resp"]
        self.nvalids = options['nvalids']
        self.respPriors = None
        self.leap = "setup.leap"
        if options["resp"] != 0:
            self.doResp = True
        else:
            self.doResp = False
        self.maxCycles = options["maxCycles"]
        self.mol2 = None
        self.frcmod = None
        # test setting this
        self.initialTarget = None
        with open(os.path.join(self.optdir, "opt_0.in"), "r") as f:
            for line in f.readlines():
                splitLine = line.split()
                if len(splitLine) > 1:
                    if splitLine[0] == "forcefield":
                        for i in range(1, 3):
                            if ".mol2" in splitLine[i]:
                                self.mol2 = splitLine[i]
                            elif ".frcmod" in splitLine[i]:
                                self.frcmod = splitLine[i]
                    if splitLine[0] == "name":
                        self.initialTarget = splitLine[1]
        if self.mol2 == None:
            raise RuntimeError("No mol2 file specified for optimization in opt_0.in")
        if self.frcmod == None:
            raise RuntimeError("No frcmod file specified for optimization in opt_0.in")
        if not os.path.isfile(os.path.join(self.optdir, self.mol2)):
            raise RuntimeError(
                f"Mol2 {self.mol2} specified in opt_0.in is not in {self.optdir}"
            )
        if not os.path.isfile(os.path.join(self.optdir, self.frcmod)):
            raise RuntimeError(
                f"Frcmod {self.frcmod} specified in opt_0.in is not in {self.optdir}"
            )
        os.chdir(self.optdir)
        os.system(f"tleap -f {self.leap} > leap.out")
        self.prmtop = None
        for f in os.listdir():
            if f.endswith(".prmtop"):
                self.prmtop = f
        if self.prmtop == None:
            raise RuntimeError(f"Leap script did not produce a .prmtop file!")
        os.chdir(self.home)

        # Initialize RESP priors
        if options["respPriors"] != 0:
            respOptions = {}
            respOptions["sampledir"] = options["sampledir"]
            respOptions["mol2"] = os.path.join(self.optdir, self.mol2)
            respOptions["mode"] = options["respPriors"]
            respOptions["prmtop"] = os.path.join(self.optdir, self.prmtop)
            self.respPriors = resp_prior.RespPriors(respOptions)

        self.train = []
        self.valid = []
        self.validInitial = []
        self.validPrevious = []
        if options["restart"]:
            self.restartCycle = self.determineRestart()
        else:
            self.restartCycle = -1

        self.targetLines = []
        self.validTargetLines = []
        self.validInitialTargetLines = []
        inTarget = False
        with open(os.path.join(self.optdir, "opt_0.in"), "r") as f:
            for line in f.readlines():
                if len(line.split()) == 0:
                    continue
                if "$target" in line:
                    inTarget = True
                if inTarget:
                    if line.split()[0] == "resp" or line.split()[0] == "w_resp":
                        continue
                    self.targetLines.append(line)
                    self.validTargetLines.append(line)
                    if line.split()[0] == "amber_leapcmd":
                        line = line.replace(line.split()[1], "setup_valid_initial.leap")
                    self.validInitialTargetLines.append(line)
        if self.doResp:
            self.targetLines.insert(1, f"w_resp {str(self.resp)}\n")
            self.targetLines.insert(1, "resp 1\n")
        with open(os.path.join(self.optdir, "setup.leap"), "r") as leapRead:
            with open(
                os.path.join(self.optdir, "setup_valid_initial.leap"), "w"
            ) as leapWrite:
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
        if not os.path.isdir(os.path.join(self.optdir, "forcefield")):
            os.mkdir(os.path.join(self.optdir, "forcefield"))
        if self.restartCycle == -1:
            copyfile(
                os.path.join(self.optdir, self.frcmod),
                os.path.join(self.optdir, "forcefield", self.frcmod),
            )
            copyfile(
                os.path.join(self.optdir, self.frcmod),
                os.path.join(self.optdir, "forcefield", f"initial_{self.frcmod}"),
            )
            copyfile(
                os.path.join(self.optdir, self.mol2),
                os.path.join(self.optdir, "forcefield", self.mol2),
            )
            copyfile(
                os.path.join(self.optdir, self.mol2),
                os.path.join(self.optdir, "forcefield", f"initial_{self.mol2}"),
            )
        # Make validation input for initial MM parameters
        with open(os.path.join(self.optdir, "valid_0.in"), "r") as srcValid:
            with open(
                os.path.join(self.optdir, "valid_0_initial.in"), "w"
            ) as destValid:
                for line in srcValid.readlines():
                    if "$target" in line:
                        break
                    destValid.write(
                        line.replace(self.frcmod, f"initial_{self.frcmod}").replace(
                            self.mol2, f"initial_{self.mol2}"
                        )
                    )

    def readOpt(self, filename):
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

        params = np.asarray(params, dtype=np.float32)
        initialParams = np.asarray(initialParams, dtype=np.float32)
        results["params"] = params
        results["labels"] = labels
        results["initialParams"] = initialParams
        return status, results

    # Unused and untested
    def changeParameter(self, inputFile, prmName, prmValue):
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

        with open("temp.txt", "w") as f:
            for line in lines:
                f.write(line)
        os.system(f"mv temp.txt {inputFile}")

    def addTargetLines(self, inputFile, targetLines, newTarget):
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
                    f.write(line.replace(self.initialTarget, newTarget))

    # Unused and untested
    def determineAdaptiveDamping(
        self, testFile, upperThreshold=0.3, lowerThreshold=0.01, adaptiveDamping=0.5
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
            if np.argwhere(diff > upperThreshold).shape[0] > 0:
                adaptiveDamping *= 2
                changeParameter(testFile, "adaptive_damping", str(adaptiveDamping))
            elif np.argwhere(diff > lowerThreshold).shape[0] == 0:
                adaptiveDamping *= 0.75
                changeParameter(testFile, "adaptive_damping", str(adaptiveDamping))
            else:
                return adaptiveDamping

    def readValid(self, filename):
        with open(filename, "r") as f:
            for line in f.readlines():
                if "Objective Function Single Point" in line:
                    return float(line.split()[6])
            raise RuntimeError(
                f"ForceBalance single-point evaluation of {filename} did not converge"
            )

    def setupInputFiles(self, i):
        # Copy previous validation and optimization FB input files to current ones
        copyfile(f"valid_{str(i - 1)}.in", f"valid_{str(i)}.in")
        copyfile(f"valid_{str(i - 1)}_initial.in", f"valid_{str(i)}_initial.in")
        copyfile(f"opt_{str(i - 1)}.in", f"opt_{str(i)}.in")
        for j in range(1, self.nvalids):
            copyfile(f"valid_{str(i - 1)}_{str(j)}.in", f"valid_{str(i)}_{str(j)}.in")
            copyfile(f"valid_{str(i - 1)}_{str(j)}_initial.in", f"valid_{str(i)}_{str(j)}_initial.in")


        # Add new targets section to each FB input file
        self.addTargetLines(
            f"opt_{str(i)}.in",
            self.targetLines,
            f"train_{str(i)}",
        )
        self.addTargetLines(
            f"valid_{str(i)}.in", self.validTargetLines, f"valid_{str(i)}"
        )
        self.addTargetLines(
            f"valid_{str(i)}_initial.in",
            self.validInitialTargetLines,
            f"valid_{str(i)}",
        )
        for j in range(1,self.nvalids):
            self.addTargetLines(
                f"valid_{str(i)}_{str(j)}.in", self.validTargetLines, f"valid_{str(i)}_{str(j)}"
            )
            self.addTargetLines(
                f"valid_{str(i)}_{str(j)}_initial.in",
                self.validInitialTargetLines,
                f"valid_{str(i)}_{str(j)}",
            )
            

    def graphResults(self):
        # Graph results so far
        cycles = len(self.valid)
        x = np.arange(cycles) + 1
        x0 = np.arange(cycles + 1)
        tickInterval = max(int(cycles / 12), 1)
        xticks = np.arange(0, cycles + 1, tickInterval)
        fig, ax = plt.subplots(figsize=(9, 6))
        ax.plot(x, self.valid, label="Validation, current parameters", marker="o")
        ax.plot(
            x, self.validPrevious, label="Validation, previous parameters", marker="o"
        )
        ax.plot(
            x, self.validInitial, label="Validation, initial parameters", marker="o"
        )
        ax.plot(x0, self.train, label="Training", marker="o")
        ax.set_xlabel("Optimization cycle", size=17)
        ax.set_ylabel("Objective function", size=17)
        ax.set_xticks(xticks)
        fig.set_dpi(200)
        plt.legend(fontsize=14)
        plt.savefig(os.path.join("..", "ObjectiveFunction.png"), bbox_inches="tight")
        plt.close()

        types = [
            "BONDSK",
            "BONDSB",
            "ANGLESK",
            "ANGLESB",
            "DIHS",
            "VDWS",
            "VDWT",
            "COUL",
        ]
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
        adds = []
        sortedParams = []
        for k in range(len(types)):
            adds.append(False)
            sortedParams.append([])
        fig, ax = plt.subplots(figsize=(9, 6))
        for k in range(len(self.labels)):
            for j in range(len(types)):
                if types[j] in self.labels[k]:
                    sortedParams[j].append(self.params[:, k])
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
            # normalizedDiff = np.zeros(i)
            mrc = np.zeros(cycles + 1)
            for k in range(cycles + 1):
                # normalizedDiff[j] = np.sqrt(np.dot(diff[:,j],diff[:,j]) / np.dot(sortedParams[i][:,j],sortedParams[i][:,j]))
                mrc[k] = np.mean(np.abs(diff[:, k]) / weights[:, k]) * 100
            # plt.plot(range(1,i+1),normalizedDiff,label=aliases[i],marker='o')
            plt.plot(range(cycles + 1), mrc, label=aliases[j], marker="o")
        plt.legend(fontsize=14)
        ax.tick_params(labelsize=14)
        ax.set_xlabel("Optimization Cycle", size=17)
        # ax.set_ylabel('Normalized RMS parameter change',size=17)
        ax.set_ylabel("Mean relative parameter change / %", size=17)
        plt.savefig(os.path.join("..", "ParameterChange.png"), bbox_inches="tight")
        plt.close()

    def sortParams(self, results, i):
        self.params[i + 1, :] = self.params[i, :]
        for j in range(len(results["labels"])):
            if self.labels[j] == results["labels"][j]:
                self.params[i + 1, j] = results["params"][j]
            else:
                for k in range(len(self.labels)):
                    if self.labels[k] == results["labels"][j]:
                        self.params[i + 1, k] = results["params"][j]
                        break

    # TODO: make each optimization occur in a separate directory
    # We assume that optimizeForcefield and all the functions it calls run in the args.optdir directory
    def optimizeForcefield(self, i):
        if i > 0:
            self.setupInputFiles(i)
            # If we're just restarting, skip if this calculation finished
            if len(self.validPrevious) <= i:
                if i > 1:
                    copyfile(
                        os.path.join("result", f"opt_{str(i - 1)}", self.frcmod),
                        os.path.join("forcefield", self.frcmod),
                    )
                    copyfile(
                        os.path.join("result", f"opt_{str(i - 1)}", self.mol2),
                        os.path.join("forcefield", self.mol2),
                    )
                os.system(
                    f"ForceBalance.py valid_{str(i)}.in > valid_{str(i)}_previous.out"
                )
                self.validPrevious.append(
                    self.readValid(f"valid_{str(i)}_previous.out")
                )
                for j in range(1, self.nvalids):
                    os.system(
                        f"ForceBalance.py valid_{str(i)}_{str(j)}.in > valid_{str(i)}_{str(j)}_previous.out"
                    )

            if len(self.train) <= i:
                if self.respPriors is not None:
                    self.respPriors.updateRespPriors(
                        i, os.path.join("forcefield", self.mol2)
                    )
                os.system(f"ForceBalance.py opt_{str(i)}.in > opt_{str(i)}.out")
                status, results = self.readOpt(f"opt_{str(i)}.out")
                if status == -1:
                    raise RuntimeError(
                        f"ForceBalance optimization of {os.path.join(self.optdir, f'opt_{str(i)}.in')} failed"
                    )
                if status == 1:
                    print("WARNING: large change in one of the parameters")
                    print(
                        "Ethan should implement adaptive changing of adaptive_damping"
                    )
                copyfile(
                    os.path.join("result", f"opt_{str(i)}", self.frcmod),
                    os.path.join("forcefield", self.frcmod),
                )
                copyfile(
                    os.path.join("result", f"opt_{str(i)}", self.mol2),
                    os.path.join("forcefield", self.mol2),
                )
                self.train.append(results["obj"])
                self.sortParams(results, i)
            if len(self.valid) <= i:
                os.system(f"ForceBalance.py valid_{str(i)}.in > valid_{str(i)}.out")
                self.valid.append(self.readValid(f"valid_{str(i)}.out"))
                for j in range(1, self.nvalids):
                    os.system(
                        f"ForceBalance.py valid_{str(i)}_{str(j)}.in > valid_{str(i)}_{str(j)}.out"
                    )
            if len(self.validInitial) <= i:
                os.system(
                    f"ForceBalance.py valid_{str(i)}_initial.in > valid_{str(i)}_initial.out"
                )
                self.validInitial.append(self.readValid(f"valid_{str(i)}_initial.out"))
            self.graphResults()
        else:
            os.system("ForceBalance.py opt_0.in > opt_0.out")
            copyfile(
                os.path.join("result", "opt_0", self.frcmod),
                os.path.join("forcefield", self.frcmod),
            )
            copyfile(
                os.path.join("result", "opt_0", self.mol2),
                os.path.join("forcefield", self.mol2),
            )
            status, results = self.readOpt("opt_0.out")
            if status != 0:
                raise RuntimeError("ForceBalance optimization of opt_0.in failed")
            self.train.append(results["obj"])
            self.labels = results["labels"]
            self.params = np.zeros((self.maxCycles + 2, len(self.labels)))
            self.params[0, :] = np.asarray(results["initialParams"])
            self.sortParams(results, i)

    def determineRestart(self):
        # Determine cycle for restart, set restart variables
        restartCycle = -1
        for i in range(self.maxCycles + 2):
            optOutput = os.path.join(self.optdir, "opt_" + str(i) + ".out")
            if os.path.isfile(optOutput):
                status, results = self.readOpt(optOutput)
                if status == 0:
                    if i == 0:
                        self.params = np.zeros(
                            (self.maxCycles + 2, len(results["labels"]))
                        )
                        self.labels = results["labels"]
                        self.params[0, :] = results["initialParams"]
                    else:
                        try:
                            v = self.readValid(
                                os.path.join(self.optdir, f"valid_{str(i)}.out")
                            )
                            vPrev = self.readValid(
                                os.path.join(
                                    self.optdir, f"valid_{str(i)}_previous.out"
                                )
                            )
                            vInitial = self.readValid(
                                os.path.join(self.optdir, f"valid_{str(i)}_initial.out")
                            )
                        except:
                            break
                        self.valid.append(v)
                        self.validPrevious.append(vPrev)
                        self.validInitial.append(vInitial)
                        if self.respPriors is not None:
                            self.respPriors.getCharges(i)
                    self.train.append(results["obj"])
                    self.sortParams(results, i)
                else:
                    break
            else:
                break
        restartCycle = i - 1
        print("Restarting optimization at cycle " + str(restartCycle + 1))
        return restartCycle
