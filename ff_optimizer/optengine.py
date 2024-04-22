import os
from pathlib import Path
from shutil import copyfile

import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np

from . import resp_prior

mpl.use("Agg")


class OptEngine:
    def setVariables(self, inp):
        self.converged = False
        self.home = os.getcwd()
        self.inp = inp
        self.optdir = Path(inp.optdir)
        self.resp = inp.resp
        self.nvalids = inp.nvalids
        self.respPriors = None
        self.leap = "setup.leap"
        if inp.resp != 0:
            self.doResp = True
        else:
            self.doResp = False
        self.maxCycles = inp.maxcycles
        self.mol2 = None
        self.frcmod = None
        # test setting this
        self.initialTarget = None
        self.targetLines = []
        self.validTargetLines = []
        self.validInitialTargetLines = []

    def readFileNames(self):
        with open(self.optdir / "opt_0.in", "r") as f:
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

    def checkFileNames(self):
        if self.mol2 == None:
            raise RuntimeError("No mol2 file specified for optimization in opt_0.in")
        if self.frcmod == None:
            raise RuntimeError("No frcmod file specified for optimization in opt_0.in")
        if not (self.optdir / self.mol2).is_file():
            raise RuntimeError(
                f"Mol2 {self.mol2} specified in opt_0.in is not in {self.optdir}"
            )
        if not (self.optdir / self.frcmod).is_file():
            raise RuntimeError(
                f"Frcmod {self.frcmod} specified in opt_0.in is not in {self.optdir}"
            )

    def testTleap(self):
        os.chdir(self.optdir)
        os.system(f"tleap -f {self.leap} > leap.out")
        self.prmtop = None
        for f in os.listdir():
            if f.endswith(".prmtop"):
                self.prmtop = f
        if self.prmtop == None:
            raise RuntimeError(f"Leap script did not produce a .prmtop file!")
        os.chdir(self.home)

    def initializeRespPriors(self):
        if self.inp.resppriors != 0:
            mol2 = self.optdir / self.mol2
            self.respPriors = resp_prior.RespPriors(self.inp, mol2, self.prmtop)

    def readTargetLines(self):
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

    def writeValidInitialLeap(self):
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

    def copyToFF(self, f):
        copyfile(self.optdir / f, self.optdir / "forcefield" / f)
        copyfile(self.optdir / f, self.optdir / "forcefield" / f"initial_{f}")

    def copyFiles(self):
        (self.optdir / "forcefield").mkdir(exist_ok=True)
        if self.restartCycle == -1:
            for f in [self.frcmod, self.mol2]:
                self.copyToFF(f)

    def makeValidIn(self):
        with open(self.optdir / "valid_0.in", "r") as srcValid:
            with open(self.optdir / "temp.txt", "w") as destValid:
                for line in srcValid.readlines():
                    if "$target" in line:
                        break
                    destValid.write(line)
        os.rename(self.optdir / "temp.txt", self.optdir / "valid_0.in")

    def makeInitialValidIn(self):
        with open(self.optdir / "valid_0.in", "r") as srcValid:
            with open(self.optdir / "valid_0_initial.in", "w") as destValid:
                for line in srcValid.readlines():
                    if "$target" in line:
                        break
                    destValid.write(
                        line.replace(self.frcmod, f"initial_{self.frcmod}").replace(
                            self.mol2, f"initial_{self.mol2}"
                        )
                    )

    def copyValids(self):
        for j in range(1, self.nvalids):
            copyfile(self.optdir / "valid_0.in", self.optdir / f"valid_0_{j}.in")
            copyfile(
                self.optdir / "valid_0_initial.in",
                self.optdir / f"valid_0_{j}_initial.in",
            )

    # We assume __init__ and all the functions it calls run from the top directory in the optimization
    # All other functions are called from within self.optdir
    def __init__(self, inp):
        # set some internal variables based on the input
        self.setVariables(inp)
        # read .frcmod and .mol2 file names from forcebalance input
        self.readFileNames()
        # check that those files exist
        self.checkFileNames()
        # check that tleap creates a prmtop
        self.testTleap()
        # initialize RESP priors
        self.initializeRespPriors()
        # determine which cycle to restart at, if applicable
        self.determineRestart()
        # read target-specific lines from forcebalance input
        self.readTargetLines()
        # write tleap file for evaluating initial parameters on validation set
        self.writeValidInitialLeap()
        # copy forcefield files to forcefield dir for forcebalance
        self.copyFiles()
        # make validation set inputs for forcebalance
        self.makeValidIn()
        self.makeInitialValidIn()
        # duplicate inputs if using multiple validation sets
        self.copyValids()

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
        oldFiles = [
            f"opt_{i - 1}.in",
            f"valid_{i - 1}.in",
            f"valid_{i - 1}_initial.in",
            f"opt_{i - 1}.in",
        ]
        newFiles = [
            f"opt_{i}.in",
            f"valid_{i}.in",
            f"valid_{i}_initial.in",
            f"opt_{i}.in",
        ]
        for j in range(1, self.nvalids):
            oldFiles.append(f"valid_{i - 1}_{j}.in")
            newFiles.append(f"valid_{i}_{j}.in")
            oldFiles.append(f"valid_{i - 1}_{j}_initial.in")
            newFiles.append(f"valid_{i}_{j}_initial.in")

        for oldFile, newFile in zip(oldFiles, newFiles):
            copyfile(oldFile, newFile)
            # Add new targets section to each FB input file
            if "opt" in newFile:
                self.addTargetLines(newFile, self.targetLines, f"train_{i}")
            elif "initial" in newFile:
                self.addTargetLines(newFile, self.validInitialTargetLines, f"valid_{i}")
            else:
                self.addTargetLines(newFile, self.validTargetLines, f"valid_{i}")

    def computeXTicks(self):
        cycles = len(self.valid)
        x = np.arange(cycles) + 1
        tickInterval = max(int(cycles / 12), 1)
        xticks = np.arange(0, cycles + 1, tickInterval)
        return cycles, x, xticks

    def graphObjectiveFunction(self):
        cycles, x, xticks = self.computeXTicks()
        fig, ax = plt.subplots(figsize=(9, 6))
        ax.plot(x, self.valid, label="Validation, current parameters", marker="o")
        ax.plot(
            x, self.validPrevious, label="Validation, previous parameters", marker="o"
        )
        ax.plot(
            x, self.validInitial, label="Validation, initial parameters", marker="o"
        )
        x0 = np.arange(cycles + 1)
        ax.plot(x0, self.train, label="Training", marker="o")
        ax.set_xlabel("Optimization cycle", size=17)
        ax.set_ylabel("Objective function", size=17)
        ax.set_xticks(xticks)
        fig.set_dpi(200)
        plt.legend(fontsize=14)
        plt.savefig(self.optdir / "ObjectiveFunction.png", bbox_inches="tight")
        plt.close()

    def computeSortedParams(self):
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
        adds = []
        sortedParams = []
        for k in range(len(types)):
            adds.append(False)
            sortedParams.append([])
        # fig, ax = plt.subplots(figsize=(9, 6))
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
        # plt.close()
        return sortedParams

    def computeMRC(self, sortedParamsType):
        cycles = len(self.valid)
        sortedParams = np.asarray(sortedParamsType, dtype=np.float32)
        diff = (sortedParams - np.roll(sortedParams, 1, axis=1))[:, 1:]
        weights = np.maximum(
            np.abs(sortedParams), np.roll(np.abs(sortedParams), 1, axis=1)
        )[:, 1:]
        # normalizedDiff = np.zeros(i)
        mrc = np.zeros(cycles + 1)
        for k in range(cycles + 1):
            # normalizedDiff[j] = np.sqrt(np.dot(diff[:,j],diff[:,j]) / np.dot(sortedParams[i][:,j],sortedParams[i][:,j]))
            mrc[k] = np.mean(np.abs(diff[:, k]) / weights[:, k]) * 100
        return mrc

    def graphMRC(self):
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

        sortedParams = self.computeSortedParams()
        fig, ax = plt.subplots(figsize=(9, 6))
        cycles, x, xticks = self.computeXTicks()
        ax.set_xticks(xticks)
        for j in range(len(aliases)):
            if len(sortedParams[j]) == 0:
                continue
            mrc = self.computeMRC(sortedParams[j])
            plt.plot(range(cycles + 1), mrc, label=aliases[j], marker="o")
        plt.legend(fontsize=14)
        ax.tick_params(labelsize=14)
        ax.set_xlabel("Optimization Cycle", size=17)
        # ax.set_ylabel('Normalized RMS parameter change',size=17)
        ax.set_ylabel("Mean relative parameter change / %", size=17)
        plt.savefig(self.optdir / "ParameterChange.png", bbox_inches="tight")
        plt.close()

    def graphResults(self):
        # Graph results so far
        self.graphObjectiveFunction()
        self.graphMRC()

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

    def copyResults(self, i):
        copyfile(
            os.path.join("result", f"opt_{i}", self.frcmod),
            os.path.join("forcefield", self.frcmod),
        )
        copyfile(
            os.path.join("result", f"opt_{i}", self.mol2),
            os.path.join("forcefield", self.mol2),
        )

    def runValidPrevious(self, i):
        # If we're just restarting, skip if this calculation finished
        if len(self.validPrevious) < i:
            os.system(f"ForceBalance.py valid_{i}.in > valid_{i}_previous.out")
            for j in range(1, self.nvalids):
                os.system(
                    f"ForceBalance.py valid_{i}_{j}.in > valid_{i}_{j}_previous.out"
                )
            self.validPrevious.append(self.readValid(f"valid_{i}_previous.out"))

    def runTraining(self, i):
        if len(self.train) <= i:
            if self.respPriors is not None:
                self.respPriors.updateRespPriors(
                    i, os.path.join("forcefield", self.mol2)
                )
            os.system(f"ForceBalance.py opt_{i}.in > opt_{i}.out")
            status, results = self.readOpt(f"opt_{i}.out")
            if status == -1:
                raise RuntimeError(
                    f"ForceBalance optimization of {os.path.join(self.optdir, f'opt_{i}.in')} failed"
                )
            if status == 1:
                print("WARNING: large change in one of the parameters")
            self.train.append(results["obj"])
            self.sortParams(results, i)

    def runValid(self, i):
        if len(self.valid) < i:
            os.system(f"ForceBalance.py valid_{i}.in > valid_{i}.out")
            for j in range(1, self.nvalids):
                os.system(f"ForceBalance.py valid_{i}_{j}.in > valid_{i}_{j}.out")
            self.valid.append(self.readValid(f"valid_{i}.out"))

    def runValidInitial(self, i):
        if len(self.validInitial) < i:
            os.system(f"ForceBalance.py valid_{i}_initial.in > valid_{i}_initial.out")
            for j in range(1, self.nvalids):
                os.system(
                    f"ForceBalance.py valid_{i}_{j}_initial.in > valid_{i}_{j}_initial.out"
                )
            self.validInitial.append(self.readValid(f"valid_{i}_initial.out"))

    def runInitialTraining(self):
        os.system("ForceBalance.py opt_0.in > opt_0.out")
        self.copyResults(0)
        status, results = self.readOpt("opt_0.out")
        if status != 0:
            raise RuntimeError("ForceBalance optimization of opt_0.in failed")
        self.train.append(results["obj"])
        self.labels = results["labels"]
        self.params = np.zeros((self.maxCycles + 2, len(self.labels)))
        self.params[0, :] = np.asarray(results["initialParams"])
        self.sortParams(results, 0)

    # We assume that optimizeForcefield and all the functions it calls run in the args.optdir directory
    def optimizeForcefield(self, i):
        if i > 0:
            # copy ff files from previous cycle (important if restarting)
            self.copyResults(i - 1)
            # create forcebalance input files
            self.setupInputFiles(i)
            # evaluate new validation set with previous parameters
            self.runValidPrevious(i)
            # optimize forcefield parameters
            self.runTraining(i)
            # copy ff files from optimization
            self.copyResults(i)
            # evaluate new validation set with new parameters
            self.runValid(i)
            # evaluate new validation set with initial parameters
            self.runValidInitial(i)
            # make some pretty graphs
            self.graphResults()
            # check if iterative optimization is done yet
            self.checkConvergence()
        else:
            # special setup for first-time training
            self.runInitialTraining()

    def checkValids(self, i, suffix=""):
        v = self.readValid(self.optdir / f"valid_{i}{suffix}.out")
        for j in range(1, self.nvalids):
            self.readValid(self.optdir / f"valid_{i}_{j}{suffix}.out")
        return v

    def computeValidDiff(self):
        vDiff = []
        for i in range(len(self.valid)):
            vDiff.append((self.valid[i] - self.validPrevious[i]) / self.valid[i] * 100)
        return vDiff

    def checkConvergence(self):
        patience = 5
        inPatience = False
        cutoff = -1 # Cutoff is 1% change in performance
        lastCycle = -1
        validDiff = self.computeValidDiff()
        for j in range(len(self.valid)):
            if not inPatience and validDiff[j] > cutoff:
                inPatience = True
                patienceCycle = j
            if inPatience and validDiff[j] < cutoff:
                inPatience = False
            if inPatience and j - patienceCycle >= patience:
                lastCycle = j
                self.converged = True
                break
        return lastCycle

    def checkOpt(self, i):
        status, results = self.readOpt(self.optdir / f"opt_{i}.out")
        if status == 0:
            if i == 0:
                self.params = np.zeros((self.maxCycles + 2, len(results["labels"])))
                self.labels = results["labels"]
                self.params[0, :] = results["initialParams"]
            self.train.append(results["obj"])
            self.sortParams(results, i)
        else:
            raise RuntimeError("ForceBalance failed!")

    def determineRestart(self):
        self.train = []
        self.valid = []
        self.validInitial = []
        self.validPrevious = []
        if not self.inp.restart:
            self.restartCycle = -1
            return
        # Determine cycle for restart, set restart variables
        for i in range(self.maxCycles + 2):
            if i > 0:
                # check if valid previous finished
                try:
                    vPrev = self.checkValids(i, "_previous")
                except:
                    break
                self.validPrevious.append(vPrev)
            # check if parameter optimization finished
            try:
                self.checkOpt(i)
            except:
                break
            if i > 0:
                # check if validation finished
                try:
                    v = self.checkValids(i)
                except:
                    break
                self.valid.append(v)
                # check if validation with initial params finished
                try:
                    vInitial = self.checkValids(i, "_initial")
                except:
                    break
                self.validInitial.append(vInitial)
                if self.respPriors is not None:
                    self.respPriors.getCharges(i)
        self.restartCycle = i

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
