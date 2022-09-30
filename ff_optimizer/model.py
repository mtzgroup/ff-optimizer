import os
from shutil import copyfile, rmtree

from . import mmengine, optengine, qmengine
from .utils import convertTCtoFB


# Template class for ff models used by ff_optimizer
# All these functions must be implemented and
# all these variables set properly for ff_optimizer to
# use the inherited class
class AbstractModel:
    def __init__(self):
        self.restartCycle = -1

    def initialCycle(self):
        pass

    def doMMSampling(self, i):
        pass

    def doQMCalculations(self, i):
        pass

    def doParameterOptimization(self, i):
        pass


class Model(AbstractModel):
    def __init__(self, args):
        self.args = args
        # Set some miscellaneous variables
        self.home = os.getcwd()
        self.optdir = args.optdir
        self.sampledir = args.sampledir
        os.rename(
            os.path.join(args.optdir, args.valid0),
            os.path.join(args.optdir, "valid_0.in"),
        )
        os.rename(
            os.path.join(args.optdir, args.opt0), os.path.join(args.optdir, "opt_0.in")
        )
        self.mdFiles = []
        self.heatCounter = 0
        for f in os.listdir(args.sampledir):
            if os.path.isfile(os.path.join(args.sampledir, f)):
                self.mdFiles.append(f)
                if f.startswith("heat"):
                    self.heatCounter += 1
        if args.resp != 0 or args.respPriors != 0:
            self.doResp = True
        else:
            self.doResp = False
        # Set up each engine
        self.optEngine = self.initializeOptEngine(args)
        self.qmEngine = self.initializeQMEngine(args)
        self.mmEngine = self.initializeMMEngine(args)
        self.restartCycle = self.optEngine.restartCycle

    def initializeOptEngine(self, args):
        optOptions = {}
        optOptions["optdir"] = args.optdir
        optOptions["sampledir"] = args.sampledir
        optOptions["respPriors"] = args.respPriors
        optOptions["resp"] = args.resp
        optOptions["maxCycles"] = args.maxcycles
        optOptions["restart"] = args.restart
        optOptions["nvalids"] = args.nvalids
        optEngine = optengine.OptEngine(optOptions)
        return optEngine

    def initializeQMEngine(self, args):
        if args.qmengine == "debug":
            qmEngine = qmengine.DebugEngine(
                os.path.join(args.sampledir, args.tctemplate),
                os.path.join(args.sampledir, args.tctemplate_long),
                doResp=self.doResp,
            )
        elif args.qmengine == "queue":
            qmEngine = qmengine.SbatchEngine(
                os.path.join(args.sampledir, args.tctemplate),
                os.path.join(args.sampledir, args.tctemplate_long),
                os.path.join(args.sampledir, args.sbatch),
                os.getenv("USER"),
                doResp=self.doResp,
            )
        elif args.qmengine == "chemcloud":
            qmEngine = qmengine.CCCloudEngine(
                os.path.join(args.sampledir, args.tctemplate),
                os.path.join(args.sampledir, args.tctemplate_long),
                doResp=self.doResp,
            )
        return qmEngine

    def initializeMMEngine(self, args):
        mmOptions = {}
        if args.conformers is None:
            mmOptions["start"] = args.start
            mmOptions["end"] = args.end
            mmOptions["split"] = args.split
            mmOptions["coordPath"] = os.path.join(args.dynamicsdir, args.coors)
        else:
            mmOptions["start"] = None
            mmOptions["end"] = None
            mmOptions["split"] = None
            mmOptions["coordPath"] = os.path.join(args.dynamicsdir, args.conformers)
        mmOptions["conformers"] = args.conformersPerSet
        mmOptions["nvalids"] = args.nvalids
        mmOptions["trainMdin"] = args.trainMdin
        mmOptions["validMdin"] = args.validMdin
        mmOptions["leap"] = "setup.leap"
        mmOptions["heatCounter"] = self.heatCounter
        if args.mmengine == "amber":
            mmEngine = mmengine.ExternalAmberEngine(mmOptions)
        return mmEngine

    def initialCycle(self):
        # Create initial target data from dynamics
        with open(os.path.join(self.optdir, self.args.opt0)) as f:
            for line in f.readlines():
                splitLine = line.split()
                if len(splitLine) > 1:
                    if splitLine[0] == "name":
                        initialTarget = splitLine[1]
        if not os.path.isdir(os.path.join(self.optdir, "targets")):
            os.mkdir(os.path.join(self.optdir, "targets"))
        path = os.path.join(self.optdir, "targets", initialTarget)
        if not os.path.isdir(path):
            os.mkdir(path)
        l = convertTCtoFB(
            os.path.join(self.args.dynamicsdir, self.args.tcout),
            os.path.join(self.args.dynamicsdir, self.args.coors),
            self.args.stride,
            self.args.start,
            self.args.end,
            os.path.join(path, "qdata.txt"),
            os.path.join(path, "all.mdcrd"),
        )
        for f in ["setup.leap", "conf.pdb", "setup_valid_initial.leap"]:
            copyfile(os.path.join(self.optdir, f), os.path.join(path, f))

        os.chdir(self.optdir)
        self.optEngine.optimizeForcefield(0)
        os.chdir(self.home)

    def doMMSampling(self, i):
        # Make sampling directory and copy files into it
        sampleName = f"{str(i)}_cycle_{str(i)}"
        samplePath = os.path.join(self.sampledir, sampleName)
        if not os.path.isdir(samplePath):
            os.mkdir(samplePath)
        elif self.restartCycle == -1:
            rmtree(samplePath)
            os.mkdir(samplePath)
        for f in os.listdir(os.path.join(self.optdir, "result", f"opt_{str(i-1)}")):
            copyfile(
                os.path.join(self.optdir, "result", "opt_" + str(i - 1), f),
                os.path.join(samplePath, f),
            )
        for f in ["conf.pdb", "setup.leap"]:
            copyfile(os.path.join(self.optdir, f), os.path.join(samplePath, f))
        for f in self.mdFiles:
            copyfile(os.path.join(self.sampledir, f), os.path.join(samplePath, f))
        os.chdir(samplePath)
        # Do MM sampling
        if i == self.restartCycle + 1:
            self.mmEngine.restart()
        else:
            self.mmEngine.getMMSamples()
        os.chdir(self.home)

    def doQMCalculations(self, i):
        # Run QM calculations for each sampling trajectory
        os.chdir(os.path.join(self.sampledir, f"{str(i)}_cycle_{str(i)}"))
        for f in os.listdir():
            if (f.startswith("train") or f.startswith("valid")) and os.path.isdir(f):
                os.chdir(f)
                if i == self.restartCycle + 1:
                    self.qmEngine.restart(".")
                else:
                    pdbs = []
                    for g in os.listdir():
                        if g.endswith(".pdb"):
                            pdbs.append(g)
                    self.qmEngine.getQMRefData(pdbs, ".")
                os.chdir("..")
        os.chdir(self.home)

    def doParameterOptimization(self, i):
        # Copy new QM data into appropriate folders
        trainFolder = os.path.join(self.optdir, "targets", f"train_{str(i)}")
        validFolder = os.path.join(self.optdir, "targets", f"valid_{str(i)}")
        if not os.path.isdir(trainFolder):
            os.mkdir(trainFolder)
        if not os.path.isdir(validFolder):
            os.mkdir(validFolder)
        for f in ["setup.leap", "conf.pdb", "setup_valid_initial.leap"]:
            copyfile(os.path.join(self.optdir, f), os.path.join(trainFolder, f))
            copyfile(os.path.join(self.optdir, f), os.path.join(validFolder, f))

        valids = []
        for f in os.listdir(os.path.join(self.sampledir, f"{str(i)}_cycle_{str(i)}")):
            if f.startswith("train") and os.path.isdir(
                os.path.join(self.sampledir, f"{str(i)}_cycle_{str(i)}", f)
            ):
                mmTrainFolder = os.path.join(
                    self.sampledir, f"{str(i)}_cycle_{str(i)}", f
                )
            elif f.startswith("valid") and os.path.isdir(
                os.path.join(self.sampledir, f"{str(i)}_cycle_{str(i)}", f)
            ):
                valids.append(
                    os.path.join(self.sampledir, f"{str(i)}_cycle_{str(i)}", f)
                )

        for f in ["all.mdcrd", "qdata.txt"]:
            copyfile(
                os.path.join(mmTrainFolder, f),
                os.path.join(self.optdir, "targets", f"train_{str(i)}", f),
            )
            copyfile(
                os.path.join(valids[0], f),
                os.path.join(validFolder, f),
            )
            for j in range(1, self.args.nvalids):
                validFolderJ = os.path.join(
                    self.optdir, "targets", f"valid_{str(i)}_{str(j)}"
                )
                if not os.path.isdir(validFolderJ):
                    os.mkdir(validFolderJ)
                copyfile(os.path.join(valids[j], f), os.path.join(validFolderJ, f))

        # Run ForceBalance on each input
        os.chdir(self.optdir)
        self.optEngine.optimizeForcefield(i)
        os.chdir(self.home)
        optResults = []
        optResults.append(self.optEngine.valid[-1])
        optResults.append(self.optEngine.valid[-1] / self.optEngine.validInitial[-1])
        optResults.append(self.optEngine.valid[-1] - self.optEngine.validPrevious[-1])
        if i > 1:
            try:
                optResults.append(self.optEngine.valid[-1] - self.optEngine.valid[-2])
            except:
                import pdb

                pdb.set_trace()

        return optResults
