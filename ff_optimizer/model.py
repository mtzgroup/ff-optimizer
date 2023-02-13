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


# Functions in this class assume that they are operating in the home directory
# where the job was started. 
class Model(AbstractModel):
    def __init__(self, args):
        self.setArgs(args)
        self.getMdFiles()
        # Set up each engine
        self.optEngine = self.initializeOptEngine(args)
        self.qmEngine = self.initializeQMEngine(args)
        self.mmEngine = self.initializeMMEngine(args)
        self.restartCycle = self.optEngine.restartCycle

    def getMdFiles(self):
        self.mdFiles = []
        self.heatCounter = 0
        for f in os.listdir(self.sampledir):
            if os.path.isfile(os.path.join(self.sampledir, f)):
                self.mdFiles.append(f)
                if f.startswith("heat"):
                    self.heatCounter += 1

    def setArgs(self, args):
        self.args = args
        # Set some miscellaneous variables
        self.home = os.getcwd()
        self.optdir = args.optdir
        self.sampledir = args.sampledir
        self.nvalids = args.nvalids
        os.rename(
            os.path.join(args.optdir, args.valid0),
            os.path.join(args.optdir, "valid_0.in"),
        )
        os.rename(
            os.path.join(args.optdir, args.opt0), os.path.join(args.optdir, "opt_0.in")
        )
        if args.resp != 0 or args.respPriors != 0:
            self.doResp = True
        else:
            self.doResp = False

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
                os.path.join(args.sampledir, args.tctemplate_backup),
                doResp=self.doResp,
            )
        elif args.qmengine == "queue":
            qmEngine = qmengine.SbatchEngine(
                os.path.join(args.sampledir, args.tctemplate),
                os.path.join(args.sampledir, args.tctemplate_backup),
                os.path.join(args.sampledir, args.sbatch),
                os.getenv("USER"),
                doResp=self.doResp,
            )
        elif args.qmengine == "chemcloud":
            qmEngine = qmengine.CCCloudEngine(
                os.path.join(args.sampledir, args.tctemplate),
                os.path.join(args.sampledir, args.tctemplate_backup),
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
        # Prepare initial target data
        path = self.createTCData()
        self.copyLeapFiles(path)

        # Do initial optimization
        os.chdir(self.optdir)
        self.optEngine.optimizeForcefield(0)
        os.chdir(self.home)

    def createTCData(self):
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
        return path

    def copyLeapFiles(self, dest, validInitial=True):
        files = ["setup.leap", "conf.pdb"]
        if validInitial:
            files += ["setup_valid_initial.leap"]
        for f in files:
            copyfile(os.path.join(self.optdir, f), os.path.join(dest, f))

    def doMMSampling(self, i):
        # Make sampling directory and copy files into it
        samplePath = self.makeSampleDir(i)
        self.copySamplingFiles(i, samplePath)
        # Do MM sampling
        os.chdir(samplePath)
        if i == self.restartCycle:
            self.mmEngine.restart()
        else:
            self.mmEngine.getMMSamples()
        os.chdir(self.home)

    def makeSampleDir(self, i):
        sampleName = f"{str(i)}_cycle_{str(i)}"
        samplePath = os.path.join(self.sampledir, sampleName)
        if not os.path.isdir(samplePath):
            os.mkdir(samplePath)
        elif self.restartCycle == -1:
            rmtree(samplePath)
            os.mkdir(samplePath)
        return samplePath

    def copyFFFiles(self, i, dest):
        for f in os.listdir(os.path.join(self.optdir, "result", f"opt_{str(i)}")):
            copyfile(
                os.path.join(self.optdir, "result", "opt_" + str(i), f),
                os.path.join(dest, f),
            )

    def copySamplingFiles(self, i, samplePath):
        self.copyFFFiles(i-1, samplePath)
        self.copyLeapFiles(samplePath, validInitial=False)
        for f in self.mdFiles:
            copyfile(os.path.join(self.sampledir, f), os.path.join(samplePath, f))

    def doQMCalculations(self, i):
        # Run QM calculations for each sampling trajectory
        os.chdir(os.path.join(self.sampledir, f"{str(i)}_cycle_{str(i)}"))
        for f in os.listdir():
            if (f.startswith("train") or f.startswith("valid")) and os.path.isdir(f):
                os.chdir(f)
                if i == self.restartCycle:
                    self.qmEngine.restart(".")
                else:
                    xyzs = self.getXYZs(".")
                    self.qmEngine.getQMRefData(xyzs, ".")
                os.chdir("..")
        os.chdir(self.home)

    def getXYZs(self, folder):
        xyzs = []
        for f in os.listdir(folder):
            if f.endswith(".xyz"):
                xyzs.append(f)
        return xyzs

    def makeFBTargets(self, i):
        if not os.path.isdir(os.path.join(self.optdir, "targets")):
            os.mkdir(os.path.join(self.optdir, "targets"))
        folders = [os.path.join(self.optdir, "targets", f"train_{str(i)}"), os.path.join(self.optdir, "targets", f"valid_{str(i)}")]
        for j in range(1, self.nvalids):
            folders.append(
                os.path.join(self.optdir, "targets", f"valid_{str(i)}_{str(j)}")
            )
        for f in folders:
            if not os.path.isdir(f):
                os.mkdir(f)
        return folders


    def doParameterOptimization(self, i):
        # Copy new QM data into appropriate folders
        targetFolders = self.makeFBTargets(i)
        for f in targetFolders:
            self.copyLeapFiles(f) 
        sampleFolders = self.getSampleFolders(i)
        self.copyQMResults(sampleFolders, targetFolders)

        # Run ForceBalance on each input
        os.chdir(self.optdir)
        self.optEngine.optimizeForcefield(i)
        os.chdir(self.home)
        self.getOptResults()

    # sampleFolders and targetFolders should be ordered such that the train
    # folder is index 0, with the valid folders in order behind it
    def copyQMResults(self, sampleFolders, targetFolders):
        for f in ["all.mdcrd", "qdata.txt"]:
            for i in range(len(sampleFolders)):
                copyfile(os.path.join(sampleFolders[i], f), os.path.join(targetFolders[i], f))
        
    def getSampleFolders(self, i):
        sampleFolders = [os.path.join(self.sampledir, f"{str(i)}_cycle_{str(i)}", "train")]
        for i in range(1, self.nvalids+1):
            sampleFolders.append(os.path.join(self.sampledir, f"{str(i)}_cycle_{str(i)}", f"valid_{i}"))
        for f in sampleFolders:
            if not os.path.isdir(f):
                raise RuntimeError(f"QM results folder {f} could not be found")
        return sampleFolders

    def getOptResults(self):
        self.optResults = []
        self.optResults.append(self.optEngine.valid[-1])
        self.optResults.append(
            self.optEngine.valid[-1] / self.optEngine.validInitial[-1]
        )
        self.optResults.append(
            self.optEngine.valid[-1] - self.optEngine.validPrevious[-1]
        )
        if len(self.optEngine.valid) > 1:
            self.optResults.append(self.optEngine.valid[-1] - self.optEngine.valid[-2])
