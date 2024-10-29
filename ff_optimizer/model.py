import os
from pathlib import Path
from shutil import copyfile, rmtree

from . import mmengine, optengine, qmengine
from .inputs import Input
from .utils import convertTCtoFB, getXYZs


# Template class for ff models used by ff_optimizer
# All these functions must be implemented and
# all these variables set properly for ff_optimizer to
# use the inherited class
class AbstractModel:
    def __init__(self):
        """
        Initialize the AbstractModel.
        """
        # self.restartCycle = -1
        # raise NotImplementedError
        # commented in case using super becomes necessary

    def initialCycle(self):
        """
        Perform the initial cycle of the model.
        This method should be implemented by subclasses.
        """
        raise NotImplementedError()

    def doMMSampling(self, i: int):
        """
        Perform MM sampling for a given iteration.
        This method should be implemented by subclasses.

        Args:
            i (int): The iteration number.
        """
        raise NotImplementedError()

    def doQMCalculations(self, i: int):
        """
        Perform QM calculations for a given iteration.
        This method should be implemented by subclasses.

        Args:
            i (int): The iteration number.
        """
        raise NotImplementedError()

    def doParameterOptimization(self, i: int):
        """
        Perform parameter optimization for a given iteration.
        This method should be implemented by subclasses.

        Args:
            i (int): The iteration number.
        """
        raise NotImplementedError()


# Functions in this class assume that they are operating in the home directory
# where the job was started.
class Model(AbstractModel):
    def __init__(self, inp: Input):
        """
        Initialize the Model.

        Args:
            inp (Input): Input object containing configuration parameters.
        """
        self.setParams(inp)
        self.getMDFiles(inp)
        # Set up each engine
        self.optEngine = self.initializeOptEngine(inp)
        self.qmEngine = self.initializeQMEngine(inp)
        self.mmEngine = self.initializeMMEngine(inp)
        self.restartCycle = self.optEngine.restartCycle

    def getMDFiles(self, inp: Input):
        """
        Get MD files and count heat files.

        Args:
            inp (Input): Input object containing configuration parameters.
        """
        self.mdFiles = []
        self.heatCounter = 0
        for f in os.listdir(self.sampledir):
            if os.path.isfile(os.path.join(self.sampledir, f)):
                self.mdFiles.append(f)
                if f.startswith("heat"):
                    self.heatCounter += 1
        # Easiest way to pass this to MMEngine is through the Inputs class
        inp.heatCounter = self.heatCounter

    def setParams(self, inp: Input):
        """
        Set parameters for the model.

        Args:
            inp (Input): Input object containing configuration parameters.
        """
        # Set some miscellaneous variables
        self.home = Path(".").cwd().absolute()
        self.optdir = inp.optdir
        self.sampledir = inp.sampledir
        self.dynamicsdir = inp.dynamicsdir
        self.inp = inp
        self.converged = False

    def initializeOptEngine(self, inp: Input) -> optengine.OptEngine:
        """
        Initialize the optimization engine.

        Args:
            inp (Input): Input object containing configuration parameters.

        Returns:
            OptEngine: The initialized optimization engine.
        """
        optEngine = optengine.OptEngine(inp)
        return optEngine

    def initializeQMEngine(self, inp: Input) -> qmengine.QMEngine:
        """
        Initialize the QM engine.

        Args:
            inp (Input): Input object containing configuration parameters.

        Returns:
            QMEngine: The initialized QM engine.
        """
        if inp.qmengine == "debug":
            qmEngine = qmengine.DebugEngine(inp)
        elif inp.qmengine == "slurm":
            qmEngine = qmengine.SlurmEngine(inp)
        elif inp.qmengine == "chemcloud":
            qmEngine = qmengine.ChemcloudEngine(inp)
        return qmEngine

    def initializeMMEngine(self, inp: Input) -> mmengine.MMEngine:
        """
        Initialize the MM engine.

        Args:
            inp (Input): Input object containing configuration parameters.

        Returns:
            MMEngine: The initialized MM engine.
        """
        if inp.mmengine == "amber":
            mmEngine = mmengine.ExternalAmberEngine(inp)
        if inp.mmengine == "openmm":
            mmEngine = mmengine.ExternalOpenMMEngine(inp)
        return mmEngine

    def initialCycle(self):
        """
        Perform the initial cycle of the model.
        """
        if self.inp.initialtraining:
            # Prepare initial target data
            path = self.createTCData()
            self.copyLeapFiles(path)
            # Do initial optimization
            os.chdir(self.optdir)
            self.optEngine.optimizeForcefield(0)
            os.chdir(self.home)
        else:
            # just put ff parameters where the next iteration expects them to be
            path = self.optdir / "result" / "opt_0"
            path.mkdir(exist_ok=True, parents=True)
            copyfile(self.optdir / self.optEngine.mol2, path / self.optEngine.mol2)
            copyfile(self.optdir / self.optEngine.frcmod, path / self.optEngine.frcmod)

    def createTCData(self):
        """
        Create initial target data from dynamics.

        Returns:
            Path: Path to the created target data.
        """
        with open(self.optdir / "opt_0.in", "r") as f:
            for line in f.readlines():
                splitLine = line.split()
                if len(splitLine) > 1:
                    if splitLine[0] == "name":
                        initialTarget = splitLine[1]
        path = self.inp.optdir / "targets" / initialTarget
        path.mkdir(parents=True, exist_ok=True)
        l = convertTCtoFB(
            self.dynamicsdir / self.inp.tcout,
            self.dynamicsdir / self.inp.coors,
            self.inp.stride,
            self.inp.start,
            self.inp.end,
            path / "qdata.txt",
            path / "all.mdcrd",
        )
        return path

    def copyLeapFiles(self, dest, validInitial=True):
        """
        Copy leap files to the destination.

        Args:
            dest (Path): Destination path for the leap files.
            validInitial (bool): Whether to copy the valid initial leap file.
        """
        files = ["setup.leap", "conf.pdb"]
        # default is to copy the valid initial leap if the input wants to
        if validInitial and self.inp.validinitial:
            files += ["setup_valid_initial.leap"]
        for f in files:
            copyfile(self.optdir / f, dest / f)

    def doMMSampling(self, i):
        """
        Perform MM sampling for a given iteration.

        Args:
            i (int): The iteration number.
        """
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
        """
        Create a sampling directory for a given iteration.

        Args:
            i (int): The iteration number.

        Returns:
            Path: Path to the created sampling directory.
        """
        sampleName = f"{str(i)}_cycle_{str(i)}"
        samplePath = self.sampledir / sampleName
        if not samplePath.exists():
            samplePath.mkdir()
        elif self.restartCycle == -1:
            rmtree(samplePath, ignore_errors=True)
            rmtree(samplePath, ignore_errors=True)
            samplePath.mkdir(exist_ok=True)
        return samplePath

    def copyFFFiles(self, i: int, dest: Path):
        """
        Copy force field files to the sampling directory.

        Args:
            i (int): The iteration number.
            dest (Path): Path to the destination directory.
        """
        resultPath = self.optdir / "result" / f"opt_{i}"
        for f in resultPath.iterdir():
            copyfile(f, dest / f.name)

    def copySamplingFiles(self, i: int, samplePath: Path):
        """
        Copy sampling files for a given iteration.

        Args:
            i (int): The iteration number.
            samplePath (Path): Path to the sampling directory.
        """
        self.copyFFFiles(i - 1, samplePath)
        self.copyLeapFiles(samplePath, validInitial=False)
        for f in self.mdFiles:
            copyfile(self.sampledir / f, samplePath / f)

    def doQMCalculations(self, i: int):
        """
        Perform QM calculations for a given iteration.

        Args:
            i (int): The iteration number.
        """
        # Run QM calculations for each sampling trajectory
        for f in (self.sampledir / f"{str(i)}_cycle_{str(i)}").iterdir():
            if (
                f.name.startswith("train") or f.name.startswith("valid")
            ) and f.is_dir():
                os.chdir(f)
                if i == self.restartCycle:
                    self.qmEngine.restart()
                else:
                    xyzs = getXYZs()
                    self.qmEngine.getQMRefData(xyzs)
                os.chdir("..")
            os.chdir(self.home)

    def makeFBTargets(self, i: int) -> list:
        """
        Create ForceBalance target folders for a given iteration.

        Args:
            i (int): The iteration number.

        Returns:
            list: List of created target folder paths.
        """
        targets = self.optdir / "targets"
        if not targets.exists():
            targets.mkdir()
        folders = [targets / f"train_{i}", targets / f"valid_{i}"]
        for j in range(1, self.inp.nvalids):
            folders.append(targets / f"valid_{i}_{j}")
        for f in folders:
            if not f.is_dir():
                f.mkdir()
        return folders

    def doParameterOptimization(self, i: int):
        """
        Perform parameter optimization for a given iteration.

        Args:
            i (int): The iteration number.
        """
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
    def copyQMResults(self, sampleFolders: list, targetFolders: list):
        """
        Copy QM results from sample folders to target folders.

        Args:
            sampleFolders (list): List of sample folder paths.
            targetFolders (list): List of target folder paths.
        """
        for f in ["all.mdcrd", "qdata.txt"]:
            for i in range(len(sampleFolders)):
                copyfile(sampleFolders[i] / f, targetFolders[i] / f)

    def getSampleFolders(self, i: int) -> list:
        """
        Get sample folders for a given iteration.

        Args:
            i (int): The iteration number.

        Returns:
            list: List of sample folder paths.

        Raises:
            RuntimeError: If a QM results folder cannot be found.
        """
        sampleFolders = [self.sampledir / f"{str(i)}_cycle_{str(i)}" / "train"]
        for j in range(1, self.inp.nvalids + 1):
            sampleFolders.append(
                self.sampledir / f"{str(i)}_cycle_{str(i)}" / f"valid_{j}"
            )
        for f in sampleFolders:
            if not f.is_dir():
                raise RuntimeError(f"QM results folder {f} could not be found")
        return sampleFolders

    def getOptResults(self):
        """
        Get optimization results and update the model's state.
        """
        self.converged = self.optEngine.converged
        self.optResults = []
        self.optResults.append(self.optEngine.train[-1])
        self.optResults.append(self.optEngine.valid[-1])
        self.optResults.append(
            (self.optEngine.validPrevious[-1] - self.optEngine.valid[-1])
            / self.optEngine.valid[-1]
            * 100
        )
        if self.inp.validinitial:
            self.optResults.append(
                self.optEngine.valid[-1] / self.optEngine.validInitial[-1]
            )
        # if len(self.optEngine.valid) > 1:
        #    self.optResults.append(self.optEngine.valid[-1] - self.optEngine.valid[-2])
