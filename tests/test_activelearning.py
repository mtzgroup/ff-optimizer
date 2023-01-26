import os
import random
from shutil import copyfile, rmtree

import numpy as np
import pytest

from ff_optimizer import active_learning, model, optengine, qmengine, utils

from . import checkUtils


class FakeModel:
    def __init__(self, args):
        self.sampledir = args.sampledir
        self.mmEngine = FakeMMEngine()


class FakeMMEngine:
    def __init__(self, options={}):
        self.prmtop = "amber.prmtop"


class FakeArgs:
    def __init__(self):
        self.activeLearning = 3
        self.sampledir = "2_sampling"
        self.optdir = "1_opt"
        self.dynamicsdir = "0_dynamics"
        self.restart = False
        self.valid0 = "valid_0.in"
        self.opt0 = "opt_0.in"
        self.resp = 0
        self.respPriors = 0
        self.maxcycles = 100
        self.qmengine = "chemcloud"
        self.tctemplate = "tc.in"
        self.tctemplate_long = "tc.in"


def monkeyInit(self, args):
    self.nmodels = args.activeLearning
    self.models = [FakeModel(args) for i in range(self.nmodels)]
    self.symbols = None


@pytest.mark.amber
def test_computeEnergyForceSP(monkeypatch):
    os.chdir(os.path.join(os.path.dirname(__file__), "active_learning"))
    monkeypatch.setattr(active_learning.ActiveLearningModel, "__init__", monkeyInit)
    args = FakeArgs()
    model = active_learning.ActiveLearningModel(args)
    geometry = np.loadtxt(
        "coords.xyz", usecols=(1, 2, 3), skiprows=2, max_rows=84
    ).flatten()
    force = np.loadtxt(
        "MMforce.xyz", usecols=(1, 2, 3), skiprows=2, max_rows=84
    ).flatten()
    results = model.computeEnergyForce([geometry], "amber.prmtop")
    checkUtils.checkArrays(results[0][1], force)


def readXYZTraj(filename):
    frame = []
    frames = []
    i = 1
    natoms = -1
    with open(filename, "r") as f:
        for line in f.readlines():
            if natoms == -1:
                natoms = int(line.split()[0])
            if i > natoms + 2:
                frames.append(np.asarray(frame, dtype=np.float32))
                frame = []
                i = 1
            if i > 2:
                splitLine = line.split()
                frame.append(splitLine[1])
                frame.append(splitLine[2])
                frame.append(splitLine[3])
            i += 1
    return frames


@pytest.mark.amber
def test_computeEnergyForceAll(monkeypatch):
    os.chdir(os.path.join(os.path.dirname(__file__), "active_learning"))
    monkeypatch.setattr(active_learning.ActiveLearningModel, "__init__", monkeyInit)
    args = FakeArgs()
    model = active_learning.ActiveLearningModel(args)
    os.chdir(os.path.join(os.path.dirname(__file__), "active_learning"))
    frames = readXYZTraj("coords.xyz")
    results = model.computeEnergyForce(frames, "amber.prmtop")
    forces = readXYZTraj("MMforce.xyz")
    for i in range(len(frames)):
        checkUtils.checkArrays(results[i][1], forces[i])


def test_collectGeometries(monkeypatch):
    os.chdir(
        os.path.join(os.path.dirname(__file__), "active_learning", "collectGeometries")
    )
    monkeypatch.setattr(active_learning.ActiveLearningModel, "__init__", monkeyInit)
    args = FakeArgs()
    model = active_learning.ActiveLearningModel(args)

    geometries = model.collectGeometries(7, 0)
    assert len(geometries) == 30
    geometries = model.collectGeometries(7, 1)
    assert len(geometries) == 24
    geometries = model.collectGeometries(7, 2)
    assert len(geometries) == 27


@pytest.mark.amber
def test_computeAll(monkeypatch):
    monkeypatch.setattr(active_learning.ActiveLearningModel, "__init__", monkeyInit)
    os.chdir(
        os.path.join(os.path.dirname(__file__), "active_learning", "collectGeometries")
    )
    args = FakeArgs()
    model = active_learning.ActiveLearningModel(args)
    geometries = model.collectGeometries(7, 1)
    os.chdir(os.path.join("..", "computeAll"))
    prmtops = [f"prm{str(i)}.prmtop" for i in range(1, 4)]
    energies, forces = model.computeAll(geometries, prmtops)
    assert energies.shape == (3, 24)
    assert forces.shape == (3, 24, 18)


def test_chooseGeometries3Models(monkeypatch):
    monkeypatch.setattr(active_learning.ActiveLearningModel, "__init__", monkeyInit)
    os.chdir(
        os.path.join(os.path.dirname(__file__), "active_learning", "chooseGeometries")
    )
    args = FakeArgs()
    model = active_learning.ActiveLearningModel(args)
    energies = []
    forces = []
    for i in range(3):
        energies.append(np.loadtxt(f"{i}_energy.txt"))
        forces.append(np.loadtxt(f"{i}_force.txt"))
    energies = np.asarray(energies, dtype=np.float32)
    forces = np.asarray(forces, dtype=np.float32)
    random.seed(404)
    newGeoms = model.chooseGeometries(energies, forces)
    print(newGeoms)
    assert newGeoms == [6, 14, 18, 20, 7, 12, 15, 16]


def test_chooseGeometriesOdd(monkeypatch):
    monkeypatch.setattr(active_learning.ActiveLearningModel, "__init__", monkeyInit)
    os.chdir(
        os.path.join(os.path.dirname(__file__), "active_learning", "chooseGeometries")
    )
    args = FakeArgs()
    model = active_learning.ActiveLearningModel(args)
    energies = []
    forces = []
    for i in range(3):
        energies.append(np.loadtxt(f"{i}_energy.txt")[:21])
        forces.append(np.loadtxt(f"{i}_force.txt")[:21, :])
    energies = np.asarray(energies, dtype=np.float32)
    forces = np.asarray(forces, dtype=np.float32)
    random.seed(404)
    newGeoms = model.chooseGeometries(energies, forces)
    print(newGeoms)
    assert newGeoms == [6, 14, 18, 20, 7, 13, 5]


def test_chooseGeometries2Models(monkeypatch):
    monkeypatch.setattr(active_learning.ActiveLearningModel, "__init__", monkeyInit)
    os.chdir(
        os.path.join(os.path.dirname(__file__), "active_learning", "chooseGeometries")
    )
    args = FakeArgs()
    args.activeLearning = 2
    model = active_learning.ActiveLearningModel(args)
    energies = []
    forces = []
    for i in range(2):
        energies.append(np.loadtxt(f"{i}_energy.txt"))
        forces.append(np.loadtxt(f"{i}_force.txt"))
    energies = np.asarray(energies, dtype=np.float32)
    forces = np.asarray(forces, dtype=np.float32)
    random.seed(404)
    newGeoms = model.chooseGeometries(energies, forces)
    assert newGeoms == [6, 14, 18, 20, 3, 11, 7, 12, 15, 21, 19, 2]


def monkeyComputeAll(self, geometries, prmtops):
    for prmtop in prmtops:
        if not os.path.isfile(prmtop):
            raise RuntimeError(f"{prmtop} is missing!")
    return [0], []


def monkeyChooseGeometries(self, energies, forces):
    return energies


def test_doActiveLearning(monkeypatch):
    monkeypatch.setattr(active_learning.ActiveLearningModel, "__init__", monkeyInit)
    monkeypatch.setattr(
        active_learning.ActiveLearningModel, "computeAll", monkeyComputeAll
    )
    monkeypatch.setattr(
        active_learning.ActiveLearningModel, "chooseGeometries", monkeyChooseGeometries
    )
    os.chdir(
        os.path.join(os.path.dirname(__file__), "active_learning", "doActiveLearning")
    )
    args = FakeArgs()
    model = active_learning.ActiveLearningModel(args)
    model.prmtop = "amber.prmtop"
    folders = ["train", "valid_1", "valid_2"]
    for i in range(1, 4):
        for j in range(1, 4):
            copyfile(
                os.path.join("ref", f"{i}.xyz"),
                os.path.join(
                    f"model_{j}", "2_sampling", "7_cycle_7", folders[i - 1], f"1.xyz"
                ),
            )
    model.doActiveLearning(7)

    refGeom = utils.readXYZ(os.path.join("ref", "1.xyz"))
    testGeom = utils.readXYZ(
        os.path.join("model_1", "2_sampling", "7_cycle_7", "train", "1.xyz")
    )
    assert checkUtils.checkArrays(testGeom, refGeom)
    testGeom = utils.readXYZ(
        os.path.join("model_2", "2_sampling", "7_cycle_7", "valid_1", "1.xyz")
    )
    assert checkUtils.checkArrays(testGeom, refGeom)
    testGeom = utils.readXYZ(
        os.path.join("model_3", "2_sampling", "7_cycle_7", "valid_2", "1.xyz")
    )
    assert checkUtils.checkArrays(testGeom, refGeom)

    refGeom = utils.readXYZ(os.path.join("ref", "3.xyz"))
    testGeom = utils.readXYZ(
        os.path.join("model_2", "2_sampling", "7_cycle_7", "train", "1.xyz")
    )
    assert checkUtils.checkArrays(testGeom, refGeom)
    testGeom = utils.readXYZ(
        os.path.join("model_3", "2_sampling", "7_cycle_7", "valid_1", "1.xyz")
    )
    assert checkUtils.checkArrays(testGeom, refGeom)
    testGeom = utils.readXYZ(
        os.path.join("model_1", "2_sampling", "7_cycle_7", "valid_2", "1.xyz")
    )
    assert checkUtils.checkArrays(testGeom, refGeom)

    refGeom = utils.readXYZ(os.path.join("ref", "2.xyz"))
    testGeom = utils.readXYZ(
        os.path.join("model_3", "2_sampling", "7_cycle_7", "train", "1.xyz")
    )
    assert checkUtils.checkArrays(testGeom, refGeom)
    testGeom = utils.readXYZ(
        os.path.join("model_1", "2_sampling", "7_cycle_7", "valid_1", "1.xyz")
    )
    assert checkUtils.checkArrays(testGeom, refGeom)
    testGeom = utils.readXYZ(
        os.path.join("model_2", "2_sampling", "7_cycle_7", "valid_2", "1.xyz")
    )
    assert checkUtils.checkArrays(testGeom, refGeom)


def monkeyInitModel(self, args):
    self.restartCycle = random.randint(1, 10)
    options = {"nvalids": 1, "restart": self.restartCycle}

    class monkeyOpt:
        def __init__(self, options):
            self.options = options
            self.restartCycle = options["restart"]

        def determineRestart(self):
            pass

    self.optEngine = monkeyOpt(options)


def test_init(monkeypatch):
    monkeypatch.setattr(model.Model, "__init__", monkeyInitModel)
    random.seed(1015)
    args = FakeArgs()
    os.chdir(os.path.join(os.path.dirname(__file__), "active_learning", "init"))
    for i in range(1, 4):
        if os.path.isdir(f"model_{i}"):
            rmtree(f"model_{i}")
    mod = active_learning.ActiveLearningModel(args)
    for i in range(1, 4):
        os.remove(os.path.join(f"model_{i}", "1_opt", "leap.out"))
    os.remove(os.path.join("1_opt", "leap.out"))
    assert mod.restartCycle == 2
    for i in range(1, 4):
        assert os.path.isfile(os.path.join(f"model_{i}", "1_opt", "opt_0.in"))
        assert os.path.isfile(os.path.join(f"model_{i}", "2_sampling", "md.in"))


def monkeyInitOpt(self, args):
    self.train = []
    self.valid = []
    self.nvalids = args["nvalids"]
    self.validPrevious = []
    self.validInitial = []
    self.maxCycles = args["maxCycles"]
    self.optdir = args["optdir"]
    self.restartCycle = self.determineRestart()
    self.respPriors = None


def monkeyInitMM(self, args):
    pass


def monkeyInitQM(self, arg1, arg2, doResp):
    pass


def monkeyRename(a, b):
    pass


def test_restart(monkeypatch):
    monkeypatch.setattr(optengine.OptEngine, "__init__", monkeyInitOpt)
    monkeypatch.setattr(qmengine.CCCloudEngine, "__init__", monkeyInitQM)
    monkeypatch.setattr(model.Model, "initializeMMEngine", monkeyInitMM)
    monkeypatch.setattr(os, "rename", monkeyRename)

    args = FakeArgs()
    args.restart = True
    os.chdir(os.path.join(os.path.dirname(__file__), "active_learning", "restart"))
    mod = active_learning.ActiveLearningModel(args)
    os.remove(os.path.join("1_opt", "leap.out"))
    assert mod.restartCycle == 2


def monkeySetupFiles(self, i):
    pass


def monkeySortParams(self, results, i):
    pass


def monkeyGraphResults(self):
    pass


def monkeyOptInit(self, args):
    self.nvalids = args["nvalids"]
    for f in os.listdir(args["optdir"]):
        if f.endswith(".mol2"):
            self.mol2 = f
        elif f.endswith(".frcmod"):
            self.frcmod = f
    self.validPrevious = []
    self.train = []
    self.valid = []
    self.validInitial = []
    self.restartCycle = 0
    self.respPriors = None


def monkeyForceBalance(command):
    out = command.split()[3]
    out.split("_")[1].split(".")[0]
    copyfile(os.path.join("reference", out), out)


def monkeyInitialize(self, args):
    pass


def clean():
    for d in ["model_1", "model_2"]:
        os.chdir(os.path.join(d, "1_opt"))
        for f in os.listdir():
            if f.endswith(".out"):
                os.remove(f)
        os.chdir("targets")
        for f in os.listdir():
            rmtree(f)
        os.chdir(os.path.join("..", "..", ".."))


def monkeyALInit(self, args):
    self.home = os.getcwd()
    self.nmodels = args.activeLearning
    args.nvalids = args.activeLearning
    self.models = []
    for i in range(1, self.nmodels + 1):
        folder = f"model_{str(i)}"
        os.chdir(folder)
        self.models.append(model.Model(args))
        os.chdir("..")
        self.models[-1].optEngine.nvalids = 1
    self.restartCycle = 0
    self.nthreads = min(os.cpu_count(), self.nmodels)


def test_doParameterOptimization(monkeypatch):
    args = FakeArgs()
    args.activeLearning = 2
    args.restart = True
    os.chdir(
        os.path.join(
            os.path.dirname(__file__), "active_learning", "doParameterOptimization"
        )
    )
    clean()
    monkeypatch.setattr(model.Model, "initializeQMEngine", monkeyInitialize)
    monkeypatch.setattr(model.Model, "initializeMMEngine", monkeyInitialize)
    monkeypatch.setattr(optengine.OptEngine, "setupInputFiles", monkeySetupFiles)
    monkeypatch.setattr(optengine.OptEngine, "__init__", monkeyOptInit)
    monkeypatch.setattr(os, "system", monkeyForceBalance)
    monkeypatch.setattr(optengine.OptEngine, "sortParams", monkeySortParams)
    monkeypatch.setattr(optengine.OptEngine, "graphResults", monkeyGraphResults)
    monkeypatch.setattr(active_learning.ActiveLearningModel, "__init__", monkeyALInit)

    m = active_learning.ActiveLearningModel(args)
    m.doParameterOptimization(1)
    assert len(m.optResults) == 3
    m.doParameterOptimization(2)
    assert len(m.optResults) == 4
    clean()

#def test_restart(monkeypatch):
#    monkeypatch.setattr(optengine.OptEngine, "__init__", monkeyInitOpt)
#    monkeypatch.setattr(qmengine.CCCloudEngine, "__init__", monkeyInitQM)
#    monkeypatch.setattr(model.Model, "initializeMMEngine", monkeyInitMM)
#    monkeypatch.setattr(os, "rename", monkeyRename)
#
#    args = FakeArgs()
#    args.restart = True
#    os.chdir("/home/curtie/7_FF_fitting/0_ff-optimizer/4_alanine_tetrapeptide/5_40_ppc/6_active_learning_3_models")
#    mod = active_learning.ActiveLearningModel(args)
#    print(mod.restartCycle)
#    assert False
#
