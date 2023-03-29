import os
from pathlib import Path
from shutil import copyfile, rmtree

from ff_optimizer import model, optengine

from . import checkUtils

home = Path(__name__).parent.absolute()


class FakeModel:
    def __init__(self, args):
        self.sampledir = Path(args.sampledir)
        self.mmEngine = FakeMMEngine()


class FakeMMEngine:
    def __init__(self, options={}):
        self.prmtop = "amber.prmtop"


class FakeArgs:
    def __init__(self):
        self.activeLearning = 1
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
        self.nvalids = 2
        self.tcout = "tc.out"
        self.coors = "coors.xyz"
        self.stride = 100
        self.start = 0
        self.end = -1


def monkeyQMRestart():
    pass


def monkeyInit(self, args):
    pass


class MonkeyQMEngine:
    def __init__(self, args):
        self.xyzs = []
        self.dirs = []

    def getQMRefData(self, xyzs):
        self.xyzs.append(xyzs)
        self.dirs.append(os.getcwd().split("/")[-1])

    def restart(self):
        self.dirs.append(os.getcwd().split("/")[-1])


def monkeyInitQM(self, args):
    return MonkeyQMEngine(args)


def monkeyInitModel(self, args):
    self.home = os.getcwd()
    self.sampledir = Path(args.sampledir)
    self.qmEngine = self.initializeQMEngine(args)
    self.restartCycle = 1


def test_doQMCalculations(monkeypatch):
    args = FakeArgs()
    os.chdir(home / "model" / "test1")
    monkeypatch.setattr(model.Model, "initializeQMEngine", monkeyInitQM)
    monkeypatch.setattr(model.Model, "__init__", monkeyInitModel)

    m = model.Model(args)
    m.doQMCalculations(3)
    testXyzs = {}
    testXyzs["train"] = ["1.xyz", "2.xyz"]
    testXyzs["valid"] = ["3.xyz", "4.xyz"]
    testXyzs["valid_2"] = ["5.xyz", "6.xyz"]
    assert len(m.qmEngine.dirs) == 3
    for i in range(len(m.qmEngine.dirs)):
        for xyz in m.qmEngine.xyzs[i]:
            assert xyz.name in testXyzs[m.qmEngine.dirs[i]]


def test_doQMCalculationsRestart(monkeypatch):
    args = FakeArgs()
    os.chdir(home / "model" / "test1")
    monkeypatch.setattr(model.Model, "initializeQMEngine", monkeyInitQM)
    monkeypatch.setattr(model.Model, "__init__", monkeyInitModel)

    m = model.Model(args)
    m.restartCycle = 2
    m.doQMCalculations(3)
    assert len(m.qmEngine.dirs) == 3


class MonkeyMMEngine:
    def __init__(self, args):
        self.didRestart = False

    def restart(self):
        self.didRestart = True

    def getMMSamples(self):
        pass


class MonkeyOptEngine:
    def __init__(self, args):
        self.restartCycle = 2


def monkeyInitOpt(self, args):
    return MonkeyOptEngine(args)


def monkeyInitMM(self, args):
    return MonkeyMMEngine(args)


def test_doMMsampling(monkeypatch):
    args = FakeArgs()
    os.chdir(home / "model" / "test2")
    path = os.path.join("2_sampling", "3_cycle_3")
    if os.path.isdir(path):
        rmtree(path)
    monkeypatch.setattr(model.Model, "initializeOptEngine", monkeyInitOpt)
    monkeypatch.setattr(model.Model, "initializeQMEngine", monkeyInit)
    monkeypatch.setattr(model.Model, "initializeMMEngine", monkeyInitMM)
    m = model.Model(args)
    m.restartCycle = 3
    m.doMMSampling(3)
    isDir = os.path.isdir(path)
    isFile = True
    for f in ["conf.pdb", "setup.leap", "md.in", "heat1.in", "sol.mol2", "sol.frcmod"]:
        if not os.path.isfile(os.path.join(path, f)):
            isFile = False
    with open(os.path.join(path, "sol.frcmod"), "r") as f:
        testLines = f.readlines()
    with open(os.path.join("1_opt", "result", "opt_2", "sol.frcmod"), "r") as f:
        refLines = f.readlines()
    rmtree(path)
    assert isFile
    assert isDir
    assert m.mmEngine.didRestart
    assert checkUtils.checkLists(testLines, refLines)


def test_doMMsamplingNoRestart(monkeypatch):
    args = FakeArgs()
    os.chdir(home / "model" / "test2")
    path = os.path.join("2_sampling", "3_cycle_3")
    if os.path.isdir(path):
        rmtree(path)
    monkeypatch.setattr(model.Model, "initializeOptEngine", monkeyInitOpt)
    monkeypatch.setattr(model.Model, "initializeQMEngine", monkeyInit)
    monkeypatch.setattr(model.Model, "initializeMMEngine", monkeyInitMM)
    m = model.Model(args)
    m.restartCycle = -1
    m.doMMSampling(3)
    rmtree(path)
    assert not m.mmEngine.didRestart


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


def clean():
    os.chdir("1_optimization")
    for f in os.listdir():
        if f.endswith(".out"):
            os.remove(f)
    if os.path.isdir("targets"):
        os.chdir("targets")
        for f in os.listdir():
            rmtree(f)
        os.chdir(os.path.join("..", ".."))
    else:
        os.chdir("..")


def test_doParameterOptimization(monkeypatch):
    args = FakeArgs()
    args.optdir = "1_optimization"
    args.sampledir = "2_mm_sampling"
    args.nvalids = 2
    os.chdir(home / "model" / "test3")
    clean()
    monkeypatch.setattr(model.Model, "initializeQMEngine", monkeyInit)
    monkeypatch.setattr(model.Model, "initializeMMEngine", monkeyInit)
    monkeypatch.setattr(optengine.OptEngine, "setupInputFiles", monkeySetupFiles)
    monkeypatch.setattr(optengine.OptEngine, "__init__", monkeyOptInit)
    monkeypatch.setattr(os, "system", monkeyForceBalance)
    monkeypatch.setattr(optengine.OptEngine, "sortParams", monkeySortParams)
    monkeypatch.setattr(optengine.OptEngine, "graphResults", monkeyGraphResults)

    m = model.Model(args)
    m.doParameterOptimization(1)
    assert len(m.optResults) == 3
    m.doParameterOptimization(2)
    assert len(m.optResults) == 4
    files = [
        "setup.leap",
        "setup_valid_initial.leap",
        "conf.pdb",
        "all.mdcrd",
        "qdata.txt",
    ]
    copiedFiles = True
    for f in files:
        for target in ["train_1", "valid_1", "valid_1_1"]:
            if not os.path.isfile(os.path.join("1_optimization", "targets", target, f)):
                copiedFiles = False
    assert copiedFiles
    testDir = os.path.join("1_optimization", "targets", "train_1")
    refDir = os.path.join("2_mm_sampling", "1_cycle_1", "train")
    assert checkUtils.checkFiles(
        os.path.join(testDir, "all.mdcrd"), os.path.join(refDir, "all.mdcrd")
    )
    assert checkUtils.checkFiles(
        os.path.join(testDir, "qdata.txt"), os.path.join(refDir, "qdata.txt")
    )
    testDir = os.path.join("1_optimization", "targets", "valid_1")
    refDir = os.path.join("2_mm_sampling", "1_cycle_1", "valid_1")
    assert checkUtils.checkFiles(
        os.path.join(testDir, "all.mdcrd"), os.path.join(refDir, "all.mdcrd")
    )
    assert checkUtils.checkFiles(
        os.path.join(testDir, "qdata.txt"), os.path.join(refDir, "qdata.txt")
    )
    testDir = os.path.join("1_optimization", "targets", "valid_1_1")
    refDir = os.path.join("2_mm_sampling", "1_cycle_1", "valid_2")
    assert checkUtils.checkFiles(
        os.path.join(testDir, "all.mdcrd"), os.path.join(refDir, "all.mdcrd")
    )
    assert checkUtils.checkFiles(
        os.path.join(testDir, "qdata.txt"), os.path.join(refDir, "qdata.txt")
    )
    clean()


def test_getMdFiles(monkeypatch):
    args = FakeArgs()
    os.chdir(home / "model" / "test2")
    monkeypatch.setattr(model.Model, "initializeOptEngine", monkeyInitOpt)
    monkeypatch.setattr(model.Model, "initializeQMEngine", monkeyInit)
    monkeypatch.setattr(model.Model, "initializeMMEngine", monkeyInitMM)
    m = model.Model(args)
    m.getMDFiles()
    assert "md.in" in m.mdFiles
    assert "heat1.in" in m.mdFiles
    assert m.heatCounter == 1


def monkeyInitSetArgsOnly(self, args):
    self.setArgs(args)


def test_setArgs(monkeypatch):
    args = FakeArgs()
    args.respPriors = 1
    args.valid0 = "v0.in"
    args.opt0 = "o0.in"
    os.chdir(home / "model" / "test4")
    v0 = Path("1_opt") / "valid_0.in"
    o0 = Path("1_opt") / "opt_0.in"
    if v0.is_file():
        v0.rename(Path("1_opt") / "v0.in")
    if o0.is_file():
        o0.rename(Path("1_opt") / "o0.in")
    monkeypatch.setattr(model.Model, "__init__", monkeyInitSetArgsOnly)
    m = model.Model(args)
    moved = True
    if v0.is_file():
        v0.rename(Path("1_opt") / "v0.in")
    else:
        moved = False
    if o0.is_file():
        o0.rename(Path("1_opt") / "o0.in")
    else:
        moved = False
    assert moved
    assert m.doResp


def monkeyConvert(a, b, c, d, e, qdata, mdcrd):
    ref = Path("ref")
    copyfile(ref / "qdata.txt", qdata)
    copyfile(ref / "all.mdcrd", mdcrd)
    return 1


def test_createTCData(monkeypatch):
    args = FakeArgs()
    os.chdir(home / "model" / "test5")
    monkeypatch.setattr(model.Model, "initializeOptEngine", monkeyInitOpt)
    monkeypatch.setattr(model.Model, "initializeQMEngine", monkeyInit)
    monkeypatch.setattr(model.Model, "initializeMMEngine", monkeyInitMM)
    monkeypatch.setattr(model, "convertTCtoFB", monkeyConvert)
    m = model.Model(args)
    m.createTCData()
    target = Path("1_opt") / "targets" / "initial"
    hasFiles = True
    if not (target / "qdata.txt").is_file():
        hasFiles = False
    if not (target / "all.mdcrd").is_file():
        hasFiles = False
    rmtree(target.parent)
    assert hasFiles


def test_copyLeapFiles(monkeypatch):
    args = FakeArgs()
    os.chdir(home / "model" / "test5")
    monkeypatch.setattr(model.Model, "initializeOptEngine", monkeyInitOpt)
    monkeypatch.setattr(model.Model, "initializeQMEngine", monkeyInit)
    monkeypatch.setattr(model.Model, "initializeMMEngine", monkeyInitMM)
    m = model.Model(args)
    dest = Path("test1")
    dest.mkdir()
    m.copyLeapFiles(dest)
    test1 = True
    for f in ["setup.leap", "conf.pdb", "setup_valid_initial.leap"]:
        if not (dest / f).is_file():
            test1 = False
    rmtree(dest)
    dest = Path("test2")
    dest.mkdir()
    m.copyLeapFiles(dest, False)
    test2 = True
    for f in ["setup.leap", "conf.pdb"]:
        if not (dest / f).is_file():
            test2 = False
    rmtree(dest)

    assert test1
    assert test2


def test_makeSampleDir(monkeypatch):
    args = FakeArgs()
    os.chdir(home / "model" / "test5")
    monkeypatch.setattr(model.Model, "initializeOptEngine", monkeyInitOpt)
    monkeypatch.setattr(model.Model, "initializeQMEngine", monkeyInit)
    monkeypatch.setattr(model.Model, "initializeMMEngine", monkeyInitMM)
    m = model.Model(args)
    m.restartCycle = 3
    m.makeSampleDir(7)
    samplePath = m.sampledir / "7_cycle_7"
    test1 = samplePath.is_dir()
    with open(samplePath / "hi.txt", "w") as f:
        f.write("hi")
    m.makeSampleDir(7)
    test2 = (samplePath / "hi.txt").is_file()
    m.restartCycle = -1
    m.makeSampleDir(7)
    test3 = not (samplePath / "hi.txt").is_file()
    assert test1
    assert test2
    assert test3


def test_copyFFFiles(monkeypatch):
    args = FakeArgs()
    os.chdir(home / "model" / "test5")
    monkeypatch.setattr(model.Model, "initializeOptEngine", monkeyInitOpt)
    monkeypatch.setattr(model.Model, "initializeQMEngine", monkeyInit)
    monkeypatch.setattr(model.Model, "initializeMMEngine", monkeyInitMM)
    m = model.Model(args)
    test = Path("test")
    test.mkdir()
    m.copyFFFiles(1, test)
    copied = True
    for f in ["test.frcmod", "test.mol2"]:
        if not (test / f).is_file():
            copied = False
    rmtree(test)
    assert copied


def test_copySamplingFiles(monkeypatch):
    args = FakeArgs()
    os.chdir(home / "model" / "test5")
    monkeypatch.setattr(model.Model, "initializeOptEngine", monkeyInitOpt)
    monkeypatch.setattr(model.Model, "initializeQMEngine", monkeyInit)
    monkeypatch.setattr(model.Model, "initializeMMEngine", monkeyInitMM)
    m = model.Model(args)
    test = Path("test")
    test.mkdir()
    m.copySamplingFiles(2, test)
    copied = True
    for f in [
        "test.frcmod",
        "test.mol2",
        "setup.leap",
        "conf.pdb",
        "md.in",
        "heat1.in",
    ]:
        if not (test / f).is_file():
            copied = False
    rmtree(test)
    assert copied


def test_makeFBTargets(monkeypatch):
    args = FakeArgs()
    args.nvalids = 3
    os.chdir(home / "model" / "test5")
    monkeypatch.setattr(model.Model, "initializeOptEngine", monkeyInitOpt)
    monkeypatch.setattr(model.Model, "initializeQMEngine", monkeyInit)
    monkeypatch.setattr(model.Model, "initializeMMEngine", monkeyInitMM)
    m = model.Model(args)
    m.makeFBTargets(3)
    testFolders = ["train_3", "valid_3", "valid_3_1", "valid_3_2"]
    isMade = True
    path = m.optdir / "targets"
    for f in testFolders:
        if not (path / f).is_dir():
            isMade = False
    rmtree(path)
    assert isMade


def test_copyQMResults(monkeypatch):
    args = FakeArgs()
    args.nvalids = 3
    os.chdir(home / "model" / "test6")
    monkeypatch.setattr(model.Model, "initializeOptEngine", monkeyInitOpt)
    monkeypatch.setattr(model.Model, "initializeQMEngine", monkeyInit)
    monkeypatch.setattr(model.Model, "initializeMMEngine", monkeyInitMM)
    m = model.Model(args)
    dests = m.makeFBTargets(7)
    srcs = m.getSampleFolders(7)
    m.copyQMResults(srcs, dests)
    copiedCorrectly = True
    for d in dests:
        for f in ["all.mdcrd", "qdata.txt"]:
            with open(d / f, "r") as g:
                name = g.readline().split()[0]
            if name != d.name:
                copiedCorrectly = False
    rmtree(m.optdir / "targets")
    assert copiedCorrectly


def test_getSampleFolders(monkeypatch):
    args = FakeArgs()
    args.nvalids = 3
    os.chdir(home / "model" / "test6")
    monkeypatch.setattr(model.Model, "initializeOptEngine", monkeyInitOpt)
    monkeypatch.setattr(model.Model, "initializeQMEngine", monkeyInit)
    monkeypatch.setattr(model.Model, "initializeMMEngine", monkeyInitMM)
    m = model.Model(args)
    sampleFolders = m.getSampleFolders(7)
    path = m.sampledir / "7_cycle_7"
    refFolders = [path / "train", path / "valid_1", path / "valid_2", path / "valid_3"]
    try:
        m.getSampleFolders(6)
        failed = False
    except:
        failed = True
    assert failed
    assert sampleFolders == refFolders


# Gonna add tests with feature to change initialCycle
def test_initialCycle():
    assert False
