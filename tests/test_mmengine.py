import os
from pathlib import Path
from shutil import copyfile, rmtree

from ff_optimizer import mmengine, utils

from . import checkUtils
from .test_inputs import getDefaults


def monkeyGetIndices(self):
    return 1, 2, 3


def test_getIndices():
    pass


def test_checkForTCFormatting():
    options = getDefaults()
    os.chdir(os.path.join(os.path.dirname(__file__), "mmengine"))
    options.dynamicsdir = Path(".").absolute()
    options.coors = "coors.xyz"
    mmEngine = mmengine.MMEngine(options)
    test = mmEngine.checkForTCFormatting()
    assert test


def test_checkForTCFormatting2():
    options = getDefaults()
    os.chdir(os.path.join(os.path.dirname(__file__), "mmengine"))
    options.dynamicsdir = Path(".").absolute()
    options.coors = "coors2.xyz"
    mmEngine = mmengine.MMEngine(options)
    test = mmEngine.checkForTCFormatting()
    assert not test


def test_countFrames():
    options = getDefaults()
    os.chdir(os.path.join(os.path.dirname(__file__), "mmengine"))
    options.dynamicsdir = Path(".").absolute()
    options.coors = "coors.xyz"
    mmEngine = mmengine.MMEngine(options)
    test = mmEngine.countFrames()
    assert test == 1000


def test_countFrames2():
    options = getDefaults()
    os.chdir(os.path.join(os.path.dirname(__file__), "mmengine"))
    options.dynamicsdir = Path(".").absolute()
    options.coors = "coors2.xyz"
    mmEngine = mmengine.MMEngine(options)
    test = mmEngine.countFrames()
    assert test == 3


def test_getIndices():
    options = getDefaults()
    os.chdir(os.path.join(os.path.dirname(__file__), "mmengine"))
    options.dynamicsdir = Path(".").absolute()
    options.coors = "coors2.xyz"
    options.start = 0
    options.split = 1
    options.end = 7
    mmEngine = mmengine.MMEngine(options)
    start, end, split = mmEngine.getIndices()
    assert start == 0
    assert end == 2
    assert split == 1


def test_getIndices2():
    options = getDefaults()
    os.chdir(os.path.join(os.path.dirname(__file__), "mmengine"))
    options.dynamicsdir = Path(".").absolute()
    options.coors = "coors2.xyz"
    mmEngine = mmengine.MMEngine(options)
    start, end, split = mmEngine.getIndices()
    assert start == 0
    assert end == 2


def test_getIndices3():
    options = getDefaults()
    os.chdir(os.path.join(os.path.dirname(__file__), "mmengine"))
    options.dynamicsdir = Path(".").absolute()
    options.coors = "coors2.xyz"
    options.split = 3
    try:
        mmengine.MMEngine(options)
        test = False
    except:
        test = True
    assert test


def test_getIndices4():
    options = getDefaults()
    os.chdir(os.path.join(os.path.dirname(__file__), "mmengine"))
    options.dynamicsdir = Path(".").absolute()
    options.coors = "coors2.xyz"
    options.start = 3
    try:
        mmengine.MMEngine(options)
        test = False
    except:
        test = True
    assert test


def monkeyGetFrame(self, frame, dest):
    folder = dest.parent
    with open(os.path.join(folder, "frames.txt"), "a") as f:
        f.write(str(frame) + "\n")


def monkeyGetFrames(self):
    return [i for i in range((self.inp.nvalids + 1) * self.inp.conformersperset)]


def monkeySetup(self):
    pass


def test_getFrames1(monkeypatch):
    options = getDefaults()
    options.split = True
    os.chdir(os.path.join(os.path.dirname(__file__), "mmengine", "restart1"))
    options.dynamicsdir = Path("..").absolute()
    # monkeypatch.setattr(mmengine.MMEngine, "getIndices", monkeyGetIndices)
    mmEngine = mmengine.MMEngine(options)
    frames = mmEngine.getFrames()
    assert len(frames) == 2
    assert mmEngine.splitIndex > frames[0]
    assert mmEngine.splitIndex <= frames[1]


def test_getMMsamples(monkeypatch):
    options = getDefaults()
    os.chdir(os.path.join(os.path.dirname(__file__), "mmengine"))
    options.dynamicsdir = Path(".").absolute()
    if os.path.isdir("testSetup"):
        rmtree("testSetup")
    os.mkdir("testSetup")
    os.chdir("testSetup")

    monkeypatch.setattr(mmengine.MMEngine, "getFrame", monkeyGetFrame)
    monkeypatch.setattr(mmengine.MMEngine, "getFrames", monkeyGetFrames)
    # monkeypatch.setattr(mmengine.MMEngine, "getIndices", monkeyGetIndices)
    monkeypatch.setattr(mmengine.MMEngine, "setup", monkeySetup)
    mmEngine = mmengine.MMEngine(options)
    mmEngine.getMMSamples()

    with open(os.path.join("train", "frames.txt"), "r") as f:
        trainCrds = f.readlines()

    with open(os.path.join("valid_1", "frames.txt"), "r") as f:
        validCrds = f.readlines()

    os.chdir("..")
    # rmtree("testSetup")
    assert trainCrds == ["0\n"]
    assert validCrds == ["1\n"]


def test_getFrames2(monkeypatch):
    options = getDefaults()
    options.split = True
    os.chdir(os.path.join(os.path.dirname(__file__), "mmengine", "restart1"))
    options.dynamicsdir = Path("..").absolute()
    # monkeypatch.setattr(mmengine.MMEngine, "getIndices", monkeyGetIndices)
    mmEngine = mmengine.MMEngine(options)
    frames = mmEngine.getFrames()
    print(frames)
    assert len(frames) == 2
    assert mmEngine.splitIndex > frames[0]
    assert mmEngine.splitIndex <= frames[1]


def test_getMMsamples(monkeypatch):
    options = getDefaults()
    os.chdir(os.path.join(os.path.dirname(__file__), "mmengine"))
    options.dynamicsdir = Path(".").absolute()
    if os.path.isdir("testSetup"):
        rmtree("testSetup")
    os.mkdir("testSetup")
    os.chdir("testSetup")

    monkeypatch.setattr(mmengine.MMEngine, "getFrame", monkeyGetFrame)
    monkeypatch.setattr(mmengine.MMEngine, "getFrames", monkeyGetFrames)
    # monkeypatch.setattr(mmengine.MMEngine, "getIndices", monkeyGetIndices)
    monkeypatch.setattr(mmengine.MMEngine, "setup", monkeySetup)
    mmEngine = mmengine.MMEngine(options)
    mmEngine.getMMSamples()

    with open(os.path.join("train", "frames.txt"), "r") as f:
        trainCrds = f.readlines()

    with open(os.path.join("valid_1", "frames.txt"), "r") as f:
        validCrds = f.readlines()

    os.chdir("..")
    # rmtree("testSetup")
    assert trainCrds == ["0\n"]
    assert validCrds == ["1\n"]


# Also tests utils.writeRst()
def test_getFrame(monkeypatch):
    options = getDefaults()
    os.chdir(os.path.join(os.path.dirname(__file__), "mmengine"))
    options.dynamicsdir = Path(".").absolute()
    mmEngine = mmengine.MMEngine(options)

    mmEngine.getFrame(23, "temp.rst7")
    testCoors = []
    with open("temp.rst7", "r") as f:
        for line in f.readlines()[2:]:
            testCoors.append(line.split())
    os.remove("temp.rst7")
    refCoors = []
    with open("23.rst7", "r") as f:
        for line in f.readlines()[2:]:
            refCoors.append(line.split())
    checkUtils.checkArrays(testCoors, refCoors)


def monkeyGetSamples(self):
    return


def monkeySample(self, frames, mdin):
    with open("sample.txt", "w") as f:
        for frame in frames:
            f.write(f"{str(frame)}.rst7\n")
    return


def monkeySetup(self):
    return


# Tests restarting a completed directory
def test_restart1(monkeypatch):
    options = getDefaults()
    os.chdir(os.path.join(os.path.dirname(__file__), "mmengine", "restart1"))
    options.dynamicsdir = Path("..").absolute()
    monkeypatch.setattr(mmengine.MMEngine, "getIndices", monkeyGetIndices)
    monkeypatch.setattr(mmengine.MMEngine, "setup", monkeySetup)
    monkeypatch.setattr(mmengine.MMEngine, "getMMSamples", monkeyGetSamples)
    monkeypatch.setattr(mmengine.MMEngine, "sample", monkeySample)
    mmEngine = mmengine.MMEngine(options)

    mmEngine.restart()

    passTest = True
    for f in os.listdir():
        if os.path.isdir(f):
            if os.path.isfile(os.path.join(f, "sample.txt")):
                passTest = False
                os.remove(os.path.join(f, "sample.txt"))
                print("oops, sampling")
    if os.path.isfile(os.path.join("valid_1", "MMFinished.txt")):
        os.remove(os.path.join("valid_1", "MMFinished.txt"))
    else:
        print("oops, no mmfinished")
        passTest = False
    assert passTest


# Tests restarting a directory with one IC complete
# Tests restarting a directory with a partially complete IC
def test_restart2(monkeypatch):
    options = getDefaults()
    os.chdir(os.path.join(os.path.dirname(__file__), "mmengine", "restart2"))
    options.dynamicsdir = Path("..").absolute()
    monkeypatch.setattr(mmengine.MMEngine, "getIndices", monkeyGetIndices)
    monkeypatch.setattr(mmengine.MMEngine, "setup", monkeySetup)
    monkeypatch.setattr(mmengine.MMEngine, "getMMSamples", monkeyGetSamples)
    monkeypatch.setattr(mmengine.MMEngine, "sample", monkeySample)
    mmEngine = mmengine.MMEngine(options)

    mmEngine.restart()

    passTest = True
    if os.path.isfile(os.path.join("train", "sample.txt")):
        passTest = False
        os.remove(os.path.join("train", "sample.txt"))
    if os.path.isfile(os.path.join("valid_1", "sample.txt")):
        with open(os.path.join("valid_1", "sample.txt"), "r") as f:
            if f.readlines() != ["928.rst7\n"]:
                passTest = False
        os.remove(os.path.join("valid_1", "sample.txt"))
    else:
        passTest = False

    try:
        os.mkdir("valid_1")
    except:
        pass
    copyfile(os.path.join("..", "928.rst7"), os.path.join("valid_1", "928.rst7"))

    assert passTest


# Tests restarting a directory with the wrong number of rsts
def test_restart3(monkeypatch):
    options = getDefaults()
    monkeypatch.setattr(mmengine.MMEngine, "getIndices", monkeyGetIndices)
    monkeypatch.setattr(mmengine.MMEngine, "setup", monkeySetup)
    monkeypatch.setattr(mmengine.MMEngine, "sample", monkeySample)
    monkeypatch.setattr(mmengine.MMEngine, "getFrame", monkeyGetFrame)
    monkeypatch.setattr(mmengine.MMEngine, "getFrames", monkeyGetFrames)
    options.nvalids = 2
    options.conformersperset = 2
    if not os.path.isdir(
        os.path.join(os.path.dirname(__file__), "mmengine", "restart3")
    ):
        os.mkdir(os.path.join(os.path.dirname(__file__), "mmengine", "restart3"))
    os.chdir(os.path.join(os.path.dirname(__file__), "mmengine", "restart3"))
    mmEngine = mmengine.MMEngine(options)

    for f in os.listdir():
        utils.rmrf(f)
    train = Path("train")
    train.mkdir()
    with open(train / "7.rst7", "w") as f:
        f.write("hi!")
    mmEngine.restart()

    passTest = True
    if train.is_dir():
        if not (train / "frames.txt").is_file():
            passTest = False
            print("didn't sample")
        if (train / "7.rst7").is_file():
            print("didn't delete")
            passTest = False
        rmtree(train)
    else:
        passTest = False
    valid1 = Path("valid_1")
    if valid1.is_dir():
        if not (valid1 / "frames.txt").is_file():
            passTest = False
        rmtree(valid1)
    else:
        passTest = False
    valid2 = Path("valid_2")
    if valid2.is_dir():
        if not (valid2 / "frames.txt").is_file():
            passTest = False
        rmtree(valid2)
    else:
        passTest = False
    assert passTest


def test_getFrames3(monkeypatch):
    options = getDefaults()
    options.split = True
    options.conformersperset = 3
    options.nvalids = 2
    os.chdir(os.path.join(os.path.dirname(__file__), "mmengine", "restart1"))
    options.dynamicsdir = Path("..").absolute()
    # monkeypatch.setattr(mmengine.MMEngine, "getIndices", monkeyGetIndices)
    mmEngine = mmengine.MMEngine(options)
    frames = mmEngine.getFrames()
    assert len(frames) == 9
    for i in range(3):
        assert mmEngine.splitIndex > frames[i]
    for i in range(3, 9):
        assert mmEngine.splitIndex <= frames[i]


def test_getMMsamples2(monkeypatch):
    options = getDefaults()
    options.conformersperset = 3
    options.nvalids = 2
    os.chdir(os.path.join(os.path.dirname(__file__), "mmengine"))
    options.dynamicsdir = Path(".").absolute()
    if os.path.isdir("testSetup"):
        rmtree("testSetup")
    os.mkdir("testSetup")
    os.chdir("testSetup")

    monkeypatch.setattr(mmengine.MMEngine, "getFrame", monkeyGetFrame)
    monkeypatch.setattr(mmengine.MMEngine, "getFrames", monkeyGetFrames)
    # monkeypatch.setattr(mmengine.MMEngine, "getIndices", monkeyGetIndices)
    monkeypatch.setattr(mmengine.MMEngine, "setup", monkeySetup)
    mmEngine = mmengine.MMEngine(options)
    mmEngine.getMMSamples()

    with open(os.path.join("train", "frames.txt"), "r") as f:
        trainCrds = f.readlines()
    validCrds = []
    with open(os.path.join("valid_1", "frames.txt"), "r") as f:
        validCrds.append(f.readlines())
    with open(os.path.join("valid_2", "frames.txt"), "r") as f:
        validCrds.append(f.readlines())

    os.chdir("..")
    rmtree("testSetup")
    print(trainCrds)
    print(validCrds[0])
    print(validCrds[1])
    assert trainCrds == ["0\n", "1\n", "2\n"]
    assert validCrds[0] == ["3\n", "4\n", "5\n"]
    assert validCrds[1] == ["6\n", "7\n", "8\n"]
