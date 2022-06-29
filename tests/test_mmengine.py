import os
from shutil import rmtree

from numpy import zeros

from ff_optimizer import mmengine

from . import checkUtils

options = {}
options["start"] = 33
options["end"] = 2000
options["split"] = 777
options["stride"] = 50
options["nvalids"] = 1
options["trainMdin"] = "md.in"
options["validMdin"] = "md.in"
options["leap"] = "setup.leap"
options["conformers"] = 1


def monkeyGetIndices(self):
    return 1, 2, 3


def test_getIndices():
    pass

def monkeyGetFrame(self, frame, dest):
    folder = dest.split('/')[0]
    with open(os.path.join(folder,"frames.txt"),'a') as f:
        f.write(str(frame) + "\n")

def monkeyGetFrames(self):
    return [i for i in range((options['nvalids'] + 1) * options['conformers'])] 

def monkeySetup(self):
    pass

def test_getFrames(monkeypatch):
    coordPath = os.path.join("..", "coors.xyz")
    options["coordPath"] = coordPath
    os.chdir(os.path.join(os.path.dirname(__file__), "mmengine", "restart1"))
    #monkeypatch.setattr(mmengine.MMEngine, "getIndices", monkeyGetIndices)
    mmEngine = mmengine.MMEngine(options)
    frames = mmEngine.getFrames()
    assert len(frames) == 2
    assert mmEngine.splitIndex > frames[0]
    assert mmEngine.splitIndex <= frames[1]

def test_getMMsamples(monkeypatch):
    coordPath = os.path.join("..", "coors.xyz")
    options["coordPath"] = coordPath
    os.chdir(os.path.join(os.path.dirname(__file__), "mmengine"))
    if os.path.isdir("testSetup"):
        rmtree("testSetup")
    os.mkdir("testSetup")
    os.chdir("testSetup")
    
    monkeypatch.setattr(mmengine.MMEngine, "getFrame", monkeyGetFrame)
    monkeypatch.setattr(mmengine.MMEngine, "getFrames", monkeyGetFrames)
    #monkeypatch.setattr(mmengine.MMEngine, "getIndices", monkeyGetIndices)
    monkeypatch.setattr(mmengine.MMEngine, "setup", monkeySetup)
    mmEngine = mmengine.MMEngine(options)
    mmEngine.getMMSamples()

    with open(os.path.join("train","frames.txt"),'r') as f:
        trainCrds = f.readlines()

    with open(os.path.join("valid_1","frames.txt"),'r') as f:
        validCrds = f.readlines()

    os.chdir("..")
    #rmtree("testSetup")
    assert trainCrds == ["0\n"]
    assert validCrds == ["1\n"]

# Also tests utils.writeRst()
def test_getFrame(monkeypatch):
    monkeypatch.setattr(mmengine.MMEngine, "getIndices", monkeyGetIndices)
    options["coordPath"] = os.path.join(
        "ff-optimizer", "tests", "mmengine", "coors.xyz"
    )
    os.chdir(os.path.dirname(__file__))
    options["start"] = None
    mmEngine = mmengine.MMEngine(options)

    mmEngine.getFrame(23, "23.rst7")
    testCoors = []
    with open("23.rst7", "r") as f:
        for line in f.readlines()[2:]:
            testCoors.append(line.split())
    os.remove("23.rst7")
    refCoors = []
    with open(os.path.join("mmengine", "23.rst7"), "r") as f:
        for line in f.readlines()[2:]:
            refCoors.append(line.split())
    checkUtils.checkArray(testCoors, refCoors)


def monkeyGetSamples(self):
    self.test[0] = -1
    return


def monkeyCpptraj(cmd):
    return


def monkeySample(self, rst, mdin):
    self.test[self.counter] = len(mdin)
    self.counter += 1
    return


# Tests restarting a completed directory
def test_restart1(monkeypatch):
    monkeypatch.setattr(mmengine.MMEngine, "getIndices", monkeyGetIndices)
    coordPath = os.path.join("mmengine", "coors.xyz")
    options["coordPath"] = coordPath
    options["leap"] = "setup.leap"
    os.chdir(os.path.join(os.path.dirname(__file__), "mmengine", "restart1"))
    mmEngine = mmengine.MMEngine(options)
    mmEngine.prmtop = "water.prmtop"

    monkeypatch.setattr(mmengine.MMEngine, "getMMSamples", monkeyGetSamples)
    monkeypatch.setattr(mmengine.MMEngine, "sample", monkeySample)
    monkeypatch.setattr(os, "system", monkeyCpptraj)

    mmEngine.test = zeros(mmEngine.options["nvalids"] + 1)
    mmEngine.counter = 0
    mmEngine.restart()

    for f in os.listdir():
        if os.path.isdir(f):
            assert os.path.isfile(os.path.join(f, "cpptraj.in"))
            os.remove(os.path.join(f, "cpptraj.in"))
    for i in range(mmEngine.options["nvalids"] + 1):
        assert mmEngine.test[i] == 0


# Tests restarting a directory with one IC complete
def test_restart2(monkeypatch):
    monkeypatch.setattr(mmengine.MMEngine, "getIndices", monkeyGetIndices)
    os.chdir(os.path.join(os.path.dirname(__file__), "mmengine", "restart2"))
    mmEngine = mmengine.MMEngine(options)
    mmEngine.prmtop = "water.prmtop"

    monkeypatch.setattr(mmengine.MMEngine, "getMMSamples", monkeyGetSamples)
    monkeypatch.setattr(mmengine.MMEngine, "sample", monkeySample)
    monkeypatch.setattr(os, "system", monkeyCpptraj)

    mmEngine.test = zeros(mmEngine.options["nvalids"] + 1)
    mmEngine.counter = 0
    mmEngine.restart()

    assert os.path.isfile(os.path.join("train_928", "cpptraj.in"))
    assert os.path.isdir("valid_828")
    os.remove(os.path.join("train_928", "cpptraj.in"))
    rmtree("valid_828")
    assert mmEngine.test[0] == 5
    assert mmEngine.test[1] == 0


# Tests restarting a directory with a partially complete IC
def test_restart3(monkeypatch):
    monkeypatch.setattr(mmengine.MMEngine, "getIndices", monkeyGetIndices)
    options["trainMdin"] = "train.in"
    options["validMdin"] = "md.in"
    os.chdir(os.path.join(os.path.dirname(__file__), "mmengine", "restart3"))
    mmEngine = mmengine.MMEngine(options)
    mmEngine.prmtop = "water.prmtop"

    monkeypatch.setattr(mmengine.MMEngine, "getMMSamples", monkeyGetSamples)
    monkeypatch.setattr(mmengine.MMEngine, "sample", monkeySample)
    monkeypatch.setattr(os, "system", monkeyCpptraj)

    mmEngine.test = zeros(mmEngine.options["nvalids"] + 1)
    mmEngine.counter = 0
    mmEngine.restart()

    assert os.path.isdir("valid_828")
    rmtree("valid_828")
    assert mmEngine.test[0] == 5
    assert mmEngine.test[1] == 8


# Tests restarting a directory with the wrong number of rsts
def test_restart4(monkeypatch):
    monkeypatch.setattr(mmengine.MMEngine, "getIndices", monkeyGetIndices)
    options["trainMdin"] = "train.in"
    options["validMdin"] = "md.in"
    os.chdir(os.path.join(os.path.dirname(__file__), "mmengine", "restart4"))
    mmEngine = mmengine.MMEngine(options)
    mmEngine.prmtop = "water.prmtop"

    monkeypatch.setattr(mmengine.MMEngine, "getMMSamples", monkeyGetSamples)
    monkeypatch.setattr(mmengine.MMEngine, "sample", monkeySample)
    monkeypatch.setattr(os, "system", monkeyCpptraj)

    mmEngine.test = zeros(mmEngine.options["nvalids"] + 1)
    mmEngine.counter = 0
    mmEngine.restart()

    assert not os.path.isdir("valid_828")
    os.mkdir("valid_828")
    os.system(f"cp ../restart1/828.rst7 .")

    assert mmEngine.test[0] == -1

def test_getFrames(monkeypatch):
    options["conformers"] = 3
    options["nvalids"] = 2
    coordPath = os.path.join("..", "coors.xyz")
    options["coordPath"] = coordPath
    os.chdir(os.path.join(os.path.dirname(__file__), "mmengine", "restart1"))
    #monkeypatch.setattr(mmengine.MMEngine, "getIndices", monkeyGetIndices)
    mmEngine = mmengine.MMEngine(options)
    frames = mmEngine.getFrames()
    assert len(frames) == 9
    for i in range(3):
        assert mmEngine.splitIndex > frames[i]
    for i in range(3,9):
        assert mmEngine.splitIndex <= frames[i]

def test_getMMsamples(monkeypatch):
    coordPath = os.path.join("..", "coors.xyz")
    options["coordPath"] = coordPath
    os.chdir(os.path.join(os.path.dirname(__file__), "mmengine"))
    if os.path.isdir("testSetup"):
        rmtree("testSetup")
    os.mkdir("testSetup")
    os.chdir("testSetup")
    
    monkeypatch.setattr(mmengine.MMEngine, "getFrame", monkeyGetFrame)
    monkeypatch.setattr(mmengine.MMEngine, "getFrames", monkeyGetFrames)
    #monkeypatch.setattr(mmengine.MMEngine, "getIndices", monkeyGetIndices)
    monkeypatch.setattr(mmengine.MMEngine, "setup", monkeySetup)
    mmEngine = mmengine.MMEngine(options)
    mmEngine.getMMSamples()

    with open(os.path.join("train","frames.txt"),'r') as f:
        trainCrds = f.readlines()
    validCrds = []
    with open(os.path.join("valid_1","frames.txt"),'r') as f:
        validCrds.append(f.readlines())
    with open(os.path.join("valid_2","frames.txt"),'r') as f:
        validCrds.append(f.readlines())

    os.chdir("..")
    rmtree("testSetup")
    assert trainCrds == ["0\n", "1\n", "2\n"]
    assert validCrds[0] == ["3\n", "4\n", "5\n"]
    assert validCrds[1] == ["6\n", "7\n", "8\n"]
