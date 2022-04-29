import pytest
from ff_optimizer import mmengine
import os
from . import checkUtils
from numpy import loadtxt
import GPUtil
from shutil import copyfile, rmtree

options = {}
options['start'] = 33
options['end'] = 2000
options['split'] = 777
options['stride'] = 50
options['coordPath'] = os.path.join("ff-optimizer","tests","mmengine","coors.xyz")
options['heatCounter'] = 8

def monkeyGetIndices(self):
    return 1,2,3

def test_AmberInit(monkeypatch):
    def monkeyFail(maxLoad=0):
        raise RuntimeError("Oops")

    def monkeyID(maxLoad=0):
        return [0]

    os.chdir(os.path.dirname(__file__))
    monkeypatch.setattr(mmengine.MMEngine,"getIndices",monkeyGetIndices)
    monkeypatch.setattr(GPUtil, "getAvailable", monkeyFail)
    mmEngine = mmengine.ExternalAmberEngine(options)
    assert mmEngine.amberExe == "pmemd"
    monkeypatch.setattr(GPUtil, "getAvailable", monkeyID)
    mmEngine = mmengine.ExternalAmberEngine(options)
    assert mmEngine.amberExe == "pmemd.cuda"

@pytest.mark.external
def test_runSander(monkeypatch):
    monkeypatch.setattr(mmengine.MMEngine,"getIndices",monkeyGetIndices)
    os.chdir(os.path.join(os.path.dirname(__file__),"mmengine"))
    options['coordPath'] = os.path.join("tests","mmengine","coors.xyz")
    mmEngine = mmengine.ExternalAmberEngine(options)
    mmEngine.amberExe = "pmemd"
    mmEngine.runSander("water.prmtop","md.in","md.out","1.rst7","md.mdcrd","md.rst7")
    testRst = loadtxt("md.rst7",skiprows=2)
    os.remove("md.rst7")
    os.remove("mdinfo")
    os.remove("md.out")
    refRst = loadtxt(os.path.join("ref","md.rst7"),skiprows=2)
    checkUtils.checkArray(testRst, refRst)

@pytest.mark.gpu
@pytest.mark.external
def test_runSanderCUDA(monkeypatch):
    monkeypatch.setattr(mmengine.MMEngine,"getIndices",monkeyGetIndices)
    os.chdir(os.path.join(os.path.dirname(__file__),"mmengine"))
    options['coordPath'] = os.path.join("tests","mmengine","coors.xyz")
    mmEngine = mmengine.ExternalAmberEngine(options)
    mmEngine.amberExe = "pmemd.cuda"
    mmEngine.runSander("water.prmtop","md.in","md.out","1.rst7","md.mdcrd","md.rst7")
    testRst = loadtxt("md.rst7",skiprows=2)
    os.remove("md.rst7")
    os.remove("mdinfo")
    os.remove("md.out")
    refRst = loadtxt(os.path.join("ref","md.rst7"),skiprows=2)
    checkUtils.checkArray(testRst, refRst)

@pytest.mark.external
def test_sample(monkeypatch):
    monkeypatch.setattr(mmengine.MMEngine,"getIndices",monkeyGetIndices)
    os.chdir(os.path.join(os.path.dirname(__file__),"mmengine","sample"))
    if os.path.isdir("928"):
        rmtree("928")
    os.mkdir("928")
    os.chdir("928")
    options['coordPath'] = "coors.xyz"
    mmEngine = mmengine.ExternalAmberEngine(options)
    mmEngine.prmtop = "water.prmtop"

    def monkeySander(self, prmtop, mdin, mdout, mdcrd, mdtraj, restart, mdvels=None):
        if not os.path.isfile(prmtop):
            raise RuntimeError("No prmtop!")
        if not os.path.isfile(mdin):
            raise RuntimeError("No mdin!")
        if not os.path.isfile(mdcrd):
            raise RuntimeError("No input crd!")
        copyfile(os.path.join("..","ref",restart),restart)
        try:
            copyfile(os.path.join("..","ref",mdtraj),mdtraj)
        except:
            1+1
        return

    monkeypatch.setattr(mmengine.ExternalAmberEngine,"runSander",monkeySander)
    mmEngine.sample("928","md.in")
    os.chdir("..")
    for i in range(1,11):
        with open(os.path.join("928",f"{str(i)}.pdb"),'r') as f:
            testLines = f.readlines()
        with open(os.path.join("ref",f"{str(i)}.pdb"),'r') as f:
            refLines = f.readlines()
        assert len(testLines) == len(refLines)
        for i in range(len(testLines)):
            assert refLines[i] == testLines[i]
    rmtree("928")
