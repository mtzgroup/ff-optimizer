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
options['split'] = 1000
options['stride'] = 50
options['coordsDir'] = "mmengine"
options['coords'] = "coors.xyz"
options['heatCounter'] = 8

def test_AmberInit(monkeypatch):
    def monkeyFail():
        raise RuntimeError("Oops")

    def monkeyID():
        return [0]

    os.chdir(os.path.dirname(__file__))
    monkeypatch.setattr(GPUtil, "getAvailable", monkeyFail)
    mmEngine = mmengine.ExternalAmberEngine(options)
    assert mmEngine.amberExe == "pmemd"
    monkeypatch.setattr(GPUtil, "getAvailable", monkeyID)
    mmEngine = mmengine.ExternalAmberEngine(options)
    assert mmEngine.amberExe == "pmemd.cuda"

@pytest.mark.external
def test_runSander():
    os.chdir(os.path.join(os.path.dirname(__file__),"mmengine"))
    options['coordsDir'] = "."
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
def test_runSanderCUDA():
    os.chdir(os.path.join(os.path.dirname(__file__),"mmengine"))
    options['coordsDir'] = "."
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
    os.chdir(os.path.join(os.path.dirname(__file__),"mmengine","sample"))
    options['coordsDir'] = ".."
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
    mmEngine.sample("928.rst7","md.in")
    for i in range(1,11):
        with open(os.path.join("928",f"{str(i)}.pdb"),'r') as f:
            testLines = f.readlines()
        with open(os.path.join("ref",f"{str(i)}.pdb"),'r') as f:
            refLines = f.readlines()
        assert len(testLines) == len(refLines)
        for i in range(len(testLines)):
            assert refLines[i] == testLines[i]
    rmtree("928")

def test_restart():
    pass
