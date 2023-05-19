import os
from pathlib import Path
from shutil import copyfile, rmtree

import GPUtil
import pytest
from numpy import loadtxt

from ff_optimizer import mmengine, utils

from . import checkUtils
from .test_inputs import getDefaults


def monkeyGetIndices(self):
    return 1, 2, 3


def monkeySander(self, prmtop, mdin, mdout, mdcrd, mdtraj, restart, mdvels=None):
    if not os.path.isfile(prmtop):
        raise RuntimeError("No prmtop!")
    if not os.path.isfile(mdin):
        raise RuntimeError("No mdin!")
    if not os.path.isfile(mdcrd):
        raise RuntimeError("No input crd!")
    copyfile(os.path.join("..", "ref", restart), restart)
    try:
        copyfile(os.path.join("..", "ref", mdtraj), mdtraj)
    except:
        pass
    return


def test_AmberInit(monkeypatch):
    options = getDefaults()
    options.heatCounter = 1

    def monkeyFail(maxLoad=0):
        raise RuntimeError("Oops")

    def monkeyID(maxLoad=0):
        return [0]

    os.chdir(os.path.dirname(__file__))
    monkeypatch.setattr(mmengine.MMEngine, "getIndices", monkeyGetIndices)
    monkeypatch.setattr(GPUtil, "getAvailable", monkeyFail)
    mmEngine = mmengine.ExternalAmberEngine(options)
    assert mmEngine.amberExe == "pmemd"
    monkeypatch.setattr(GPUtil, "getAvailable", monkeyID)
    mmEngine = mmengine.ExternalAmberEngine(options)
    assert mmEngine.amberExe == "pmemd.cuda"


@pytest.mark.amber
def test_runSander(monkeypatch):
    inp = getDefaults()
    inp.heatCounter = 1
    inp.dynamicsdir = Path("tests") / "mmengine"
    monkeypatch.setattr(mmengine.MMEngine, "getIndices", monkeyGetIndices)
    os.chdir(os.path.join(os.path.dirname(__file__), "mmengine"))
    mmEngine = mmengine.ExternalAmberEngine(inp)
    mmEngine.amberExe = "pmemd"
    mmEngine.runSander(
        "water.prmtop", "md.in", "md.out", "1.rst7", "md.mdcrd", "md.rst7"
    )
    testRst = loadtxt("md.rst7", skiprows=2)
    os.remove("md.rst7")
    os.remove("mdinfo")
    os.remove("md.out")
    refRst = loadtxt(os.path.join("ref", "md.rst7"), skiprows=2)
    checkUtils.checkArrays(testRst, refRst)


@pytest.mark.gpu
@pytest.mark.amber
def test_runSanderCUDA(monkeypatch):
    inp = getDefaults()
    inp.heatCounter = 1
    inp.dynamicsdir = Path("tests") / "mmengine"
    monkeypatch.setattr(mmengine.MMEngine, "getIndices", monkeyGetIndices)
    os.chdir(os.path.join(os.path.dirname(__file__), "mmengine"))
    mmEngine = mmengine.ExternalAmberEngine(inp)
    mmEngine.amberExe = "pmemd.cuda"
    mmEngine.runSander(
        "water.prmtop", "md.in", "md.out", "1.rst7", "md.mdcrd", "md.rst7"
    )
    testRst = loadtxt("md.rst7", skiprows=2)
    os.remove("md.rst7")
    os.remove("mdinfo")
    os.remove("md.out")
    refRst = loadtxt(os.path.join("ref", "md.rst7"), skiprows=2)
    checkUtils.checkArrays(testRst, refRst)


@pytest.mark.amber
def test_sample(monkeypatch):
    inp = getDefaults()
    inp.heatCounter = 1
    monkeypatch.setattr(mmengine.MMEngine, "getIndices", monkeyGetIndices)
    monkeypatch.setattr(mmengine.ExternalAmberEngine, "runSander", monkeySander)
    os.chdir(os.path.join(os.path.dirname(__file__), "mmengine", "sample"))
    inp.dynamicsdir = Path("..").absolute()

    if os.path.isdir("train"):
        rmtree("train")
    os.mkdir("train")

    copyfile("928.rst7", os.path.join("train", "928.rst7"))
    mmEngine = mmengine.ExternalAmberEngine(inp)
    mmEngine.prmtop = "water.prmtop"
    mmEngine.symbols = utils.getSymbolsFromPrmtop("water.prmtop")

    os.chdir("train")
    mmEngine.sample([928], "md.in")
    os.chdir("..")
    for i in range(1, 11):
        testXyz = loadtxt(
            os.path.join("train", f"{str(i)}.xyz"), skiprows=2, usecols=(1, 2, 3)
        )
        refXyz = loadtxt(
            os.path.join("ref", f"{str(i)}.xyz"), skiprows=2, usecols=(1, 2, 3)
        )
        assert checkUtils.checkArrays(testXyz, refXyz)
    # rmtree("train")


@pytest.mark.amber
def test_sample_conformers(monkeypatch):
    inp = getDefaults()
    inp.heatCounter = 1
    monkeypatch.setattr(mmengine.MMEngine, "getIndices", monkeyGetIndices)
    monkeypatch.setattr(mmengine.ExternalAmberEngine, "runSander", monkeySander)
    os.chdir(os.path.join(os.path.dirname(__file__), "mmengine", "sample_conformers"))
    inp.dynamicsdir = Path("..").absolute()
    mmEngine = mmengine.ExternalAmberEngine(inp)
    mmEngine.prmtop = "water.prmtop"
    mmEngine.symbols = utils.getSymbolsFromPrmtop("water.prmtop")
    if os.path.isdir("valid_1"):
        rmtree("valid_1")
    os.mkdir("valid_1")
    copyfile("928.rst7", os.path.join("valid_1", "928.rst7"))
    copyfile("135.rst7", os.path.join("valid_1", "135.rst7"))
    copyfile("253.rst7", os.path.join("valid_1", "253.rst7"))

    os.chdir("valid_1")
    mmEngine.sample([928, 135, 253], "md.in")
    os.chdir("..")
    for i in range(1, 31):
        testXyz = loadtxt(
            os.path.join("valid_1", f"{str(i)}.xyz"), skiprows=2, usecols=(1, 2, 3)
        )
        refXyz = loadtxt(
            os.path.join("ref", f"{str(i)}.xyz"), skiprows=2, usecols=(1, 2, 3)
        )
        assert checkUtils.checkArrays(testXyz, refXyz)
    # rmtree("valid_1")
