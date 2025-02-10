import os
import pytest

from numpy import loadtxt

from ff_optimizer import mmengine, utils

from . import checkUtils
from .test_inputs import getDefaults


def monkeyGetIndices(self):
    return 1, 2, 3


def clean():
    for f in os.listdir():
        if f.endswith(".xyz") or f.endswith(".nc") or f.endswith(".out"):
            os.remove(f)


def test_sample(monkeypatch):
    options = getDefaults()
    monkeypatch.setattr(mmengine.MMEngine, "getIndices", monkeyGetIndices)
    os.chdir(os.path.join(os.path.dirname(__file__), "mmengine"))
    mmEngine = mmengine.ExternalOpenMMEngine(options)
    mmEngine.prmtop = "bicarb_cluster.prmtop"
    os.chdir("sample_openmm")
    clean()
    mmEngine.symbols = utils.getSymbolsFromPrmtop(mmEngine.prmtop)
    mmEngine.sample([12345], "md.py")
    passTest = True
    for i in range(1, 9):
        testXYZ = loadtxt(f"{i}.xyz", usecols=(1, 2, 3), skiprows=2)
        refXYZ = loadtxt(os.path.join("ref", f"{i}.xyz"), usecols=(1, 2, 3), skiprows=2)
        if not checkUtils.checkArrays(testXYZ, refXYZ, 1e-2):
            passTest = False
    clean()
    assert passTest
