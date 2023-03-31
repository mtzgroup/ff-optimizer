import os

from numpy import loadtxt

from ff_optimizer import mmengine, utils

from . import checkUtils

options = {}
options["start"] = 33
options["end"] = 2000
options["split"] = 777
options["stride"] = 50
options["coordPath"] = os.path.join("ff-optimizer", "tests", "mmengine", "coors.xyz")
options["heatCounter"] = 8


def monkeyGetIndices(self):
    return 1, 2, 3


def clean():
    for f in os.listdir():
        if f.endswith(".xyz") or f.endswith(".nc") or f.endswith(".out"):
            os.remove(f)


def test_sample(monkeypatch):
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
        if not checkUtils.checkArrays(testXYZ, refXYZ, 1e-3):
            passTest = False
    clean()
    assert passTest
