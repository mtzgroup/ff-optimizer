import os
from pathlib import Path

import numpy as np
import pytest
from chemcloud.models import Molecule

from ff_optimizer import utils

from . import checkUtils

home = Path(__file__).parent.absolute()

# check that grads and energies are being read in correctly
def test_readGradFromTCout():
    os.chdir(home)
    testEnergy, testGrads = utils.readGradFromTCout(
        os.path.join("qmengine", "test.out")
    )
    energy = np.loadtxt(os.path.join("qmengine", "energy.txt"))
    grads = np.loadtxt(os.path.join("qmengine", "grads.txt")).flatten()
    assert checkUtils.checkFloat(energy, testEnergy)
    assert checkUtils.checkArrays(grads, testGrads)


def test_readGradFromTCout_TCCloud():
    os.chdir(home)
    testEnergy, testGrads = utils.readGradFromTCout(
        os.path.join("qmengine", "result.txt")
    )
    energy = np.loadtxt(os.path.join("qmengine", "energy.txt"))
    grads = np.loadtxt(os.path.join("qmengine", "grads.txt")).flatten()
    assert checkUtils.checkFloat(energy, testEnergy)
    assert checkUtils.checkArrays(grads, testGrads)


def test_convertTCtoFB():
    os.chdir(home / "utils" / "dynamics")
    utils.convertTCtoFB("tc.out", "coors.xyz", 100)
    mdcrdCorrect = checkUtils.checkFileFloatsNoWhitespace("refAll1.mdcrd", "all.mdcrd")
    qdataCorrect = checkUtils.checkFileFloatsNoWhitespace("refQdata1.txt", "qdata.txt")
    os.remove("all.mdcrd")
    os.remove("qdata.txt")
    assert mdcrdCorrect
    assert qdataCorrect


def test_convertTCtoFB2():
    os.chdir(home / "utils" / "dynamics")
    utils.convertTCtoFB("tc.out", "coors.xyz", 100, 100, 700)
    mdcrdCorrect = checkUtils.checkFileFloatsNoWhitespace("refAll2.mdcrd", "all.mdcrd")
    checkUtils.checkFileFloatsNoWhitespace("refQdata2.txt", "qdata.txt")
    os.remove("all.mdcrd")
    os.remove("qdata.txt")
    assert mdcrdCorrect
    # assert qdataCorrect


def test_readEsp():
    os.chdir(home)
    testEspXYZ, testEsp = utils.readEsp(os.path.join("utils", "esp.xyz"))
    refEspXYZ = np.loadtxt(os.path.join("utils", "espXYZ.txt")).flatten()
    refEsp = np.loadtxt(os.path.join("utils", "esp.txt"))
    assert checkUtils.checkArrays(testEspXYZ, refEspXYZ, 0.00001)
    assert checkUtils.checkArrays(testEsp, refEsp, 0.00001)


def test_convertPDBtoMolecule():
    os.chdir(home)
    testMol = utils.convertPDBtoMolecule(os.path.join("utils", "test.pdb"))
    refMol = Molecule.from_file(os.path.join("utils", "test.xyz"))
    assert checkUtils.checkArrays(testMol.geometry, refMol.geometry)
    assert len(testMol.symbols) == len(refMol.symbols)
    for i in range(len(testMol.symbols)):
        assert testMol.symbols[i] == refMol.symbols[i]

def test_writeRst():
    os.chdir(home)
    coords = np.loadtxt(
        os.path.join("mmengine", "test.xyz"), skiprows=2, usecols=(1, 2, 3)
    )
    utils.writeRst(list(coords), coords.shape[0], "test.rst7")
    testCoors = []
    with open("test.rst7", "r") as f:
        for line in f.readlines()[2:]:
            for token in line.split():
                testCoors.append(token)
    os.remove("test.rst7")
    refCoors = []
    with open(os.path.join("mmengine", "test.rst7"), "r") as f:
        for line in f.readlines()[2:]:
            for token in line.split():
                refCoors.append(token)
    checkUtils.checkArrays(testCoors, refCoors)


# check that pdbs are read in correctly
def test_readPDB1():
    os.chdir(home)
    testCoords = utils.readPDB("qmengine/test.pdb")
    coords = np.loadtxt("qmengine/coords.txt").flatten()
    assert checkUtils.checkArrays(coords, testCoords)


def test_readPDB2():
    os.chdir(os.path.join(home, "utils"))
    coors = utils.readPDB("1.pdb")
    ref = np.loadtxt("1.txt").flatten()
    checkUtils.checkArrays(coors, ref)


def test_writePDB():
    os.chdir(os.path.join(home, "utils"))
    coords = utils.readPDB("3.pdb")
    utils.writePDB(coords, "3_test.pdb", "2.pdb")
    with open("3_test.pdb", "r") as f:
        testLines = f.readlines()
    with open("3.pdb", "r") as f:
        refLines = f.readlines()
    assert testLines == refLines


def test_readXYZ1():
    os.chdir(os.path.join(home, "utils"))
    coords = utils.readXYZ("test.xyz")
    assert checkUtils.checkFloat(coords[0], 0.081)
    assert checkUtils.checkFloat(coords[2], -0.195)
    assert checkUtils.checkFloat(coords[-1], -2.544)


def test_readXYZ2():
    os.chdir(os.path.join(home, "utils"))
    coords, symbols = utils.readXYZ("test.xyz", readSymbols=True)
    assert symbols[0] == "O"
    assert symbols[-1] == "H"


def test_writeXYZ():
    os.chdir(os.path.join(home, "utils"))
    refCoords, refSymbols = utils.readXYZ("test.xyz", readSymbols=True)
    utils.writeXYZ(refCoords, refSymbols, "temp.xyz")
    testCoords, testSymbols = utils.readXYZ("temp.xyz", readSymbols=True)
    os.remove("temp.xyz")
    assert checkUtils.checkArrays(refCoords, testCoords)
    assert checkUtils.checkLists(list(refSymbols), list(testSymbols))


def test_convertNCtoXYZs():
    os.chdir(os.path.join(home, "utils", "test_nc"))
    _, symbols = utils.readXYZ("ref_1.xyz", readSymbols=True)
    offset = 10
    numXYZs = utils.convertNCtoXYZs("test.nc", symbols, offset=offset)
    passTest = True
    for i in range(1, 11):
        refCoords = utils.readXYZ(f"ref_{i}.xyz")
        testCoords = utils.readXYZ(f"{i+offset}.xyz")
        if not checkUtils.checkArrays(refCoords, testCoords):
            passTest = False
        os.remove(f"{i+offset}.xyz")
    assert passTest


def test_convertNCtoXYZs2():
    os.chdir(os.path.join(home, "utils", "test_nc"))
    numXYZs = utils.convertNCtoXYZs("test_vel.nc", [], offset=0)
    assert numXYZs == 0


def test_getSymbolsFromPrmtop():
    os.chdir(os.path.join(home, "utils"))
    _, refSymbols = utils.readXYZ("prmtop_test.xyz", readSymbols=True)
    testSymbols = utils.getSymbolsFromPrmtop("test.prmtop")
    assert checkUtils.checkLists(list(refSymbols), testSymbols)
