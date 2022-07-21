import os

import numpy as np
from chemcloud.models import Molecule

from ff_optimizer import utils

from . import checkUtils


# check that grads and energies are being read in correctly
def test_readGradFromTCout():
    os.chdir(os.path.dirname(__file__))
    testEnergy, testGrads = utils.readGradFromTCout(
        os.path.join("qmengine", "test.out")
    )
    energy = np.loadtxt(os.path.join("qmengine", "energy.txt"))
    grads = np.loadtxt(os.path.join("qmengine", "grads.txt")).flatten()
    assert checkUtils.checkFloat(energy, testEnergy)
    assert checkUtils.checkArray(grads, testGrads)


def test_readGradFromTCout_TCCloud():
    os.chdir(os.path.dirname(__file__))
    testEnergy, testGrads = utils.readGradFromTCout(
        os.path.join("qmengine", "result.txt")
    )
    energy = np.loadtxt(os.path.join("qmengine", "energy.txt"))
    grads = np.loadtxt(os.path.join("qmengine", "grads.txt")).flatten()
    assert checkUtils.checkFloat(energy, testEnergy)
    assert checkUtils.checkArray(grads, testGrads)


def test_readOpt():
    pass


def test_readEsp():
    os.chdir(os.path.dirname(__file__))
    testEspXYZ, testEsp = utils.readEsp(os.path.join("utils", "esp.xyz"))
    refEspXYZ = np.loadtxt(os.path.join("utils", "espXYZ.txt")).flatten()
    refEsp = np.loadtxt(os.path.join("utils", "esp.txt"))
    assert checkUtils.checkArray(testEspXYZ, refEspXYZ, 0.00001)
    assert checkUtils.checkArray(testEsp, refEsp, 0.00001)


def test_convertPDBtoMolecule():
    os.chdir(os.path.dirname(__file__))
    testMol = utils.convertPDBtoMolecule(os.path.join("utils", "test.pdb"))
    refMol = Molecule.from_file(os.path.join("utils", "test.xyz"))
    assert checkUtils.checkArray(testMol.geometry, refMol.geometry)
    assert len(testMol.symbols) == len(refMol.symbols)
    for i in range(len(testMol.symbols)):
        assert testMol.symbols[i] == refMol.symbols[i]


def test_writeRst():
    os.chdir(os.path.dirname(__file__))
    coords = np.loadtxt(
        os.path.join("mmengine", "23.xyz"), skiprows=2, usecols=(1, 2, 3)
    )
    utils.writeRst(list(coords), coords.shape[0], "23.rst7")
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


# check that pdbs are read in correctly
def test_readPDB1():
    os.chdir(os.path.dirname(__file__))
    testCoords = utils.readPDB("qmengine/test.pdb")
    coords = np.loadtxt("qmengine/coords.txt").flatten()
    assert checkUtils.checkArray(coords, testCoords)


def test_readPDB2():
    os.chdir(os.path.join(os.path.dirname(__file__), "utils"))
    coors = utils.readPDB("1.pdb")
    ref = np.loadtxt("1.txt").flatten()
    checkUtils.checkArray(coors, ref)

def test_writePDB():    
    os.chdir(os.path.join(os.path.dirname(__file__), "utils"))
    coords = utils.readPDB("3.pdb")
    utils.writePDB(coords, "3_test.pdb", "2.pdb")
    with open("3_test.pdb", 'r') as f:
        testLines = f.readlines()
    with open("3.pdb", 'r') as f:
        refLines = f.readlines()
    assert testLines == refLines
