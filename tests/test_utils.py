from ff_optimizer import utils
import os
import numpy as np
from . import checkUtils
from tccloud.models import Molecule

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
    coords = np.loadtxt(os.path.join("mmengine","23.xyz"),skiprows=2, usecols=(1,2,3))
    utils.writeRst(list(coords),coords.shape[0],"23.rst7")
    testCoors = []
    with open("23.rst7",'r') as f:
        for line in f.readlines()[2:]:
            testCoors.append(line.split())
    os.remove("23.rst7")
    refCoors = []
    with open(os.path.join("mmengine","23.rst7"),'r') as f:
        for line in f.readlines()[2:]:
            refCoors.append(line.split())
    checkUtils.checkArray(testCoors, refCoors)
    

