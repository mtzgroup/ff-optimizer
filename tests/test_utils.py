from ff_optimizer import utils
import os
import numpy as np
from . import checkUtils

def test_convertPDBtoXYZ():
    pass

# check that grads and energies are being read in correctly
def test_readGradFromTCout():
    os.chdir(os.path.dirname(__file__))
    testEnergy, testGrads = utils.readGradFromTCout(os.path.join("qmengine","test.out"))
    energy = np.loadtxt(os.path.join("qmengine","energy.txt"))
    grads = np.loadtxt(os.path.join("qmengine","grads.txt")).flatten()
    assert checkUtils.checkFloat(energy, testEnergy)
    assert checkUtils.checkArray(grads, testGrads)

def test_readGradFromTCout_TCCloud():
    os.chdir(os.path.dirname(__file__))
    testEnergy, testGrads = utils.readGradFromTCout(os.path.join("qmengine","result.txt"))
    energy = np.loadtxt(os.path.join("qmengine","energy.txt"))
    grads = np.loadtxt(os.path.join("qmengine","grads.txt")).flatten()
    assert checkUtils.checkFloat(energy, testEnergy)
    assert checkUtils.checkArray(grads, testGrads)

def test_readOpt():
    pass
