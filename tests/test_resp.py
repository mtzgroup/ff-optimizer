import pytest
import os
from ff_optimizer import resp_prior
from numpy import loadtxt
from . import checkUtils
from shutil import copyfile

def test_findRepeatIndex():
    os.chdir(os.path.join(os.path.dirname(__file__),"resp"))
    options = {}
    options['sampleDir'] = "sample"
    options['mol2'] = "acn.mol2"
    options['mode'] = 1
    respPriors = resp_prior.RespPriors(options)
    j = respPriors.findRepeatIndex(3)
    assert j == 0
    j = respPriors.findRepeatIndex(4)
    assert j == 5

def test_readCharges():
    os.chdir(os.path.join(os.path.dirname(__file__),"resp"))
    options = {}
    options['sampleDir'] = "sample"
    options['mol2'] = "dasa.mol2"
    options['mode'] = 1
    respPriors = resp_prior.RespPriors(options)
    with open("resp.out",'r') as f:
        lines = f.readlines()
    esp, resp = respPriors.readCharges(lines)
    assert esp[0] == 0.837268
    assert esp[-1] == 0.530720
    assert resp[0] == 0.651138
    assert resp[-1] == 0.493376

def test_getCharges(monkeypatch):
    os.chdir(os.path.join(os.path.dirname(__file__),"resp"))
    options = {}
    options['sampleDir'] = "sample"
    options['mol2'] = "dasa.mol2"
    options['mode'] = 1
    respPriors = resp_prior.RespPriors(options)
    respPriors.getCharges(1)
    assert respPriors.allResp[0][0] == 0.527414
    assert respPriors.allEsp[0][-1] == 0.444999
    assert respPriors.allResp[1][-1] == 0.415671
    assert respPriors.allEsp[1][0] == 0.754497

# incidentally also tests findRepeatIndex
def test_setMol2Charges():
    os.chdir(os.path.join(os.path.dirname(__file__),"resp"))
    options = {}
    options['sampleDir'] = "sample"
    options['mol2'] = "dasa.mol2"
    options['mode'] = 1
    respPriors = resp_prior.RespPriors(options)
    charges = loadtxt("charges.txt")
    respPriors.setMol2Charges(charges, "test.mol2")
    passed = checkUtils.checkFileFloatsNoWhitespace("test.mol2","ref.mol2")
    copyfile("dasa.mol2","test.mol2")
    assert passed

def test_setPriors1():
    os.chdir(os.path.join(os.path.dirname(__file__),"resp"))
    options = {}
    options['sampleDir'] = "sample"
    options['mol2'] = "dasa.mol2"
    options['mode'] = 1
    respPriors = resp_prior.RespPriors(options)
    priors = loadtxt("priors.txt")
    copyfile("inp1.in", "test1.in")
    respPriors.setPriors(priors, "test1.in")
    passed = checkUtils.checkFileFloatsNoWhitespace("ref1.in", "test1.in")
    os.remove("test1.in")
    assert passed

def test_setPriors2():
    os.chdir(os.path.join(os.path.dirname(__file__),"resp"))
    options = {}
    options['sampleDir'] = "sample"
    options['mol2'] = "dasa.mol2"
    options['mode'] = 1
    respPriors = resp_prior.RespPriors(options)
    priors = loadtxt("priors.txt")
    copyfile("inp2.in", "test2.in")
    respPriors.setPriors(priors, "test2.in")
    passed = checkUtils.checkFileFloatsNoWhitespace("ref2.in", "test2.in")
    #os.remove("test2.in")
    assert passed

def test_getRepeats1():
    os.chdir(os.path.join(os.path.dirname(__file__),"resp"))
    options = {}
    options['sampleDir'] = "sample"
    options['mol2'] = "dasa.mol2"
    options['mode'] = 1
    respPriors = resp_prior.RespPriors(options)
    print(respPriors.repeats)
    refRepeats = [[0], [1], [2], [3], [4], [5, 6, 7], [], [], [8], [9, 10, 11], [], [], [12], [13], [14], [15], [16], [17], [18], [19], [20], [21], [22, 23], [], [24], [25, 26, 27], [], [], [28], [29, 30], [], [31], [32, 33, 34], [], [], [35], [36], [37], [38], [39], [40], [41]]
    assert checkUtils.checkLists(refRepeats, respPriors.repeats)

def test_getRepeats2():
    os.chdir(os.path.join(os.path.dirname(__file__),"resp"))
    options = {}
    options['sampleDir'] = "sample"
    options['mol2'] = "acn.mol2"
    options['mode'] = 1
    respPriors = resp_prior.RespPriors(options)
    refRepeats = [[0, 3], [1], [], [], [], [4, 5]]
    assert checkUtils.checkLists(refRepeats, respPriors.repeats)

def test_computeChargeDistribution():
    os.chdir(os.path.join(os.path.dirname(__file__),"resp"))
    options = {}
    options['sampleDir'] = "sample"
    options['mol2'] = "dasa.mol2"
    options['mode'] = 1
    respPriors = resp_prior.RespPriors(options)
    respPriors.getCharges(1)
    respPriors.computeChargeDistributions()
    assert checkUtils.checkFloat(respPriors.espMeans[0], 0.7691842)
    assert checkUtils.checkFloat(respPriors.espStdevs[0], 0.071253054)
    assert checkUtils.checkFloat(respPriors.espMeans[5], 0.13353156)
    assert checkUtils.checkFloat(respPriors.espStdevs[5], 0.01762259)

def test_computePriors():
    os.chdir(os.path.join(os.path.dirname(__file__),"resp"))
    options = {}
    options['sampleDir'] = "sample"
    options['mol2'] = "dasa.mol2"
    options['mode'] = 2
    respPriors = resp_prior.RespPriors(options)
    respPriors.getCharges(1)
    respPriors.computeChargeDistributions()
    priors = respPriors.computePriors()
    assert checkUtils.checkFloat(priors[5],0.05675933895)
    assert priors[6] == -1
    assert checkUtils.checkFloat(priors[0],0.1956745629)

# All individual functions are tested, so we just need to make sure this runs
def test_updateRespPriors():
    os.chdir(os.path.join(os.path.dirname(__file__),"resp"))
    options = {}
    options['sampleDir'] = "sample"
    options['mol2'] = "dasa.mol2"
    options['mode'] = 1
    respPriors = resp_prior.RespPriors(options)
    copyfile("dasa.mol2","update.mol2")
    copyfile("inp1.in", "update.in")
    respPriors.updateRespPriors(1, "update.mol2", "update.in")
    os.remove("update.mol2")
    os.remove("update.in")
