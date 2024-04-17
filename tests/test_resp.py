import os
import pytest
from pathlib import Path
from shutil import copyfile

from numpy import loadtxt

from ff_optimizer import resp_prior

from . import checkUtils
from .test_inputs import getDefaults


def test_findRepeatIndex():
    os.chdir(os.path.join(os.path.dirname(__file__), "resp"))
    options = getDefaults()
    options.resppriors = 1
    options.sampledir = Path("sample")
    mol2 = "acn.mol2"
    prmtop = "dasa.prmtop"
    respPriors = resp_prior.RespPriors(options, mol2, prmtop)
    j = respPriors.findRepeatIndex(3)
    assert j == 0
    j = respPriors.findRepeatIndex(4)
    assert j == 5


def test_readCharges():
    os.chdir(os.path.join(os.path.dirname(__file__), "resp"))
    options = getDefaults()
    options.resppriors = 1
    options.sampledir = Path("sample")
    mol2 = "dasa.mol2"
    prmtop = "dasa.prmtop"
    respPriors = resp_prior.RespPriors(options, mol2, prmtop)
    with open("resp.out", "r") as f:
        lines = f.readlines()
    esp, resp = respPriors.readCharges(lines)
    assert esp[0] == 0.837268
    assert esp[-1] == 0.530720
    assert resp[0] == 0.651138
    assert resp[-1] == 0.493376


def test_getCharges(monkeypatch):
    os.chdir(os.path.join(os.path.dirname(__file__), "resp"))
    options = getDefaults()
    options.resppriors = 1
    options.sampledir = Path("sample")
    mol2 = "dasa.mol2"
    prmtop = "dasa.prmtop"
    respPriors = resp_prior.RespPriors(options, mol2, prmtop)
    respPriors.getCharges(1)
    assert respPriors.allResp[1][0] == 0.527414
    assert respPriors.allEsp[1][-1] == 0.444999
    assert respPriors.allResp[2][-1] == 0.415671
    assert respPriors.allEsp[2][0] == 0.754497


# incidentally also tests findRepeatIndex
def test_setMol2Charges():
    os.chdir(os.path.join(os.path.dirname(__file__), "resp"))
    options = getDefaults()
    options.resppriors = 1
    options.sampledir = Path("sample")
    mol2 = "dasa.mol2"
    prmtop = "dasa.prmtop"
    respPriors = resp_prior.RespPriors(options, mol2, prmtop)
    charges = loadtxt("charges.txt")
    respPriors.setMol2Charges(charges, "test.mol2")
    passed = checkUtils.checkFileFloatsNoWhitespace("test.mol2", "ref.mol2")
    copyfile("dasa.mol2", "test.mol2")
    assert passed


def test_setPriors1():
    os.chdir(os.path.join(os.path.dirname(__file__), "resp"))
    options = getDefaults()
    options.resppriors = 1
    options.sampledir = Path("sample")
    mol2 = "dasa.mol2"
    prmtop = "dasa.prmtop"
    respPriors = resp_prior.RespPriors(options, mol2, prmtop)
    priors = loadtxt("priors.txt")
    copyfile("inp1.in", "test1.in")
    respPriors.setPriors(priors, "test1.in")
    passed = checkUtils.checkFileFloatsNoWhitespace("ref1.in", "test1.in")
    os.remove("test1.in")
    assert passed


def test_setPriors2():
    os.chdir(os.path.join(os.path.dirname(__file__), "resp"))
    options = getDefaults()
    options.resppriors = 1
    options.sampledir = Path("sample")
    mol2 = "dasa.mol2"
    prmtop = "dasa.prmtop"
    respPriors = resp_prior.RespPriors(options, mol2, prmtop)
    priors = loadtxt("priors.txt")
    copyfile("inp2.in", "test2.in")
    respPriors.setPriors(priors, "test2.in")
    passed = checkUtils.checkFileFloatsNoWhitespace("ref2.in", "test2.in")
    # os.remove("test2.in")
    assert passed


def test_getRepeats1():
    os.chdir(os.path.join(os.path.dirname(__file__), "resp"))
    options = getDefaults()
    options.resppriors = 1
    options.sampledir = Path("sample")
    mol2 = "dasa.mol2"
    prmtop = "dasa.prmtop"
    respPriors = resp_prior.RespPriors(options, mol2, prmtop)
    refRepeats = [
        [0],
        [1],
        [2],
        [3],
        [4],
        [5, 6, 7],
        [],
        [],
        [8],
        [9, 10, 11],
        [],
        [],
        [12],
        [13],
        [14],
        [15],
        [16],
        [17],
        [18],
        [19],
        [20],
        [21],
        [22, 23],
        [],
        [24],
        [25, 26, 27],
        [],
        [],
        [28],
        [29, 30],
        [],
        [31],
        [32, 33, 34],
        [],
        [],
        [35],
        [36],
        [37],
        [38],
        [39],
        [40],
        [41],
    ]
    assert checkUtils.checkLists(refRepeats, respPriors.repeats)


def test_getRepeats2():
    os.chdir(os.path.join(os.path.dirname(__file__), "resp"))
    options = getDefaults()
    options.resppriors = 1
    options.sampledir = Path("sample")
    mol2 = "acn.mol2"
    prmtop = "dasa.prmtop"
    respPriors = resp_prior.RespPriors(options, mol2, prmtop)
    refRepeats = [[0, 3], [1], [], [], [], [4, 5]]
    assert checkUtils.checkLists(refRepeats, respPriors.repeats)


def test_computeChargeDistribution():
    os.chdir(os.path.join(os.path.dirname(__file__), "resp"))
    options = getDefaults()
    options.resppriors = 1
    options.sampledir = Path("sample")
    mol2 = "dasa.mol2"
    prmtop = "dasa.prmtop"
    respPriors = resp_prior.RespPriors(options, mol2, prmtop)
    respPriors.getCharges(1)
    respPriors.computeChargeDistributions()
    assert checkUtils.checkFloats(respPriors.espMeans[0], 0.7691842)
    assert checkUtils.checkFloats(respPriors.espStdevs[0], 0.071253054)
    assert checkUtils.checkFloats(respPriors.espMeans[5], 0.13353156)
    assert checkUtils.checkFloats(respPriors.espStdevs[5], 0.01762259)

def test_computePriors():
    os.chdir(os.path.join(os.path.dirname(__file__), "resp"))
    options = getDefaults()
    options.resppriors = 2
    options.sampledir = Path("sample")
    mol2 = "dasa.mol2"
    prmtop = "dasa.prmtop"
    respPriors = resp_prior.RespPriors(options, mol2, prmtop)
    respPriors.getCharges(1)
    respPriors.computeChargeDistributions()
    priors = respPriors.computePriors()
    print(priors)
    assert checkUtils.checkFloats(priors[5], 0.05675933895)
    assert priors[6] == -1
    assert checkUtils.checkFloats(priors[0], 0.1956745629)


# All individual functions are tested, so we just need to make sure this runs
# ^ Pretty sus
def test_updateRespPriors():
    os.chdir(os.path.join(os.path.dirname(__file__), "resp"))
    options = getDefaults()
    options.resppriors = 1
    options.sampledir = Path("sample")
    mol2 = "dasa.mol2"
    prmtop = "dasa.prmtop"
    respPriors = resp_prior.RespPriors(options, mol2, prmtop)
    copyfile("dasa.mol2", "update.mol2")
    copyfile("inp1.in", "opt_1.in")
    respPriors.updateRespPriors(1, "update.mol2")
    os.remove("update.mol2")
    os.remove("opt_1.in")
