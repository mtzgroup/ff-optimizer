import os
import random
from shutil import copyfile, rmtree

import numpy as np
import pytest

from ff_optimizer import active_learning, model, optengine, qmengine, utils

from . import checkUtils


class FakeModel:
    def __init__(self, args):
        self.sampledir = args.sampledir
        self.mmEngine = FakeMMEngine()


class FakeMMEngine:
    def __init__(self, options={}):
        self.prmtop = "amber.prmtop"


class FakeArgs:
    def __init__(self):
        self.activeLearning = 1
        self.sampledir = "2_sampling"
        self.optdir = "1_opt"
        self.dynamicsdir = "0_dynamics"
        self.restart = False
        self.valid0 = "valid_0.in"
        self.opt0 = "opt_0.in"
        self.resp = 0
        self.respPriors = 0
        self.maxcycles = 100
        self.qmengine = "chemcloud"
        self.tctemplate = "tc.in"
        self.tctemplate_long = "tc.in"
        self.nvalids = 2


def monkeyQMRestart(directory):
    pass

def monkeyInit(self, args):
    pass

class MonkeyQMEngine():
    def __init__(self, args):
        self.pdbs = []
        self.dirs= []

    def getQMRefData(self, pdbs, directory):
        self.pdbs.append(pdbs)
        self.dirs.append(os.getcwd().split('/')[-1])

    def restart(self, directory):
        self.dirs.append(os.getcwd().split('/')[-1])

def monkeyInitQM(self, args):
    return MonkeyQMEngine(args)

def monkeyInitModel(self, args):
    self.home = os.getcwd()
    self.sampledir = args.sampledir
    self.qmEngine = self.initializeQMEngine(args)
    self.restartCycle = 1

def test_doQMCalculations(monkeypatch):
    args = FakeArgs()
    os.chdir(os.path.join(os.path.dirname(__file__), "model", "test1"))
    monkeypatch.setattr(model.Model,"initializeQMEngine",monkeyInitQM)
    monkeypatch.setattr(model.Model,"__init__",monkeyInitModel)

    m = model.Model(args)
    m.doQMCalculations(3)
    testPdbs = {}
    testPdbs['train'] = ['1.pdb', '2.pdb']
    testPdbs['valid'] = ['3.pdb', '4.pdb']
    testPdbs['valid_2'] = ['5.pdb', '6.pdb']
    assert len(m.qmEngine.dirs) == 3
    for i in range(len(m.qmEngine.dirs)):
        for pdb in m.qmEngine.pdbs[i]:
            assert pdb in testPdbs[m.qmEngine.dirs[i]]

def test_doQMCalculationsRestart(monkeypatch):
    args = FakeArgs()
    os.chdir(os.path.join(os.path.dirname(__file__), "model", "test1"))
    monkeypatch.setattr(model.Model,"initializeQMEngine",monkeyInitQM)
    monkeypatch.setattr(model.Model,"__init__",monkeyInitModel)

    m = model.Model(args)
    m.restartCycle = 2
    m.doQMCalculations(3)
    assert len(m.qmEngine.dirs) == 3
    
class MonkeyMMEngine():
    def __init__(self, args):
        self.didRestart = False

    def restart(self):
        self.didRestart = True

    def getMMSamples(self):
        pass

class MonkeyOptEngine():
    def __init__(self, args):
        self.restartCycle = 2


def monkeyInitOpt(self, args):
    return MonkeyOptEngine(args)

def monkeyInitMM(self, args):
    return MonkeyMMEngine(args)

def test_doMMsampling(monkeypatch):
    args = FakeArgs()
    os.chdir(os.path.join(os.path.dirname(__file__), "model", "test2"))
    path = os.path.join("2_sampling", "3_cycle_3")
    if os.path.isdir(path):
        rmtree(path)
    monkeypatch.setattr(model.Model,"initializeOptEngine",monkeyInitOpt)
    monkeypatch.setattr(model.Model,"initializeQMEngine",monkeyInit)
    monkeypatch.setattr(model.Model,"initializeMMEngine",monkeyInitMM)
    m = model.Model(args)
    m.restartCycle = 2
    m.doMMSampling(3)
    isDir = os.path.isdir(path)
    isFile = True
    for f in ["conf.pdb", "setup.leap", "md.in", "heat1.in", "sol.mol2", "sol.frcmod"]:
        if not os.path.isfile(os.path.join(path, f)):
            isFile = False
    with open(os.path.join(path, "sol.frcmod"), 'r') as f:
        testLines = f.readlines()
    with open(os.path.join("1_opt","result","opt_2","sol.frcmod"), 'r') as f:
        refLines = f.readlines()
    rmtree(path)
    assert isFile
    assert isDir
    assert m.mmEngine.didRestart
    assert checkUtils.checkLists(testLines, refLines)

def test_doMMsamplingNoRestart(monkeypatch):
    args = FakeArgs()
    os.chdir(os.path.join(os.path.dirname(__file__), "model", "test2"))
    path = os.path.join("2_sampling", "3_cycle_3")
    if os.path.isdir(path):
        rmtree(path)
    monkeypatch.setattr(model.Model,"initializeOptEngine",monkeyInitOpt)
    monkeypatch.setattr(model.Model,"initializeQMEngine",monkeyInit)
    monkeypatch.setattr(model.Model,"initializeMMEngine",monkeyInitMM)
    m = model.Model(args)
    m.restartCycle = -1
    m.doMMSampling(3)
    rmtree(path)
    assert not m.mmEngine.didRestart

def monkeySetupFiles(self, i):
    pass

def monkeySortParams(self, results, i):
    pass

def monkeyGraphResults(self):
    pass

def monkeyOptInit(self, args):
    self.nvalids = args['nvalids']
    for f in os.listdir(args['optdir']):
        if f.endswith(".mol2"):
            self.mol2 = f
        elif f.endswith(".frcmod"):
            self.frcmod = f
    self.validPrevious = []
    self.train = []
    self.valid = []
    self.validInitial = []
    self.restartCycle = 0
    self.respPriors = None
        
def monkeyForceBalance(command):
    out = command.split()[3]
    cycle = out.split("_")[1].split(".")[0]
    copyfile(os.path.join("reference",out), out)

def clean():
    os.chdir("1_optimization")
    for f in os.listdir():
        if f.endswith(".out"):
            os.remove(f)
    os.chdir("targets")
    for f in os.listdir():
        rmtree(f)
    os.chdir(os.path.join("..",".."))

def test_doParameterOptimization(monkeypatch):
    args = FakeArgs()
    args.optdir = "1_optimization"
    args.sampledir = "2_mm_sampling"
    args.nvalids = 2
    os.chdir(os.path.join(os.path.dirname(__file__), "model", "test3"))
    clean()
    monkeypatch.setattr(model.Model,"initializeQMEngine",monkeyInit)
    monkeypatch.setattr(model.Model,"initializeMMEngine",monkeyInit)
    monkeypatch.setattr(optengine.OptEngine,"setupInputFiles",monkeySetupFiles)
    monkeypatch.setattr(optengine.OptEngine,"__init__",monkeyOptInit)
    monkeypatch.setattr(os,"system",monkeyForceBalance)
    monkeypatch.setattr(optengine.OptEngine,"sortParams",monkeySortParams)
    monkeypatch.setattr(optengine.OptEngine,"graphResults",monkeyGraphResults)
    
    m = model.Model(args)
    m.doParameterOptimization(1)
    assert len(m.optResults) == 3
    optResults2 = m.doParameterOptimization(2)
    assert len(m.optResults) == 4
    files = ["setup.leap", "setup_valid_initial.leap", "conf.pdb", "all.mdcrd", "qdata.txt"]
    copiedFiles = True
    for f in files:
        for target in ["train_1", "valid_1", "valid_1_1"]:
            if not os.path.isfile(os.path.join("1_optimization", "targets", target, f)):
                copiedFiles = False
    assert copiedFiles
    testDir = os.path.join("1_optimization", "targets", "train_1")
    refDir = os.path.join("2_mm_sampling","1_cycle_1","train")
    assert checkUtils.checkFiles(os.path.join(testDir, "all.mdcrd"), os.path.join(refDir, "all.mdcrd"))
    assert checkUtils.checkFiles(os.path.join(testDir, "qdata.txt"), os.path.join(refDir, "qdata.txt"))
    testDir = os.path.join("1_optimization", "targets", "valid_1")
    refDir = os.path.join("2_mm_sampling","1_cycle_1","valid")
    assert checkUtils.checkFiles(os.path.join(testDir, "all.mdcrd"), os.path.join(refDir, "all.mdcrd"))
    assert checkUtils.checkFiles(os.path.join(testDir, "qdata.txt"), os.path.join(refDir, "qdata.txt"))
    testDir = os.path.join("1_optimization", "targets", "valid_1_1")
    refDir = os.path.join("2_mm_sampling","1_cycle_1","valid_1")
    assert checkUtils.checkFiles(os.path.join(testDir, "all.mdcrd"), os.path.join(refDir, "all.mdcrd"))
    assert checkUtils.checkFiles(os.path.join(testDir, "qdata.txt"), os.path.join(refDir, "qdata.txt"))
    clean()

def test_initialCycle():
    pass

def test_init():
    pass

