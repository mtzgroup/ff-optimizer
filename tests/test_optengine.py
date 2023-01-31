import errno
import pytest
import os
from shutil import copyfile, rmtree

import numpy as np

from ff_optimizer import optengine

from . import checkUtils


def monkeyGraph():
    pass


def cleanOptDir(optdir, removeResult=False):
    try:
        rmtree(os.path.join(optdir, "forcefield"))
    except:
        pass
    if os.path.isfile(os.path.join(optdir, "valid_0_initial.in")):
        os.remove(os.path.join(optdir, "valid_0_initial.in"))
    if os.path.isfile(os.path.join(optdir, "setup_valid_initial.leap")):
        os.remove(os.path.join(optdir, "setup_valid_initial.leap"))
    if removeResult:
        try:
            rmtree(os.path.join(optdir, "result"))
        except:
            pass
    for f in os.listdir(optdir):
        if (
            f.startswith("prev")
            or f.endswith(".inpcrd")
            or f.endswith(".prmtop")
            or f == "leap.out"
            or f == "leap.log"
        ):
            os.remove(os.path.join(optdir, f))


def test_init1():
    os.chdir(os.path.dirname(__file__))
    options = {}
    options["optdir"] = os.path.join("optengine", "opt1")
    cleanOptDir(options["optdir"])
    options["resp"] = 0
    options["respPriors"] = 0
    options["restart"] = False
    options["maxCycles"] = 10
    options["nvalids"] = 1
    optEngine = optengine.OptEngine(options)

    assert not optEngine.doResp
    assert optEngine.restartCycle == -1
    refTargetLines = [
        "$target\n",
        "type                abinitio_amber\n",
        "name                dynamics\n",
        "amber_leapcmd       setup.leap\n",
        "$end\n",
    ]
    assert checkUtils.checkLists(refTargetLines, optEngine.targetLines)

    refTargetLines = [
        "$target\n",
        "type                abinitio_amber\n",
        "name                dynamics\n",
        "amber_leapcmd       setup.leap\n",
        "$end\n",
    ]
    assert checkUtils.checkLists(refTargetLines, optEngine.validTargetLines)
    refTargetLines

    refTargetLines = [
        "$target\n",
        "type                abinitio_amber\n",
        "name                dynamics\n",
        "amber_leapcmd       setup_valid_initial.leap\n",
        "$end\n",
    ]
    assert checkUtils.checkLists(refTargetLines, optEngine.validInitialTargetLines)

    files = sorted(os.listdir(os.path.join(options["optdir"], "forcefield")))
    refFiles = ["dasa.frcmod", "dasa.mol2", "initial_dasa.frcmod", "initial_dasa.mol2"]
    assert checkUtils.checkLists(files, refFiles)
    rmtree(os.path.join(options["optdir"], "forcefield"))

    assert optEngine.mol2 == "dasa.mol2"
    assert optEngine.frcmod == "dasa.frcmod"
    assert optEngine.initialTarget == "dynamics"

    with open(os.path.join(options["optdir"], "valid_0_initial.in"), "r") as f:
        lines = f.readlines()

    refLines = [
        "$options\n",
        "jobtype             single\n",
        "forcefield          initial_dasa.frcmod initial_dasa.mol2\n",
        "backup              false\n",
        "$end\n",
        "\n",
    ]
    assert checkUtils.checkLists(lines, refLines)
    cleanOptDir(options["optdir"])


def test_init2():
    os.chdir(os.path.dirname(__file__))
    options = {}
    options["optdir"] = os.path.join("optengine", "opt2")
    cleanOptDir(options["optdir"])
    options["resp"] = 0
    options["restart"] = False
    options["maxCycles"] = 10
    options["respPriors"] = 0
    options["nvalids"] = 1
    optEngine = optengine.OptEngine(options)

    assert not optEngine.doResp
    assert optEngine.restartCycle == -1
    refTargetLines = [
        "$target\n",
        "type                  abinitio_amber        # The target type; fitting ab initio data using Amber\n",
        'name stuff                       # Also the subdirectory containing data within "targets"\n',
        "weight                1.0                   # Weight of this target relative to others\n",
        "writelevel            1                     # Amount of data printed to temp directory\n",
        "amber_leapcmd         setup.leap\n",
        "$end\n",
    ]
    checkUtils.checkLists(refTargetLines, optEngine.targetLines)
    assert checkUtils.checkLists(refTargetLines, optEngine.targetLines)

    refTargetLines = [
        "$target\n",
        "type                  abinitio_amber        # The target type; fitting ab initio data using Amber\n",
        'name stuff                       # Also the subdirectory containing data within "targets"\n',
        "weight                1.0                   # Weight of this target relative to others\n",
        "writelevel            1                     # Amount of data printed to temp directory\n",
        "amber_leapcmd         setup.leap\n",
        "$end\n",
    ]
    assert checkUtils.checkLists(refTargetLines, optEngine.validTargetLines)

    refTargetLines = [
        "$target\n",
        "type                  abinitio_amber        # The target type; fitting ab initio data using Amber\n",
        'name stuff                       # Also the subdirectory containing data within "targets"\n',
        "weight                1.0                   # Weight of this target relative to others\n",
        "writelevel            1                     # Amount of data printed to temp directory\n",
        "amber_leapcmd         setup_valid_initial.leap\n",
        "$end\n",
    ]
    assert checkUtils.checkLists(refTargetLines, optEngine.validInitialTargetLines)

    files = sorted(os.listdir(os.path.join(options["optdir"], "forcefield")))
    refFiles = ["initial_sol.frcmod", "initial_sol.mol2", "sol.frcmod", "sol.mol2"]
    assert checkUtils.checkLists(files, refFiles)
    rmtree(os.path.join(options["optdir"], "forcefield"))

    assert optEngine.mol2 == "sol.mol2"
    assert optEngine.frcmod == "sol.frcmod"
    assert optEngine.initialTarget == "stuff"

    with open(os.path.join(options["optdir"], "valid_0_initial.in"), "r") as f:
        lines = f.readlines()

    refLines = [
        "$options\n",
        "jobtype             single\n",
        "forcefield          initial_sol.frcmod initial_sol.mol2\n",
        "backup              false\n",
        "$end\n",
        "\n",
    ]
    assert checkUtils.checkLists(lines, refLines)
    cleanOptDir(options["optdir"])


def test_init1_RESP():
    os.chdir(os.path.dirname(__file__))
    options = {}
    options["optdir"] = os.path.join("optengine", "opt1")
    cleanOptDir(options["optdir"])
    options["resp"] = 1
    options["restart"] = False
    options["maxCycles"] = 10
    options["respPriors"] = 0
    options["nvalids"] = 1
    optEngine = optengine.OptEngine(options)

    assert optEngine.doResp
    assert optEngine.restartCycle == -1
    refTargetLines = [
        "$target\n",
        "resp 1\n",
        "w_resp 1\n",
        "type                abinitio_amber\n",
        "name                dynamics\n",
        "amber_leapcmd       setup.leap\n",
        "$end\n",
    ]
    assert checkUtils.checkLists(refTargetLines, optEngine.targetLines)

    refTargetLines = [
        "$target\n",
        "type                abinitio_amber\n",
        "name                dynamics\n",
        "amber_leapcmd       setup.leap\n",
        "$end\n",
    ]
    assert checkUtils.checkLists(refTargetLines, optEngine.validTargetLines)
    refTargetLines

    refTargetLines = [
        "$target\n",
        "type                abinitio_amber\n",
        "name                dynamics\n",
        "amber_leapcmd       setup_valid_initial.leap\n",
        "$end\n",
    ]
    assert checkUtils.checkLists(refTargetLines, optEngine.validInitialTargetLines)

    files = sorted(os.listdir(os.path.join(options["optdir"], "forcefield")))
    refFiles = ["dasa.frcmod", "dasa.mol2", "initial_dasa.frcmod", "initial_dasa.mol2"]
    assert checkUtils.checkLists(files, refFiles)
    rmtree(os.path.join(options["optdir"], "forcefield"))

    assert optEngine.mol2 == "dasa.mol2"
    assert optEngine.frcmod == "dasa.frcmod"
    assert optEngine.initialTarget == "dynamics"

    with open(os.path.join(options["optdir"], "valid_0_initial.in"), "r") as f:
        lines = f.readlines()

    refLines = [
        "$options\n",
        "jobtype             single\n",
        "forcefield          initial_dasa.frcmod initial_dasa.mol2\n",
        "backup              false\n",
        "$end\n",
        "\n",
    ]
    assert checkUtils.checkLists(lines, refLines)
    cleanOptDir(options["optdir"])


def test_readOpt1():
    os.chdir(os.path.dirname(__file__))
    options = {}
    options["optdir"] = os.path.join("optengine", "opt1")
    options["resp"] = 1
    options["restart"] = False
    options["maxCycles"] = 10
    options["respPriors"] = 0
    options["nvalids"] = 1
    optEngine = optengine.OptEngine(options)
    cleanOptDir(options["optdir"])

    status, results = optEngine.readOpt(os.path.join("optengine", "maxsteps.out"))
    assert status == 2


def test_readOpt2():
    os.chdir(os.path.dirname(__file__))
    options = {}
    options["optdir"] = os.path.join("optengine", "opt1")
    options["resp"] = 1
    options["restart"] = False
    options["maxCycles"] = 10
    options["respPriors"] = 0
    options["nvalids"] = 1
    optEngine = optengine.OptEngine(options)
    cleanOptDir(options["optdir"])

    status, results = optEngine.readOpt(os.path.join("optengine", "fail.out"))
    assert status == -1


def test_readOpt3():
    os.chdir(os.path.dirname(__file__))
    options = {}
    options["optdir"] = os.path.join("optengine", "opt1")
    options["resp"] = 1
    options["restart"] = False
    options["maxCycles"] = 10
    options["respPriors"] = 0
    options["nvalids"] = 1
    optEngine = optengine.OptEngine(options)
    cleanOptDir(options["optdir"])

    status, results = optEngine.readOpt(os.path.join("optengine", "success.out"))
    assert status == 0
    refParams = [-6.5120e-01, 3.2560e-01, 1.0000e02, 1.0452e02, 1.7564e00, 2.9777e-01]
    assert checkUtils.checkArrays(refParams, results["params"])
    refInitialParams = [
        -9.4824e-01,
        4.7412e-01,
        1.0000e02,
        1.0452e02,
        1.7564e00,
        1.7912e-01,
    ]
    assert checkUtils.checkArrays(refInitialParams, results["initialParams"])
    refLabels = [
        "COUL:SOL-1",
        "COUL:SOL-2",
        "ANGLESK/HW-OW-HW",
        "ANGLESB/HW-OW-HW",
        "VDWS/OW",
        "VDWT/OW",
    ]
    assert checkUtils.checkLists(refLabels, results["labels"])


def test_addTargetLines1():
    os.chdir(os.path.join(os.path.dirname(__file__), "optengine"))
    options = {}
    options["optdir"] = "opt1"
    options["resp"] = 1
    options["restart"] = False
    options["maxCycles"] = 10
    options["respPriors"] = 0
    options["nvalids"] = 1
    optEngine = optengine.OptEngine(options)
    cleanOptDir(options["optdir"])
    copyfile("opt_4.in", "opt_4_2.in")
    optEngine.addTargetLines("opt_4_2.in", optEngine.targetLines, "train_4")
    with open("opt_4.in", "r") as f:
        refLines = f.readlines()
    with open("opt_4_2.in", "r") as f:
        testLines = f.readlines()
    os.remove("opt_4_2.in")
    assert checkUtils.checkLists(refLines, testLines)


def test_addTargetLines2():
    os.chdir(os.path.join(os.path.dirname(__file__), "optengine"))
    options = {}
    options["optdir"] = "opt1"
    options["resp"] = 0
    options["restart"] = False
    options["maxCycles"] = 10
    options["respPriors"] = 0
    options["nvalids"] = 1
    optEngine = optengine.OptEngine(options)
    cleanOptDir(options["optdir"])
    copyfile("opt_4.in", "opt_5.in")
    optEngine.addTargetLines("opt_5.in", optEngine.targetLines, "train_5")
    with open("opt_4.in", "r") as f:
        refLines = f.readlines()
    refLines += [
        "\n",
        "$target\n",
        "type                abinitio_amber\n",
        "name                train_5\n",
        "amber_leapcmd       setup.leap\n",
        "$end\n",
    ]
    with open("opt_5.in", "r") as f:
        testLines = f.readlines()
    os.remove("opt_5.in")
    assert checkUtils.checkLists(refLines, testLines)


def test_readValid1():
    os.chdir(os.path.join(os.path.dirname(__file__), "optengine"))
    options = {}
    options["optdir"] = "opt1"
    options["resp"] = 0
    options["restart"] = False
    options["maxCycles"] = 10
    options["respPriors"] = 0
    options["nvalids"] = 1
    optEngine = optengine.OptEngine(options)
    cleanOptDir(options["optdir"])
    test = optEngine.readValid("valid.out")
    assert checkUtils.checkFloat(test, 1.39869922, 0.000000001)


def test_readValid2():
    os.chdir(os.path.join(os.path.dirname(__file__), "optengine"))
    options = {}
    options["optdir"] = "opt1"
    options["resp"] = 0
    options["restart"] = False
    options["maxCycles"] = 10
    options["respPriors"] = 0
    options["nvalids"] = 1
    optEngine = optengine.OptEngine(options)
    cleanOptDir(options["optdir"])
    try:
        optEngine.readValid("validFail.out")
        assert False
    except RuntimeError:
        pass


def test_setupInputFiles():
    os.chdir(os.path.join(os.path.dirname(__file__), "optengine"))
    options = {}
    options["optdir"] = "opt1"
    options["resp"] = 0.5
    options["restart"] = False
    options["maxCycles"] = 10
    options["respPriors"] = 0
    options["nvalids"] = 1
    optEngine = optengine.OptEngine(options)
    cleanOptDir(options["optdir"])
    optEngine.optdir = "inputs"
    os.chdir("inputs")
    optEngine.setupInputFiles(8)
    with open("opt_8.in", "r") as f:
        testLinesOpt = f.readlines()
    with open(os.path.join("..", "ref_inputs", "opt_8.in")) as f:
        refLinesOpt = f.readlines()
    with open("valid_8.in", "r") as f:
        testLinesValid = f.readlines()
    with open(os.path.join("..", "ref_inputs", "valid_8.in")) as f:
        refLinesValid = f.readlines()
    with open("valid_8_initial.in", "r") as f:
        testLinesValidInitial = f.readlines()
    with open(os.path.join("..", "ref_inputs", "valid_8_initial.in")) as f:
        refLinesValidInitial = f.readlines()
    os.remove("opt_8.in")
    os.remove("valid_8.in")
    os.remove("valid_8_initial.in")
    os.chdir("..")
    checkUtils.checkLists(testLinesOpt, refLinesOpt)
    checkUtils.checkLists(testLinesValid, refLinesValid)
    checkUtils.checkLists(testLinesValidInitial, refLinesValidInitial)


def monkeyForceBalance(command):
    split = command.split()
    inp = split[1]
    out = split[3]
    name = inp.split("_")[1].replace(".in", "")
    with open(inp, "r") as f:
        for line in f.readlines():
            splitLine = line.split()
            if len(splitLine) > 0:
                if splitLine[0] == "forcefield":
                    ffFiles = splitLine[1:]
    if not os.path.isfile(inp):
        raise FileNotFoundError(errno.ENOENT, os.strerror(errno.ENOENT), inp)
    if not os.path.isfile(os.path.join("ref", out)):
        raise FileNotFoundError(
            errno.ENOENT, os.strerror(errno.ENOENT), os.path.join("ref", out)
        )
    copyfile(os.path.join("ref", out), out)
    if inp.split("_")[0] == "opt":
        if not os.path.isdir("result"):
            os.mkdir("result")
        if not os.path.isdir(os.path.join("result", f"opt_{name}")):
            os.mkdir(os.path.join("result", f"opt_{name}"))
        for f in ffFiles:
            copyfile(os.path.join("ref", f), os.path.join("result", f"opt_{name}", f))
            if f.endswith(".mol2"):
                copyfile(os.path.join("forcefield", f), f"prev_{f}")


def test_optimizeForcefield0(monkeypatch):
    os.chdir(os.path.join(os.path.dirname(__file__), "optengine"))
    options = {}
    options["optdir"] = "opt1"
    cleanOptDir(options["optdir"])
    options["resp"] = 0
    options["restart"] = False
    options["maxCycles"] = 10
    options["respPriors"] = 0
    options["nvalids"] = 1
    optEngine = optengine.OptEngine(options)
    monkeypatch.setattr(os, "system", monkeyForceBalance)
    os.chdir("opt1")
    optEngine.optimizeForcefield(0)
    filesFound = True
    if not os.path.isfile(os.path.join("result", "opt_0", "dasa.mol2")):
        filesFound = False
    if not os.path.isfile(os.path.join("result", "opt_0", "dasa.frcmod")):
        filesFound = False
    if not os.path.isfile("opt_0.out"):
        filesFound = False
    with open(os.path.join("forcefield", "dasa.mol2"), "r") as f:
        refMol2Lines = f.readlines()
    with open(os.path.join("result", "opt_0", "dasa.mol2"), "r") as f:
        testMol2Lines = f.readlines()
    with open(os.path.join("forcefield", "dasa.frcmod"), "r") as f:
        refLines = f.readlines()
    with open(os.path.join("result", "opt_0", "dasa.frcmod"), "r") as f:
        testLines = f.readlines()
    os.remove("opt_0.out")
    os.chdir("..")
    cleanOptDir(options["optdir"], removeResult=True)
    assert filesFound
    assert checkUtils.checkLists(testMol2Lines, refMol2Lines)
    assert checkUtils.checkLists(testLines, refLines)


def test_optimizeForcefield1(monkeypatch):
    os.chdir(os.path.join(os.path.dirname(__file__), "optengine"))
    options = {}
    options["optdir"] = "restart1"
    cleanOptDir(options["optdir"])
    options["resp"] = 0
    options["restart"] = True
    options["maxCycles"] = 10
    options["respPriors"] = 0
    options["nvalids"] = 1
    optEngine = optengine.OptEngine(options)
    monkeypatch.setattr(os, "system", monkeyForceBalance)
    monkeypatch.setattr(optEngine, "graphResults", monkeyGraph)
    os.chdir("restart1")
    if not os.path.isdir("result"):
        os.mkdir("result")
    if not os.path.isdir(os.path.join("result", "opt_3")):
        os.mkdir(os.path.join("result", "opt_3"))
    copyfile(
        os.path.join("ref", "dasa_prev.mol2"),
        os.path.join("result", "opt_3", "dasa.mol2"),
    )
    copyfile(
        os.path.join("ref", "dasa_prev.frcmod"),
        os.path.join("result", "opt_3", "dasa.frcmod"),
    )
    optEngine.optimizeForcefield(4)
    assert os.path.isfile(os.path.join("result", "opt_4", "dasa.mol2"))
    assert os.path.isfile(os.path.join("result", "opt_4", "dasa.frcmod"))
    assert os.path.isfile("opt_4.out")
    with open(os.path.join("forcefield", "dasa.mol2"), "r") as f:
        refLines = f.readlines()
    with open(os.path.join("result", "opt_4", "dasa.mol2"), "r") as f:
        testLines = f.readlines()
    assert checkUtils.checkLists(testLines, refLines)
    with open(os.path.join("forcefield", "dasa.frcmod"), "r") as f:
        refLines = f.readlines()
    with open(os.path.join("result", "opt_4", "dasa.frcmod"), "r") as f:
        testLines = f.readlines()
    assert checkUtils.checkLists(testLines, refLines)
    os.chdir("..")
    cleanOptDir(options["optdir"], removeResult=True)


# Finished FB cycle
def test_determineRestart1():
    os.chdir(os.path.join(os.path.dirname(__file__), "optengine"))
    for f in os.listdir("restart1"):
        if "4" in f:
            os.remove(os.path.join("restart1", f))
    options = {}
    options["optdir"] = "restart1"
    options["resp"] = 0
    options["restart"] = True
    options["maxCycles"] = 10
    options["respPriors"] = 0
    options["nvalids"] = 1
    optEngine = optengine.OptEngine(options)
    cleanOptDir(options["optdir"])
    assert checkUtils.checkFloat(optEngine.params[0, 0], 266.5)
    assert checkUtils.checkFloat(optEngine.params[0, -1], -6.0860)
    assert checkUtils.checkFloat(optEngine.params[4, 0], 235.8)
    assert checkUtils.checkFloat(optEngine.params[4, -1], -4.4187)
    assert len(optEngine.validPrevious) == 3
    assert len(optEngine.train) == 4
    assert len(optEngine.validInitial) == 3
    assert len(optEngine.valid) == 3
    assert optEngine.restartCycle == 4


# Failed in opt
def test_determineRestart2():
    os.chdir(os.path.join(os.path.dirname(__file__), "optengine"))
    options = {}
    options["optdir"] = "restart2"
    options["resp"] = 0
    options["restart"] = True
    options["maxCycles"] = 10
    options["respPriors"] = 0
    options["nvalids"] = 1
    optEngine = optengine.OptEngine(options)
    cleanOptDir(options["optdir"])
    assert len(optEngine.validPrevious) == 3
    assert len(optEngine.train) == 3
    assert len(optEngine.valid) == 2
    assert len(optEngine.validInitial) == 2
    assert optEngine.restartCycle == 3


# Failed in valid_initial
def test_determineRestart3():
    os.chdir(os.path.join(os.path.dirname(__file__), "optengine"))
    options = {}
    options["optdir"] = "restart3"
    options["resp"] = 0
    options["restart"] = True
    options["maxCycles"] = 10
    options["respPriors"] = 0
    options["nvalids"] = 1
    optEngine = optengine.OptEngine(options)
    cleanOptDir(options["optdir"])
    print(optEngine.params[:6, 0])
    assert checkUtils.checkFloat(optEngine.params[0, 0], 297.1)
    assert checkUtils.checkFloat(optEngine.params[1, 0], 223.05)
    assert checkUtils.checkFloat(optEngine.params[2, 0], 216.97)
    assert checkUtils.checkFloat(optEngine.params[3, 0], 245.75)
    assert checkUtils.checkFloat(optEngine.params[4, 0], 245.28)
    assert checkUtils.checkFloat(optEngine.params[5, 0], 0)
    assert len(optEngine.train) == 4
    assert len(optEngine.valid) == 3
    assert len(optEngine.validPrevious) == 3
    assert len(optEngine.validInitial) == 2
    assert optEngine.restartCycle == 3


def test_sortParams1():
    os.chdir(os.path.join(os.path.dirname(__file__), "optengine"))
    for f in os.listdir("restart1"):
        if "4" in f:
            os.remove(os.path.join("restart1", f))
    options = {}
    options["optdir"] = "restart1"
    options["resp"] = 0
    options["restart"] = True
    options["maxCycles"] = 10
    options["respPriors"] = 0
    options["nvalids"] = 1
    optEngine = optengine.OptEngine(options)
    cleanOptDir(options["optdir"])
    assert optEngine.labels[0] == "BONDSK/cy-c3"
    assert optEngine.labels[-1] == "IDIHSK/c2-h4-ce-n3"
    refParams = [2.6650e02, 1.6724e02, 1.9135e02, 2.1897e02, 2.3580e02]
    assert checkUtils.checkArrays(refParams, optEngine.params[:5, 0])
    refParams = [-6.0860e00, -5.3713e00, -4.3450e00, -4.5274e00, -4.4187e00]
    checkUtils.checkArrays(refParams, optEngine.params[:5, -1])


def test_sortParams2():
    os.chdir(os.path.join(os.path.dirname(__file__), "optengine"))
    options = {}
    options["optdir"] = "restart3"
    options["resp"] = 0
    options["restart"] = True
    options["maxCycles"] = 10
    options["respPriors"] = 0
    options["nvalids"] = 1
    optEngine = optengine.OptEngine(options)
    cleanOptDir(options["optdir"])
    refParams = [1.6289e-01, 1.6289e-01, 1.4744e-01, 1.4744e-01, 1.4744e-01]
    assert checkUtils.checkArrays(refParams, optEngine.params[:5, 137])


def test_respPriors(monkeypatch):
    def monkeyGraph():
        pass

    os.chdir(os.path.join(os.path.dirname(__file__), "optengine"))
    options = {}
    options["sampledir"] = os.path.join("..", "resp", "sample")
    options["respPriors"] = 1
    options["optdir"] = "restart1"
    cleanOptDir(options["optdir"])
    options["resp"] = 0
    options["restart"] = True
    options["maxCycles"] = 10
    options["nvalids"] = 1
    optEngine = optengine.OptEngine(options)
    optEngine.respPriors.allResp
    monkeypatch.setattr(os, "system", monkeyForceBalance)
    monkeypatch.setattr(optEngine, "graphResults", monkeyGraph)
    os.chdir("restart1")
    if not os.path.isdir("result"):
        os.mkdir("result")
    if not os.path.isdir(os.path.join("result", "opt_3")):
        os.mkdir(os.path.join("result", "opt_3"))
    copyfile(
        os.path.join("ref", "dasa_prev.mol2"),
        os.path.join("result", "opt_3", "dasa.mol2"),
    )
    copyfile(
        os.path.join("ref", "dasa_prev.frcmod"),
        os.path.join("result", "opt_3", "dasa.frcmod"),
    )
    optEngine.optimizeForcefield(4)
    copyfile("opt_4.in", "priors.in")
    assert checkUtils.checkFileFloatsNoWhitespace(
        "opt_4.in", os.path.join("ref", "opt_4_priors.in")
    )

    assert checkUtils.checkFileFloatsNoWhitespace(
        "prev_dasa.mol2", os.path.join("ref", "dasa_priors.mol2")
    )
    os.chdir("..")
    cleanOptDir(options["optdir"], removeResult=True)


def test_restartResp():
    os.chdir(os.path.join(os.path.dirname(__file__), "optengine"))
    for f in os.listdir("restart1"):
        if "4" in f:
            os.remove(os.path.join("restart1", f))
    options = {}
    options["optdir"] = "restart1"
    options["resp"] = 0
    options["restart"] = True
    options["maxCycles"] = 10
    options["respPriors"] = 1
    options["sampledir"] = os.path.join("..", "resp", "sample")
    options["nvalids"] = 1
    optEngine = optengine.OptEngine(options)
    cleanOptDir(options["optdir"])
    arr = np.asarray(optEngine.respPriors.allResp, dtype=np.float32)
    assert checkUtils.checkFloat(arr[5, 0], 2.0)
    assert checkUtils.checkFloat(arr[8, 0], 3.0)


def test_setupInputFiles_multipleValids():
    os.chdir(os.path.join(os.path.dirname(__file__), "optengine"))
    options = {}
    options["optdir"] = "opt1"
    options["resp"] = 0.5
    options["restart"] = False
    options["maxCycles"] = 10
    options["respPriors"] = 0
    options["nvalids"] = 3
    optEngine = optengine.OptEngine(options)
    cleanOptDir(options["optdir"])
    optEngine.optdir = "inputs"
    os.chdir("inputs")
    optEngine.setupInputFiles(8)
    with open("valid_8_1.in", "r") as f:
        testLinesValid1 = f.readlines()
    with open("valid_8_2.in", "r") as f:
        testLinesValid2 = f.readlines()
    with open(os.path.join("..", "ref_inputs", "valid_8_1.in")) as f:
        refLinesValid1 = f.readlines()
    with open(os.path.join("..", "ref_inputs", "valid_8_2.in")) as f:
        refLinesValid2 = f.readlines()
    with open("valid_8_1_initial.in", "r") as f:
        testLinesValidInitial1 = f.readlines()
    with open(os.path.join("..", "ref_inputs", "valid_8_1_initial.in")) as f:
        refLinesValidInitial1 = f.readlines()
    with open("valid_8_2_initial.in", "r") as f:
        testLinesValidInitial2 = f.readlines()
    with open(os.path.join("..", "ref_inputs", "valid_8_2_initial.in")) as f:
        refLinesValidInitial2 = f.readlines()
    os.remove("opt_8.in")
    os.remove("valid_8.in")
    os.remove("valid_8_initial.in")
    os.remove("valid_8_1.in")
    os.remove("valid_8_2.in")
    os.remove("valid_8_1_initial.in")
    os.remove("valid_8_2_initial.in")
    checkUtils.checkLists(testLinesValid1, refLinesValid1)
    checkUtils.checkLists(testLinesValidInitial1, refLinesValidInitial1)
    checkUtils.checkLists(testLinesValid2, refLinesValid2)
    checkUtils.checkLists(testLinesValidInitial2, refLinesValidInitial2)


def test_optimizeForcefield_multipleValids(monkeypatch):
    def monkeyForceBalance(command):
        split = command.split()
        inp = split[1]
        out = split[3]
        if not os.path.isfile(inp):
            raise FileNotFoundError(errno.ENOENT, os.strerror(errno.ENOENT), inp)
        if not os.path.isfile(os.path.join("ref", out)):
            raise FileNotFoundError(
                errno.ENOENT, os.strerror(errno.ENOENT), os.path.join("ref", out)
            )
        copyfile(os.path.join("ref", out), out)

    def monkeySort(results, i):
        return

    os.chdir(os.path.join(os.path.dirname(__file__), "optengine"))
    options = {}
    options["optdir"] = "opt3"
    cleanOptDir(options["optdir"])
    options["resp"] = 0
    options["restart"] = False
    options["maxCycles"] = 10
    options["respPriors"] = 0
    options["nvalids"] = 3
    optEngine = optengine.OptEngine(options)
    monkeypatch.setattr(os, "system", monkeyForceBalance)
    monkeypatch.setattr(optEngine, "graphResults", monkeyGraph)
    monkeypatch.setattr(optEngine, "sortParams", monkeySort)
    os.chdir("opt3")
    os.mkdir("result")
    os.mkdir(os.path.join("result", "opt_1"))
    for f in ["dasa.frcmod", "dasa.mol2"]:
        copyfile(f, os.path.join("result", "opt_1", f))
    os.mkdir(os.path.join("result", "opt_0"))
    for f in ["dasa.frcmod", "dasa.mol2"]:
        copyfile(f, os.path.join("result", "opt_0", f))
    optEngine.optimizeForcefield(1)
    files = os.listdir()
    valid1 = 0
    valid2 = 0
    for f in files:
        if "1_1" in f and f.endswith(".out"):
            valid1 += 1
            os.remove(f)
        elif "1_2" in f and f.endswith(".out"):
            valid2 += 1
            os.remove(f)
    os.chdir("..")
    cleanOptDir(options["optdir"], removeResult=True)
    assert valid1 == 3
    assert valid2 == 3


def test_determineRestart_multipleValids():
    os.chdir(os.path.join(os.path.dirname(__file__), "optengine"))
    for f in os.listdir("restart4"):
        if "4" in f:
            os.remove(os.path.join("restart4", f))
    options = {}
    options["optdir"] = "restart4"
    options["resp"] = 0
    options["restart"] = True
    options["maxCycles"] = 10
    options["respPriors"] = 0
    options["nvalids"] = 3
    optEngine = optengine.OptEngine(options)
    cleanOptDir(options["optdir"])
    assert optEngine.restartCycle == 2


def monkeyForcebalance(command):
    with open("fb.log", "a") as f:
        f.write(command + "\n")
    inFile = command.split()[1]
    outFileSplit = command.split()[3].replace(".out", "").split("_")
    if outFileSplit[0] == "opt":
        fbType = "opt"
    else:
        if len(outFileSplit) == 2:
            fbType = "valid"
        else:
            fbType = outFileSplit[2]
    with open(inFile, "r") as f:
        lines = f.readlines()
    for line in lines:
        if len(line.split()) > 0:
            if line.split()[0] == "forcefield":
                for token in line.split():
                    if "frcmod" in token:
                        frcmod = token
    with open(os.path.join("forcefield", frcmod), "r") as f:
        commentLine = f.readline()
    if fbType not in commentLine:
        print(commentLine)
        print(command)
        raise RuntimeError("Ran FB on wrong parameters")


def monkeySortParams(self, results, i):
    pass


def monkeyGraphResults(self):
    pass


def monkeySetupInputFiles(self, i):
    pass


def monkeyReadValid(self, filename):
    return 1


def monkeyReadOpt(self, filename):
    result = {}
    result["obj"] = 1
    return 0, result


def setupFFdir(optdir):
    os.chdir(optdir)
    if os.path.isdir("forcefield"):
        rmtree("forcefield")
    os.mkdir("forcefield")
    copyfile("dasa.frcmod", os.path.join("forcefield", "dasa.frcmod"))
    copyfile("dasa.mol2", os.path.join("forcefield", "dasa.mol2"))
    copyfile("dasa.frcmod", os.path.join("forcefield", "initial_dasa.frcmod"))
    copyfile("dasa.mol2", os.path.join("forcefield", "initial_dasa.mol2"))
    os.chdir("..")


# Restart from new FB cycle
def test_restart1Params(monkeypatch):
    monkeypatch.setattr(optengine.OptEngine, "sortParams", monkeySortParams)
    monkeypatch.setattr(optengine.OptEngine, "graphResults", monkeyGraphResults)
    monkeypatch.setattr(optengine.OptEngine, "setupInputFiles", monkeySetupInputFiles)

    os.chdir(os.path.join(os.path.dirname(__file__), "optengine"))
    options = {}
    options["optdir"] = "params1"
    options["resp"] = 0
    options["restart"] = True
    options["maxCycles"] = 10
    options["respPriors"] = 0
    options["nvalids"] = 1
    optEngine = optengine.OptEngine(options)
    monkeypatch.setattr(os, "system", monkeyForcebalance)
    monkeypatch.setattr(optengine.OptEngine, "readValid", monkeyReadValid)
    monkeypatch.setattr(optengine.OptEngine, "readOpt", monkeyReadOpt)
    setupFFdir("params1")
    os.chdir("params1")
    if os.path.isfile("fb.log"):
        os.remove("fb.log")
    optEngine.optimizeForcefield(optEngine.restartCycle)
    with open("fb.log", "r") as f:
        lines = f.readlines()
    os.remove("fb.log")
    rmtree("forcefield")
    cleanOptDir(".")
    os.chdir("..")
    assert len(lines) == 4


# Restart from failed param opt
def test_restart2Params(monkeypatch):
    monkeypatch.setattr(optengine.OptEngine, "sortParams", monkeySortParams)
    monkeypatch.setattr(optengine.OptEngine, "graphResults", monkeyGraphResults)
    monkeypatch.setattr(optengine.OptEngine, "setupInputFiles", monkeySetupInputFiles)
    testDir = "params2"

    os.chdir(os.path.join(os.path.dirname(__file__), "optengine"))
    options = {}
    options["optdir"] = testDir
    options["resp"] = 0
    options["restart"] = True
    options["maxCycles"] = 10
    options["respPriors"] = 0
    options["nvalids"] = 1
    optEngine = optengine.OptEngine(options)
    monkeypatch.setattr(os, "system", monkeyForcebalance)
    monkeypatch.setattr(optengine.OptEngine, "readValid", monkeyReadValid)
    monkeypatch.setattr(optengine.OptEngine, "readOpt", monkeyReadOpt)
    setupFFdir(testDir)
    os.chdir(testDir)
    if os.path.isfile("fb.log"):
        os.remove("fb.log")
    optEngine.optimizeForcefield(optEngine.restartCycle)
    with open("fb.log", "r") as f:
        lines = f.readlines()
    os.remove("fb.log")
    rmtree("forcefield")
    cleanOptDir(".")
    os.chdir("..")
    assert len(lines) == 3


# Restart from failed validation after param opt
def test_restart3Params(monkeypatch):
    monkeypatch.setattr(optengine.OptEngine, "sortParams", monkeySortParams)
    monkeypatch.setattr(optengine.OptEngine, "graphResults", monkeyGraphResults)
    monkeypatch.setattr(optengine.OptEngine, "setupInputFiles", monkeySetupInputFiles)
    testDir = "params3"

    os.chdir(os.path.join(os.path.dirname(__file__), "optengine"))
    options = {}
    options["optdir"] = testDir
    options["resp"] = 0
    options["restart"] = True
    options["maxCycles"] = 10
    options["respPriors"] = 0
    options["nvalids"] = 1
    optEngine = optengine.OptEngine(options)
    monkeypatch.setattr(os, "system", monkeyForcebalance)
    monkeypatch.setattr(optengine.OptEngine, "readValid", monkeyReadValid)
    monkeypatch.setattr(optengine.OptEngine, "readOpt", monkeyReadOpt)
    setupFFdir(testDir)
    os.chdir(testDir)
    if os.path.isfile("fb.log"):
        os.remove("fb.log")
    optEngine.optimizeForcefield(optEngine.restartCycle)
    with open("fb.log", "r") as f:
        lines = f.readlines()
    os.remove("fb.log")
    rmtree("forcefield")
    cleanOptDir(".")
    os.chdir("..")
    assert len(lines) == 2


# Restart from failed validation on initial params
def test_restart3Params(monkeypatch):
    monkeypatch.setattr(optengine.OptEngine, "sortParams", monkeySortParams)
    monkeypatch.setattr(optengine.OptEngine, "graphResults", monkeyGraphResults)
    monkeypatch.setattr(optengine.OptEngine, "setupInputFiles", monkeySetupInputFiles)
    testDir = "params4"

    os.chdir(os.path.join(os.path.dirname(__file__), "optengine"))
    options = {}
    options["optdir"] = testDir
    options["resp"] = 0
    options["restart"] = True
    options["maxCycles"] = 10
    options["respPriors"] = 0
    options["nvalids"] = 1
    optEngine = optengine.OptEngine(options)
    monkeypatch.setattr(os, "system", monkeyForcebalance)
    monkeypatch.setattr(optengine.OptEngine, "readValid", monkeyReadValid)
    monkeypatch.setattr(optengine.OptEngine, "readOpt", monkeyReadOpt)
    setupFFdir(testDir)
    os.chdir(testDir)
    if os.path.isfile("fb.log"):
        os.remove("fb.log")
    optEngine.optimizeForcefield(optEngine.restartCycle)
    with open("fb.log", "r") as f:
        lines = f.readlines()
    os.remove("fb.log")
    rmtree("forcefield")
    cleanOptDir(".")
    os.chdir("..")
    assert len(lines) == 1
