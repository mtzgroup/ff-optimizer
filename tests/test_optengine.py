import errno
import os
from shutil import copyfile, rmtree

import numpy as np

from ff_optimizer import optengine

from . import checkUtils


def cleanOptDir(optdir):
    try:
        rmtree(os.path.join(optdir, "forcefield"))
    except:
        pass
    if os.path.isfile(os.path.join(optdir, "valid_0_initial.in")):
        os.remove(os.path.join(optdir, "valid_0_initial.in"))
    if os.path.isfile(os.path.join(optdir, "setup_valid_initial.leap")):
        os.remove(os.path.join(optdir, "setup_valid_initial.leap"))
    try:
        rmtree(os.path.join(optdir, "result"))
    except:
        pass
    for f in os.listdir(optdir):
        if f.startswith("prev") or f.endswith(".inpcrd") or f.endswith(".prmtop") or f == "leap.out":
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
    assert checkUtils.checkArray(refParams, results["params"])
    refInitialParams = [
        -9.4824e-01,
        4.7412e-01,
        1.0000e02,
        1.0452e02,
        1.7564e00,
        1.7912e-01,
    ]
    assert checkUtils.checkArray(refInitialParams, results["initialParams"])
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
    optEngine.optdir = "."
    optEngine.setupInputFiles(8)
    with open("opt_8.in", "r") as f:
        testLinesOpt = f.readlines()
    with open(os.path.join("inputs", "opt_8.in")) as f:
        refLinesOpt = f.readlines()
    with open("valid_8.in", "r") as f:
        testLinesValid = f.readlines()
    with open(os.path.join("inputs", "valid_8.in")) as f:
        refLinesValid = f.readlines()
    with open("valid_8_initial.in", "r") as f:
        testLinesValidInitial = f.readlines()
    with open(os.path.join("inputs", "valid_8_initial.in")) as f:
        refLinesValidInitial = f.readlines()
    os.remove("opt_8.in")
    os.remove("valid_8.in")
    os.remove("valid_8_initial.in")
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
    assert os.path.isfile(os.path.join("result", "opt_0", "dasa.mol2"))
    assert os.path.isfile(os.path.join("result", "opt_0", "dasa.frcmod"))
    assert os.path.isfile("opt_0.out")
    with open(os.path.join("forcefield", "dasa.mol2"), "r") as f:
        refLines = f.readlines()
    with open(os.path.join("result", "opt_0", "dasa.mol2"), "r") as f:
        testLines = f.readlines()
    assert checkUtils.checkLists(testLines, refLines)
    with open(os.path.join("forcefield", "dasa.frcmod"), "r") as f:
        refLines = f.readlines()
    with open(os.path.join("result", "opt_0", "dasa.frcmod"), "r") as f:
        testLines = f.readlines()
    assert checkUtils.checkLists(testLines, refLines)
    os.chdir("..")
    cleanOptDir(options["optdir"])


def test_optimizeForcefield1(monkeypatch):
    def monkeyGraph():
        pass

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
    cleanOptDir(options["optdir"])


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
    assert optEngine.restartCycle == 3


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
    assert optEngine.restartCycle == 2


# Failed in a validation
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
    assert optEngine.restartCycle == 2


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
    checkUtils.checkArray(refParams, optEngine.params[:5, 0])
    refParams = [-6.0860e00, -5.3713e00, -4.3450e00, -4.5274e00, -4.4187e00]
    checkUtils.checkArray(refParams, optEngine.params[:5, -1])


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
    checkUtils.checkArray(refParams, optEngine.params[:5, 137])


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
    cleanOptDir(options["optdir"])


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
