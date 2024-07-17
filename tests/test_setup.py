import os
from pathlib import Path
from shutil import copyfile, rmtree

from ff_optimizer import inputs, model, setup

from . import checkUtils
from .test_inputs import getDefaults


def test_getCharge():
    os.chdir(os.path.dirname(__file__))
    sampledir = Path("setup")
    charge = setup.getCharge(sampledir)
    assert charge == 0


def test_writeTCFiles():
    os.chdir(os.path.dirname(__file__))
    os.chdir("setup")
    inp = getDefaults()
    setup.writeTCFiles(inp, 0)
    test1 = checkUtils.checkFiles("tc_template.in", "ref_tc_template.in")
    test2 = checkUtils.checkFiles("tc_template_backup.in", "ref_tc_template_backup.in")
    os.remove("tc_template.in")
    os.remove("tc_template_backup.in")
    assert test1 and test2


def test_writeMDFiles():
    os.chdir(os.path.dirname(__file__))
    os.chdir("setup")
    getDefaults()
    setup.writeMDFiles()
    test = True
    for i in range(1, 9):
        test = test and checkUtils.checkFiles(f"heat{i}.in", f"ref_heat{i}.in")
        os.remove(f"heat{i}.in")
    test = test and checkUtils.checkFiles(f"md.in", "ref_md.in")
    os.remove("md.in")
    assert test


def test_setupFF():
    os.chdir(os.path.dirname(__file__))
    optdir = Path("setup").absolute()
    os.chdir(optdir / "setupFF")
    setup.setupFF(Path("water.xyz"), optdir)
    for f in os.listdir():
        if not f.endswith(".xyz"):
            os.remove(f)
    os.chdir(optdir)
    test = checkUtils.checkFiles("water.mol2", "ref_water.mol2")
    test = test and checkUtils.checkFiles("water.frcmod", "ref_water.frcmod")
    os.remove("water.mol2")
    os.remove("water.frcmod")
    assert test


def test_setupForceBalance():
    os.chdir(os.path.dirname(__file__))
    os.chdir("setup")
    inp = getDefaults()
    frcmod = "ref_water.frcmod"
    mol2 = "ref_water.mol2"
    resname = "WAT"
    setup.setupForceBalance(inp, frcmod, mol2, resname)
    test = checkUtils.checkFiles("opt_0.in", "ref_opt_0.in")
    test = test and checkUtils.checkFiles("valid_0.in", "ref_valid_0.in")
    test = test and checkUtils.checkFiles("setup.leap", "ref_setup.leap")
    os.remove("opt_0.in")
    os.remove("valid_0.in")
    os.remove("setup.leap")
    assert test


def test_editFrcmod():
    os.chdir(os.path.dirname(__file__))
    os.chdir("setup")
    copyfile("dasa.frcmod", "test.frcmod")
    setup.editFrcmod("test.frcmod")
    test = checkUtils.checkFiles("test.frcmod", "ref_dasa.frcmod")
    os.remove("test.frcmod")
    assert test

def test_setup():
    os.chdir(os.path.dirname(__file__))
    os.chdir("setup")
    os.chdir("test")
    # if the new input file and model get initialized and pass all their
    # internal checks, then we pass this test as well
    newInp = setup.setup("wat.xyz")
    rmtree("1_opt")
    rmtree("2_sampling")
    os.remove("new_input.yaml")
