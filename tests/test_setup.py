import os
import pytest
from pathlib import Path
from shutil import copyfile
import numpy as np

from ff_optimizer import setup, utils

from . import checkUtils


def getSetup():
    xyz = Path("setup") / "test" / "wat.xyz"
    st = setup.Setup(xyz)
    return st


def test_getCharge():
    os.chdir(os.path.dirname(__file__))
    st = getSetup()
    os.chdir("setup")
    charge = st.getCharge("ref_water.mol2")
    assert charge == 0


def test_writeTCFiles():
    os.chdir(os.path.dirname(__file__))
    st = getSetup()
    os.chdir("setup")
    st.writeTCFiles()
    test1 = checkUtils.checkFiles("tc_template.in", "ref_tc_template.in")
    test2 = checkUtils.checkFiles("tc_template_backup.in", "ref_tc_template_backup.in")
    os.remove("tc_template.in")
    os.remove("tc_template_backup.in")
    assert test1 and test2


def test_writeMDFiles():
    os.chdir(os.path.dirname(__file__))
    st = getSetup()
    os.chdir("setup")
    st.writeMDFiles()
    test = True
    for i in range(1, 9):
        test = test and checkUtils.checkFiles(f"heat{i}.in", f"ref_heat{i}.in")
        os.remove(f"heat{i}.in")
    test = test and checkUtils.checkFiles(f"md.in", "ref_md.in")
    os.remove("md.in")
    assert test


def test_setupFF():
    os.chdir(os.path.dirname(__file__))
    st = getSetup()
    optdir = Path("setup").absolute()
    st.inp.optdir = optdir
    st.xyz = Path("water.xyz")
    os.chdir(optdir / "setupFF")
    print(st.inp.optdir)
    st.setupFF()
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
    st = getSetup()
    os.chdir("setup")
    frcmod = "ref_water.frcmod"
    mol2 = "ref_water.mol2"
    resname = "WAT"
    st.setupForceBalance(frcmod, mol2, resname)
    test = checkUtils.checkFiles("opt_0.in", "ref_opt_0.in")
    test = test and checkUtils.checkFiles("valid_0.in", "ref_valid_0.in")
    test = test and checkUtils.checkFiles("setup.leap", "ref_setup.leap")
    os.remove("opt_0.in")
    os.remove("valid_0.in")
    os.remove("setup.leap")
    assert test


def test_editFrcmod():
    os.chdir(os.path.dirname(__file__))
    st = getSetup()
    os.chdir("setup")
    copyfile("dasa.frcmod", "test.frcmod")
    st.editFrcmod("test.frcmod")
    test = checkUtils.checkFiles("test.frcmod", "ref_dasa.frcmod")
    os.remove("test.frcmod")
    assert test


# need to run all of setup because we get the electron count from antechamber
# def test_spinmult():
#    os.chdir(os.path.dirname(__file__))
#    os.chdir("setup")
#    os.chdir("test")
#    st = setup.Setup("wat.xyz", charge=1)
#    inp = st.setup()
#    with open(inp.sampledir / inp.tctemplate, "r") as f:
#        for line in f.readlines():
#            if "spinmult" in line:
#                spinmult = int(line.split()[1])
#    utils.rmrf("1_opt")
#    utils.rmrf("2_sampling")
#    os.remove("input.yaml")
#    assert spinmult == 2


def test_charges():
    os.chdir(Path(os.path.dirname(__file__)))
    st = getSetup()
    os.chdir("setup")
    os.chdir("charges")
    os.system("ff-opt setup --charge -1 bicarb.xyz")
    print(os.listdir("1_opt"))
    mol2charge = st.getCharge(Path("1_opt") / "bicarb.mol2")
    os.chdir(Path(os.path.dirname(__file__)) / "setup" / "charges")
    with open(Path("2_sampling") / "tc_template.in") as f:
        for line in f.readlines():
            if "charge" in line:
                tcCharge = int(line.split()[1])
    with open(Path("2_sampling") / "tc_template_backup.in") as f:
        for line in f.readlines():
            if "charge" in line:
                tcChargeBackup = int(line.split()[1])
    utils.rmrf("1_opt")
    utils.rmrf("2_sampling")
    os.remove("input.yaml")
    assert mol2charge == -1
    assert tcCharge == -1
    assert tcChargeBackup == -1

def test_getSpinMult():
    os.chdir(Path(os.path.dirname(__file__)))
    os.chdir("setup")
    st = setup.Setup("bcla.xyz", 0)
    mult1 = st.spinMult
    st = setup.Setup("bcla.xyz", -7)
    mult2 = st.spinMult
    assert mult1 == 1
    assert mult2 == 2

def test_readCharges(): 
    os.chdir(Path(os.path.dirname(__file__)))
    os.chdir("setup")
    st = setup.Setup("bcla.xyz", 0)
    charges = list(np.loadtxt("charges.txt", dtype=str))
    with open("tc.out", "r") as f:
        lines = list(f.readlines())
    testCharges = st.readCharges(lines)
    assert checkUtils.checkListsFloats(charges, testCharges)

def test_changeCharges():
    os.chdir(Path(os.path.dirname(__file__)))
    os.chdir("setup")
    copyfile("bclg_wrong_charges.mol2", "bclg.mol2")
    charges = list(np.loadtxt("charges.txt", dtype=str))
    st = setup.Setup("bcla.xyz", 0)
    st.changeCharges(charges, "bclg.mol2")
    test = checkUtils.checkFiles("bclg.mol2", "ref_bclg.mol2")
    os.remove("bclg.mol2")
    assert test
    

