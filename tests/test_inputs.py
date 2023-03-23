import os
from pathlib import Path
from shutil import copyfile, rmtree

from ff_optimizer import inputs

from . import checkUtils
import pytest

home = Path(__name__).parent.absolute()

def getDefaults():
    os.chdir(home / "inputs" / "default") 
    inp = inputs.Inputs.fromYaml("nothing.yaml")
    return inp

def test_defaults():
    os.chdir(home / "inputs" / "default") 
    inp = inputs.Inputs.fromYaml("nothing.yaml")
    assert type(inp) is inputs.Inputs
    assert type(inp.generalInput) is inputs.GeneralInput
    assert type(inp.mmInput) is inputs.MMInput
    assert type(inp.qmInput) is inputs.QMInput

def test_checkForFile1():
    os.chdir(home / "inputs")
    try:
        inputs.checkForFile("default", "1_opt")
        failed = False
    except FileNotFoundError:
        failed = True
    assert failed

def test_checkForFile2():
    os.chdir(home / "inputs")
    try:
        inputs.checkForFile("default", "1_opt", isFile=False)
        found = True
    except FileNotFoundError:
        found = False
    assert found
    
# Should be more tests here?
