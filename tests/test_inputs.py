import os
from pathlib import Path

from ff_optimizer import inputs

home = Path(__name__).parent.absolute()


def fakePostInit(self):
    self.pathify()


def getDefaults():
    setattr(inputs.Input, "__post_init__", fakePostInit)
    with open("nothing.yaml", "w") as f:
        f.write(" ")
    inp = inputs.Input.fromYaml("nothing.yaml")
    os.remove("nothing.yaml")
    return inp


def test_defaults():
    inp = getDefaults()
    assert type(inp) is inputs.Input


def test_checkForFile1():
    os.chdir(home / "inputs")
    try:
        inputs.checkForFile(Path("default") / "1_opt")
        failed = False
    except FileNotFoundError:
        failed = True
    assert failed


def test_checkForFile2():
    os.chdir(home / "inputs")
    try:
        inputs.checkForFile(Path("default") / "1_opt", isFile=False)
        found = True
    except FileNotFoundError:
        found = False
    assert found


# Should be more tests here?
