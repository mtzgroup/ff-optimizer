import os
from pathlib import Path

from ff_optimizer import inputs

home = Path(__file__).parent.absolute()


def fakePostInit(self):
    self.pathify()


def getDefaults():
    postInit = inputs.Input.__post_init__
    setattr(inputs.Input, "__post_init__", fakePostInit)
    with open("nothing.yaml", "w") as f:
        f.write(" ")
    inp = inputs.Input.fromYaml("nothing.yaml")
    setattr(inputs.Input, "__post_init__", postInit)
    os.remove("nothing.yaml")
    return inp


def test_defaults():
    inp = getDefaults()
    assert isinstance(inp, inputs.Input)
    assert isinstance(inp.optdir, Path)


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


# @pytest.mark.debug
def test_post_init():
    os.chdir(home / "inputs")
    inp = inputs.Input.fromYaml("input.yaml")
    assert inp.maxcycles == -1


# Should be more tests here?
