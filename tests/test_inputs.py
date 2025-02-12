import os
from pathlib import Path

from ff_optimizer import inputs

home = Path(__file__).parent.absolute()


def fakePostInit(self):
    self.dynamicsdir = Path("0_dynamics")
    self.pathify()


def fakePostInit2(self):
    self.pathify()


def getDefaults(fake=fakePostInit):
    postInit = inputs.Input.__post_init__
    setattr(inputs.Input, "__post_init__", fake)
    with open("nothing.yaml", "w") as f:
        f.write(" ")
    inp = inputs.Input.fromYaml(Path("nothing.yaml"))
    setattr(inputs.Input, "__post_init__", postInit)
    os.remove("nothing.yaml")
    return inp


def test_defaults():
    inp = getDefaults()
    assert isinstance(inp, inputs.Input)
    assert isinstance(inp.optdir, Path)
    assert isinstance(inp.dynamicsdir, Path)


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


def test_post_init():
    os.chdir(home / "inputs")
    inp = inputs.Input.fromYaml("input.yaml")
    assert inp.maxcycles == -1


def test_setupDynamicsFolder():
    os.chdir(home / "inputs")
    inp = getDefaults(fakePostInit2)
    inp.setupDynamicsFolder()
    xyz = Path(inp.optdir) / "conf.xyz"
    madeXYZ = xyz.is_file()
    assert madeXYZ
    os.remove(xyz)
    assert inp.coors == "conf.xyz"
    assert inp.optdir == inp.dynamicsdir


def test_checkFiles1(monkeypatch):
    os.chdir(home / "inputs" / "test1")
    monkeypatch.setattr(inputs.Input, "__post_init__", fakePostInit2)
    inp = inputs.Input.fromYaml("input.yaml")
    test1 = True
    try:
        inp.checkFiles()
    except:
        test1 = False

    os.rename("1_opt", "oops")
    test2 = False
    try:
        inp.checkFiles()
    except:
        test2 = True
    os.rename("oops", "1_opt")

    os.rename("wat.xyz", "oops")
    test3 = False
    try:
        inp.checkFiles()
    except:
        test3 = True
    os.rename("oops", "wat.xyz")
    assert test1
    assert test2
    assert test3


def test_checkFiles2(monkeypatch):
    os.chdir(home / "inputs" / "test2")
    monkeypatch.setattr(inputs.Input, "__post_init__", fakePostInit2)
    inp = inputs.Input.fromYaml("input.yaml")
    # base test should pass checkFiles
    test1 = True
    inp.checkFiles()
    try:
        inp.checkFiles()
    except:
        test1 = False

    os.chdir("1_opt")
    os.remove("conf.xyz")

    # renaming things should not pass checkFiles
    test2 = []
    for i, file in enumerate(["conf.pdb", "setup.leap", inp.opt0, inp.valid0]):
        os.rename(file, "oops")
        temp = False
        try:
            inp.checkFiles()
        except Exception as e:
            print(e)
            temp = True
        test2.append(temp)
        os.rename("oops", file)
    assert test1
    for test in test2:
        assert test


def test_checkFiles3(monkeypatch):
    os.chdir(home / "inputs" / "test3")
    monkeypatch.setattr(inputs.Input, "__post_init__", fakePostInit)
    inp = inputs.Input.fromYaml("input.yaml")
    # base test should pass checkFiles
    test1 = True
    inp.checkFiles()
    try:
        inp.checkFiles()
    except:
        test1 = False

    test2 = []
    os.chdir("0_dynamics")
    for i, file in enumerate([inp.coors, inp.tcout, inp.conformers]):
        os.rename(file, "oops")
        temp = False
        try:
            inp.checkFiles()
        except Exception as e:
            print(e)
            temp = True
        test2.append(temp)
        os.rename("oops", file)

    os.chdir(home / "inputs" / "test3" / "2_sampling")
    test3 = False
    os.rename(inp.sbatchtemplate, "oops")
    try:
        inp.checkFiles()
    except:
        test3 = True
    os.rename("oops", inp.sbatchtemplate)
    assert test1
    for test in test2:
        assert test
    assert test3
