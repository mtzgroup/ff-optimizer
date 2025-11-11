import os
from pathlib import Path
from shutil import copyfile, copytree, rmtree

from ff_optimizer import qmengine

from . import checkUtils
from .test_inputs import getDefaults
import pytest

def clean():
    for f in Path(".").iterdir():
        if f.name.startswith("scr"):
            rmtree(f)
        if f.name.endswith("out") or f.name.endswith("in"):
            if f.name != "tc.in" and f.name != "tc_backup.in":
                os.remove(f)
        if f.name == "all.mdcrd" or f.name == "qdata.txt":
            os.remove(f)
    for i in range(1, 26):
        if Path(f"sge_{i}.sh").exists():
            os.remove(f"sge_{i}.sh")

class TestSGEEngine:

    def test_writeSGEScript(self):
        os.chdir(os.path.dirname(__file__))
        inp = getDefaults()
        inp.sampledir = Path("sgeengine")
        inp.tctemplate = "tc.in"
        inp.tctemplate_backup = "tc_backup.in"
        inp.sgetemplate = "sge_template.sh"
        sgeEngine = qmengine.SGEEngine(inp)
        sgeEngine.writeSGEScript("23", "sge_23.sh")
        test = checkUtils.checkFiles("sge_23.sh", "sgeengine/sge_ref.sh")
        os.remove("sge_23.sh")
        assert test

    def test_getQMRefData(self, monkeypatch):
        os.chdir(os.path.dirname(__file__))
        inp = getDefaults()
        inp.sampledir = Path("sgeengine")
        inp.tctemplate = "tc.in"
        inp.tctemplate_backup = "tc_backup.in"
        inp.sgetemplate = "sge_template.sh"
        sgeEngine = qmengine.SGEEngine(inp)

        def monkeySGECommand(self, command):
            return "this is a thing"

        def monkeyWaitForJobs(self, jobIDs):
            for i in range(1, len(jobIDs) + 1):
                copytree(f"ref/scr.{i}", f"scr.{i}")
                copyfile(f"ref/tc_{i}.out", f"tc_{i}.out")

        monkeypatch.setattr(qmengine.SGEEngine, "waitForJobs", monkeyWaitForJobs)
        monkeypatch.setattr(qmengine.SGEEngine, "sgeCommand", monkeySGECommand)

        xyzs = [Path(f"{i}.xyz") for i in range(1, 26)]
        os.chdir("sgeengine")
        clean()
        sgeEngine.getQMRefData(xyzs)
        removed = True
        for f in os.listdir():
            if f.startswith("scr"):
                removed = False
        checkAllMdcrd = checkUtils.checkFiles("all.mdcrd", "ref/all.mdcrd")
        checkQdata = checkUtils.checkFiles("qdata.txt", "ref/qdata.txt")
        clean()
        assert removed
        assert checkAllMdcrd
        assert checkQdata