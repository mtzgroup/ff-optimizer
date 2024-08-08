import os
import pytest
from pathlib import Path
from shutil import copyfile, copytree, rmtree

from ff_optimizer import qmengine

from . import checkUtils
from .test_inputs import getDefaults


def clean():
    for f in Path(".").iterdir():
        if f.name.startswith("scr"):
            rmtree(f)
        if f.name.endswith("out") or f.name.endswith("in"):
            if f.name != "tc.in" and f.name != "tc_backup.in":
                os.remove(f)
        if f.name == "all.mdcrd" or f.name == "qdata.txt":
            os.remove(f)


class TestSbatchEngine:

    # check that sbatch files are being read in correctly
    # def test_readSbatchFile(self):
    #    os.chdir(os.path.dirname(__file__))
    #    sbatchEngine = qmengine.SbatchEngine("qmengine/tc.in","qmengine/tc_backup.in","qmengine/sbatch.sh","test")
    #    for line in sbatchEngine.sbatchLines:
    #        print(line)
    #    assert sbatchEngine.tcVersion == "TeraChem/2021.02-intel-2017.8.262-CUDA-9.0.176"
    #    assert sbatchEngine.sbatchOptions == ['#SBATCH -t 2:00:00\n','#SBATCH --gres gpu:2\n','#SBATCH --mem=16GB\n','#SBATCH --exclude=fire-11-01\n']

    # check that sbatch files are being written correctly
    def test_writeSbatchFile(self):
        os.chdir(os.path.dirname(__file__))
        inp = getDefaults()
        inp.sampledir = Path("slurmengine")
        inp.tctemplate = "tc.in"
        inp.tctemplate_backup = "tc_backup.in"
        inp.sbatchtemplate = "sbatch_template.sh"
        sbatchEngine = qmengine.SlurmEngine(inp)
        sbatchEngine.writeSbatchFile("23", "sbatch_23.sh")
        test = checkUtils.checkFiles("sbatch_23.sh", "slurmengine/sbatch_ref.sh")
        os.remove("sbatch_23.sh")
        assert test

    def test_getQMRefData(self, monkeypatch):
        os.chdir(os.path.dirname(__file__))
        inp = getDefaults()
        inp.sampledir = Path("slurmengine")
        inp.tctemplate = "tc.in"
        inp.tctemplate_backup = "tc_backup.in"
        slurmEngine = qmengine.SlurmEngine(inp)

        def monkeySlurmCommand(self, command):
            return "this is a thing"

        def monkeyWaitForJobs(self, jobIDs):
            for i in range(1, len(jobIDs) + 1):
                copytree(f"ref/scr.{i}", f"scr.{i}")
                copyfile(f"ref/tc_{i}.out", f"tc_{i}.out")

        monkeypatch.setattr(qmengine.SlurmEngine, "waitForJobs", monkeyWaitForJobs)
        monkeypatch.setattr(qmengine.SlurmEngine, "slurmCommand", monkeySlurmCommand)

        xyzs = [Path(f"{i}.xyz") for i in range(1, 26)]
        os.chdir("slurmengine")
        clean()
        slurmEngine.getQMRefData(xyzs)
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

    def test_getQMRefDataResp(self, monkeypatch):
        os.chdir(os.path.dirname(__file__))
        inp = getDefaults()
        inp.sampledir = Path("slurmengine")
        inp.tctemplate = "tc.in"
        inp.tctemplate_backup = "tc_backup.in"
        inp.resp = 1.0
        slurmEngine = qmengine.SlurmEngine(inp)

        def monkeySlurmCommand(self, command):
            return "this is a thing"

        def monkeyWaitForJobs(self, jobIDs):
            for i in range(1, len(jobIDs) + 1):
                copytree(f"ref/scr.{i}", f"scr.{i}")
                copyfile(f"ref/tc_{i}.out", f"tc_{i}.out")

        monkeypatch.setattr(qmengine.SlurmEngine, "waitForJobs", monkeyWaitForJobs)
        monkeypatch.setattr(qmengine.SlurmEngine, "slurmCommand", monkeySlurmCommand)

        xyzs = [Path(f"{i}.xyz") for i in range(1, 26)]
        os.chdir("slurmengine")
        clean()
        slurmEngine.getQMRefData(xyzs)
        removed = True
        for f in os.listdir():
            if f.startswith("scr"):
                removed = False

        copied = True
        for i in range(1, 26):
            if not Path(f"esp_{i}.xyz").exists():
                copied = False

        checkAllMdcrd = checkUtils.checkFiles("all.mdcrd", "ref/all.mdcrd")
        checkQdata = checkUtils.checkFiles("qdata.txt", "ref/qdataResp.txt")
        clean()
        assert removed
        assert copied
        assert checkAllMdcrd
        assert checkQdata
