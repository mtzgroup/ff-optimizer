import os
import pytest
from pathlib import Path

from qcio import ProgramInput, ProgramOutput, Provenance, Structure
from qcparse import parse

from ff_optimizer import qmengine

from . import checkUtils
from . import chemcloud_mocks as ccm
from .test_inputs import getDefaults


def test_init():
    os.chdir(os.path.dirname(__file__))
    inp = getDefaults()
    inp.tctemplate = "qmengine/tc.in"
    inp.tctemplate_backup = "qmengine/tc_backup.in"
    inp.sampledir = Path("")
    chemcloudEngine = qmengine.ChemcloudEngine(inp)
    assert chemcloudEngine.specialKeywords["method"] == "b3lyp"
    assert chemcloudEngine.specialKeywords["basis"] == "6-31gss"
    assert chemcloudEngine.keywords["dftd"] == "d3"
    assert chemcloudEngine.specialKeywords["charge"] == 0
    assert chemcloudEngine.specialKeywords["spinmult"] == 1
    assert chemcloudEngine.backupKeywords["threall"] == "1.0e-14"
    assert chemcloudEngine.backupKeywords["diismaxvecs"] == "40"
    assert chemcloudEngine.backupKeywords["maxit"] == "200"
    assert "method" not in chemcloudEngine.keywords.keys()
    assert "method" not in chemcloudEngine.backupKeywords.keys()


def test_loadStructureFromXYZ1():
    os.chdir(os.path.dirname(__file__))
    inp = getDefaults()
    inp.tctemplate = "qmengine/tc.in"
    inp.tctemplate_backup = "qmengine/tc_backup.in"
    inp.sampledir = Path("")
    chemcloudEngine = qmengine.ChemcloudEngine(inp)
    chemcloudEngine.specialKeywords["spinmult"] = 3
    chemcloudEngine.specialKeywords["charge"] = 0
    mol = chemcloudEngine.loadStructureFromXYZ("qmengine/1.xyz")
    assert mol.multiplicity == 3
    assert mol.charge == 0


# this test checks if spinmult is provided in the tc.in files
# which is not necessary anymore
"""
def test_loadStructureFromXYZ2():
    os.chdir(os.path.dirname(__file__))
    inp = getDefaults()
    inp.tctemplate = "qmengine/tc_wrong.in"
    inp.tctemplate_backup = "qmengine/tc_backup.in"
    inp.sampledir = Path("")
    passTest = True
    try:
        qmengine.ChemcloudEngine(inp)
        passTest = False
    except:
        pass
    assert passTest
"""


def test_loadStructureFromXYZ3():
    os.chdir(os.path.dirname(__file__))
    inp = getDefaults()
    inp.tctemplate = "qmengine/tc.in"
    inp.tctemplate_backup = "qmengine/tc_backup.in"
    inp.sampledir = Path("")
    chemcloudEngine = qmengine.ChemcloudEngine(inp)
    chemcloudEngine.specialKeywords = {}
    mol = chemcloudEngine.loadStructureFromXYZ("qmengine/1.xyz")
    assert mol.multiplicity == 1
    assert mol.charge == 0


def test_loadStructureFromXYZ4():
    os.chdir(os.path.dirname(__file__))
    inp = getDefaults()
    inp.tctemplate = "qmengine/tc.in"
    inp.tctemplate_backup = "qmengine/tc_backup.in"
    inp.sampledir = Path("")
    chemcloudEngine = qmengine.ChemcloudEngine(inp)
    chemcloudEngine.specialKeywords = {}
    chemcloudEngine.specialKeywords["spinmult"] = 3
    mol = chemcloudEngine.loadStructureFromXYZ("qmengine/1.xyz")
    assert mol.multiplicity == 3
    assert mol.charge == 0


def test_initResp():
    os.chdir(os.path.dirname(__file__))
    inp = getDefaults()
    inp.tctemplate = "qmengine/tc.in"
    inp.tctemplate_backup = "qmengine/tc_backup.in"
    inp.resp = 1.0
    inp.sampledir = Path("")
    chemcloudEngine = qmengine.ChemcloudEngine(inp)
    assert chemcloudEngine.keywords["resp"] == "yes"
    assert chemcloudEngine.backupKeywords["resp"] == "yes"


def test_getQMRefData(monkeypatch):
    os.chdir(os.path.dirname(__file__))
    inp = getDefaults()
    inp.tctemplate = "qmengine/tc.in"
    inp.tctemplate_backup = "qmengine/tc_backup.in"
    inp.sampledir = Path("")
    chemcloudEngine = qmengine.ChemcloudEngine(inp)
    calcDir = Path("chemcloudengine")
    xyzs = [f"{i}.xyz" for i in range(1, 11)]

    def monkeyRunJobs(self, xyzs, useBackup=False):
        refXyzs = ["3.xyz", "6.xyz", "9.xyz"]
        if "1.xyz" in xyzs:
            assert len(xyzs) == 10
            assert not useBackup
            return refXyzs
        else:
            assert checkUtils.checkLists(refXyzs, xyzs)
            assert useBackup

    def monkeyRead(*args):
        return [], [], [], [], []

    def monkeyWrite(*args):
        pass

    monkeypatch.setattr(qmengine.ChemcloudEngine, "runJobs", monkeyRunJobs)
    monkeypatch.setattr(qmengine.QMEngine, "readQMRefData", monkeyRead)
    monkeypatch.setattr(qmengine.QMEngine, "writeFBdata", monkeyWrite)
    os.chdir(calcDir)
    chemcloudEngine.getQMRefData(xyzs)


def test_createProgramInputsResp():
    os.chdir(os.path.dirname(__file__))
    inp = getDefaults()
    inp.tctemplate = "qmengine/tc.in"
    inp.tctemplate_backup = "qmengine/tc_backup.in"
    inp.sampledir = Path("")
    inp.resp = 1.0
    chemcloudEngine = qmengine.ChemcloudEngine(inp)
    xyzs = ["qmengine/test.xyz"]
    chemcloudEngine.createProgramInputs(xyzs)
    assert chemcloudEngine.doResp
    resp = ["resp", "yes"]
    assert resp in chemcloudEngine.inputSettings
    assert resp in chemcloudEngine.backupInputSettings


def test_writeResultResp(monkeypatch):
    os.chdir(os.path.dirname(__file__))
    inp = getDefaults()
    inp.tctemplate = "qmengine/tc.in"
    inp.tctemplate_backup = "qmengine/tc_backup.in"
    inp.sampledir = Path("")
    inp.resp = 1.90
    chemcloudEngine = qmengine.ChemcloudEngine(inp)
    prov = Provenance(program="terachem")
    mod = {"method": "hf", "basis": "sto-3g"}
    mol = Structure.open("qmengine/test/1.xyz")
    sp = ProgramInput(model=mod, structure=mol, calctype="energy", extras={"id": 999})
    res = parse("qmengine/test/tc_1.out", "terachem")
    with open("qmengine/test/tc_1.out", "r") as f:
        stdoutLines = f.readlines()
    stdout = "".join(stdoutLines)
    res.add_file("qmengine/test/esp.xyz")
    out = ProgramOutput(
        input_data=sp,
        results=res,
        provenance=prov,
        # files=files,
        stdout=stdout,
        success=True,
    )
    print(res.files.keys())
    chemcloudEngine.writeResult(out)
    wroteTcout = Path("tc_999.out").exists()
    wroteEsp = Path("esp_999.xyz").exists()
    os.remove("tc_999.out")
    os.remove("esp_999.xyz")
    assert wroteTcout
    assert wroteEsp


def test_dumpFailedJobs(monkeypatch):
    os.chdir(os.path.dirname(__file__))
    inp = getDefaults()
    inp.tctemplate = "qmengine/tc.in"
    inp.tctemplate_backup = "qmengine/tc_backup.in"
    inp.sampledir = Path("")
    chemcloudEngine = qmengine.ChemcloudEngine(inp)

    os.chdir("chemcloudengine")
    xyzs = ["3.xyz", "6.xyz"]
    inputs = chemcloudEngine.createProgramInputs(xyzs)
    outputs = []
    prov = Provenance(program="terachem")
    mod = {"method": "hf", "basis": "sto-3g"}
    mol = Structure.open("3.xyz")
    res = parse("tc_1.out", "terachem")
    # res.Files = {}
    # res.OptimizationResults = {}
    sp = ProgramInput(model=mod, structure=mol, calctype="energy", extras={"id": 3})
    out = ProgramOutput(
        input_data=sp,
        results=res,
        provenance=prov,
        success=False,
        traceback="Oops",
    )
    outputs.append(out)
    sp = ProgramInput(model=mod, structure=mol, calctype="energy", extras={"id": 6})
    out = ProgramOutput(
        input_data=sp,
        results=res,
        provenance=prov,
        success=False,
        traceback="Oops",
    )
    outputs.append(out)
    qmengine.dumpFailedJobs(inputs, outputs)
    pass3Input = checkUtils.checkFiles("3_input.yaml", "ref_3_input.yaml")
    pass3Output = checkUtils.checkFiles("3_output.yaml", "ref_3_output.yaml")
    pass6Input = checkUtils.checkFiles("6_input.yaml", "ref_6_input.yaml")
    pass6Output = checkUtils.checkFiles("6_output.yaml", "ref_6_output.yaml")
    os.remove("3_input.yaml")
    os.remove("3_output.yaml")
    os.remove("6_input.yaml")
    os.remove("6_output.yaml")
    assert pass3Input
    assert pass3Output
    assert pass6Input
    assert pass6Output


def monkeyWrite(self, output):
    pass


def test_getFailedJobs(monkeypatch):
    os.chdir(os.path.dirname(__file__))
    inp = getDefaults()
    inp.tctemplate = "qmengine/tc.in"
    inp.tctemplate_backup = "qmengine/tc_backup.in"
    inp.sampledir = Path("")
    ccEngine = qmengine.ChemcloudEngine(inp)
    os.chdir("chemcloudengine")

    monkeypatch.setattr(qmengine.ChemcloudEngine, "writeResult", monkeyWrite)
    inputs = []
    for i in range(1, 11):
        if i % 3 == 0:
            inputs.append(ccm.inputFromTCout(f"tc_{i}_success.out", extras={"id": i}))
        else:
            inputs.append(ccm.inputFromTCout(f"tc_{i}.out", extras={"id": i}))
    outputs = [
        ccm.programOutputFromTCout(f"tc_{i}.out", inp=inputs[i - 1])
        for i in range(1, 11)
    ]
    retryXyzs = ccEngine.getFailedJobs(outputs)
    assert retryXyzs == ["3.xyz", "6.xyz", "9.xyz"]


# haha
def monkeyDump(inputs, outputs):
    pass


def test_runJobs(monkeypatch):
    os.chdir(os.path.dirname(__file__))
    inp = getDefaults()
    inp.tctemplate = "qmengine/tc.in"
    inp.tctemplate_backup = "qmengine/tc_backup.in"
    inp.sampledir = Path("")
    monkeypatch.setattr(qmengine, "dumpFailedJobs", monkeyDump)
    monkeypatch.setattr(qmengine.ChemcloudEngine, "writeResult", monkeyWrite)
    ccEngine = qmengine.ChemcloudEngine(inp)
    os.chdir("chemcloudengine")

    xyzs = [f"{i}.xyz" for i in range(1, 11)]
    inputs = []
    for i in range(1, 11):
        if i % 3 == 0:
            inputs.append(ccm.inputFromTCout(f"tc_{i}_success.out", extras={"id": i}))
        else:
            inputs.append(ccm.inputFromTCout(f"tc_{i}.out", extras={"id": i}))
    outputs = [
        ccm.programOutputFromTCout(f"tc_{i}.out", inp=inputs[i - 1])
        for i in range(1, 11)
    ]
    ccm.patcher(monkeypatch, qmengine, outputs)
    retryXyzs = ccEngine.runJobs(xyzs)
    assert retryXyzs == ["3.xyz", "6.xyz", "9.xyz"]


def test_runJobsFail(monkeypatch):
    os.chdir(os.path.dirname(__file__))
    inp = getDefaults()
    inp.tctemplate = "qmengine/tc.in"
    inp.tctemplate_backup = "qmengine/tc_backup.in"
    inp.sampledir = Path("")
    monkeypatch.setattr(qmengine, "dumpFailedJobs", monkeyDump)
    monkeypatch.setattr(qmengine.ChemcloudEngine, "writeResult", monkeyWrite)
    ccEngine = qmengine.ChemcloudEngine(inp)
    os.chdir("chemcloudengine")

    xyzs = [f"{i}.xyz" for i in range(1, 11)]
    inputs = []
    for i in range(1, 11):
        if i % 3 == 0:
            inputs.append(ccm.inputFromTCout(f"tc_{i}_success.out", extras={"id": i}))
        else:
            inputs.append(ccm.inputFromTCout(f"tc_{i}.out", extras={"id": i}))
    outputs = [
        ccm.programOutputFromTCout(f"tc_{i}.out", inp=inputs[i - 1])
        for i in range(1, 8)
    ]
    ccm.patcher(monkeypatch, qmengine, outputs)
    crash = ""
    try:
        ccEngine.runJobs(xyzs)
    except Exception as e:
        crash = e
    assert str(crash) == "ChemCloud did not return the same number of outputs as inputs"


def test_runJobsBackup(monkeypatch):
    os.chdir(os.path.dirname(__file__))
    inp = getDefaults()
    inp.tctemplate = "qmengine/tc.in"
    inp.tctemplate_backup = "qmengine/tc_backup.in"
    inp.sampledir = Path("")

    def monkeyDumpKeywords(inputs, outputs):
        with open("temp.txt", "w") as f:
            for key, val in inputs[0].keywords.items():
                f.write(f"{key} {val}\n")

    monkeypatch.setattr(qmengine, "dumpFailedJobs", monkeyDumpKeywords)
    monkeypatch.setattr(qmengine.ChemcloudEngine, "writeResult", monkeyWrite)
    ccEngine = qmengine.ChemcloudEngine(inp)
    os.chdir("chemcloudengine")
    xyzs = [f"{i}.xyz" for i in range(1, 11)]
    inputs = []
    for i in range(1, 11):
        if i % 3 == 0:
            inputs.append(ccm.inputFromTCout(f"tc_{i}_success.out", extras={"id": i}))
        else:
            inputs.append(ccm.inputFromTCout(f"tc_{i}.out", extras={"id": i}))
    outputs = [
        ccm.programOutputFromTCout(f"tc_{i}.out", inp=inputs[i - 1])
        for i in range(1, 11)
    ]
    ccm.patcher(monkeypatch, qmengine, outputs)
    crash = ""
    try:
        retryXyzs = ccEngine.runJobs(xyzs, useBackup=True)
    except Exception as e:
        crash = str(e)
        print(crash)
    keywords = {}
    with open("temp.txt", "r") as f:
        for line in f.readlines():
            split = line.split()
            keywords[split[0]] = split[1]
    os.remove("temp.txt")
    assert crash.startswith("Job ids ['3', '6', '9']")
    assert keywords["threall"] == "1.0e-14"


def test_runJobs1(monkeypatch):
    os.chdir(os.path.dirname(__file__))
    inp = getDefaults()
    inp.tctemplate = "qmengine/tc.in"
    inp.tctemplate_backup = "qmengine/tc_backup.in"
    inp.sampledir = Path("")
    monkeypatch.setattr(qmengine, "dumpFailedJobs", monkeyDump)
    monkeypatch.setattr(qmengine.ChemcloudEngine, "writeResult", monkeyWrite)
    ccEngine = qmengine.ChemcloudEngine(inp)
    os.chdir("chemcloudengine")

    outputs = [ccm.programOutputFromTCout("tc_1.out")]
    ccm.patcher(monkeypatch, qmengine, outputs)

    def monkeyGetFailedJobs(self, outputs):
        assert isinstance(outputs, list)
        assert len(outputs) == 1
        return []

    monkeypatch.setattr(qmengine.ChemcloudEngine, "getFailedJobs", monkeyGetFailedJobs)
    ccEngine.runJobs(["1.xyz"])
