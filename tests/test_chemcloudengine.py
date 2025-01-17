import os
from pathlib import Path

from qcio import ProgramInput, ProgramOutput, Provenance, Structure
from qcparse import parse

from ff_optimizer import qmengine

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
    xyzs = []
    outputs = []
    prov = Provenance(program="terachem")
    mod = {"method": "hf", "basis": "sto-3g"}
    for i in range(1, 11):
        xyzs.append(f"{i}.xyz")
        mol = Structure.open(calcDir / f"{i}.xyz")
        sp = ProgramInput(model=mod, structure=mol, calctype="energy", extras={"id": i})
        try:
            res = parse(calcDir / f"tc_{i}.out", "terachem")
            out = ProgramOutput(
                input_data=sp, results=res, provenance=prov, success=True
            )
            # out = SinglePointOutput(input_data=sp, results=res, provenance=prov)
        except Exception:
            out = ProgramOutput(
                input_data=sp,
                results=res,
                provenance=prov,
                success=False,
                traceback="Oops",
            )
            # out = ProgramFailure(input_data=sp, results={}, provenance=prov)
        outputs.append(out)

    refIds = ["3", "6", "9"]
    retryOutputs = []
    for i in refIds:
        mol = Structure.open(calcDir / f"{i}.xyz")
        sp = ProgramInput(model=mod, structure=mol, calctype="energy", extras={"id": i})
        res = parse(calcDir / f"tc_{i}_success.out", "terachem")
        out = ProgramOutput(input_data=sp, results=res, provenance=prov, success=True)
        retryOutputs.append(out)

    def monkeyComputeBatch(programInputs):
        if len(programInputs) == 10:
            return 0, outputs
        else:
            print(len(programInputs))
            ids = sorted([inp.extras["id"] for inp in programInputs])
            assert len(programInputs) == 3
            assert programInputs[0].keywords["threall"] == "1.0e-14"
            assert programInputs[0].keywords["diismaxvecs"] == "40"
            for i in range(len(ids)):
                assert ids[i] == refIds[i]
            return 0, retryOutputs

    def monkeyRead(*args):
        return [], [], [], [], []

    def monkeyWrite(*args):
        pass

    monkeypatch.setattr(chemcloudEngine, "computeBatch", monkeyComputeBatch)
    monkeypatch.setattr(chemcloudEngine, "writeResult", monkeyWrite)
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


def test_computeBatch(monkeypatch):
    os.chdir(os.path.dirname(__file__))
    inp = getDefaults()
    inp.tctemplate = "qmengine/tc.in"
    inp.tctemplate_backup = "qmengine/tc_backup.in"
    inp.sampledir = Path("")
    chemcloudEngine = qmengine.ChemcloudEngine(inp)

    class MonkeyFutureResult:
        def __init__(self, programInputs):
            self.ids = []
            self.id = 123
            for programInput in programInputs:
                self.ids.append(programInput.extras["id"])

        def get(self):
            results = []
            for id in self.ids:
                try:
                    result = parse(f"tc_{id}.out", "terachem")
                except:
                    result = None
                results.append(result)
            return results

    def monkeyCompute(engine, programInputs, collect_files=True):
        return MonkeyFutureResult(programInputs)

    os.chdir("chemcloudengine")
    xyzs = []
    for f in os.listdir():
        if f.endswith("xyz"):
            xyzs.append(f)
    inputs = chemcloudEngine.createProgramInputs(xyzs)
    monkeypatch.setattr(chemcloudEngine.client, "compute", monkeyCompute)

    # Check default batch size
    status, results = chemcloudEngine.computeBatch(inputs)
    if os.path.isfile(os.path.join("chemcloudengine", "jobs.txt")):
        os.remove(os.path.join("chemcloudengine", "jobs.txt"))
    assert status == 0
    assert len(results) == 10

    # Check small batch size
    chemcloudEngine.batchsize = 3
    status, results = chemcloudEngine.computeBatch(inputs)
    if os.path.isfile(os.path.join("chemcloudengine", "jobs.txt")):
        os.remove(os.path.join("chemcloudengine", "jobs.txt"))
    assert status == 0
    assert len(results) == 10

    # Check large batch size
    chemcloudEngine.batchsize = 77
    status, results = chemcloudEngine.computeBatch(inputs)
    if os.path.isfile(os.path.join("chemcloudengine", "jobs.txt")):
        os.remove(os.path.join("chemcloudengine", "jobs.txt"))
    assert status == 0
    assert len(results) == 10

    # Check that batch resizing works
    class MonkeyResultFail:
        def __init__(self, programInputs):
            pass

        def get(self):
            raise RuntimeError("oops")

    def monkeyComputeFailOnce(engine, programInputs, collect_files=True):
        if len(programInputs) > 5:
            return MonkeyResultFail(programInputs)
        else:
            return MonkeyFutureResult(programInputs)

    monkeypatch.setattr(chemcloudEngine.client, "compute", monkeyComputeFailOnce)
    chemcloudEngine.batchSize = 6
    status, results = chemcloudEngine.computeBatch(inputs)
    if os.path.isfile(os.path.join("chemcloudengine", "jobs.txt")):
        os.remove(os.path.join("chemcloudengine", "jobs.txt"))
    assert status == 0
    assert len(results) == 10
    assert chemcloudEngine.batchSize == 3

    chemcloudEngine.batchSize = 12
    status, results = chemcloudEngine.computeBatch(inputs)
    if os.path.isfile(os.path.join("chemcloudengine", "jobs.txt")):
        os.remove(os.path.join("chemcloudengine", "jobs.txt"))
    assert status == 0
    assert len(results) == 10
    assert chemcloudEngine.batchSize == 5

    # Check failure
    def monkeyComputeFail(engine, programInputs, collect_files=True):
        return MonkeyResultFail(programInputs)

    monkeypatch.setattr(chemcloudEngine.client, "compute", monkeyComputeFail)
    status, results = chemcloudEngine.computeBatch(inputs)
    if os.path.isfile(os.path.join("chemcloudengine", "jobs.txt")):
        os.remove(os.path.join("chemcloudengine", "jobs.txt"))
    assert status == -1
    assert len(results) == 0

    os.chdir("..")

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
        #files=files,
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
