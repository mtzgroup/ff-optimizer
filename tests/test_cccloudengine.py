import os

from chemcloud.models import AtomicResult, FailedOperation
from qcelemental.util.serialization import json_loads

from ff_optimizer import qmengine


def test_init():
    os.chdir(os.path.dirname(__file__))
    chemcloudEngine = qmengine.CCCloudEngine("qmengine/tc.in", "qmengine/tc_backup.in")
    assert chemcloudEngine.specialKeywords["method"] == "b3lyp"
    assert chemcloudEngine.specialKeywords["basis"] == "6-31gss"
    assert chemcloudEngine.keywords["dftd"] == "d3"
    assert chemcloudEngine.specialKeywords["charge"] == 0
    assert chemcloudEngine.specialKeywords["spinmult"] == 1
    assert chemcloudEngine.backupKeywords["threall"] == "1.0e-14"
    assert chemcloudEngine.backupKeywords["diismaxvecs"] == "40"
    assert chemcloudEngine.backupKeywords["maxit"] == "200"


def test_loadMoleculeFromXYZ1():
    os.chdir(os.path.dirname(__file__))
    ccEngine = qmengine.CCCloudEngine("qmengine/tc.in", "qmengine/tc_backup.in")
    ccEngine.specialKeywords["spinmult"] = 3
    ccEngine.specialKeywords["charge"] = 0
    mol = ccEngine.loadMoleculeFromXYZ("qmengine/1.xyz")
    assert mol.molecular_multiplicity == 3
    assert mol.molecular_charge == 0


def test_loadMoleculeFromXYZ2():
    os.chdir(os.path.dirname(__file__))
    ccEngine = qmengine.CCCloudEngine("qmengine/tc.in", "qmengine/tc_backup.in")
    ccEngine.specialKeywords = {}
    ccEngine.specialKeywords["charge"] = 1
    mol = ccEngine.loadMoleculeFromXYZ("qmengine/1.xyz")
    assert mol.molecular_multiplicity == 2
    assert mol.molecular_charge == 1


def test_loadMoleculeFromXYZ3():
    os.chdir(os.path.dirname(__file__))
    ccEngine = qmengine.CCCloudEngine("qmengine/tc.in", "qmengine/tc_backup.in")
    ccEngine.specialKeywords = {}
    mol = ccEngine.loadMoleculeFromXYZ("qmengine/1.xyz")
    assert mol.molecular_multiplicity == 1
    assert mol.molecular_charge == 0


def test_loadMoleculeFromXYZ4():
    os.chdir(os.path.dirname(__file__))
    ccEngine = qmengine.CCCloudEngine("qmengine/tc.in", "qmengine/tc_backup.in")
    ccEngine.specialKeywords = {}
    ccEngine.specialKeywords["spinmult"] = 3
    mol = ccEngine.loadMoleculeFromXYZ("qmengine/1.xyz")
    assert mol.molecular_multiplicity == 3
    assert mol.molecular_charge == 0


def test_initResp():
    os.chdir(os.path.dirname(__file__))
    chemcloudEngine = qmengine.CCCloudEngine(
        "qmengine/tc.in", "qmengine/tc_backup.in", doResp=True
    )
    assert chemcloudEngine.keywords["resp"] == "yes"
    assert chemcloudEngine.backupKeywords["resp"] == "yes"


def test_getQMRefData(monkeypatch):
    os.chdir(os.path.dirname(__file__))
    chemcloudEngine = qmengine.CCCloudEngine("qmengine/tc.in", "qmengine/tc_backup.in")
    calcDir = "chemcloudengine"
    xyzs = []
    jsons = []
    for f in os.listdir(calcDir):
        if f.endswith("xyz"):
            xyzs.append(f)
        if f.endswith("json") and "success" not in f:
            jsons.append(f)
    results = []
    for json in jsons:
        with open(os.path.join(calcDir, json), "r") as f:
            data = f.read()
            try:
                results.append(AtomicResult(**json_loads(data)))
            except:
                results.append(FailedOperation(**json_loads(data)))
    refIds = ["3", "6", "9"]
    retryResults = []
    for i in refIds:
        with open(os.path.join(calcDir, f"tc_{i}_success.json"), "r") as f:
            retryResults.append(AtomicResult(**json_loads(f.read())))

    def monkeyCompute(atomicInputs):
        print(len(atomicInputs))
        if len(atomicInputs) == 10:
            return 0, results
        else:
            ids = sorted([input.id for input in atomicInputs])
            assert len(atomicInputs) == 3
            assert atomicInputs[0].keywords["threall"] == "1.0e-14"
            assert atomicInputs[0].keywords["diismaxvecs"] == "40"
            for i in range(len(ids)):
                assert ids[i] == refIds[i]
            print(retryResults)
            return 0, retryResults

    def monkeyRead(*args):
        return [], [], [], [], []

    def monkeyWrite(*args):
        pass

    monkeypatch.setattr(chemcloudEngine, "computeBatch", monkeyCompute)
    monkeypatch.setattr(chemcloudEngine, "writeResult", monkeyWrite)
    monkeypatch.setattr(qmengine.QMEngine, "readQMRefData", monkeyRead)
    monkeypatch.setattr(qmengine.QMEngine, "writeFBdata", monkeyWrite)
    os.chdir(calcDir)
    chemcloudEngine.getQMRefData(xyzs)


def test_createAtomicInputsResp():
    os.chdir(os.path.dirname(__file__))
    chemcloudEngine = qmengine.CCCloudEngine(
        "qmengine/tc.in", "qmengine/tc_backup.in", doResp=True
    )
    xyzs = ["qmengine/test.xyz"]
    atomicInputs = chemcloudEngine.createAtomicInputs(xyzs)
    ainput = atomicInputs[0]
    assert ainput.protocols.native_files.all == "all"
    assert ainput.extras["tcfe:keywords"]["native_files"] == ["esp.xyz"]


def test_computeBatch(monkeypatch):
    os.chdir(os.path.dirname(__file__))
    chemcloudEngine = qmengine.CCCloudEngine(
        os.path.join("qmengine", "tc.in"), os.path.join("qmengine", "tc_backup.in")
    )

    class MonkeyFutureResult:
        def __init__(self, atomicInputs):
            self.ids = []
            self.id = 123
            for atomicInput in atomicInputs:
                self.ids.append(atomicInput.id)

        def get(self):
            results = []
            for id in self.ids:
                with open(f"tc_{str(id)}.json", "r") as f:
                    data = f.read()
                    try:
                        result = AtomicResult(**json_loads(data))
                    except:
                        result = FailedOperation(**json_loads(data))
                results.append(result)
            return results

    def monkeyCompute(atomicInputs, engine):
        return MonkeyFutureResult(atomicInputs)

    os.chdir("chemcloudengine")
    xyzs = []
    for f in os.listdir():
        if f.endswith("xyz"):
            xyzs.append(f)
    inputs = chemcloudEngine.createAtomicInputs(xyzs)
    print(inputs)
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
        def __init__(self, atomicInputs):
            pass

        def get(self):
            raise RuntimeError("oops")

    def monkeyComputeFailOnce(atomicInputs, engine):
        if len(atomicInputs) > 5:
            return MonkeyResultFail(atomicInputs)
        else:
            return MonkeyFutureResult(atomicInputs)

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
    def monkeyComputeFail(atomicInputs, engine):
        return MonkeyResultFail(atomicInputs)

    monkeypatch.setattr(chemcloudEngine.client, "compute", monkeyComputeFail)
    status, results = chemcloudEngine.computeBatch(inputs)
    if os.path.isfile(os.path.join("chemcloudengine", "jobs.txt")):
        os.remove(os.path.join("chemcloudengine", "jobs.txt"))
    assert status == -1
    assert len(results) == 0

    os.chdir("..")
