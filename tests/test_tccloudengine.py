import os

from qcelemental.util.serialization import json_loads
from qccloud.models import AtomicResult, FailedOperation

from ff_optimizer import qmengine


def test_init():
    os.chdir(os.path.dirname(__file__))
    qccloudEngine = qmengine.TCCloudEngine("qmengine/tc.in", "qmengine/tc_backup.in")
    assert qccloudEngine.method == "b3lyp"
    assert qccloudEngine.basis == "6-31gss"
    assert qccloudEngine.keywords["dftd"] == "d3"
    assert qccloudEngine.keywords["charge"] == "0"
    assert qccloudEngine.backupKeywords["threall"] == "1.0e-14"
    assert qccloudEngine.backupKeywords["diismaxvecs"] == "40"
    assert qccloudEngine.backupKeywords["maxit"] == "200"


def test_initResp():
    os.chdir(os.path.dirname(__file__))
    qccloudEngine = qmengine.TCCloudEngine(
        "qmengine/tc.in", "qmengine/tc_backup.in", doResp=True
    )
    assert qccloudEngine.keywords["resp"] == "yes"
    assert qccloudEngine.backupKeywords["resp"] == "yes"


def test_getQMRefData(monkeypatch):
    os.chdir(os.path.dirname(__file__))
    qccloudEngine = qmengine.TCCloudEngine("qmengine/tc.in", "qmengine/tc_backup.in")
    calcDir = "qccloudengine"
    pdbs = []
    jsons = []
    for f in os.listdir(calcDir):
        if f.endswith("pdb"):
            pdbs.append(f)
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

    monkeypatch.setattr(qccloudEngine, "computeBatch", monkeyCompute)
    monkeypatch.setattr(qccloudEngine, "writeResult", monkeyWrite)
    monkeypatch.setattr(qmengine.QMEngine, "readQMRefData", monkeyRead)
    monkeypatch.setattr(qmengine.QMEngine, "writeFBdata", monkeyWrite)
    qccloudEngine.getQMRefData(pdbs, calcDir)


def test_createAtomicInputsResp():
    os.chdir(os.path.dirname(__file__))
    qccloudEngine = qmengine.TCCloudEngine(
        "qmengine/tc.in", "qmengine/tc_backup.in", doResp=True
    )
    pdbs = ["qmengine/test.pdb"]
    atomicInputs = qccloudEngine.createAtomicInputs(pdbs)
    ainput = atomicInputs[0]
    assert ainput.protocols.native_files.all == "all"
    assert ainput.extras["tcfe:keywords"]["native_files"] == ["esp.xyz"]


def test_computeBatch(monkeypatch):
    os.chdir(os.path.dirname(__file__))
    qccloudEngine = qmengine.TCCloudEngine(os.path.join("qmengine","tc.in"), os.path.join("qmengine","tc_backup.in"))

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

    os.chdir("qccloudengine")
    pdbs = []
    for f in os.listdir():
        if f.endswith("pdb"):
            pdbs.append(f)
    inputs = qccloudEngine.createAtomicInputs(pdbs)
    print(inputs)
    monkeypatch.setattr(qccloudEngine.client, "compute", monkeyCompute)

    # Check default batch size
    status, results = qccloudEngine.computeBatch(inputs)
    if os.path.isfile(os.path.join("qccloudengine","jobs.txt")):
        os.remove(os.path.join("qccloudengine","jobs.txt"))
    assert status == 0
    assert len(results) == 10

    # Check small batch size
    qccloudEngine.batchsize = 3
    status, results = qccloudEngine.computeBatch(inputs)
    if os.path.isfile(os.path.join("qccloudengine","jobs.txt")):
        os.remove(os.path.join("qccloudengine","jobs.txt"))
    assert status == 0
    assert len(results) == 10

    # Check large batch size
    qccloudEngine.batchsize = 77
    status, results = qccloudEngine.computeBatch(inputs)
    if os.path.isfile(os.path.join("qccloudengine","jobs.txt")):
        os.remove(os.path.join("qccloudengine","jobs.txt"))
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

    monkeypatch.setattr(qccloudEngine.client, "compute", monkeyComputeFailOnce)
    qccloudEngine.batchSize = 6
    status, results = qccloudEngine.computeBatch(inputs)
    if os.path.isfile(os.path.join("qccloudengine","jobs.txt")):
        os.remove(os.path.join("qccloudengine","jobs.txt"))
    assert status == 0
    assert len(results) == 10
    assert qccloudEngine.batchSize == 3

    qccloudEngine.batchSize = 12
    status, results = qccloudEngine.computeBatch(inputs)
    if os.path.isfile(os.path.join("qccloudengine","jobs.txt")):
        os.remove(os.path.join("qccloudengine","jobs.txt"))
    assert status == 0
    assert len(results) == 10
    assert qccloudEngine.batchSize == 5

    # Check failure
    def monkeyComputeFail(atomicInputs, engine):
        return MonkeyResultFail(atomicInputs)

    monkeypatch.setattr(qccloudEngine.client, "compute", monkeyComputeFail)
    status, results = qccloudEngine.computeBatch(inputs)
    if os.path.isfile(os.path.join("qccloudengine","jobs.txt")):
        os.remove(os.path.join("qccloudengine","jobs.txt"))
    assert status == -1
    assert len(results) == 0

    os.chdir("..")
