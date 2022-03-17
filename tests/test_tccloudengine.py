from ff_optimizer import qmengine
import os
from qcelemental.util.serialization import json_loads
from tccloud.models import AtomicResult, FailedOperation
import sys

def test_init():
    os.chdir(os.path.dirname(__file__))
    tccloudEngine = qmengine.TCCloudEngine("qmengine/tc.in","qmengine/tc_backup.in")
    assert tccloudEngine.method == "b3lyp"
    assert tccloudEngine.basis == "6-31gss"
    assert tccloudEngine.keywords['dftd'] == 'd3'
    assert tccloudEngine.keywords['charge'] == '0' 
    assert tccloudEngine.backupKeywords['threall'] == '1.0e-14'
    assert tccloudEngine.backupKeywords['diismaxvecs'] == '40'
    assert tccloudEngine.backupKeywords['maxit'] == '200'

def testGetQMRefData(monkeypatch):
    os.chdir(os.path.dirname(__file__))
    tccloudEngine = qmengine.TCCloudEngine("qmengine/tc.in","qmengine/tc_backup.in")
    calcDir = "tccloudengine"
    pdbs = []
    jsons = []
    for f in os.listdir(calcDir):
        if f.endswith("pdb"):
            pdbs.append(f)
        if f.endswith("json"):
            jsons.append(f)
    results = []
    for json in jsons:
        with open(os.path.join(calcDir,json),'r') as f:
            data = f.read()
            try:
                results.append(AtomicResult(**json_loads(data)))
            except:
                results.append(FailedOperation(**json_loads(data)))
    refIds = ['3','6','9']
    retryResults = []
    for i in refIds:
        with open(os.path.join(calcDir,f"tc_{i}_success.json"),'r') as f:
            retryResults.append(AtomicResult(**json_loads(f.read())))

    def monkeyCompute(atomicInputs):
        if len(atomicInputs) == 10:
            return 0, results
        else:
            ids = sorted([input.id for input in atomicInputs])
            assert len(atomicInputs) == 3
            assert atomicInputs[0].keywords['threall'] == '1.0e-14'
            assert atomicInputs[0].keywords['diismaxvecs'] == '40'
            for i in range(len(ids)):
                assert ids[i] == refIds[i]
            return 0, retryResults

    def monkeyRead(*args):
        return [], [], [], [], []

    def monkeyWrite(*args):
        pass

    monkeypatch.setattr(tccloudEngine,"computeBatch",monkeyCompute)
    monkeypatch.setattr(tccloudEngine,"writeResult",monkeyWrite)
    monkeypatch.setattr(qmengine.QMEngine,"readQMRefData",monkeyRead)
    monkeypatch.setattr(qmengine.QMEngine,"writeFBdata",monkeyWrite)
    tccloudEngine.getQMRefData(pdbs,calcDir)
