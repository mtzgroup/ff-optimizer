from ff_optimizer import qmengine
import numpy as np
from . import checkUtils
import os
import sys

class TestSbatchEngine:
    
    # check that sbatch files are being read in correctly
    #def test_readSbatchFile(self):  
    #    os.chdir(os.path.dirname(__file__))
    #    sbatchEngine = qmengine.SbatchEngine("qmengine/tc.in","qmengine/tc_backup.in","qmengine/sbatch.sh","test")
    #    for line in sbatchEngine.sbatchLines:
    #        print(line)
    #    assert sbatchEngine.tcVersion == "TeraChem/2021.02-intel-2017.8.262-CUDA-9.0.176"
    #    assert sbatchEngine.sbatchOptions == ['#SBATCH -t 2:00:00\n','#SBATCH --gres gpu:2\n','#SBATCH --mem=16GB\n','#SBATCH --exclude=fire-11-01\n']
        
    # check that sbatch files are being written correctly
    def test_writeSbatchFile(self):
        os.chdir(os.path.dirname(__file__))
        sbatchEngine = qmengine.SbatchEngine("qmengine/tc.in","qmengine/tc_backup.in","qmengine/sbatch.sh","test")
        sbatchEngine.writeSbatchFile('23',"sbatch_23.sh")
        refLines = []
        with open("qmengine/sbatch_23.sh",'r') as refF:
            for line in refF.readlines():
                refLines.append(line)
        testLines = []
        with open("sbatch_23.sh",'r') as testF:
            for line in testF.readlines():
                testLines.append(line)
        assert len(refLines) == len(testLines)
        for i in range(len(refLines)):
            assert refLines[i] == testLines[i]
        os.remove("sbatch_23.sh")

    def test_writeSbatchFileResp(self):
        os.chdir(os.path.dirname(__file__))
        sbatchEngine = qmengine.SbatchEngine("qmengine/tc.in","qmengine/tc_backup.in","qmengine/sbatch.sh","test",doResp=True)
        sbatchEngine.writeSbatchFile('23',"sbatch_23.sh")
        refLines = []
        with open("qmengine/sbatch_23_resp.sh",'r') as refF:
            for line in refF.readlines():
                refLines.append(line)
        testLines = []
        with open("sbatch_23.sh",'r') as testF:
            for line in testF.readlines():
                testLines.append(line)
        assert len(refLines) == len(testLines)
        for i in range(len(refLines)):
            assert refLines[i] == testLines[i]
        os.remove("sbatch_23.sh")
    # check that job submission is being done correctly
    #def test_getQMRefData(self, monkeypatch):
    #    os.chdir(os.path.dirname(__file__))
    #    sbatchEngine = sbatchEngine.QMEngine("qmengine/tc.in","qmengine/tc_backup.in","qmengine/sbatch.sh")

