from ff_optimizer import qmengine, utils
import numpy as np
from . import checkUtils
import os
import sys
from tccloud.models import AtomicResult

class MonkeyProperties():
    def __init__(self, energy):
        self.return_energy = energy

class MonkeyResult():
    
    def __init__(self, id, energy, grad):
        self.id = id
        self.return_result = grad
        self.properties = MonkeyProperties(energy)
        self.success = True

class TestQMEngine:
        
    # check that input settings are read in correctly
    def test_readInput(self):
        os.chdir(os.path.dirname(__file__))
        qmEngine = qmengine.QMEngine("qmengine/tc.in","qmengine/tc_backup.in")
        inputSettings = [['basis','6-31gss'],['method','b3lyp'],['charge','0'],['dftd','d3']]
        backupInputSettings = [['basis','6-31gss'],['method','b3lyp'],['charge','0'],['dftd','d3'],['threall','1.0e-14'],['diismaxvecs','40'],['maxit','200']] 
        assert qmEngine.inputSettings == inputSettings
        assert qmEngine.backupInputSettings == backupInputSettings
        
    # check that settings are written out correctly
    def test_writeInputFile(self):
        os.chdir(os.path.dirname(__file__))
        qmEngine = qmengine.QMEngine("qmengine/tc.in","qmengine/tc_backup.in")
        inputSettings = [['basis','6-31gss'],['method','b3lyp'],['charge','0'],['dftd','d3']]
        testPath = "qmengine/test.in"
        qmEngine.writeInputFile(inputSettings, "coors.xyz", testPath)
        with open(testPath,'r') as f:
            assert "coordinates coors.xyz" in f.readline()
            assert "run gradient" in f.readline()
        readSettings = qmEngine.readInputFile(testPath)
        os.remove(testPath)
        assert readSettings == inputSettings

    # check that pdbs are read in correctly
    def test_readPDB(self):
        os.chdir(os.path.dirname(__file__))
        qmEngine = qmengine.QMEngine("qmengine/tc.in","qmengine/tc_backup.in")
        testCoords = utils.readPDB("qmengine/test.pdb")
        coords = np.loadtxt("qmengine/coords.txt").flatten()
        assert checkUtils.checkArray(coords, testCoords)
    
    # check that FB data is written correctly
    def test_writeFBdata(self):
        os.chdir(os.path.dirname(__file__))
        qmEngine = qmengine.QMEngine("qmengine/tc.in","qmengine/tc_backup.in")

        coords = []
        energies = []
        grads = []
        for i in range(1,26):
            coords.append(utils.readPDB(f"qmengine/test/{str(i)}.pdb"))
            energy, grad = utils.readGradFromTCout(f"qmengine/test/tc_{str(i)}.out")
            energies.append(energy)
            grads.append(grad)
        qmEngine.writeFBdata(energies,grads,coords)
        refLines = []
        with open("qmengine/test/all.mdcrd",'r') as refF:
            for line in refF.readlines():
                if len(line.split()) > 0:
                    refLines.append(line)
        testLines = []
        with open("all.mdcrd",'r') as testF:
            for line in testF.readlines():
                if len(line.split()) > 0:
                    testLines.append(line)
        assert len(refLines) == len(testLines)
        for i in range(1,len(refLines)):
            assert refLines[i] == testLines[i]
        refLines = []
        with open("qmengine/test/qdata.txt",'r') as refF:
            for line in refF.readlines():
                if len(line.split()) > 0:
                    refLines.append(line)
        testLines = []
        with open("qdata.txt",'r') as testF:
            for line in testF.readlines():
                if len(line.split()) > 0:
                    testLines.append(line)
        assert len(refLines) == len(testLines)
        for i in range(len(refLines)):
            refLine = refLines[i].split()
            testLine = testLines[i].split()
            assert len(refLine) == len(testLine)
            for j in range(len(refLine)):
                try:
                    assert checkUtils.checkFloat(refLine[j],testLine[j],0.0001)
                except:
                    assert refLine[j] == testLine[j]
        os.remove("qdata.txt")
        os.remove("all.mdcrd")

    def test_writeResult(self):
        os.chdir(os.path.dirname(__file__))
        tccloudEngine = qmengine.TCCloudEngine("qmengine/tc.in","qmengine/tc_backup.in")
        energy = np.loadtxt(os.path.join("qmengine","energy.txt"))
        grads = np.loadtxt(os.path.join("qmengine","grads.txt"))
        result = MonkeyResult(10,energy,grads)
        tccloudEngine.writeResult(result)
        refLines = []
        testLines = []
        with open(os.path.join("qmengine","result.txt"),'r') as refF:
            for line in refF.readlines():
                refLines.append(line)
        with open("tc_10.out",'r') as testF:
            for line in testF.readlines():
                testLines.append(line)
        assert len(refLines) == len(testLines)
        for i in range(len(refLines)):
            refLine = refLines[i].split()
            testLine = testLines[i].split()
            assert len(refLine) == len(testLine)
            for j in range(len(refLine)):
                try:
                    check = checkUtils.checkFloat(refLine[j],testLine[j])
                except:
                    check = (refLine[j] == testLine[j])
                assert check
        #os.remove("tc_10.out")

