from ff_optimizer import qmengine
import numpy as np
from . import utils
import os
import sys

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

    # check that grads and energies are being read in correctly
    def test_readGrads(self):
        os.chdir(os.path.dirname(__file__))
        qmEngine = qmengine.QMEngine("qmengine/tc.in","qmengine/tc_backup.in")
        testEnergy, testGrads = qmEngine.readGradFromTCout("qmengine/test.out")
        energy = np.loadtxt("qmengine/energy.txt")
        grads = np.loadtxt("qmengine/grads.txt").flatten()
        assert utils.checkFloat(energy, testEnergy)
        assert utils.checkArray(grads, testGrads)

    # check that pdbs are read in correctly
    def test_readPDB(self):
        os.chdir(os.path.dirname(__file__))
        qmEngine = qmengine.QMEngine("qmengine/tc.in","qmengine/tc_backup.in")
        testCoords = qmEngine.readPDB("qmengine/test.pdb")
        coords = np.loadtxt("qmengine/coords.txt").flatten()
        assert utils.checkArray(coords, testCoords)
    
    # check that FB data is written correctly
    def test_writeFBdata(self):
        os.chdir(os.path.dirname(__file__))
        qmEngine = qmengine.QMEngine("qmengine/tc.in","qmengine/tc_backup.in")

        coords = []
        energies = []
        grads = []
        for i in range(1,26):
            coords.append(qmEngine.readPDB(f"qmengine/test/{str(i)}.pdb"))
            energy, grad = qmEngine.readGradFromTCout(f"qmengine/test/tc_{str(i)}.out")
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
            assert refLines[i] == testLines[i]
        os.remove("qdata.txt")
        os.remove("all.mdcrd")

