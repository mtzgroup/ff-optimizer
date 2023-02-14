import os
from pathlib import Path

import numpy as np
from chemcloud.models import AtomicResult
from qcelemental.util.serialization import json_loads

from ff_optimizer import qmengine, utils

from . import checkUtils


class MonkeyProperties:
    def __init__(self, energy):
        self.return_energy = energy


class MonkeyResult:
    def __init__(self, id, energy, grad):
        self.id = id
        self.return_result = grad
        self.properties = MonkeyProperties(energy)
        self.success = True

    def json(self):
        return "test"


class TestQMEngine:

    # check that input settings are read in correctly
    def test_readInput(self):
        os.chdir(os.path.dirname(__file__))
        qmEngine = qmengine.QMEngine("qmengine/tc.in", "qmengine/tc_backup.in")
        inputSettings = [
            ["basis", "6-31gss"],
            ["method", "b3lyp"],
            ["charge", "0"],
            ["dftd", "d3"],
        ]
        backupInputSettings = [
            ["basis", "6-31gss"],
            ["method", "b3lyp"],
            ["charge", "0"],
            ["dftd", "d3"],
            ["threall", "1.0e-14"],
            ["diismaxvecs", "40"],
            ["maxit", "200"],
        ]
        assert qmEngine.inputSettings == inputSettings
        assert qmEngine.backupInputSettings == backupInputSettings

    def test_readInputResp(self):
        os.chdir(os.path.dirname(__file__))
        qmEngine = qmengine.QMEngine(
            "qmengine/tcResp.in", "qmengine/tc_backup.in", doResp=True
        )
        inputSettings = [
            ["basis", "6-31gss"],
            ["method", "b3lyp"],
            ["charge", "0"],
            ["dftd", "d3"],
        ]
        backupInputSettings = [
            ["basis", "6-31gss"],
            ["method", "b3lyp"],
            ["charge", "0"],
            ["dftd", "d3"],
            ["threall", "1.0e-14"],
            ["diismaxvecs", "40"],
            ["maxit", "200"],
        ]
        assert qmEngine.inputSettings == inputSettings
        assert qmEngine.backupInputSettings == backupInputSettings

    # check that settings are written out correctly
    def test_writeInputFile(self):
        os.chdir(os.path.dirname(__file__))
        qmEngine = qmengine.QMEngine("qmengine/tc.in", "qmengine/tc_backup.in")
        inputSettings = [
            ["basis", "6-31gss"],
            ["method", "b3lyp"],
            ["charge", "0"],
            ["dftd", "d3"],
        ]
        testPath = "qmengine/test.in"
        qmEngine.writeInputFile(inputSettings, "coors.xyz", testPath)
        with open(testPath, "r") as f:
            assert "coordinates coors.xyz" in f.readline()
            assert "run gradient" in f.readline()
        readSettings = qmEngine.readInputFile(testPath)
        os.remove(testPath)
        assert readSettings == inputSettings

    def test_writeInputFileResp(self):
        os.chdir(os.path.dirname(__file__))
        qmEngine = qmengine.QMEngine(
            "qmengine/tc.in", "qmengine/tc_backup.in", doResp=True
        )
        inputSettings = [
            ["basis", "6-31gss"],
            ["method", "b3lyp"],
            ["charge", "0"],
            ["dftd", "d3"],
        ]
        testPath = "qmengine/testResp.in"
        qmEngine.writeInputFile(inputSettings, "coors.xyz", testPath)
        with open(testPath, "r") as f:
            assert "coordinates coors.xyz" in f.readline()
            assert "run gradient" in f.readline()
            assert "resp yes" in f.readline()
        readSettings = qmEngine.readInputFile(testPath)
        os.remove(testPath)
        assert readSettings == inputSettings

    # check that settings are written out correctly
    def test_writeInputFile(self):
        os.chdir(os.path.dirname(__file__))
        qmEngine = qmengine.QMEngine("qmengine/tc.in", "qmengine/tc_backup.in")
        inputSettings = [
            ["basis", "6-31gss"],
            ["method", "b3lyp"],
            ["charge", "0"],
            ["dftd", "d3"],
        ]
        testPath = "qmengine/test.in"
        qmEngine.writeInputFile(inputSettings, "coors.xyz", testPath)
        with open(testPath, "r") as f:
            assert "coordinates coors.xyz" in f.readline()
            assert "run gradient" in f.readline()
        readSettings = qmEngine.readInputFile(testPath)
        os.remove(testPath)
        assert readSettings == inputSettings

    def test_writeInputFileResp(self):
        os.chdir(os.path.dirname(__file__))
        qmEngine = qmengine.QMEngine(
            "qmengine/tc.in", "qmengine/tc_backup.in", doResp=True
        )
        inputSettings = [
            ["basis", "6-31gss"],
            ["method", "b3lyp"],
            ["charge", "0"],
            ["dftd", "d3"],
        ]
        testPath = "qmengine/testResp.in"
        qmEngine.writeInputFile(inputSettings, "coors.xyz", testPath)
        with open(testPath, "r") as f:
            assert "coordinates coors.xyz" in f.readline()
            assert "run gradient" in f.readline()
            assert "resp yes" in f.readline()
        readSettings = qmEngine.readInputFile(testPath)
        os.remove(testPath)
        assert readSettings == inputSettings

    # check that pdbs are read in correctly
    def test_readPDB(self):
        os.chdir(os.path.dirname(__file__))
        qmengine.QMEngine("qmengine/tc.in", "qmengine/tc_backup.in")
        testCoords = utils.readPDB("qmengine/test.pdb")
        coords = np.loadtxt("qmengine/coords.txt").flatten()
        assert checkUtils.checkArrays(coords, testCoords)

    # check that FB data is written correctly
    def test_writeFBdata(self):
        os.chdir(os.path.dirname(__file__))
        qmEngine = qmengine.QMEngine("qmengine/tc.in", "qmengine/tc_backup.in")

        coords = []
        energies = []
        grads = []
        for i in range(1, 26):
            coords.append(utils.readPDB(f"qmengine/test/{str(i)}.pdb"))
            energy, grad = utils.readGradFromTCout(f"qmengine/test/tc_{str(i)}.out")
            energies.append(energy)
            grads.append(grad)
        qmEngine.writeFBdata(energies, grads, coords)
        refLines = []
        with open("qmengine/test/all.mdcrd", "r") as refF:
            for line in refF.readlines():
                if len(line.split()) > 0:
                    refLines.append(line)
        testLines = []
        with open("all.mdcrd", "r") as testF:
            for line in testF.readlines():
                if len(line.split()) > 0:
                    testLines.append(line)
        assert len(refLines) == len(testLines)
        for i in range(1, len(refLines)):
            assert refLines[i] == testLines[i]
        refLines = []
        with open("qmengine/test/qdata.txt", "r") as refF:
            for line in refF.readlines():
                if len(line.split()) > 0:
                    refLines.append(line)
        testLines = []
        with open("qdata.txt", "r") as testF:
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
                    assert checkUtils.checkFloat(refLine[j], testLine[j], 0.0001)
                except:
                    assert refLine[j] == testLine[j]
        os.remove("qdata.txt")
        os.remove("all.mdcrd")

    def test_writeFBdataResp(self):
        os.chdir(os.path.dirname(__file__))
        qmEngine = qmengine.QMEngine(
            "qmengine/tc.in", "qmengine/tc_backup.in", doResp=True
        )

        coords = []
        energies = []
        grads = []
        espXYZs = []
        esps = []
        for i in range(1, 26):
            coords.append(utils.readPDB(f"qmengine/test/{str(i)}.pdb"))
            energy, grad = utils.readGradFromTCout(f"qmengine/test/tc_{str(i)}.out")
            energies.append(energy)
            grads.append(grad)
            espXYZ, esp = utils.readEsp(f"qmengine/test/esp_{str(i)}.xyz")
            espXYZs.append(espXYZ)
            esps.append(esp)
        qmEngine.writeFBdata(energies, grads, coords, espXYZs, esps)
        refLines = []
        with open("qmengine/qdataResp.txt", "r") as refF:
            for line in refF.readlines():
                if len(line.split()) > 0:
                    refLines.append(line)
        testLines = []
        with open("qdata.txt", "r") as testF:
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
                    assert checkUtils.checkFloat(refLine[j], testLine[j], 0.0001)
                except:
                    assert refLine[j] == testLine[j]
        os.remove("qdata.txt")
        os.remove("all.mdcrd")

    def test_writeResult(self):
        qmEngine = qmengine.QMEngine("qmengine/tc.in", "qmengine/tc_backup.in")
        os.chdir(os.path.dirname(__file__))
        with open(os.path.join("qmengine", "tc_1_ref.json"), "r") as f:
            refResult = AtomicResult(**json_loads(f.read()))
        qmEngine.writeResult(
            os.path.join("qmengine", "tc_1.out"), os.path.join("qmengine", "1.xyz")
        )
        with open(os.path.join("qmengine", "tc_1.json"), "r") as f:
            testResult = AtomicResult(**json_loads(f.read()))
        os.remove(os.path.join("qmengine", "tc_1.json"))
        assert checkUtils.checkFloat(
            refResult.properties.return_energy,
            testResult.properties.return_energy,
            0.00001,
        )
        assert checkUtils.checkArrays(
            refResult.return_result, testResult.return_result, 0.00001
        )

    def test_restart(self, monkeypatch):
        os.chdir(Path(__file__).parent / "qmengine")
        qmEngine = qmengine.QMEngine("tc.in", "tc_backup.in")

        def monkeyCompute(xyzs, folder):
            xyzs = sorted(xyzs)
            with open("xyzs.txt", "w") as f:
                for xyz in xyzs:
                    f.write(str(xyz) + "\n")
            return xyzs

        monkeypatch.setattr(qmEngine, "getQMRefData", monkeyCompute)
        os.chdir("restart")
        qmEngine.restart()
        testXyzs = []
        with open("xyzs.txt", "r") as f:
            for line in f.readlines():
                testXyzs.append(line.replace("\n", ""))
        os.remove("xyzs.txt")
        refXyzs = ["3.xyz", "6.xyz", "9.xyz"]
        assert len(testXyzs) == len(refXyzs)
        for i in range(len(testXyzs)):
            assert testXyzs[i] == refXyzs[i]
