import os
from pathlib import Path

import chemcloud as cc
import numpy as np
import pytest
from qcio import Files, ProgramOutput
from qcparse import parse

from . import checkUtils
from . import chemcloud_mocks as ccm


def test_inputFromTCin():
    os.chdir(os.path.dirname(__file__))
    os.chdir("mockchemcloud")
    inp = ccm.inputFromTCin("test.in")
    inp.save("test.yaml")
    test = checkUtils.checkFiles("test.yaml", "ref.yaml")
    os.remove("test.yaml")
    assert test


def test_inputFromTCin2():
    os.chdir(os.path.dirname(__file__))
    os.chdir("mockchemcloud")
    inp = ccm.inputFromTCin("resp.in", extras={"id": 3})
    inp.save("test.yaml")
    test = checkUtils.checkFiles("test.yaml", "ref_resp.yaml")
    os.remove("test.yaml")
    assert test


@pytest.mark.chemcloud
def test_RunInputFromTCin():
    os.chdir(os.path.dirname(__file__))
    os.chdir("mockchemcloud")
    inp = ccm.inputFromTCin("test.in")
    out = cc.compute("terachem", inp)
    ccVib = out.results.freqs_wavenumber
    ref = parse("ref.out", "terachem")
    tcVib = ref.freqs_wavenumber
    test = checkUtils.checkListsFloats(ccVib, tcVib, thresh=0.001)


def test_inputFromTCout():
    os.chdir(os.path.dirname(__file__))
    os.chdir("mockchemcloud")
    inp = ccm.inputFromTCout("ref.out")
    inp.save("test.yaml")
    test = checkUtils.checkFiles("test.yaml", "ref_from_out.yaml")
    os.remove("test.yaml")
    assert test


def test_collectFiles1():
    os.chdir(os.path.dirname(__file__))
    os.chdir("mockchemcloud")
    inp = ccm.inputFromTCin("resp.in")
    files = ccm.collectFiles("scr.wat", inp)
    assert len(files.keys()) == 11
    for key in files.keys():
        if key.startswith("scr"):
            file = key.replace("geometry", "wat", 2)
            lines = ccm.readLines(file)
            assert lines == files[key]


def test_collectFiles2():
    os.chdir(os.path.dirname(__file__))
    os.chdir("mockchemcloud")
    inp = ccm.inputFromTCin("wat1.in")
    files = ccm.collectFiles("scr.oops", inp)
    assert len(files.keys()) == 2
    lines = ccm.readLines("ccinput.in")
    assert files["tc.in"] == lines
    refCoords = np.loadtxt("wat1.xyz", skiprows=2, usecols=(1, 2, 3))
    coordLines = files["geometry.xyz"].split("\n")[2:]
    testCoords = [line.split()[1:] for line in coordLines][:3]
    testCoords = np.asarray(testCoords)
    assert checkUtils.checkArrays(refCoords, testCoords, 1e-8)


def test_programOutputFromTCout1():
    os.chdir(os.path.dirname(__file__))
    os.chdir("mockchemcloud")
    output = ccm.programOutputFromTCout("resp.out", getFiles=True)
    charges = readCharges(output.stdout)
    assert checkUtils.checkFloats(charges[0], -1.153575)
    assert output.input_data.model.method == "Hartree-Fock"
    assert len(output.results.files.keys()) == 11
    assert output.success


def test_programOutputFromTCout2():
    os.chdir(os.path.dirname(__file__))
    os.chdir("mockchemcloud")
    output = ccm.programOutputFromTCout("resp.out", getFiles=False)
    charges = readCharges(output.stdout)
    assert checkUtils.checkFloats(charges[0], -1.153575)
    assert output.input_data.model.method == "Hartree-Fock"
    assert len(output.results.files.keys()) == 0
    assert output.success


def test_programOutputFromTCout3():
    os.chdir(os.path.dirname(__file__))
    os.chdir("mockchemcloud")
    inp = ccm.inputFromTCin("test.in")
    output = ccm.programOutputFromTCout("ref.out", inp=inp)
    assert output.input_data.model.method == "ub3lyp"
    assert len(output.results.files.keys()) == 0
    assert output.success


def test_programOutputFromTCout4():
    os.chdir(os.path.dirname(__file__))
    os.chdir("mockchemcloud")
    inp = ccm.inputFromTCin("fail1.in")
    output = ccm.programOutputFromTCout("fail1.out", inp=inp)
    assert not output.success
    assert isinstance(output.results, Files)
    assert not output.results.files
    assert len(output.traceback) > 0


def test_tclinesFromInput():
    os.chdir(os.path.dirname(__file__))
    os.chdir("mockchemcloud")
    inp = ccm.inputFromTCin("wat1.in")
    lines = ccm.tclinesFromInput(inp)
    with open("temp.in", "w") as f:
        f.write(lines)
    test = checkUtils.checkFiles("temp.in", "ccinput.in")
    os.remove("temp.in")
    assert test


def readCharges(stdout):
    lines = stdout.split("\n")
    charges = []
    i = 0
    for line in lines:
        if "ESP restrained charges:" in line:
            break
        i = i + 1
    for line in lines[i + 3 :]:
        if "-----------------------------------------------------------" in line:
            break
        splitLine = line.split()
        charges.append(splitLine[4])
    return charges


def computeWithFiles():
    os.chdir(os.path.dirname(__file__))
    os.chdir("mockchemcloud")
    inputs = [ccm.inputFromTCin("wat1.in"), ccm.inputFromTCin("wat2.in")]
    outputs = cc.compute("terachem", inputs, collect_files=True)
    assert len(outputs) == 2
    assert isinstance(outputs[0], ProgramOutput)
    result = outputs[0].results
    assert checkUtils.checkFloats(result.energy, -75.5788691372, 0.00000001)
    assert checkUtils.checkFloats(result.gradient[0, 0], -0.0006380621, 0.000001)
    charges = readCharges(outputs[0].stdout)
    assert checkUtils.checkFloats(charges[0], -0.260148)
    espLines = result.files["scr.geometry/esp.xyz"].split("\n")
    esp = espLines[2].split()[4]
    assert checkUtils.checkFloats(esp, 0.260066444823)
    with open(Path("scr.wat1") / "wat1.geometry", "r") as f:
        lines = list(f.readlines())
    lines = "".join(lines).split()
    assert checkUtils.checkListsFloats(
        lines, result.files["scr.geometry/geometry.geometry"].split()
    )


@pytest.mark.chemcloud
def test_computeWithFiles(monkeypatch):
    computeWithFiles()
    os.chdir(os.path.dirname(__file__))
    os.chdir("mockchemcloud")
    outputs = ["wat1.out", "wat2.out"]
    ccm.patcher(monkeypatch, cc, outputs)
    computeWithFiles()


def computeNoFiles():
    os.chdir(os.path.dirname(__file__))
    os.chdir("mockchemcloud")
    inputs = [ccm.inputFromTCin("wat1.in"), ccm.inputFromTCin("wat2.in")]
    outputs = cc.compute("terachem", inputs)
    assert len(outputs) == 2
    assert isinstance(outputs[1], ProgramOutput)
    result = outputs[1].results
    assert checkUtils.checkFloats(result.energy, -75.5769406865, 0.00000001)
    assert checkUtils.checkFloats(result.gradient[0, 0], 0.0069999056, 0.000001)
    assert len(result.files.keys()) == 0


@pytest.mark.chemcloud
def test_computeNoFiles(monkeypatch):
    computeNoFiles()
    os.chdir(os.path.dirname(__file__))
    os.chdir("mockchemcloud")
    outputs = ["wat1.out", "wat2.out"]
    ccm.patcher(monkeypatch, cc, outputs)
    computeNoFiles()


def compute1list():
    os.chdir(os.path.dirname(__file__))
    os.chdir("mockchemcloud")
    inputs = [ccm.inputFromTCin("wat1.in", extras={"bloop": True})]
    outputs = cc.compute("terachem", inputs)
    assert isinstance(outputs, list)
    result = outputs[0].results
    assert checkUtils.checkFloats(result.energy, -75.5788691372, 0.00000001)
    assert checkUtils.checkFloats(result.gradient[0, 0], -0.0006380621, 0.000001)
    assert outputs[0].input_data.extras["bloop"]


def compute1input():
    os.chdir(os.path.dirname(__file__))
    os.chdir("mockchemcloud")
    inputs = ccm.inputFromTCin("wat1.in", extras={"bloop": True})
    outputs = cc.compute("terachem", inputs)
    assert isinstance(outputs, ProgramOutput)
    result = outputs.results
    assert checkUtils.checkFloats(result.energy, -75.5788691372, 0.00000001)
    assert checkUtils.checkFloats(result.gradient[0, 0], -0.0006380621, 0.000001)
    assert outputs.input_data.extras["bloop"]


@pytest.mark.chemcloud
def test_compute1(monkeypatch):
    compute1list()
    compute1input()
    os.chdir(os.path.dirname(__file__))
    os.chdir("mockchemcloud")
    outputs = ["wat1.out"]
    ccm.patcher(monkeypatch, cc, outputs)
    compute1list()
    compute1input()


def computeFuture():
    os.chdir(os.path.dirname(__file__))
    os.chdir("mockchemcloud")
    inputs = [ccm.inputFromTCin("wat1.in"), ccm.inputFromTCin("wat2.in")]
    future = cc.compute("terachem", inputs, return_future=True)
    outputs = future.get()
    assert len(outputs) == 2
    assert isinstance(outputs[0], ProgramOutput)


@pytest.mark.chemcloud
def test_computeFuture(monkeypatch):
    computeFuture()
    os.chdir(os.path.dirname(__file__))
    os.chdir("mockchemcloud")
    outputs = ["wat1.out", "wat2.out"]
    ccm.patcher(monkeypatch, cc, outputs)
    computeFuture()


def computeFuture1list():
    os.chdir(os.path.dirname(__file__))
    os.chdir("mockchemcloud")
    inputs = [ccm.inputFromTCin("wat1.in")]
    future = cc.compute("terachem", inputs, return_future=True)
    outputs = future.get()
    assert isinstance(outputs, list)


def computeFuture1input():
    os.chdir(os.path.dirname(__file__))
    os.chdir("mockchemcloud")
    inputs = ccm.inputFromTCin("wat1.in")
    future = cc.compute("terachem", inputs, return_future=True)
    outputs = future.get()
    assert isinstance(outputs, ProgramOutput)


@pytest.mark.chemcloud
def test_computeFuture1(monkeypatch):
    computeFuture1list()
    computeFuture1input()
    os.chdir(os.path.dirname(__file__))
    os.chdir("mockchemcloud")
    outputs = ["wat1.out"]
    ccm.patcher(monkeypatch, cc, outputs)
    computeFuture1list()
    computeFuture1input()


def computeFailure():
    os.chdir(os.path.dirname(__file__))
    os.chdir("mockchemcloud")
    inputs = [ccm.inputFromTCin("fail1.in"), ccm.inputFromTCin("fail2.in")]
    outputs = cc.compute("terachem", inputs)
    assert len(outputs) == 2
    assert isinstance(outputs[0], ProgramOutput)
    assert isinstance(outputs[0].results, Files)
    assert not outputs[1].results.files
    assert not outputs[0].success
    assert not outputs[1].success
    assert len(outputs[0].traceback) > 0


@pytest.mark.chemcloud
def test_computeFailure(monkeypatch):
    computeFailure()
    os.chdir(os.path.dirname(__file__))
    os.chdir("mockchemcloud")
    outputs = ["fail1.out", "fail2.out"]
    ccm.patcher(monkeypatch, cc, outputs)
    computeFailure()
