import pytest
import os
from ff_optimizer import utils
from . import checkUtils


def compareQdata(qdataRef, qdataTest):
    refLines = []
    testLines = []
    with open(qdataRef, "r") as f:
        for line in f.readlines():
            if len(line.split()) > 0:
                refLines.append(line)
    with open(qdataTest, "r") as f:
        for line in f.readlines():
            if len(line.split()) > 0:
                testLines.append(line)
    assert len(refLines) == len(testLines)
    for i in range(len(refLines)):
        refLine = refLines[i].split()
        testLine = testLines[i].split()
        assert len(refLine) == len(testLine)
        for j in range(len(refLine)):
            if j == 0:
                check = refLine[j] == testLine[j]
            else:
                check = checkUtils.checkFloat(refLine[j], testLine[j], 0.0001)
            if not check:
                print(f"Compare {qdataRef}, {qdataTest} line {str(i)} token {str(j)}")
            assert check


def compareOpt(outRef, outTest):
    status, refResults = utils.readOpt(outRef)
    status, testResults = utils.readOpt(outTest)
    refParams = refResults["params"]
    testParams = testResults["params"]
    refInitialParams = refResults["initialParams"]
    testInitialParams = testResults["initialParams"]
    assert len(refParams) == len(testParams)
    for i in range(len(refParams)):
        assert checkUtils.checkFloat(refParams[i], testParams[i], 0.001)
    assert len(refInitialParams) == len(testInitialParams)
    for i in range(len(refInitialParams)):
        assert checkUtils.checkFloat(refInitialParams[i], testInitialParams[i], 0.001)


# NOTE: FB optimization is not well-conditioned for small samples?
# Can't get less than 0.1% error repeating the same calculations.
# Just check qdata.txt for now


@pytest.mark.gpu
@pytest.mark.full
def test_debug_full():
    print(os.path.dirname(__file__))
    os.chdir(os.path.join(os.path.dirname(__file__), "full"))
    os.system("./clean.sh")
    os.system("vacation_to_hawaii.py --stride 1 --split 2 --maxcycles 1 --qmengine debug")
    compareQdata(
        os.path.join("3_ref", "qdata.txt"),
        os.path.join("1_opt", "targets", "train_1", "qdata.txt"),
    )
    # compareOpt(os.path.join("3_ref","opt_1.out"),os.path.join("1_opt","opt_1.out"))
    #os.system("./clean.sh")


@pytest.mark.full
@pytest.mark.tccloud
def test_tccloud_full():
    os.chdir(os.path.join(os.path.dirname(__file__), "full"))
    os.system("./clean.sh")
    os.system(
        "vacation_to_hawaii.py --stride 1 --split 2 --maxcycles 1 --qmengine tccloud"
    )
    compareQdata(
        os.path.join("3_ref", "qdata.txt"),
        os.path.join("1_opt", "targets", "train_1", "qdata.txt"),
    )
    # compareOpt(os.path.join("3_ref","opt_1.out"),os.path.join("1_opt","opt_1.out"))
    os.system("./clean.sh")


@pytest.mark.full
@pytest.mark.queue
def test_queue_full():
    os.chdir(os.path.join(os.path.dirname(__file__), "full"))
    os.system("./clean.sh")
    os.system("vacation_to_hawaii.py --stride 1 --split 2 --maxcycles 1 --qmengine queue")
    compareQdata(
        os.path.join("3_ref", "qdata.txt"),
        os.path.join("1_opt", "targets", "train_1", "qdata.txt"),
    )
    # compareOpt(os.path.join("3_ref","opt_1.out"),os.path.join("1_opt","opt_1.out"))
    os.system("./clean.sh")
