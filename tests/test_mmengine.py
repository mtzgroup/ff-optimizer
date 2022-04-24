import pytest
from ff_optimizer import mmengine
import os
from . import checkUtils
from numpy import loadtxt
import GPUtil

options = {}
options['start'] = 33
options['end'] = 2000
options['split'] = 1000
options['stride'] = 50
options['coordsDir'] = "mmengine"
options['coords'] = "coors.xyz"

def test_getIndices():
    os.chdir(os.path.dirname(__file__))
    mmEngine = mmengine.MMEngine(options)
    testStart, testEnd, testSplit = mmEngine.getIndices()
    frames = loadtxt("mmengine/frames.txt",dtype=int)
    pass
    #assert frames[testStart] == options['start']
    #assert frames[testEnd] == options['end']
    #assert frames[testSplit] == options['split']

def test_getFrames():
    pass

# Also tests utils.writeRst()
def test_getFrame():
    os.chdir(os.path.dirname(__file__))
    options['start'] = None
    mmEngine = mmengine.MMEngine(options)
    mmEngine.getFrame(23, "23.rst7")
    testCoors = []
    with open("23.rst7",'r') as f:
        for line in f.readlines()[2:]:
            testCoors.append(line.split())
    os.remove("23.rst7")
    refCoors = []
    with open(os.path.join("mmengine","23.rst7"),'r') as f:
        for line in f.readlines()[2:]:
            refCoors.append(line.split())
    checkUtils.checkArray(testCoors, refCoors)

