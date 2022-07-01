import pytest
from ff_optimizer import active_learning
import numpy as np
import os
from . import checkUtils

def test_computeEnergyForceSP():
    os.chdir(os.path.join(os.path.dirname(__file__), "active_learning"))
    geometry = np.loadtxt("coords.xyz",usecols=(1,2,3),skiprows=2,max_rows=84).flatten()
    force = np.loadtxt("MMforce.xyz",usecols=(1,2,3),skiprows=2,max_rows=84).flatten()
    
    results = active_learning.computeEnergyForce([geometry], "amber.prmtop")
    checkUtils.checkArray(results[0][1],force)

def readXYZTraj(filename):
    frame = []
    frames = []
    i = 1
    natoms = -1
    with open(filename, 'r') as f:
        for line in f.readlines():
            if natoms == -1:
                natoms = int(line.split()[0])
            if i > natoms + 2:
                frames.append(np.asarray(frame,dtype=np.float32))
                frame = []
                i = 1
            if i > 2:
                splitLine = line.split()
                frame.append(splitLine[1])
                frame.append(splitLine[2])
                frame.append(splitLine[3])
            i += 1
    return frames

def test_computeEnergyForceAll():
    os.chdir(os.path.join(os.path.dirname(__file__), "active_learning"))
    frames = readXYZTraj("coords.xyz")
    results = active_learning.computeEnergyForce(frames, "amber.prmtop")
    forces = readXYZTraj("MMforce.xyz")
    for i in range(len(frames)):
        checkUtils.checkArray(results[i][1],forces[i])
