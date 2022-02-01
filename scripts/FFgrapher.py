#!/usr/bin/env python

import os
import argparse
import random
import matplotlib as mpl
import numpy as np
mpl.use('Agg')
import matplotlib.pyplot as plt
import subprocess
from time import sleep
from textwrap import dedent

# Some helper functions
def addTargetLines(inputFile, targetLines, initialTarget, newTarget):
    addedLines = False
    with open(inputFile,'r') as f:
        for line in f.readlines():
            splitLine = line.split()
            if len(splitLine) > 1:
                if splitLine[0] == "name" and splitLine[1] == newTarget:
                    addedLines = True
    if not addedLines:
        with open(inputFile,'a') as f:
            f.write("\n")
            for line in targetLines:
                f.write(line.replace(initialTarget,newTarget))

def slurmCommand(command):
    maxTries = 100
    i = 1
    done = False
    while i < maxTries and not done:
        try:
            output = subprocess.check_output(command)
            done = True
        except:
            sleep(2)
            i += 1
    if not done:
        raise RuntimeError("Slurm command " + command + " failed")
    return output

def convertTCtoFB(tcout, coors, stride, start=None, end=None, qdata="qdata.txt", mdcrd="all.mdcrd"):

    molSize = 0
    coords = []
    frame = []
    energies = []
    grads = []
    coordIndex = []
    gradIndex = []
    lineCounter = 0
    frameStart = -1
    coorEnergies = []
    
    with open(coors,'r') as f:
        maxFrame = -1
        for line in f.readlines():
            splitLine = line.split()
            if lineCounter == 0:
                molSize = int(splitLine[0])
            if lineCounter == 1:
                index = int(splitLine[2]) + 1 # coor files are 0-indexed
                energy = splitLine[0]
                if frameStart == -1:
                    frameStart = index
            elif lineCounter > 1 and lineCounter < molSize + 2:
                for token in splitLine[1:]:
                    frame.append(token)
            if lineCounter == molSize + 1:
                if index > maxFrame:
                    maxFrame = index
                    coords.append(frame)
                    coordIndex.append(index)
                    coorEnergies.append(energy)
                else:
                    j = len(coords) - 1
                    while coordIndex[j] > index:
                        j -= 1
                    coords[j] = frame
                    coordIndex[j] = index
                    coorEnergies[j] = energy
                frame = []
                lineCounter = -1
            lineCounter = lineCounter + 1
    
    gradCounter = molSize + 37
    index = 0
    frameStart = -1
    
    with open(tcout,'r') as f:
        maxFrame = -1
        for line in f.readlines():
            if "MD STEP" in line: 
                index = int(line.split()[4])
                if frameStart == -1:
                    if len(grads) > 0:
                        frameStart = index - 1 # first step is unnumbered "step 0"
                        gradIndex[0] = index - 1
                    else:
                        frameStart = index
            if "FINAL ENERGY" in line:
                energy = float(line.split()[2]) 
            if "Gradient units" in line:
                gradCounter = 0
                frame = []
            if gradCounter < molSize + 2 and gradCounter > 2:
                for token in line.split():
                    frame.append(token)
            if gradCounter == molSize + 2:
                for token in line.split():
                    frame.append(token)
                if index > maxFrame:
                    maxFrame = index
                    grads.append(frame)
                    gradIndex.append(index)
                    energies.append(energy)
                else:
                    j = len(grads) - 1
                    while gradIndex[j] > index:
                        j -= 1
                    grads[j] = frame
                    gradIndex[j] = index
                    energies[j] = energy
            gradCounter = gradCounter + 1
    
    if start == None:
        start = gradIndex[0]
    if end == None:
        end = gradIndex[-1]

    jobCounter = 0
    gradsCopy = []
    coordsCopy = []
    indices = []
    precision = len(coorEnergies[0].split('.')[1])
    for i in range(len(gradIndex)):
        for j in range(len(coordIndex)):
            if gradIndex[i] == coordIndex[j]:
                eFormat = "%." + str(precision) + "f"
                gradEnergy = eFormat % energies[i]
                if coorEnergies[j] !=  gradEnergy:
                    #print("Mismatched energies in step " + str(gradIndex[i]))
                    #raise RuntimeError("Mismatched energies from " + tcout + " and " + coors)
                    break
                indices.append([i,j])

    
    usedIndices = []
    lastFrame = -args.stride - 37
    for i in range(len(indices)):
        if gradIndex[indices[i][1]] >= start and gradIndex[indices[i][1]] <= end:
            if gradIndex[indices[i][1]] - lastFrame >= stride:
                lastFrame = gradIndex[indices[i][1]]
                usedIndices.append(indices[i])

    with open(qdata,'w') as f:
        for j in usedIndices:
            f.write("\n")
            f.write("JOB " + str(jobCounter) + "\n")
            coordLine = "COORDS "
            for coord in coords[j[1]]:
                coordLine = coordLine + str(coord) + " "
            gradLine = "FORCES "
            for grad in grads[j[0]]:
                gradLine = gradLine + str(grad) + " "
            f.write(coordLine + "\n")
            f.write("ENERGY " + str(energies[j[0]]) + "\n")
            f.write(gradLine + "\n")
            jobCounter = jobCounter + 1
            
    with open(mdcrd,'w') as f:
        f.write("converted from TC with convertTCMDtoFB.py\n")
        for j in usedIndices:
            tokenCounter = 1
            for coord in coords[j[1]]:
                f.write("%8.3f" % float(coord))
                if tokenCounter == 10:
                    f.write("\n")
                    tokenCounter = 1
                else:
                    tokenCounter = tokenCounter + 1
            if tokenCounter != 1:
                f.write("\n")
    return jobCounter

def changeParameter(inputFile, prmName, prmValue):
    changed = False
    lines = []
    with open(inputFile,'r') as f:
        for line in f.readlines():
            if prmName in line:
                line = line.replace(line.split()[1],prmValue)
                changed = True
            lines.append(line)
        if not changed:    
            lines.insert(1,f"{prmName} {prmValue}")

    with open("temp.txt",'w') as f:
        for line in lines:
            f.write(line)
    os.system(f"mv temp.txt {inputFile}")

def readOpt(filename):
    inInitialParams = False
    inFinalParams = False
    inFinalObj = False
    params = []
    initialParams = []
    labels = []
    results = {}
    status = -1
    with open(filename,'r') as f:
        for line in f.readlines():
            if "-------" in line:
                inFinalParams = False
                inInitialParams = False
            if inFinalParams:
                if "=======" not in line:
                    splitLine = line.split()
                    params.append(splitLine[2])
                    labels.append(splitLine[5])
            if inInitialParams:
                if "=======" not in line:
                    initialParams.append(line.split()[2])
            if inFinalObj:
                results['obj'] = float(line.split()[5])
                inFinalObj = False
            if "Final physical parameters" in line:
                inFinalParams = True
            if "Starting parameter indices" in line:
                inInitialParams = True
            if "Final objective function" in line:
                inFinalObj = True
            if "Optimization Converged" in line:
                status = 0
            if "Maximum number of optimization steps reached" in line:
                status = 2

    #if status == -1:
    #    raise RuntimeError("ForceBalance optimization of " + filename + " failed")

    params = np.asarray(params,dtype=np.float32)
    initialParams = np.asarray(initialParams,dtype=np.float32)
    results['params'] = params
    results['labels'] = labels
    results['initialParams'] = initialParams
    return status, results

def determineAdaptiveDamping(testFile,upperThreshold=0.3,lowerThreshold=0.01,adaptiveDamping=0.5):
    changeParameter(testFile,"adaptive_damping",str(adaptiveDamping))
    maxCycles = 100
    testOut = testFile.split('.')[0] + ".out"
    for j in range(maxCycles):
        os.system(f"ForceBalance.py {testFile} > {testOut}")
        status, results = readOpt(testOut)
        diff = np.abs(results['params'] - results['initialParams']) / np.maximum(results['params'], results['initialParams'])
        if np.argwhere(diff > upperThreshold).shape[0] > 0:
            adaptiveDamping *= 2
            changeParameter(testFile,"adaptive_damping",str(adaptiveDamping))
        elif np.argwhere(diff > lowerThreshold).shape[0] == 0:
            adaptiveDamping *= 0.75
            changeParameter(testFile,"adaptive_damping",str(adaptiveDamping))
        else:
            return adaptiveDamping

def readPDB(filename):
    coords = []
    with open(filename,'r') as f:
        for line in f.readlines():
            if "ATOM" in line or "HETATM" in line:
                splitLine = line.split()
                coords.append(splitLine[5])
                coords.append(splitLine[6])
                coords.append(splitLine[7])
    return coords

def readRst(filename):
    coords = []
    lineCounter = 0
    with open(filename,'r') as f:
        for line in f.readlines():
            if lineCounter > 1:
                for coord in line.split():
                    coords.append(coord)
            lineCounter = lineCounter + 1
    return coords

def readGradFromTCout(filename):
    inGradient = False
    gradCounter = 0
    molSize = 0
    grads = []
    energy = 0
    with open(filename,'r') as f:
        for line in f.readlines():
            if "Total atoms" in line:
                molSize = int(line.split()[2])
            if "FINAL ENERGY" in line:
                energy = float(line.split()[2]) 
            if "Gradient units" in line:
                inGradient = True
            if gradCounter < molSize + 3 and gradCounter > 2:
                for token in line.split():
                    grads.append(token)
            if inGradient:
                gradCounter = gradCounter + 1
    if not inGradient:
        print("File %s did not complete calculation" % filename)
        grads = -1
    return energy, grads

def writeData(coords, energy, grads):
    with open("qdata.txt",'a') as f:
        coordLine = "COORDS "
        for coord in coords:
            coordLine = coordLine + str(coord) + " "
        gradLine = "FORCES "
        for grad in grads:
            gradLine = gradLine + str(grad) + " "
        f.write(coordLine + "\n")
        f.write("ENERGY " + str(energy) + "\n")
        f.write(gradLine + "\n")
        f.write("\n")
        
    with open("all.mdcrd",'a') as f:
        tokenCounter = 1
        for coord in coords:
            f.write("%8.3f" % float(coord))
            if tokenCounter == 10:
                f.write("\n")
                tokenCounter = 1
            else:
                tokenCounter = tokenCounter + 1
        if tokenCounter != 1:
            f.write("\n")

def readValid(filename):
    with open(filename,'r') as f:
        for line in f.readlines():
            if "Objective Function Single Point" in line:
                return float(line.split()[6])
        raise RuntimeError("ForceBalance single-point evaluation of " + filename + " did not converge")

def writeRst(frame, index, molSize, dest):
    fileName = os.path.join(dest,str(index) + ".rst7")
    with open(fileName,'w') as f:
        f.write("sampled from frame " + str(index) + "\n")
        f.write(str(molSize) + "\n")
        for i in range(len(frame)):
            f.write("%12.7f%12.7f%12.7f" % (float(frame[i][0]),float(frame[i][1]), float(frame[i][2])))
            if int(i / 2) * 2 != i:
                 f.write("\n")

def getFrame(index, coordFile, dest):
    frame = []
    atoms = []
    inFrame = False
    lineCounter = 1
    frameCounter = 1
    with open(coordFile,'r') as f:
        molSize = int(f.readline())
        for line in f.readlines():
            if inFrame:
                atoms.append(line.split()[0])
                frame.append(line.split()[1:])
                lineCounter = lineCounter + 1    
            if "frame " in line:
                if frameCounter == index:
                    inFrame = True
                frameCounter += 1
            if lineCounter > molSize:
                inFrame = False
    writeRst(frame, index, molSize, dest)

# Summary stuff
summary = dedent('''\
Summary of necessary directories and files
Dynamics directory (--dynamicsdir, 0_dynamics) contains:
    -- QM dynamics coordinates (--coors, coors.xyz)
    -- TeraChem output file (--tcout, tc.out)
FB Optimization directory (--optdir, 1_opt) contains:
    -- FB optimization input file (--opt0, opt_0.in)
    -- FB validation input file (--valid0, valid_0.in)
    -- PDB conformation file (conf.pdb)
    -- Tleap setup file (setup.leap)
    -- Parameter frcmod file (*.frcmod)
    -- Parameter mol2 file (*.mol2)
MM sampling directory (--sampledir, 2_sampling) contains:
    -- Amber equilibration inputs (heat*in)
    -- Amber MM sampling input (md.in)
    -- Cpptraj input file (cpptraj.in)
    -- TeraChem input file for fire (tc_template.in)
    -- Backup TeraChem input file if job fails (tc_template_long.in)
    -- sbatch template file for fire (sbatch_template.sh)
''')
parser = argparse.ArgumentParser(epilog=summary,formatter_class=argparse.RawDescriptionHelpFormatter)
parser.add_argument('--optdir', help="Directory where ForceBalance optimization is performed, default is 1_opt", type=str,default="1_opt")
args = parser.parse_args()

# Set up the calculation
train = []
valid = []
validInitial = []
validPrevious = []
validFinalTarget = []

# Determine cycle for restart, set restart variables
totalCycles = -1
types = ['BONDSK','BONDSB','ANGLESK','ANGLESB','DIHS']
aliases = ['Bond strength', 'Bond length', 'Angle strength', 'Equilibrium angle','Dihedral strength']
colors = ['blue','green','firebrick','goldenrod','orange','purple','lightskyblue','olive']
optCounter = 0
for f in os.listdir(args.optdir):
    if "opt_" in f and ".out" in f:
        optCounter += 1

for i in range(optCounter+1):
    optOutput = os.path.join(args.optdir,"opt_" + str(i) + ".out")
    if os.path.isfile(optOutput):
        status, results = readOpt(optOutput)
        if status == 0:
            if i == 0:
                params = np.zeros((optCounter+1,len(results['labels'])))
                labels = results['labels']
                params[0,:] = results['initialParams']
            else:
                try:
                    valid.append(readValid(os.path.join(args.optdir,"valid_" + str(i) + ".out")))
                    validPrevious.append(readValid(os.path.join(args.optdir,"valid_" + str(i) + "_previous.out")))
                    validInitial.append(readValid(os.path.join(args.optdir,"valid_" + str(i) + "_initial.out")))
                except:
                    print("valids didn't complete")
                    break
            train.append(results['obj'])
            for j in range(len(results['labels'])):
                if labels[j] == results['labels'][j]:
                    params[i+1,j] = results['params'][j]
                else:
                    for k in range(j+1,len(labels)):
                        if labels[k] == results['labels'][j]:
                            params[i+1,k] = results['params'][j]
                            break
        else:
            print("opt didn't complete " + str(status))
            break
    else:
        break
    totalCycles = i
i = totalCycles

for i in range(totalCycles+1):
    if os.path.isfile(os.path.join(args.optdir,"valid_" + str(i) + "_final_target.out")):
        try:
            validFinalTarget.append(readValid(os.path.join(args.optdir,"valid_" + str(i) + "_final_target.out")))
        except:
            src = os.path.join(args.optdir,"result","opt_" + str(i),"*")
            dest = os.path.join(args.optdir,"forcefield")
            os.system(f"cp {src} {dest}")
            os.chdir(args.optdir)
            os.system("ForceBalance.py valid_" + str(totalCycles) + ".in > valid_" + str(i) + "_final_target.out")
            os.chdir("..")
            validFinalTarget.append(readValid(os.path.join(args.optdir,"valid_" + str(i) + "_final_target.out")))
    else:
        src = os.path.join(args.optdir,"result","opt_" + str(i),"*")
        dest = os.path.join(args.optdir,"forcefield")
        os.system(f"cp {src} {dest}")
        os.chdir(args.optdir)
        os.system("ForceBalance.py valid_" + str(totalCycles) + ".in > valid_" + str(i) + "_final_target.out")
        os.chdir("..")
        validFinalTarget.append(readValid(os.path.join(args.optdir,"valid_" + str(i) + "_final_target.out")))

# Graph results so far
x = range(1,i + 1)
x0 = np.arange(i + 1)
if i < 25:
    xticks = x0
else:
    xticks = np.arange(0,i + 1, 2)
fig, ax = plt.subplots(figsize=(9,6))
ax.plot(x, valid, label = 'Validation, current parameters', marker='o')
ax.plot(x, validPrevious, label = 'Validation, previous parameters', marker='o')
ax.plot(x, validInitial, label = 'Validation, initial parameters', marker='o')
ax.plot(x0, train, label = 'Training',marker='o')
ax.set_xlabel("Optimization cycle", size = 17)
ax.set_ylabel("Objective function", size = 17)
ax.set_xticks(xticks)
fig.set_dpi(200)
plt.legend(fontsize=14)
plt.savefig('ObjectiveFunction.png',bbox_inches='tight')
plt.close()

fig, ax = plt.subplots(figsize=(9,6))
ax.plot(x0, validFinalTarget, label = 'Performance on final validation set', marker='o')
ax.set_xlabel("Optimization cycle", size = 17)
ax.set_ylabel("Objective function", size = 17)
ax.set_xticks(xticks)
fig.set_dpi(200)
plt.legend(fontsize=14)
plt.savefig('ObjectiveFunctionFinalTarget.png',bbox_inches='tight')
plt.close()

types = ['BONDSK','BONDSB','ANGLESK','ANGLESB','DIHS','VDWS','VDWT','COUL']
aliases = ['Bond strength', 'Bond length', 'Angle strength', 'Equilibrium angle','Dihedral strength','LJ sigma','LJ epsilon','Atomic charge']
colors = ['blue','green','firebrick','goldenrod','orange','purple','lightskyblue','olive']
name = "opt_" + str(i) + ".out"
for j in range(len(results['labels'])):
    if labels[j] == results['labels'][j]:
        params[i+1,j] = results['params'][j]
    else:
        for k in range(j+1,len(labels)):
            if labels[k] == results['labels'][j]:
                params[i+1,k] = results['params'][j]
                break
adds = []
scs = []
legendLabels = []
sortedParams = []
for k in range(len(types)):
    adds.append(False)
    sortedParams.append([])
fig, ax = plt.subplots(figsize=(9,6))
for k in range(len(labels)):
    for j in range(len(types)):
        if types[j] in labels[k]:
            sortedParams[j].append(params[:,k])
            #sc = ax.scatter(range(1,i + 1),relError[:,i],label=types[j],c=colors[j])
            if not adds[j]:
                #scs.append(sc)
                #legendLabels.append(types[j])
                adds[j] = True
            break

#ax.set_ylim([-100,100])
#plt.legend(scs,legendLabels,fontsize=14)
#plt.savefig('params.png',bbox_inches='tight')
plt.close()
fig, ax = plt.subplots(figsize=(9,6))
ax.set_xticks(xticks)
for j in range(len(types)):
    if len(sortedParams[j]) == 0:
        continue
    sortedParams[j] = np.asarray(sortedParams[j],dtype=np.float32)
    diff = (sortedParams[j] - np.roll(sortedParams[j],1,axis=1))[:,1:]
    weights = np.maximum(np.abs(sortedParams[j]),np.roll(np.abs(sortedParams[j]),1,axis=1))[:,1:]
    normalizedDiff = np.zeros(i)
    mrc = np.zeros(i+1)
    for k in range(i+1):
        #normalizedDiff[j] = np.sqrt(np.dot(diff[:,j],diff[:,j]) / np.dot(sortedParams[i][:,j],sortedParams[i][:,j]))
        mrc[k] = np.mean(np.abs(diff[:,k]) / weights[:,k]) * 100
    #plt.plot(range(1,i+1),normalizedDiff,label=aliases[i],marker='o')
    plt.plot(range(i+1),mrc,label=aliases[j],marker='o')
plt.legend(fontsize=14)
ax.tick_params(labelsize=14)
ax.set_xlabel('Optimization Cycle',size=17)
#ax.set_ylabel('Normalized RMS parameter change',size=17)
ax.set_ylabel('Mean relative parameter change / %',size=17)
plt.savefig("ParameterChange.png",bbox_inches='tight')
plt.close()
