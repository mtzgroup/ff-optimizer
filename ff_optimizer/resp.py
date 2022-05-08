import numpy as np
import os
from scipy.stats import norm
from tccloud.models import AtomicResult
from qcelemental.util.serialization import json_loads
#import matplotlib as mpl
#import matplotlib.pyplot as plt
#mpl.use("Agg")

class RespPriors:

    def __init__(self, options):
        self.sampleDir = options['sampleDir']
        self.allEsp = []
        self.allResp = []
        self.repeats, self.resName = getRepeats(options['mol2'])
        self.mode = options['mode']

    def findRepeatIndex(self, repeats, idx):
        loc = -1
        for j, idxs in enumerate(repeats):
            for i in idxs:
                if i == idx:
                    loc = j
        return loc
    
    def readCharges(self, lines):
        esp = []
        resp = []
        inEsp = False
        inResp = False
        almostInEsp = False
        almostInResp = False
        natoms = 0
        for line in lines():
            if "Total atoms:" in line:
                natoms = int(line.split()[2])
            if inEsp:
                esp.append(float(line.split()[4]))
                if len(esp) == natoms:
                    inEsp = False
            if inResp:
                resp.append(float(line.split()[4]))
                if len(resp) == natoms:
                    inResp = False
            if almostInEsp and "--------------------------------------" in line:
                inEsp = True
                almostInEsp = False
            if almostInResp and "--------------------------------------" in line:
                inResp = True
                almostInResp = False
            if "ESP unrestraint charges:" in line:
                almostInEsp = True
            if "ESP restraint charges:" in line:
                almostInResp = True
        if len(esp) == 0 or len(resp) == 0:
            raise RuntimeError("No charges in lines")
        return esp, resp
    
    def getCharges(self, i):
        cycleDir = f"{str(i)}_cycle_{str(i)}"
        for f in os.listdir(os.path.join(self.sampleDir, cycleDir)):
            if f.startswith("train") and os.path.isdir(os.path.join(self.sampleDir, cycleDir, f)):
                trainDir = os.path.join(self.sampleDir, cycleDir, f)
        for f in os.listdir(trainDir):
            if f.startswith("tc") and f.endswith("json"):
                try:
                    with open(os.path.join(trainDir,f), "r") as json:
                        result = AtomicResult(**json_loads(json.read()))
                    lines = result['stdout']
                    esp, resp = getCharges(lines)
                except:
                    name = f.split('.')[0]
                    with open(os.path.join(trainDir, f"{name}.out"),'r') as tcout
                        lines = tcout.readlines()
                    esp, resp = getCharges(lines)
                self.allEsp.append(esp)
                self.allResp.append(resp)

    def computeChargeDistributions(self):
        espArray = np.array(self.allEsp, dtype=np.float32, copy=True)
        respArray = np.array(allResp, dtype=np.float32, copy=True)
        natoms = espArray.shape[1]
        unitLength = len(self.repeats)
        if len(esp) % units != 0:
            raise ValueError("Must have an integer number of molecules to be fitted in the RESP calculation")
        self.espMeans = np.zeros(unitLength)
        self.espStdevs = np.ones(unitLength) * -1
        self.respMeans = np.zeros(unitLength)
        self.respStdevs = np.ones(unitLength) * -1
        for i, idxs in enumerate(repeats):
            if len(idxs) == 0:
                continue
            elif len(idxs) == 1:
                self.espMeans[i], self.espStdevs[i] = norm.fit(espArray[:,i::unitLength].flatten())
                self.respMeans[i], self.respStdevs[i] = norm.fit(respArray[:,i::unitLength].flatten())
            else:
                chargeList = []
                espList = []
                respList = []
                for idx in idxs:
                    espList += list(espArray[:,i::unitLength].flatten())
                    respList += list(allResp[:,i::unitLength].flatten())
                self.espMeans[i], self.espStdevs[i] = norm.fit(np.asarray(espList,dtype=np.float32))
                self.respMeans[i], self.respStdevs[i] = norm.fit(np.asarray(respList,dtype=np.float32))
    
        #fig, ax = plt.subplots()
        #plt.hist(allCharges[:,2::unitLength].flatten(),density=True,bins=12,edgecolor="black")
        #xmin, xmax = plt.xlim()
        #x = np.linspace(xmin, xmax, 100)
        #y = norm.pdf(x, means[2], stdevs[2])
        #plt.plot(x, y)
        #fig.set_dpi(200)
        #plt.savefig("RESP_charge_dist.png",bbox_inches="tight")
    
    def computePriors(self, cdf=0.95):
        # Priors determined by distribution of RESP charges
        if self.mode == 1:
            return self.respStdevs
        # Priors determined by a Gaussian centered on RESP charges with a standard deviation such that
        # the given value of the cdf occurs at the same charge for this distribution and the ESP distribution
        elif self.mode == 2:
            ppf = norm.ppf(cdf)
            # ESP - RESP shouldn't ever be negative, but with averages over many geometries, we have no such guarantee for means
            return np.abs(self.espMeans - self.respMeans) / ppf + self.espStdev
            
    def setMol2Charges(self, charges, mol2):
        inCharges = False
        with open(mol2, 'r') as inF:
            with open("temp.txt",'w') as outF:
                for line in inF.readlines():
                    if "@<TRIPOS>BOND" in line:
                        inCharges = False
                    if inCharges:
                        idx = int(line.split()[0]) - 1
                        i = findRepeatIndex(repeats, idx)
                        if i == -1:
                            outF.write(line)
                        else:
                            line2 = line.replace(line.split()[8],str(charges[i]))
                            outF.write(line2)
                    else:
                        outF.write(line)
                    if "@<TRIPOS>ATOM" in line:
                        inCharges = True
        os.rename("temp.txt",mol2)

    def setPriors(self, priors, inputFile):
        inOptions = False
        inPriors = False
        added = False
        with open(inputFile, 'r') as inF:
            with open("temp.txt",'w') as outF:
                for line in f.readlines():
                    if inOptions:
                        if inPriors:
                            if "COUL" in line:
                                continue
                            if line.split()[0] == "/priors":
                                for i in range(priors.shape[0]);
                                    if priors[i] > 0:
                                        outF.write(f"COUL:{self.resName}-{str(i+1)} : {str(priors[i])}\n")
                                inPriors = False
                                added = True
                            outF.write(line)
                        if line.split()[0] == "priors":
                            inPriors = True
                            outF.write(line)
                        if line.split()[0] == "$end" and added = False:
                            for i in range(priors.shape[0]);
                                if priors[i] > 0:
                                    outF.write(f"COUL:{self.resName}-{str(i+1)} : {str(priors[i])}\n")
                            outF.write(line)
                    elif "$options" in line:
                        inOptions = True
                        outF.write(line)
                    else:
                        outF.write(line)
        os.rename("temp.txt",inputFile)
    
    def getRepeats(self, mol2):
        inCharges = False
        repeats = []
        with open(mol2, 'r') as f:
            f.readline()
            f.readline()
            line = f.readline()
            unitLength = int(line.split()[0])
        for i in range(unitLength):
            repeats.append([])
        with open(mol2, 'r') as f:
            for line in f.readlines():
                if "@<TRIPOS>BOND" in line:
                    inCharges = False
                if inCharges:
                    splitLine = line.split()
                    if len(splitLine) >=11:
                        if splitLine[10] == "RPT":
                            index = int(splitLine[12].split('-')[-1]) - 1
                            repeats[index].append(int(splitLine[0])-1)
                        elif splitLine[10] == "PRM":
                            index = int(splitLine[0]) - 1
                            repeats[index].append(index)
                if "@<TRIPOS>ATOM" in line:
                    inCharges = True
        return repeats
    
    def updateRespPriors(self, i, mol2, inputFile):
        self.getCharges(i)
        self.computeChargeDistributions()
        priors = self.computePriors()
        self.setPriors(priors, inputFile)
        self.setMol2Charges(self.respMeans, mol2)

