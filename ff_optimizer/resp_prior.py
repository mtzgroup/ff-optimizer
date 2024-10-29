import os
from pathlib import Path
from typing import List, Tuple

import numpy as np
from scipy.stats import norm

from .inputs import Input

# import matplotlib as mpl
# import matplotlib.pyplot as plt
# mpl.use("Agg")


class RespPriors:
    def __init__(self, inp: Input, mol2: str | Path, prmtop: str | Path):
        """
        Initialize the RespPriors class.

        Args:
            inp: Input object containing configuration parameters.
            mol2 (str): Path to the mol2 file.
            prmtop (str): Path to the prmtop file.
        """
        # We run RespPriors from within optdir but initialize in the main directory
        # getCharges will determine automatically where we are
        self.sampledir = inp.sampledir
        self.prmtop = prmtop
        self.allEsp = []
        self.allResp = []
        self.getRepeats(mol2)
        self.mode = inp.resppriors
        self.getUnits()

    def getUnits(self):
        """
        Get the number of units from the prmtop file.
        """
        self.units = 0
        inResidues = False
        almostInResidues = False
        print(os.getcwd())
        with open(self.prmtop, "r") as f:
            for line in f.readlines():
                if "%FLAG" in line:
                    inResidues = False
                if inResidues:
                    self.units += line.count(self.resName)
                if "%FORMAT" in line and almostInResidues:
                    almostInResidues = False
                    inResidues = True
                if "%FLAG RESIDUE_LABEL" in line:
                    almostInResidues = True

    def findRepeatIndex(self, idx: int) -> int:
        """
        Find the repeat index for a given atom index.

        Args:
            idx (int): Atom index.

        Returns:
            int: Repeat index, or -1 if not found.
        """
        loc = -1
        for j, idxs in enumerate(self.repeats):
            for i in idxs:
                if i == idx:
                    loc = j
        return loc

    def readCharges(self, lines: List[str]) -> Tuple[List[float], List[float]]:
        """
        Read ESP and RESP charges from TeraChem output lines.

        Args:
            lines (list): List of lines from TeraChem output.

        Returns:
            tuple: ESP charges and RESP charges.

        Raises:
            RuntimeError: If no charges are found in the lines.
        """
        esp = []
        resp = []
        inEsp = False
        inResp = False
        almostInEsp = False
        almostInResp = False
        natoms = 0
        for line in lines:
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

    def getCharges(self, i: int):
        """
        Get charges from TeraChem output files for a given iteration.

        Args:
            i (int): Iteration number.

        Raises:
            RuntimeError: If no training directory is found.
        """
        cycleDir = f"{i}_cycle_{i}"
        path = self.sampledir.absolute() / cycleDir
        trainDir = ""
        for f in path.iterdir():
            if f.name.startswith("train") and (path / f).is_dir():
                trainDir = path / f
        if trainDir == "":
            raise RuntimeError(
                f"No training directory found in {os.path.join(self.sampledir, cycleDir)}"
            )
        outs = []
        for f in os.listdir(trainDir):
            if f.startswith("tc") and f.endswith(".out"):
                outs.append(f)

        outs = sorted(outs)
        for out in outs:
            with open(trainDir / out, "r") as tcout:
                lines = tcout.readlines()
            esp, resp = self.readCharges(lines)
            self.allEsp.append(esp)
            self.allResp.append(resp)

    def computeChargeDistributions(self):
        """
        Compute charge distributions for ESP and RESP charges.

        Raises:
            ValueError: If the number of atoms is not divisible by the number of units.
        """
        espArray = np.array(self.allEsp, dtype=np.float32, copy=True)
        respArray = np.array(self.allResp, dtype=np.float32, copy=True)
        espArray.shape[1]
        unitLength = len(self.repeats)
        if espArray.shape[1] % self.units != 0:
            raise ValueError(
                "Must have an integer number of molecules to be fitted in the RESP calculation"
            )
        self.espMeans = np.zeros(unitLength)
        self.espStdevs = np.ones(unitLength) * -1
        self.respMeans = np.zeros(unitLength)
        self.respStdevs = np.ones(unitLength) * -1
        for i, idxs in enumerate(self.repeats):
            if len(idxs) == 0:
                continue
            elif len(idxs) == 1:
                self.espMeans[i], self.espStdevs[i] = norm.fit(
                    espArray[:, i::unitLength].flatten()
                )
                self.respMeans[i], self.respStdevs[i] = norm.fit(
                    respArray[:, i::unitLength].flatten()
                )
            else:
                espList = []
                respList = []
                for idx in idxs:
                    espList += list(espArray[:, idx::unitLength].flatten())
                    respList += list(respArray[:, idx::unitLength].flatten())
                self.espMeans[i], self.espStdevs[i] = norm.fit(
                    np.asarray(espList, dtype=np.float32)
                )
                self.respMeans[i], self.respStdevs[i] = norm.fit(
                    np.asarray(respList, dtype=np.float32)
                )

        # fig, ax = plt.subplots()
        # plt.hist(allCharges[:,2::unitLength].flatten(),density=True,bins=12,edgecolor="black")
        # xmin, xmax = plt.xlim()
        # x = np.linspace(xmin, xmax, 100)
        # y = norm.pdf(x, means[2], stdevs[2])
        # plt.plot(x, y)
        # fig.set_dpi(200)
        # plt.savefig("RESP_charge_dist.png",bbox_inches="tight")

    def computePriors(self, cdf: float = 0.95) -> np.ndarray:
        """
        Compute priors based on the specified mode.

        Args:
            cdf (float, optional): Cumulative distribution function value. Defaults to 0.95.

        Returns:
            numpy.ndarray: Computed priors.
        """
        # Priors determined by distribution of RESP charges
        if self.mode == 1:
            return self.respStdevs
        # Priors determined by a Gaussian centered on RESP charges with a standard deviation such that
        # the given value of the cdf occurs at the same charge for this distribution and the ESP distribution
        elif self.mode == 2:
            ppf = norm.ppf(cdf)
            return np.abs(self.espMeans - self.respMeans) / ppf + self.espStdevs

    def setMol2Charges(self, charges: List[float], mol2: str | Path):
        """
        Set charges in the mol2 file.

        Args:
            charges (list): List of charges to set.
            mol2 (str): Path to the mol2 file.
        """
        inCharges = False
        with open(mol2, "r") as inF:
            with open("temp.txt", "w") as outF:
                for line in inF.readlines():
                    if "@<TRIPOS>BOND" in line:
                        inCharges = False
                    if inCharges:
                        idx = int(line.split()[0]) - 1
                        i = self.findRepeatIndex(idx)
                        if i == -1:
                            outF.write(line)
                        else:
                            line2 = line.replace(line.split()[8], str(charges[i]))
                            outF.write(line2)
                    else:
                        outF.write(line)
                    if "@<TRIPOS>ATOM" in line:
                        inCharges = True
        os.rename("temp.txt", mol2)

    def setPriors(self, priors: List[float], inputFile: str | Path):
        """
        Set priors in the input file.

        Args:
            priors (list): List of priors to set.
            inputFile (str): Path to the input file.
        """
        inOptions = False
        inPriors = False
        added = False
        with open(inputFile, "r") as inF:
            with open("temp.txt", "w") as outF:
                for line in inF.readlines():
                    if inOptions and len(line.split()) > 0:
                        if inPriors:
                            if "COUL" in line:
                                continue
                            if line.split()[0] == "/priors":
                                for i in range(priors.shape[0]):
                                    if priors[i] > 0:
                                        outF.write(
                                            f"COUL:{self.resName}-{str(i+1)} : {str(priors[i])}\n"
                                        )
                                inPriors = False
                                added = True
                        if line.split()[0] == "priors":
                            inPriors = True
                        if line.split()[0] == "$end":
                            if not added:
                                outF.write("priors\n")
                                for i in range(priors.shape[0]):
                                    if priors[i] > 0:
                                        outF.write(
                                            f"COUL:{self.resName}-{str(i+1)} : {str(priors[i])}\n"
                                        )
                                outF.write("/priors\n")
                            inOptions = False
                        outF.write(line)
                    elif "$options" in line:
                        inOptions = True
                        outF.write(line)
                    else:
                        outF.write(line)
        os.rename("temp.txt", inputFile)

    def getRepeats(self, mol2: str | Path):
        """
        Get repeats from the mol2 file.

        Args:
            mol2 (str): Path to the mol2 file.
        """
        inCharges = False
        self.repeats = []
        with open(mol2, "r") as f:
            f.readline()
            self.resName = f.readline().split()[0]
            line = f.readline()
            unitLength = int(line.split()[0])
        for i in range(unitLength):
            self.repeats.append([])
        with open(mol2, "r") as f:
            for line in f.readlines():
                if "@<TRIPOS>BOND" in line:
                    inCharges = False
                if inCharges:
                    splitLine = line.split()
                    if len(splitLine) >= 11:
                        if splitLine[10] == "RPT":
                            index = int(splitLine[12].split("-")[-1]) - 1
                            self.repeats[index].append(int(splitLine[0]) - 1)
                        elif splitLine[10] == "PRM":
                            index = int(splitLine[0]) - 1
                            self.repeats[index].append(index)
                if "@<TRIPOS>ATOM" in line:
                    inCharges = True

    def updateRespPriors(self, i: int, mol2: str | Path):
        """
        Update RESP priors for a given iteration.

        Args:
            i (int): Iteration number.
            mol2 (str): Path to the mol2 file.
        """
        self.getCharges(i)
        self.computeChargeDistributions()
        priors = self.computePriors()
        inputFile = f"opt_{str(i)}.in"
        self.setPriors(priors, inputFile)
        self.setMol2Charges(self.respMeans, mol2)
