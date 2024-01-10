import os
from multiprocessing import Pool
from random import sample
from shutil import copytree, rmtree

import numpy as np
try:
    import sander
except:
    pass

from . import utils
from .model import AbstractModel, Model

# from time import perf_counter


def mmSampleStar(inp):
    mmSampleParallel(*inp)


def mmSampleParallel(model, folder, i):
    os.chdir(folder)
    model.doMMSampling(i)


def paramOptStar(inp):
    model = paramOptParallel(*inp)
    return model


def paramOptParallel(model, folder, i):
    os.chdir(folder)
    model.doParameterOptimization(i)
    return model


def sanderEnergyForce(geometry):
    sander.set_positions(geometry)
    energy, force = sander.energy_forces()
    return [energy.tot, force]


class ActiveLearningModel(AbstractModel):
    def __init__(self, inp):
        self.setParams(inp)
        os.chdir(inp.optdir)
        os.system(f"tleap -f setup.leap > leap.out")
        self.prmtop = None
        for f in os.listdir():
            if f.endswith(".prmtop"):
                self.prmtop = f
        if self.prmtop == None:
            raise RuntimeError("setup.leap file did not produce a .prmtop file")
        os.chdir(self.home)
        for i in range(1, self.nmodels + 1):
            folder = f"model_{str(i)}"
            if os.path.isdir(folder):
                if not inp.restart:
                    rmtree(folder)
            else:
                os.mkdir(folder)
            if not os.path.isdir(os.path.join(folder, inp.optdir)):
                copytree(inp.optdir, os.path.join(folder, inp.optdir))
            if not os.path.isdir(os.path.join(folder, inp.sampledir)):
                copytree(inp.sampledir, os.path.join(folder, inp.sampledir))
            os.chdir(folder)
            # Want to sample multiple validation sets, but don't need to evaluate them all
            # So we need to compute restart cycle manually from this level
            # There must be a better way to do this
            self.models.append(Model(inp))
            self.models[-1].optEngine.nvalids = 1
            if inp.restart:
                self.models[-1].optEngine.determineRestart()
                self.models[-1].restartCycle = self.models[-1].optEngine.restartCycle
            os.chdir(self.home)
        if inp.restart:
            self.restartCycle = min([m.restartCycle for m in self.models])
        else:
            self.restartCycle = -1
        self.symbols = None

    def setParams(self, inp):
        inp.nvalids = inp.activelearning - 1
        self.home = os.getcwd()
        self.nmodels = inp.activelearning
        self.nthreads = min(os.cpu_count(), self.nmodels)
        self.models = []
        self.trainGeometries = None
        self.validGeometries = None
        # account for the fact that models are in self.home/model_{i}
        self.dynamicsdir = os.path.join("..", inp.dynamicsdir)

    def initialCycle(self):
        for i in range(1, self.nmodels + 1):
            os.chdir(f"model_{str(i)}")
            self.models[i - 1].initialCycle()
            os.chdir(self.home)

    def doMMSampling(self, i):
        ## serial code
        # for j in range(1, self.nmodels + 1):
        #    os.chdir(f"model_{str(j)}")
        #    models[j-1].doMMSampling(i)
        #    os.chdir(home)
        tasks = [(self.models[j], f"model_{str(j+1)}", i) for j in range(self.nmodels)]
        with Pool(self.nthreads) as p:
            p.map(mmSampleStar, tasks)
        self.doActiveLearning(i)

    def doQMCalculations(self, i):
        for j in range(1, self.nmodels + 1):
            os.chdir(f"model_{str(j)}")
            self.models[j - 1].doQMCalculations(i)
            os.chdir(self.home)

    def doParameterOptimization(self, i):
        # serial code
        # optResults = []
        # for i in range(1, self.nmodels + 1):
        #    os.chdir(f"model_{str(i)}")
        #    result = models[i].doParameterOptimization(i)
        #    optResults.append(result)
        #    os.chdir(home)
        ## for now we arbitrarily print only the first model's info
        # return optResults[0]
        tasks = [(self.models[j], f"model_{str(j+1)}", i) for j in range(self.nmodels)]
        with Pool(self.nthreads) as p:
            self.models = p.map(paramOptStar, tasks)
        # for now we arbitrarily print only the first model's info
        self.optResults = self.models[0].optResults

    def doActiveLearning(self, i):
        sampleDirs = [
            os.path.join(
                f"model_{str(k+1)}",
                self.models[k].sampledir,
                f"{str(i)}_cycle_{str(i)}",
            )
            for k in range(self.nmodels)
        ]

        prmtops = [
            os.path.join(sampleDirs[k], self.prmtop) for k in range(self.nmodels)
        ]
        dirs = ["train"]
        for k in range(1, self.nmodels):
            dirs.append(f"valid_{str(k)}")
        for k in range(self.nmodels):
            geometries = self.collectAll(i, k)
            energies, forces = self.computeAll(geometries, prmtops)
            if self.trainGeometries == self.validGeometries:
                geometryIndices = self.chooseGeometries(
                    energies, forces, self.trainGeometries
                )
                newGeometries = [geometries[i] for i in geometryIndices]
                for j in range(self.nmodels):
                    self.writeGeoms(
                        newGeometries, os.path.join(sampleDirs[j], dirs[j - k])
                    )
            else:
                geometryIndices = self.chooseGeometries(
                    energies, forces, self.trainGeometries
                )
                trainGeometries = [geometries[i] for i in geometryIndices]
                geometryIndices = self.chooseGeometries(
                    energies, forces, self.validGeometries
                )
                validGeometries = [geometries[i] for i in geometryIndices]
                for j in range(self.nmodels):
                    if j == k:
                        self.writeGeoms(
                            trainGeometries, os.path.join(sampleDirs[j], dirs[j - k])
                        )
                    else:
                        self.writeGeoms(
                            validGeometries, os.path.join(sampleDirs[j], dirs[j - k])
                        )

    def computeEnergyForce(self, geometries, prmtop):

        with sander.setup(prmtop, geometries[0], None, sander.gas_input()):
            # parallel code
            # paraStart = perf_counter()
            # with multiprocessing.Pool(None) as pool:
            #    results = pool.map(sanderEnergyForce,geometries,2)
            # paraEnd = perf_counter()

            # serial code
            results = []
            for i in range(len(geometries)):
                results.append(sanderEnergyForce(geometries[i]))
            # serEnd = perf_counter()
            # paraTime = paraEnd - paraStart
            # serTime = serEnd - paraEnd
            # print(paraTime, serTime)

        return results

    def collectAll(self, i, j):
        geometries = []
        dirs = ["train"]
        for k in range(1, self.nmodels):
            dirs.append(f"valid_{str(k)}")

        for k in range(self.nmodels):
            sampleDir = os.path.join(
                f"model_{str(k+1)}",
                self.models[k].sampledir,
                f"{str(i)}_cycle_{str(i)}",
                dirs[k - j],
            )
            # If restarting, we redo the .nc processing (in case pooling
            # was already done and the original xyzs overwritten)
            if i == self.restartCycle:
                currentDir = os.getcwd()
                os.chdir(sampleDir)
                ncs = []
                prmtop = None
                for f in os.listdir(".."):
                    if f.endswith(".prmtop"):
                        prmtop = f
                for f in os.listdir():
                    if f.endswith(".nc"):
                        ncs.append(f)
                self.symbols = utils.getSymbolsFromPrmtop(os.path.join("..", prmtop))
                for nc in ncs:
                    utils.convertNCtoXYZs(nc, self.symbols)
                os.chdir(currentDir)
            tempGeoms = self.collectGeometries(sampleDir)
            if self.trainGeometries is None and k == j:
                self.trainGeometries = len(tempGeoms)
            elif self.validGeometries is None and k != j:
                self.validGeometries = len(tempGeoms)
            geometries += tempGeoms

        return geometries

    def collectGeometries(self, sampleDir):
        geometries = []
        for f in os.listdir(sampleDir):
            if ".xyz" in f:
                fPath = os.path.join(sampleDir, f)
                if self.symbols is None:
                    geometry, self.symbols = utils.readXYZ(fPath, readSymbols=True)
                    geometries.append(geometry)
                else:
                    geometries.append(utils.readXYZ(fPath))
                os.remove(fPath)
        return geometries

    def computeAll(self, geometries, prmtops):
        energies = []
        forces = []
        allEnergies = []
        allForces = []
        for prmtop in prmtops:
            results = self.computeEnergyForce(geometries, prmtop)
            for result in results:
                energies.append(result[0])
                forces.append(result[1])
            allEnergies.append(energies)
            allForces.append(forces)
            energies = []
            forces = []

        allEnergies = np.asarray(allEnergies, dtype=np.float32)
        allForces = np.asarray(allForces, dtype=np.float32)

        return allEnergies, allForces

    def chooseGeometries(self, energies, forces, geomsNeeded):
        if self.nmodels == 2:
            # energySpread = np.abs(energies[0,:] - energies[1,:])
            forceSpread = np.linalg.norm(forces[0, :, :] - forces[1, :, :], axis=1)
        else:
            # energySpread = np.std(energies, axis=1)
            forceSpread = np.linalg.norm(np.std(forces, axis=0), axis=1)
        # For now, we only consider spread in forces
        indices = list(np.argsort(forceSpread))
        # Grab half of geometries by spread criterion; other half randomly
        newGeoms = [indices[-k] for k in range(1, int((geomsNeeded + 1) / 2) + 1)]
        newGeoms += sample(
            indices[: energies.shape[1] - len(newGeoms)], geomsNeeded - len(newGeoms)
        )
        return newGeoms

    def writeGeoms(self, geometries, dest):
        for i in range(1, len(geometries) + 1):
            utils.writeXYZ(
                geometries[i - 1], self.symbols, os.path.join(dest, f"{str(i)}.xyz")
            )
