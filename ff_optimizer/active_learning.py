import sander
from .model import Model, AbstractModel
from .utils import readPDB, writePDB
from shutil import copytree
from random import sample
import os
import numpy as np
#import multiprocessing
#from time import perf_counter

def sanderEnergyForce(geometry):
    sander.set_positions(geometry)
    energy, force = sander.energy_forces()
    return [energy.tot, force]
    
class ActiveLearningModel(AbstractModel):

    def __init__(self, args):
        self.home = os.getcwd()
        self.nmodels = args.activeLearning
        self.models = []
        # account for the fact that models are in self.home/model_{i}
        args.dynamicsdir = os.path.join("..",args.dynamicsdir)
        for i in range(1,self.nmodels+1):
            folder = f"model_{str(i)}"
            if not os.path.isdir(folder):
                os.mkdir(folder)
            copytree(args.optdir, os.path.join(folder,args.optdir))
            copytree(args.sampledir, os.path.join(folder,args.sampledir))
            os.chdir(folder)
            self.models.append(Model(args))
            os.chdir(self.home)
        self.restartCycle = min([model.restartCycle for model in self.models])
        self.templatePdb = os.path.join(args.optdir, "conf.pdb")

    def initialCycle(self):
        for i in range(1,self.nmodels+1):
            os.chdir(f"model_{str(i)}")
            models[i].initialCycle()
            os.chdir(home)

    def doMMSampling(self, i):
        for i in range(1,self.nmodels+1):
            os.chdir(f"model_{str(i)}")
            models[i].doMMSampling(i)
            os.chdir(home)
        self.doActiveLearning(i)

    def doQMCalculations(self, i):
        for i in range(1,self.nmodels+1):
            os.chdir(f"model_{str(i)}")
            models[i].doQMCalculations(i)
            os.chdir(home)

    def doParameterOptimization(self, i):
        # serial code
        for i in range(1,self.nmodels+1):
            os.chdir(f"model_{str(i)}")
            models[i].doParameterOptimization(i)
            os.chdir(home)

    def doActiveLearning(self, i):
        sampleDirs = [os.path.join(f"model_{str(k+1)}", self.models[k].sampledir, f"{str(i)}_cycle_{str(i)}") for k in range(self.nmodels)]
        prmtops = [os.path.join(sampleDirs[k], self.models[k].mmEngine.prmtop) for k in range(self.nmodels)]
        dirs = ["train"]
        for k in range(1,self.nmodels):
            dirs.append(f"valid_{str(k)}")
        for k in range(self.nmodels):
            geometries = self.collectGeometries(i, k)
            energies, forces = self.computeAll(geometries, prmtops)
            geometryIndices = self.chooseGeometries(energies, forces)
            newGeometries = [geometries[i] for i in geometryIndices]
            for j in range(self.nmodels):
                self.writeGeoms(newGeometries, os.path.join(sampleDirs[j], dirs[j-k]))

    def computeEnergyForce(self, geometries, prmtop):
    
        with sander.setup(prmtop, geometries[0], None, sander.gas_input()):
            #parallel code
            #paraStart = perf_counter()
            #with multiprocessing.Pool(None) as pool:
            #    results = pool.map(sanderEnergyForce,geometries,2)
            #paraEnd = perf_counter()
    
            #serial code
            results = []
            for i in range(len(geometries)):
                results.append(sanderEnergyForce(geometries[i]))
            #serEnd = perf_counter()
            #paraTime = paraEnd - paraStart
            #serTime = serEnd - paraEnd
            #print(paraTime, serTime)
    
        return results
    
    def collectGeometries(self, i, j):
        geometries = []
        dirs = ["train"]
        for k in range(1,self.nmodels):
            dirs.append(f"valid_{str(k)}")

        for k in range(self.nmodels):
            sampleDir = os.path.join(f"model_{str(k+1)}",self.models[k].sampledir, f"{str(i)}_cycle_{str(i)}", dirs[k-j])
            for f in os.listdir(sampleDir):
                if ".pdb" in f:
                    geometries.append(readPDB(os.path.join(sampleDir,f)))
                    #os.remove(f)
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
    
    def chooseGeometries(self, energies, forces):
        geomsNeeded = int(energies.shape[1] / self.nmodels)
        print(geomsNeeded)
        if self.nmodels == 2:
            #energySpread = np.abs(energies[0,:] - energies[1,:])
            forceSpread = np.linalg.norm(forces[0,:,:] - forces[1,:,:], axis=1)
        else:
            #energySpread = np.std(energies, axis=1)
            forceSpread = np.linalg.norm(np.std(forces, axis=0), axis=1)
        # For now, we only consider spread in forces
        indices = list(np.argsort(forceSpread))
        # Grab half of geometries by spread criterion; other half randomly
        newGeoms = [indices[-k] for k in range(1,int((geomsNeeded+1)/2) + 1)]
        newGeoms += sample(indices[:energies.shape[1]-len(newGeoms)], geomsNeeded-len(newGeoms))
        return newGeoms

    def writeGeoms(self, geometries, dest):
        for i in range(1,len(geometries)+1):
            writePDB(geometries[i-1], os.path.join(dest, f"{str(i)}.pdb"), self.templatePdb) 
