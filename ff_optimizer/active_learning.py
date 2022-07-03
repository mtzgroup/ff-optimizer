import sander
from .model import Model, AbstractModel
from .utils import readPDB
#import multiprocessing
#from time import perf_counter

def sanderEnergyForce(geometry):
    sander.set_positions(geometry)
    energy, force = sander.energy_forces()
    return [energy.tot, force]
    
class ActiveLearningModel(AbstractModel):

    def __init__(self, args):
        self.home = os.getcwd()
        self.nmodels = args.nmodels
        self.models = []
        args.dynamicsdir = os.path.join("..",args.dynamicsdir)
        for i in range(self.nmodels):
            folder = f"model_{str(i)}"
            if not os.path.isdir(folder):
                os.mkdir(folder)
            os.cwd(folder)
            models.append(Model(args))
            os.cwd(home)

    def initialCycle(self):
        for i in range(self.nmodels):
            os.cwd(f"model_{str(i)}")
            models[i].initialCycle()
            os.cwd(home)

    def doMMSampling(self, i):
        for i in range(self.nmodels):
            os.cwd(f"model_{str(i)}")
            models[i].doMMSampling(i)
            os.cwd(home)

    def doQMCalculations(self, i):
        for i in range(self.nmodels):
            os.cwd(f"model_{str(i)}")
            models[i].doQMCalculations(i)
            os.cwd(home)

    def doParameterOptimization(self, i):
        pass
    
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
    
    def collectGeometries(self, i):
        geometries = []
        for i in range(self.nmodels):
            sampledir = os.path.join(f"model_{str(i)}",self.models[i].sampledir, f"{str(i)}_cycle_{str(i)}")
                #if ".pdb" in f:
                    #geometries.append(readPDB(os.path.join(directory,f)))
        return geometries
                    
    
    def computeAll(geometries, prmtop):
        energies = []
        forces = []
        for i in range(self.models):
            results = self.computeEnergyForce(geometries, prmtop)
            for result in results:
                energies.append(result[0])
                forces.append(result[1])
    
        return energies, forces
    
    #def chooseGeometries(energies, forces):
    #    if self.models == 2:
    #        energySpread = np.abs(energies[0,:] - energies[1,:])
    #        forceSpread = np.linalg.norm(forces[0,:,:] - forces[1,:,:], axis=1)
    #    else:
    #        energySpread = np.std(energies, axis=1)
    #        forceSpread = np.linalg.norm(np.std(forces, axis=2), axis=1)
    #    pass
    #
