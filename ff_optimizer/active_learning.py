import sander
#import multiprocessing
#from time import perf_counter

class ActiveLearning:

    def __init__(self, options):
        self.models = options['models']
        pass
    
    def sanderEnergyForce(geometry):
        print(type(geometry))
        sander.set_positions(geometry)
        energy, force = sander.energy_forces()
        return [energy.tot, force]
    
    def computeEnergyForce(geometries, prmtop):
    
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
    
    def collectGeometries(dirs):
        geometries = []
        for directory in dirs:
            for f in os.listdir(directory):
                if ".pdb" in f:
                    geometries.append(readPDB(os.path.join(directory,f)))
        return geometries
                    
    
    def computeAll(geometries, prmtop):
        energies = []
        forces = []
        for i in range(self.models):
            results = computeEnergyForce(geometries, prmtop)
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
