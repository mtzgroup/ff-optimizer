import os
import subprocess

class QMEngine():

    def __init__(self,inputFile:str,backupInputFile:str):
        self.inputSettings = self.readInputFile(inputFile)
        self.backupInputSettings = self.readInputFile(backupInputFile)

    def readGradFromTCout(self,outFile:str):
        inGradient = False
        gradCounter = 0
        molSize = 0
        grads = []
        energy = 0
        finished = False
        with open(outFile, "r") as f:
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
                if "Job finished" in line:
                    finished = True
        if not finished:
            print("File %s did not complete calculation" % filename)
            grads = -1
        return energy, grads

    def readInputFile(self,inputFile:str):
        settings = []
        with open(inputFile,'r') as f:
            for line in f.readlines():
                splitLine = line.split()
                # ignore commented lines, the coordinates specification, and run gradient
                if len(splitLine) > 0:
                    if splitLine[0][0] != '#' and splitLine[0].lower() != 'coordinates' and splitLine[0].lower() != 'run':
                        setting = []
                        # for options > 2 tokens long
                        for token in splitLine:
                            # ignore comments at the end of the line
                            if token[0] != '#':
                                setting.append(token)
                            else:
                                break
                        settings.append(setting)
        return settings

    def writeInputFile(self, settings:list, coordinates:str, fileName:str):
        with open(fileName,'w') as f:
            f.write(f"coordinates {coordinates}\n")
            f.write(f"run gradient\n")
            for setting in settings:
                # reassemble setting (if > 2 tokens long)
                for token in setting:
                    f.write(f"{token} ")
                f.write("\n")

    def readPDB(self,pdb:str):
        coords = []
        with open(pdb, 'r') as f:
            for line in f.readlines():
                if line.startswith("ATOM") or line.startswith("HETATM"):
                    splitLine = line.split()
                    coords.append(splitLine[5])
                    coords.append(splitLine[6])
                    coords.append(splitLine[7])
        return coords

    def writeFBdata(self,energies:list, grads:list, coords:list):
        with open("qdata.txt",'w') as f:
            for i in range(len(energies)):
                f.write(f"JOB {i+1}\n")
                coordLine = "COORDS "
                for coord in coords[i]:
                    coordLine = coordLine + str(coord) + " "
                gradLine = "FORCES "
                for grad in grads[i]:
                    gradLine = gradLine + str(grad) + " "
                f.write(coordLine + "\n")
                f.write(f"ENERGY {str(energies[i])}\n")
                f.write(gradLine + "\n\n")

        with open("all.mdcrd",'w') as f:
            f.write("Converted from pdb by QMEngine\n")
            for i in range(len(energies)):
                tokenCounter = 1
                for coord in coords[i]:
                    f.write("%8.3f" % float(coord))
                    if tokenCounter == 10:
                        f.write("\n")
                        tokenCounter = 1
                    else:
                        tokenCounter += 1
                if tokenCounter != 1:
                    f.write("\n")

    def getQMRefData(self, pdbs:list, calcDir:str):
        pass

    def restart(self, calcDir:str):
        pdbs = []
        os.chdir(calcDir)
        for f in os.listdir():
            if f.endswith(".pdb"):
                name = f.split('.')[0]
                energy, status = self.readGradFromTCout(f"tc_{name}.out")
                if status == -1:
                    pdbs.append(f)
        self.getQMRefData(pdbs, ".") 

class SbatchEngine(QMEngine):

    def __init__(self, inputFile:str, backupInputFile:str, sbatchFile:str, user:str):
        self.user = user
        self.readSbatchFile(sbatchFile)
        QMEngine.__init__(self,inputFile,backupInputFile)

    def readSbatchFile(self, sbatchFile:str):
        tcVersion = "TeraChem/2021.02-intel-2017.8.262-CUDA-9.0.176"
        sbatchOptions = []
        with open(sbatchFile,'r') as f:
            for line in f.readlines():
                if line.startswith("#SBATCH"):
                    # Grab option from #SBATCH option=value
                    option = line.split()[1].split('=')[0]
                    # exclude options we're setting ourselves
                    if option != '--fin' and option != '--fout' and option != '-J':
                        sbatchOptions.append(line)
                # get TC version from ml line
                if "TeraChem" in line and not line.startswith('#'):  
                    if "ml" in line or "module load" in line:
                        for token in line:
                            if token.startswith("TeraChem"):
                                tcVersion = token
        self.tcVersion = tcVersion
        self.sbatchOptions = sbatchOptions

    def writeSbatchFile(self, index:str, fileName:str):
        index = str(index)
        with open(fileName,'w') as f:
            f.write("#!/bin/bash\n\n")
            f.write(f"#SBATCH -J FB_ref_gradient_{index}\n")
            f.write(f"#SBATCH --fin=tc_{index}.in,{index}.pdb,tc_{index}_backup.in\n")
            f.write(f"#SBATCH --fout=tc_{index}.out\n")
            for line in self.sbatchOptions:
                f.write(line)
            f.write("cd $SCRATCH\n\n")
            f.write(f"module load {self.tcVersion}\n")
            f.write(f"terachem tc_{index}.in > tc_{index}.out\n")
            f.write(f"if [ $(grep -c \"Job finished\" tc_{index}.out) -ne 1 ]\n")
            f.write("then\n")
            f.write(f"  terachem tc_{index}_backup.in > tc_{index}.out\n")
            f.write("fi\n")

    def slurmCommand(self,command:list):
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
            raise RuntimeError(f"Slurm command {command} failed")
        return output

    def getQMRefData(self, pdbs:list, calcDir:str):
        os.chdir(calcDir)
        jobIDs = []
        for pdb in pdbs:
            name = pdb.split('.')[0]
            super().writeInputFile(self.inputSettings, pdb, f"tc_{name}.in")
            super().writeInputFile(self.backupInputSettings, pdb, f"tc_{name}.in")
            self.writeSbatchFile(name,f"sbatch_{name}.sh")
            job = self.slurmCommand(["sbatch",f"sbatch_{name}.sh"])
            jobIDs.append(job.split()[3])

        while len(jobIDs) > 0:
            runningIDs = []
            sleep(10)
            status = slurmCommand(["squeue","-o","%.12i","-u",self.user])
            for runningID in status.replace(b" ", b"").replace(b'"', b"").split(b"\n")[1:]:
                for submitID in jobIDs:
                    if runningID == submitID:
                        runningIDs.append(runningID)
            jobIDs = runningIDs

        allPdbs = []
        for f in os.listdir():
            if f.endswith(".pdb"):
                allPdbs.append(f)
        energies = []
        grads = []
        coords = []
        for pdb in allPdbs:
            name = pdb.split('.')[0]
            energy, grad = super().readGradFromTCout(f"tc_{name}.out")
            if grad == -1:
                raise RuntimeError("TC job " + os.path.join(calcDir,f"tc_{name}.out") + " failed")
            coord = super().readPDB(pdb)
            energies.append(energy)
            grads.append(grad)
            coords.append(coord)
        super().writeFBdata(energies,grads,coords)

class debugEngine(QMEngine):

    def __init__(self, inputFile, backupInputFile):
        super().__init__(inputFile, backupInputFile)

    def getQMRefData(self, pdbs:list, calcDir:str):
        os.chdir(calcDir)
        energies = []
        grads = []
        coords = []
        for pdb in pdbs:
            name = pdb.split('.')[0]
            super().writeInputFile(self.inputSettings, pdb, f"tc_{name}.in")
            os.system(f"terachem tc_{name}.in > tc_{name}.out")
            energy, grad = super().readGradFromTCout(f"tc_{name}.out")
            if grad == -1:  
                super().writeInputFile(self.backupInputSettings, pdb, f"tc_{name}_backup.in")
                os.system(f"terachem tc_{name}_backup.in > tc_{name}.out")
                energy, grad = super().readGradFromTCout(f"tc_{name}.out")
                if grad == -1:
                    raise RuntimeError("TC job " + os.path.join(calcDir,f"tc_{name}.out") + " failed")
            coord = super().readPDB(pdb)
            energies.append(energy)
            grads.append(grad)
            coords.append(coord)
        super().writeFBdata(energies,grads,coords)

class TCCloudEngine(QMEngine):
    
    def __init__(self, inputFile:str, backInputFile:str, batchSize = 10):
        self.batchSize = batchSize
        super().__init__(inputFile, backupInputFile)
        
    def computeBatch(self, pdbs:list, calcDir:str):
        pass
