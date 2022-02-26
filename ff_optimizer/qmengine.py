import os
import subprocess
from tccloud import TCClient
from tccloud.models import AtomicInput, Molecule
from . import utils, units
from math import ceil
from time import sleep
import numpy as np

class QMEngine():

    def __init__(self,inputFile:str,backupInputFile:str,doResp=False):
        self.inputSettings = self.readInputFile(inputFile)
        self.backupInputSettings = self.readInputFile(backupInputFile)
        self.doResp = doResp

    def readInputFile(self,inputFile:str):
        settings = []
        with open(inputFile,'r') as f:
            for line in f.readlines():
                splitLine = line.split()
                # ignore commented lines, the coordinates specification, and run gradient
                if len(splitLine) > 0:
                    if splitLine[0][0] != '#' and splitLine[0].lower() != 'coordinates' and splitLine[0].lower() != 'run' and splitLine[0].lower() != "resp":
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
            if self.doResp:
                f.write(f"resp yes\n")
            for setting in settings:
                # reassemble setting (if > 2 tokens long)
                for token in setting:
                    f.write(f"{token} ")
                f.write("\n")


    def writeFBdata(self, energies:list, grads:list, coords:list, espXYZs = None, esps = None):
        with open("qdata.txt",'w') as f:
            for i in range(len(energies)):
                f.write(f"JOB {i+1}\n")
                coordLine = "COORDS "
                for coord in coords[i]:
                    coordLine = coordLine + str(round(float(coord),3)) + " "
                gradLine = "FORCES "
                for grad in grads[i]:
                    gradLine = gradLine + str(round(float(grad),5)) + " "
                f.write(coordLine + "\n")
                f.write(f"ENERGY {str(round(float(energies[i]),6))}\n")
                f.write(gradLine + "\n")
                if espXYZs is not None:
                    espXYZLine = "ESPXYZ "
                    for xyz in espXYZs[i]:
                        espXYZLine = espXYZLine + str(round(float(xyz),6)) + " "
                    f.write(espXYZLine + "\n")
                    espLine = "ESPVAL "
                    for val in esps[i]:
                        espLine = espLine + str(round(float(val),6)) + " "
                    f.write(espLine + "\n\n")
                else:
                    f.write("\n")

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

    def readQMRefData(self):
        coords = []
        energies = []
        grads = []
        if self.doResp:
            espXYZs = []
            esps = []
        else:
            espXYZs = None
            esps = None
        for f in os.listdir():
            if f.endswith(".pdb"):
                coord = utils.readPDB(f)
                name = f.split('.')[0]
                tcOut = f"tc_{name}.out"
                energy, grad = utils.readGradFromTCout(tcOut)
                if type(grad) == int:
                    if grad == -1:
                        raise RuntimeError(f"Terachem job {tcOut} in {os.getcwd()} did not succeed!")
                energies.append(energy)
                grads.append(grad)
                coords.append(coord)
                if self.doResp:
                    espXYZ, esp = utils.readEsp(f"esp_{name}.xyz")
                    espXYZs.append(espXYZ)
                    esps.append(esp)
        return energies, grads, coords, espXYZs, espreadQMRefData()
        super().writeFBdata(energies,grads,coords)

class SbatchEngine(QMEngine):

    def __init__(self, inputFile:str, backupInputFile:str, sbatchFile:str, user:str, doResp=False):
        self.user = user
        self.readSbatchFile(sbatchFile)
        super().__init__(self,inputFile,backupInputFile,doResp)


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
            if self.doResp:
                f.write(f"#SBATCH --fout=tc_{index}.out,esp_{index}.xyz\n")
            else:
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
            if self.doResp:
                f.write(f"mv scr.{index}/esp.xyz esp_{index}.xyz\n")

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
            raise RuntimeError(f"Slurm command {str(command)} failed")
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
            status = self.slurmCommand(["squeue","-o","%.12i","-u",self.user])
            for runningID in status.replace(b" ", b"").replace(b'"', b"").split(b"\n")[1:]:
                for submitID in jobIDs:
                    if runningID == submitID:
                        runningIDs.append(runningID)
            jobIDs = runningIDs

        energies, grads, coords, espXYZs, esps = super().readQMRefData()
        super().writeFBdata(energies,grads,coords, espXYZs, esps)

class DebugEngine(QMEngine):

    def __init__(self, inputFile, backupInputFile, doResp=False):
        super().__init__(inputFile, backupInputFile, doResp)

    def getQMRefData(self, pdbs:list, calcDir:str):
        os.chdir(calcDir)
        energies = []
        grads = []
        coords = []
        for pdb in pdbs:
            name = pdb.split('.')[0]
            super().writeInputFile(self.inputSettings, pdb, f"tc_{name}.in")
            os.system(f"terachem tc_{name}.in > tc_{name}.out")
            energy, grad = utils.readGradFromTCout(f"tc_{name}.out")
            if grad == -1:  
                super().writeInputFile(self.backupInputSettings, pdb, f"tc_{name}_backup.in")
                os.system(f"terachem tc_{name}_backup.in > tc_{name}.out")
            if self.doResp:
                espXYZ, esps = utils.readEsp(f"scr.{name}/esp.xyz")
        energies, grads, coords, espXYZs, esps = super().readQMRefData()
        super().writeFBdata(energies,grads,coords,espXYZs,esps)

class TCCloudEngine(QMEngine):
    
    def __init__(self, inputFile:str, backupInputFile:str, batchSize = None):
        self.batchSize = batchSize
        self.client = TCClient()
        super().__init__(inputFile, backupInputFile)
        self.keywords = {}
        self.backupKeywords = {}
        for setting in self.inputSettings:
            if setting[0] == 'method':
                self.method = setting[1]
            elif setting[0] == 'basis':
                self.basis = setting[1]
            else:
                keyword = ""
                for token in setting[1:]:
                    keyword += f"{token} "
                self.keywords[setting[0]] = keyword
        for setting in self.backupInputSettings:
            if setting[0] == 'method' or setting[0] == 'basis':
                pass
            else:
                keyword = ""
                for token in setting[1:]:
                    keyword += f"{token} "
                self.backupKeywords[setting[0]] = keyword
                    
    def computeBatch(self, atomicInputs:list):
        print("Submitting jobs to TCCloud")
        status = 0
        if self.batchSize == None:
            batchSize = len(atomicInputs)
        else:
            batchSize = self.batchSize
        results = []
        for i in range(ceil(1.0*len(atomicInputs)/batchSize)):
            resultsBatch = []
            atomicInputsBatch = []
            for j in range(batchSize):
                if i * batchSize + j < len(atomicInputs):
                    atomicInputsBatch.append(atomicInputs[i * batchSize + j])
            try:
                futureResultBatch = self.client.compute(atomicInputsBatch,engine="terachem_pbs")
                resultsBatch = futureResultBatch.get()
            except:
                self.batchSize = int(batchSize/2)
                print(f"Submission failed; resubmitting with batch size {str(self.batchSize)}")
                sleep(30)
                if self.batchSize < 2:
                    status = -1
                    break
                tempStatus, resultsBatch = self.computeBatch(atomicInputsBatch)
                if tempStatus == -1:
                    status = -1
            for result in resultsBatch:
                results.append(result)
        return status, results

    def createAtomicInputs(self, pdbs:list):
        mod = {}
        mod["method"] = self.method
        mod["basis"] = self.basis
        atomicInputs = []
        for pdb in pdbs:
            jobId = pdb.split('.')[0]
            xyz = utils.convertPDBtoXYZ(pdb)
            mol = Molecule.from_file(xyz)
            atomicInput = AtomicInput(molecule=mol,model=mod,driver="gradient",keywords=self.keywords,id=jobId)
            atomicInputs.append(atomicInput)
        return atomicInputs

    def writeResult(self, result):
        np.savetxt(f"tc_{str(result.id)}.out",np.asarray(result.return_result,dtype=np.float32),header=f"TCCloud gradient output file. Energy = {str(result.properties.return_energy)}")
        
    def getQMRefData(self, pdbs:list, calcDir:str):
        atomicInputs = self.createAtomicInputs(pdbs)
        status, results  = self.computeBatch(atomicInputs)
        retryInputs = []
        failedIndices = []
        for i in range(len(results)):
            if results[i].success:
                self.writeResult(results[i]) 
            else:
                retryInput = AtomicInput(molecule=atomicInputs[i].molecule,model=atomicInputs[i].model,driver="gradient",keywords=self.backupKeywords,id=atomicInputs[i].id)
                retryInputs.append(retryInput)
                failedIndices.append(i)
        if status == -1:
            raise RuntimeError("Batch resubmission reached size 1; QM calculations incomplete")
        if len(failedIndices) > 0:
            failedIndex = -1
            status, retryResults = self.computeBatch(retryInputs)
            for result in retryResults:
                if result.success:
                    self.writeResult(result)
                else:
                    failedIndex = result.id
            if status == -1:
                raise RuntimeError("Batch resubmission reached size 1; QM calculations incomplete")
        energies, grads, coords, espXYZs, esps = super().readQMRefData()
        super().writeFBdata(energies,grads,coords,espXYZs,esps) 
