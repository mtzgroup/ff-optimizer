import os
import subprocess
from tccloud import TCClient
from tccloud.models import AtomicInput, Molecule, AtomicResult
from . import utils, units
from math import ceil
from time import sleep
import numpy as np
from qcelemental.util.serialization import json_loads
from qcelemental.models.results import AtomicResultProperties
from qcelemental.models import Provenance

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
                # ignore commented lines, the coordinates specification, run gradient, and resp keywords
                if len(splitLine) > 0:
                    if splitLine[0][0] != '#' and splitLine[0].lower() != 'coordinates' and splitLine[0].lower() != 'run' and splitLine[0].lower() != "resp" and splitLine[0].lower() != "esp_restraint_a" and splitLine[0].lower() != "esp_restraint_b":
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
                f.write("resp yes\n")
                f.write("esp_restraint_a 0\n")
                f.write("esp_restraint_b 0\n")
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

    def writeResult(self, tcOut, pdb):
        energy, grad = utils.readGradFromTCout(tcOut)
        mol = utils.convertPDBtoMolecule(pdb)
        properties = AtomicResultProperties()
        model = {"method" : "see tc.in", "basis" : "see tc.in"}
        provenance = Provenance(**{"creator" : "ffoptimizer (file-based)"})
        name = tcOut.split('.')[0]
        jobId = name.replace('tc_','')
        if type(grad) == int:
            if grad == -1:
                result = AtomicResult(**{"molecule" : mol, "driver" : "gradient", "model" : model, "provenance" : provenance, "properties" : properties, "return_result" : np.zeros(mol.geometry.shape), "success" : False, "id" : jobId})
                with open(f"{name}.json",'w') as f:
                    f.write(result.json())
                return
        properties = AtomicResultProperties(**{"return_energy" : energy})
        result = AtomicResult(**{"molecule" : mol, "driver" : "gradient", "model" : model, "provenance" : provenance, "properties" : properties, "return_result" : grad, "success" : True, "id" : jobId})
        #os.remove(tcOut)
        with open(f"{name}.json",'w') as f:
            f.write(result.json())

    def readResult(self, json):
        with open(json,'r') as f:
            try:
                result = AtomicResult(**json_loads(f.read()))
            except:
                raise RuntimeError(f"Corrupted json {json} could not be read in at {os.getcwd}")
        return result
        
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
                json = f"tc_{name}.json"
                result = self.readResult(json)
                if not result.success:
                    raise RuntimeError(f"Terachem job {json} in {os.getcwd()} did not succeed!")
                energies.append(result.properties.return_energy)
                grads.append(result.return_result.flatten())
                coords.append(coord)
                if self.doResp:
                    espXYZ, esp = utils.readEsp(f"esp_{name}.xyz")
                    espXYZs.append(espXYZ)
                    esps.append(esp)
        return energies, grads, coords, espXYZs, esps

    def restart(self, calcDir:str):
        pdbs = []
        os.chdir(calcDir)
        for f in os.listdir():
            if f.endswith(".pdb"):
                name = f.split('.')[0]
                json = f"tc_{name}.json"
                if os.path.isfile(json):
                    try:
                        with open(json, 'r') as j:
                            result = AtomicResult(**json_loads(j.read()))
                    except:
                        pdbs.append(f)
                        continue
                    if not result.success:
                        pdbs.append(f)
                else:
                    pdbs.append(f)
        self.getQMRefData(pdbs, ".") 

    # This is the function which must be implemented by inherited classes
    def getQMRefData(self, pdbs:list, calcDir:str):
        pass

class SbatchEngine(QMEngine):

    def __init__(self, inputFile:str, backupInputFile:str, sbatchFile:str, user:str, doResp=False):
        self.user = user
        self.readSbatchFile(sbatchFile)
        super().__init__(inputFile,backupInputFile,doResp)

    def readSbatchFile(self, sbatchFile:str):
        tcVersion = None #"TeraChem/2021.02-intel-2017.8.262-CUDA-9.0.176"
        sbatchLines = []
        with open(sbatchFile,'r') as f:
            for line in f.readlines():
                if line.startswith("#SBATCH"):
                    # Grab option from #SBATCH option=value
                    option = line.split()[1].split('=')[0]
                    # exclude options we're setting ourselves
                    if option != '--fin' and option != '--fout' and option != '-J':
                        sbatchLines.append(line)
                    self.optionsEnd = len(sbatchLines)
                # get TC version from ml line
                elif "TeraChem" in line and not line.startswith('#'):  
                    if "ml" in line or "module load" in line:
                        for token in line.split():
                            if token.startswith("TeraChem"):
                                tcVersion = token
                                sbatchLines.append(line)
                                break
                else:
                    sbatchLines.append(line)
        if tcVersion == None:
            sbatchLines.append("module load TeraChem/2021.02-intel-2017.8.262-CUDA-9.0.176\n")
        self.sbatchLines = sbatchLines

    def writeSbatchFile(self, index:str, fileName:str):
        sbatchLines = self.sbatchLines.copy()
        index = str(index)
        if self.doResp:
            sbatchLines.insert(self.optionsEnd,f"#SBATCH --fout=tc_{index}.out,esp_{index}.xyz\n")
        else:
            sbatchLines.insert(self.optionsEnd,f"#SBATCH --fout=tc_{index}.out\n")
        sbatchLines.insert(self.optionsEnd,f"#SBATCH --fin=tc_{index}.in,{index}.pdb,tc_{index}_backup.in\n")
        sbatchLines.insert(self.optionsEnd,f"#SBATCH -J FB_ref_gradient_{index}\n")
        with open(fileName,'w') as f:
            for line in sbatchLines:
                f.write(line)
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
            super().writeInputFile(self.backupInputSettings, pdb, f"tc_{name}_backup.in")
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

        for pdb in pdbs:
            name = pdb.split('.')[0]
            super().writeResult(f"tc_{name}.out",pdb)

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
            super().writeResult(f"tc_{name}.out",pdb)

        energies, grads, coords, espXYZs, esps = super().readQMRefData()
        super().writeFBdata(energies,grads,coords,espXYZs,esps)

class TCCloudEngine(QMEngine):
    

    def __init__(self, inputFile:str, backupInputFile:str, batchSize = None, doResp = False):
        if batchSize is None:
            self.batchSize = 100
        else:
            self.batchSize = batchSize
        self.client = TCClient()
        super().__init__(inputFile, backupInputFile, doResp)
        self.keywords = {}
        self.backupKeywords = {}
        for setting in self.inputSettings:
            if setting[0] == 'method':
                self.method = setting[1]
            elif setting[0] == 'basis':
                self.basis = setting[1]
            else:
                if len(setting) == 2:
                    self.keywords[setting[0]] = setting[1]
                else:
                    keyword = setting[1]
                    for token in setting[2:]:
                        keyword += f" {token}"
                    self.keywords[setting[0]] = keyword
        for setting in self.backupInputSettings:
            if setting[0] == 'method' or setting[0] == 'basis':
                pass
            else:
                if len(setting) == 2:
                    self.backupKeywords[setting[0]] = setting[1]
                else:
                    keyword = setting[1]
                    for token in setting[2:]:
                        keyword += f" {token}"
                    self.backupKeywords[setting[0]] = keyword
        if doResp:
            self.keywords["resp"] = "yes"
            self.keywords["esp_restraint_a"] = "0"
            self.keywords["esp_restraint_b"] = "0"
            self.backupKeywords["resp"] = "yes"
            self.backupKeywords["esp_restraint_a"] = "0"
            self.backupKeywords["esp_restraint_b"] = "0"

                    
    def computeBatch(self, atomicInputs:list):
        status = 0
        # If there are no jobs to run after restart
        if len(atomicInputs) == 0:
            return status, []
        batchSize = min(self.batchSize, len(atomicInputs))
        results = []
        stride = int(len(atomicInputs) / batchSize)
        try:
            # HOW TO RESTART IF CODE FAILS AFTER SUBMISSION
            futureResults = [self.client.compute(atomicInputs[i::stride],engine="terachem_fe") for i in range(stride)]
            resultBatches = [futureResults[i].get() for i in range(stride)]
            for batch in resultBatches:
                for result in batch:
                    results.append(result)
        except Exception as e:
            print(e)
            self.batchSize = int(batchSize/2)
            print(f"Submission failed; resubmitting with batch size {str(self.batchSize)}")
            #sleep(30)
            if self.batchSize < 2:
                status = -1
                return status, results
            tempStatus, results = self.computeBatch(atomicInputs)
            if tempStatus == -1:
                status = -1
        return status, results

    def createAtomicInputs(self, pdbs:list, useBackup=False):
        mod = {}
        mod["method"] = self.method
        mod["basis"] = self.basis
        atomicInputs = []
        if useBackup:
            keywords = self.backupKeywords
        else:
            keywords = self.keywords
        for pdb in sorted(pdbs):
            jobId = pdb.split('.')[0]
            mol = utils.convertPDBtoMolecule(pdb)
            if self.doResp:
                atomicInput = AtomicInput(molecule=mol,model=mod,driver="gradient",keywords=keywords,id=jobId,protocols={"native_files" : "all"},extras={"tcfe:keywords" : {"native_files" : ["esp.xyz"]}})
            else:
                atomicInput = AtomicInput(molecule=mol,model=mod,driver="gradient",keywords=keywords,id=jobId)
            atomicInputs.append(atomicInput)
        return atomicInputs

    def writeResult(self, result):
        with open(f"tc_{str(result.id)}.json",'w') as f:
            f.write(result.json())
        if self.doResp:
            with open(f"esp_{str(result.id)}.xyz",'w') as f:
                try:
                    f.write(result.native_files["esp.xyz"])
                except:
                    print("Job {str(result.id)} in {os.getcwd()} is missing esp file!")
                    pass
        
    def getQMRefData(self, pdbs:list, calcDir:str):
        cwd = os.getcwd()
        os.chdir(calcDir)
        atomicInputs = self.createAtomicInputs(pdbs)
        status, results  = self.computeBatch(atomicInputs)
        retryPdbs = []
        for result in results:
            if result.success:
                self.writeResult(result) 
            else:
                retryPdbs.append(f"{result.input_data['id']}.pdb")
        if status == -1:
            raise RuntimeError("Batch resubmission reached size 1; QM calculations incomplete")
        if len(retryPdbs) > 0:
            failedIndices = []
            retryInputs = self.createAtomicInputs(retryPdbs, useBackup=True)
            batchSize = self.batchSize
            status, retryResults = self.computeBatch(retryInputs)
            self.batchSize = batchSize
            for result in retryResults:
                if result.success:
                    self.writeResult(result)
                else:
                    failedIndices.append(result.input_data['id'])
                    print(f"Job id {result.input_data['id']} suffered a {result.error.error_type}")
                    print(result.error.error_message)
            if len(failedIndices) > 0:
                raise RuntimeError(f"Job ids {str(failedIndices)} in {os.getcwd()} failed twice!")
            if status == -1:
                raise RuntimeError("Batch resubmission reached size 1; QM calculations incomplete")
        energies, grads, coords, espXYZs, esps = super().readQMRefData()
        super().writeFBdata(energies,grads,coords,espXYZs,esps) 
        os.chdir(cwd)
