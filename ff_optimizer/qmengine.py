import os
import subprocess
import traceback
from shutil import copyfile, rmtree
from time import sleep

from chemcloud import CCClient
from qcio import ProgramInput, SinglePointResults, Structure
from qcparse import parse
from pathlib import Path

from . import utils


class QMEngine:
    def __init__(self, inp):
        if inp.resp == 0 and inp.resppriors == 0:
            self.doResp = False
        else:
            self.doResp = True
        self.inputSettings = self.readInputFile(inp.sampledir / inp.tctemplate)
        self.backupInputSettings = self.readInputFile(
            inp.sampledir / inp.tctemplate_backup
        )

    def readInputFile(self, inputFile: str):
        settings = []
        with open(inputFile, "r") as f:
            for line in f.readlines():
                splitLine = line.split()
                # ignore commented lines, the coordinates specification, run gradient, and resp keywords
                if len(splitLine) > 0:
                    if (
                        splitLine[0][0] != "#"
                        and splitLine[0].lower() != "coordinates"
                        and splitLine[0].lower() != "run"
                        and splitLine[0].lower() != "resp"
                        # and splitLine[0].lower() != "esp_restraint_a"
                        # and splitLine[0].lower() != "esp_restraint_b"
                    ):
                        setting = []
                        # for options > 2 tokens long
                        for token in splitLine:
                            # ignore comments at the end of the line
                            if token[0] != "#":
                                setting.append(token)
                            else:
                                break
                        settings.append(setting)
        if self.doResp:
            settings.append(["resp", "yes"])
        return settings

    def writeInputFile(self, settings: list, coordinates: str, fileName: str):
        with open(fileName, "w") as f:
            f.write(f"coordinates {coordinates}\n")
            f.write(f"run gradient\n")
            # if self.doResp:
            # f.write("resp yes\n")
            # f.write("esp_restraint_a 0\n")
            # f.write("esp_restraint_b 0\n")
            for setting in settings:
                # reassemble setting (if > 2 tokens long)
                for token in setting:
                    f.write(f"{token} ")
                f.write("\n")

    def writeFBdata(
        self, energies: list, grads: list, coords: list, espXYZs=None, esps=None
    ):
        with open("qdata.txt", "w") as f:
            for i in range(len(energies)):
                f.write(f"JOB {i+1}\n")
                coordLine = "COORDS "
                for coord in coords[i]:
                    coordLine = coordLine + str(round(float(coord), 5)) + " "
                gradLine = "FORCES "
                for grad in grads[i]:
                    gradLine = gradLine + str(round(float(grad), 8)) + " "
                f.write(coordLine + "\n")
                f.write(f"ENERGY {str(round(float(energies[i]),7))}\n")
                f.write(gradLine + "\n")
                if espXYZs is not None:
                    espXYZLine = "ESPXYZ "
                    for xyz in espXYZs[i]:
                        espXYZLine = espXYZLine + str(round(float(xyz), 6)) + " "
                    f.write(espXYZLine + "\n")
                    espLine = "ESPVAL "
                    for val in esps[i]:
                        espLine = espLine + str(round(float(val), 6)) + " "
                    f.write(espLine + "\n\n")
                else:
                    f.write("\n")

        with open("all.mdcrd", "w") as f:
            f.write("Converted from xyz by QMEngine\n")
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

    def readResult(self, js):
        result = SinglePointResults.open(js)
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
        for f in utils.getXYZs():
            coord = utils.readXYZ(f)
            name = utils.getName(f)
            result = parse(f"tc_{name}.out", "terachem")
            if result.energy is None or result.gradient is None:
                raise RuntimeError(
                    f"Terachem job tc_{name}.out in {os.getcwd()} did not succeed!"
                )
            energies.append(result.energy)
            grads.append(result.gradient.flatten())
            coords.append(coord)
            if self.doResp:
                espXYZ, esp = utils.readEsp(f"esp_{name}.xyz")
                espXYZs.append(espXYZ)
                esps.append(esp)
        return energies, grads, coords, espXYZs, esps

    def restart(self):
        xyzs = []
        for f in utils.getXYZs():
            name = utils.getName(f)
            out = f"tc_{name}.out"
            try:
                parse(out, "terachem")
            except Exception as e:
                print(f)
                print(e)
                xyzs.append(f)
        self.getQMRefData(xyzs)

    # This is the function which must be implemented by inherited classes
    def getQMRefData(self, xyzs: list):
        raise NotImplementedError()


class SlurmEngine(QMEngine):
    def __init__(
        self,
        inp,
    ):
        self.readSbatchFile(inp.sampledir / inp.sbatchtemplate)
        self.inp = inp
        super().__init__(inp)

    def readSbatchFile(self, sbatchFile: str):
        with open(sbatchFile, "r") as f:
            self.sbatchLines = f.readlines()

    def replaceVars(self, line, index):
        tcin = f"tc_{index}.in"
        tcbackup = f"tc_{index}_backup.in"
        line = line.replace("JOBID", index)
        line = line.replace("TCTEMPLATEBACKUP", tcbackup)
        line = line.replace("TCTEMPLATE", tcin)
        return line

    def writeSbatchFile(self, index: str, fileName: str):
        index = str(index)
        with open(fileName, "w") as f:
            for line in self.sbatchLines:
                f.write(self.replaceVars(line, index))

    def slurmCommand(self, command: list):
        maxTries = 10
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

    def waitForJobs(self, jobIDs):
        while len(jobIDs) > 0:
            runningIDs = []
            sleep(10)
            # status = self.slurmCommand(["squeue", "-o", "%.12i", "-u", self.user])
            status = self.slurmCommand(["squeue", "-o", "%.12i"])
            for runningID in (
                status.replace(b" ", b"").replace(b'"', b"").split(b"\n")[1:]
            ):
                for submitID in jobIDs:
                    if runningID == submitID:
                        runningIDs.append(runningID)
            jobIDs = runningIDs

    def getQMRefData(self, xyzs: list):
        jobIDs = []
        for xyz in xyzs:
            name = xyz.name.split(".")[0]
            super().writeInputFile(self.inputSettings, xyz, f"tc_{name}.in")
            super().writeInputFile(
                self.backupInputSettings, xyz, f"tc_backup_{name}.in"
            )
            self.writeSbatchFile(name, f"sbatch_{name}.sh")
            job = self.slurmCommand(["sbatch", f"sbatch_{name}.sh"])
            jobIDs.append(job.split()[3])
        self.waitForJobs(jobIDs)

        for xyz in xyzs:
            name = xyz.name.split(".")[0]
            if self.doResp:
                copyfile(f"scr.{name}/esp.xyz", f"esp_{name}.xyz")
            scr = Path(f"scr.{name}")
            if scr.is_dir():
                rmtree(scr)

        energies, grads, coords, espXYZs, esps = super().readQMRefData()
        super().writeFBdata(energies, grads, coords, espXYZs, esps)


class DebugEngine(QMEngine):
    def __init__(self, inp):
        super().__init__(inp)

    def getQMRefData(self, xyzs: list):
        energies = []
        grads = []
        coords = []
        for xyz in xyzs:
            name = xyz.name.split(".")[0]
            super().writeInputFile(self.inputSettings, xyz, f"tc_{name}.in")
            os.system(f"terachem tc_{name}.in > tc_{name}.out")
            energy, grad = utils.readGradFromTCout(f"tc_{name}.out")
            if grad == -1:
                super().writeInputFile(
                    self.backupInputSettings, xyz, f"tc_{name}_backup.in"
                )
                os.system(f"terachem tc_{name}_backup.in > tc_{name}.out")
            if self.doResp:
                espXYZ, esps = utils.readEsp(f"scr.{name}/esp.xyz")

        energies, grads, coords, espXYZs, esps = super().readQMRefData()
        super().writeFBdata(energies, grads, coords, espXYZs, esps)


class ChemcloudEngine(QMEngine):
    def __init__(self, inp):
        self.batchSize = inp.batchsize
        self.retries = inp.retries
        self.client = CCClient()
        super().__init__(inp)
        self.setKeywords(self.doResp)

    def setKeywords(self, doResp):
        self.keywords = {}
        self.backupKeywords = {}
        self.specialKeywords = {}
        self.keywords = self.sortSettings(self.inputSettings)
        self.backupKeywords = self.sortSettings(self.backupInputSettings, True)
        self.checkSpecialKeywords()
        if doResp:
            self.keywords["resp"] = "yes"
            # self.keywords["esp_restraint_a"] = "0"
            # self.keywords["esp_restraint_b"] = "0"
            self.backupKeywords["resp"] = "yes"
            # self.backupKeywords["esp_restraint_a"] = "0"
            # self.backupKeywords["esp_restraint_b"] = "0"

    def sortSettings(self, settings, backup=False):
        keywords = {}
        specialSettings = ["method", "basis", "charge", "spinmult"]
        for setting in settings:
            if setting[0] in specialSettings:
                if not backup:
                    self.specialKeywords[setting[0]] = setting[1]
            else:
                if len(setting) == 2:
                    keywords[setting[0]] = setting[1]
                else:
                    keyword = setting[1]
                    for token in setting[2:]:
                        keyword += f" {token}"
                    keywords[setting[0]] = keyword
        return keywords

    def loadStructureFromXYZ(self, xyz):
        tempMol = Structure.open(xyz)
        kwargs = tempMol.model_dump()
        if "charge" in self.specialKeywords.keys():
            kwargs["charge"] = self.specialKeywords["charge"]
        if "spinmult" in self.specialKeywords.keys():
            kwargs["multiplicity"] = self.specialKeywords["spinmult"]
        mol = Structure(**kwargs)
        return mol

    def checkSpecialKeywords(self):
        if "method" not in self.specialKeywords.keys():
            raise ValueError("ES method not specified in tc template file")
        if "basis" not in self.specialKeywords.keys():
            raise ValueError("ES basis not specified in tc template file")
        if "charge" in self.specialKeywords.keys():
            try:
                self.specialKeywords["charge"] = int(self.specialKeywords["charge"])
            except:
                raise ValueError("Molecular charge must be an int")
        else:
            raise ValueError(
                "Molecular charge ('charge') must be specified in input file"
            )
        if "spinmult" in self.specialKeywords.keys():
            try:
                self.specialKeywords["spinmult"] = int(self.specialKeywords["spinmult"])
            except:
                raise ValueError("Molecular spin multiplicity must be an integer")
        else:
            raise ValueError(
                "Spin multiplicity ('spinmult') must be specified in input file"
            )

    def computeBatch(self, programInputs: list):
        status = 0
        # If there are no jobs to run after restart
        if len(programInputs) == 0:
            return status, []
        batchSize = min(self.batchSize, len(programInputs))
        outputs = []
        stride = int(len(programInputs) / batchSize)
        # import pickle
        # with open("inputs.pickle", "wb") as f:
        #    pickle.dump(programInputs, f)
        try:
            # HOW TO RESTART IF CODE FAILS AFTER SUBMISSION?
            futureResults = [
                self.client.compute(
                    "terachem", programInputs[i::stride], collect_files=True
                )
                for i in range(stride)
            ]
            outputBatches = [futureResults[i].get() for i in range(stride)]
            for batch in outputBatches:
                for output in batch:
                    outputs.append(output)
        except Exception as e:
            traceback.print_exc()
            print(e)
            self.batchSize = int(batchSize / 2)
            print(
                f"Submission failed; resubmitting with batch size {str(self.batchSize)}"
            )
            # sleep(30)
            if self.batchSize < 2:
                status = -1
                return status, outputs
            tempStatus, outputs = self.computeBatch(programInputs)
            if tempStatus == -1:
                status = -1
        return status, outputs

    def createProgramInputs(self, xyzs: list, useBackup=False):
        mod = {}
        mod["method"] = self.specialKeywords["method"]
        mod["basis"] = self.specialKeywords["basis"]
        programInputs = []
        if useBackup:
            keywords = self.backupKeywords
        else:
            keywords = self.keywords
        for xyz in sorted(xyzs):
            jobID = utils.getName(xyz)
            mol = self.loadStructureFromXYZ(xyz)
            programInput = ProgramInput(
                structure=mol,
                model=mod,
                calctype="gradient",
                keywords=keywords,
                extras={"id": jobID},
            )
            programInputs.append(programInput)
        return programInputs

    def writeResult(self, output):
        jobID = output.input_data.extras["id"]
        out = f"tc_{jobID}.out"
        with open(out, "w") as f:
            f.write(output.stdout)
        if self.doResp and output.success:
            with open(f"esp_{jobID}.xyz", "w") as f:
                f.write(output.files["esp.xyz"])

    def getFailedJobs(self, outputs):
        retryXyzs = []
        for output in outputs:
            self.writeResult(output)
            if not output.success:
                jobID = output.input_data.extras["id"]
                retryXyzs.append(f"{jobID}.xyz")
        return retryXyzs

    def runJobs(self, xyzs: list, useBackup=False):
        programInputs = self.createProgramInputs(xyzs, useBackup=useBackup)
        status, outputs = self.computeBatch(programInputs)
        if len(outputs) != len(programInputs):
            dumpFailedJobs(programInputs, outputs)
            raise RuntimeError(
                "ChemCloud did not return the same number of outputs as inputs"
            )
        retryXyzs = self.getFailedJobs(outputs)
        if status == -1:
            dumpFailedJobs(programInputs, outputs)
            raise RuntimeError(
                "Batch resubmission reached size 1; QM calculations incomplete"
            )
        return retryXyzs

    def getQMRefData(self, xyzs: list):
        cwd = os.getcwd()
        retryXyzs = self.runJobs(xyzs)
        for _ in range(self.retries):
            if len(retryXyzs) == 0:
                break
            retryXyzs = self.runJobs(retryXyzs, useBackup=True)
        if len(retryXyzs) > 0:
            # for result in [self.readResult(f"tc_{xyz.split('.')[0]}.json") for xyz in retryXyzs]:
            #    print(
            #        f"Job id {result.input_data['id']} suffered a {result.error.error_type}"
            #    )
            #    print(result.error.error_message)
            raise RuntimeError(
                f"Job ids {sorted([xyz.split('.')[0] for xyz in retryXyzs])} in {os.getcwd()} failed {str(self.retries)} times!"
            )
        energies, grads, coords, espXYZs, esps = super().readQMRefData()
        super().writeFBdata(energies, grads, coords, espXYZs, esps)
        os.chdir(cwd)


def dumpFailedJobs(programInputs, outputs):
    for inp in programInputs:
        name = inp.extras["id"] + "_input.yaml"
        inp.save(name)
    for out in outputs:
        name = out.input_data.extras["id"] + "_output.yaml"
        out.save(name)
