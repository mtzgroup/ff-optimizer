try:
    import sander
except:
    pass
import os
from random import randint
from shutil import copyfile, rmtree
import GPUtil
from .utils import writeRst

class MMEngine:
    def __init__(self, options):
        self.options = options
        # Account for being in basedir/sampledir/x_cycle_x/
        self.coordPath = os.path.join("..", "..", options["coordPath"])
        self.startIndex, self.endIndex, self.splitIndex = self.getIndices()

    # 0-indexed
    # TODO: figure out what to do with this awful function
    # It really shouldn't be my responsibility to untangle an awful coors.xyz file
    def getIndices(self):
        start = self.options["start"]
        end = self.options["end"]
        split = self.options["split"]
        self.coordIndex = []
        terachemFormat = True
        print(os.getcwd())
        with open(self.options["coordPath"], "r") as f:
            try:
                natoms = int(f.readline())
                if f.readline().split()[1] != "frame":
                    terachemFormat = False
                    print(
                        "Warning: coordinates file {self.options['coords']} does not have TeraChem formatted comment lines"
                    )
                    print("Assuming coordinates are contiguous and 0-indexed")
            except:
                raise RuntimeError(
                    "XYZ coordinates file {self.options['coords']} is not correctly formatted"
                )
        if terachemFormat:
            splitIndex = None
            with open(self.options["coordPath"], "r") as f:
                try:
                    for line in f.readlines():
                        if "frame" in line:
                            self.coordIndex.append(int(line.split()[2]))
                except:
                    raise RuntimeError(
                        "XYZ coordinates file {self.options['coords']} is not correctly formatted"
                    )
            nframes = len(self.coordIndex)
            startIndex = 0
            if start != None:
                while self.coordIndex[startIndex] < start and startIndex < nframes - 1:
                    startIndex += 1
            if startIndex == nframes:
                raise ValueError("No frames after start index!")
            endIndex = len(self.coordIndex) - 1
            if end != None:
                while self.coordIndex[endIndex] > end and endIndex > 0:
                    endIndex -= 1
            if endIndex == 0:
                raise ValueError("No frames before end index!")
            splitIndex = 0
            if split != None:
                while self.coordIndex[splitIndex] < split and splitIndex < nframes - 1:
                    splitIndex += 1
                if splitIndex == 0:
                    raise ValueError(f"There must be frames before {str(split)}")
                elif splitIndex >= len(self.coordIndex):
                    raise ValueError(f"There must be frames after {str(split)}")
            return startIndex, endIndex, splitIndex
        else:
            with open(self.coordPath, "r") as f:
                try:
                    natoms = int(f.readline())
                    lineCounter = 0
                    for line in f.readlines():
                        lineCounter += 1
                    nframes = int(lineCounter / (natoms + 2))
                except:
                    raise RuntimeError(
                        "XYZ coordinates file {self.options['coords']} is not correctly formatted"
                    )
            if end > nframes:
                end = nframes
            if split is not None:
                if split > nframes:
                    raise ValueError("There must be frames after split")
            if start > nframes:
                raise ValueError("There must be frames after start")
            return start, end, split

    def getFrames(self, samplePath="."):
        frames = []
        if self.splitIndex is None:
            for i in range((self.options["nvalids"] + 1) * self.options["conformers"]):
                frames.append(self.coordIndex[randint(self.startIndex, self.endIndex)])
        else:
            for i in range(self.options["conformers"]):
                frames.append(
                    self.coordIndex[randint(self.startIndex, self.splitIndex - 1)]
                )
            for i in range(self.options["nvalids"] * self.options["conformers"]):
                frames.append(self.coordIndex[randint(self.splitIndex, self.endIndex)])

        return frames

    # 0-indexed
    def getFrame(self, index, dest):
        frame = []
        with open(self.coordPath, "r") as f:
            natoms = int(f.readline())
            lines = f.readlines()
            for i in range(len(lines)):
                if int(i / (natoms + 2)) == index and (i + 1) % (natoms + 2) > 1:
                    frame.append(lines[i].split()[1:])
                elif int((i + 1) / (natoms + 2)) > index:
                    break
        writeRst(frame, natoms, dest)

    def getMMSamples(self):
        self.setup()
        frames = self.getFrames()
        name = str(frames[0])
        if not os.path.isdir("train"):
            os.mkdir("train")
        trainFrames = frames[:self.options['conformers']]
        for frame in trainFrames:
            self.getFrame(frame, os.path.join("train", f"{str(frame)}.rst7"))
        os.chdir("train")
        self.sample(self.options["trainMdin"])
        with open("MMFinished.txt", "w") as f:
            f.write("MM sampling finished")
        os.chdir("..")
        for i in self.options["nvalids"]:
            validName = f"valid_{str(i)}"
            if not os.path.isdir(validName):
                os.mkdir(validName)
            validFrames = frames[(i+1) * self.options['conformers'] : (i+2) * self.options['conformers']]
            for frame in validFrames:
                self.getFrame(frame, os.path.join(validName),f"{str(frame)}.rst7")
            os.chdir(validName)
            self.sample(validFrames, self.options["validMdin"])
            with open("MMFinished.txt", "w") as f:
                f.write("MM sampling finished\n")
                f.write("Remove this file and the .nc file\n")
                f.write("If you want to force a recalculation of the MM sampling")
            os.chdir("..")

    # This is the function that has to be implemented for every MMEngine
    def sample(self, frames, mdin):
        pass

    def setup(self):
        os.system(f"tleap -f {self.options['leap']} > leap.out")
        self.prmtop = None
        for f in os.listdir():
            if f.endswith(".prmtop"):
                self.prmtop = f
        if self.prmtop == None:
            raise RuntimeError(
                f"Tleap failed to create a new .prmtop file, check {os.path.join(os.getcwd(),'leap.out')} for more information"
            )

    def restart(self):
        rsts = []
        for f in os.listdir():
            if f.endswith(".rst7"):
                rsts.append(f)
        # Check to make sure we have the right number of ICs
        if len(rsts) != self.options["nvalids"] + 1:
            for rst in rsts:
                os.remove(rst)
            for f in os.listdir():
                if os.path.isdir(f):
                    rmtree(f)
            self.getMMSamples()
            return

        rsts = sorted(rsts)
        self.setup()
        for i in range(len(rsts)):
            name = rsts[i].split(".")[0]
            folder = "."
            for f in os.listdir():
                if f.endswith(f"_{name}"):
                    folder = f
            # Skip finished jobs
            if os.path.isfile(os.path.join(folder, "MMFinished.txt")):
                continue
            if not os.path.isfile(os.path.join(folder, f"{name}.nc")):
                if folder == ".":
                    if i == 0:
                        isValid = False
                        for f in os.listdir():
                            if f.startswith("train"):
                                isValid = True
                    else:
                        isValid = True
                    if isValid:
                        folder = f"valid_{str(name)}"
                    else:
                        folder = f"train_{str(name)}"
                    os.mkdir(folder)
                elif folder.startswith("train"):
                    isValid = False
                else:
                    isValid = True
                os.chdir(folder)
                if not isValid:
                    self.sample(name, self.options["trainMdin"])
                else:
                    self.sample(name, self.options["validMdin"])
                os.chdir("..")
            # TODO: should distinguish between restarting from MD and from cpptraj
            else:
                os.chdir(folder)
                with open("cpptraj.in", "w") as f:
                    f.write(f"loadcrd {name}.nc coors\n")
                    f.write(f"crdout coors {name}.pdb multi\n")
                    f.write("exit\n")
                try:
                    os.system(f"cpptraj -p {self.prmtop} -i cpptraj.in > cpptraj.out")
                except Exception as e:
                    print(e)
                    raise RuntimeError(
                        f"Error in trajectory postprocessing in {os.getcwd()}"
                    )
                os.chdir("..")


class AmberEngine(MMEngine):
    def __init__(self, options):
        super().__init__(options)
        inputs = sander.gas_input()
        sander.setup(prmtop, options["coordinates"], None, inputs)

class ExternalAmberEngine(MMEngine):
    def __init__(self, options):
        self.options = options
        self.heatCounter = options["heatCounter"]
        try:
            deviceIDs = GPUtil.getAvailable(maxLoad=0.1)
            if len(deviceIDs) > 0:
                print("Nvidia GPUs detected; defaulting to pmemd.cuda")
                self.amberExe = "pmemd.cuda"
                os.environ["CUDA_VISIBLE_DEVICES"] = str(deviceIDs[0])
            else:
                print("No Nvidia GPUs available; defaulting to sander")
                self.amberExe = "pmemd"
        except:
            print("No Nvidia GPUs available; defaulting to sander")
            self.amberExe = "pmemd"
        super().__init__(options)

    def runSander(self, prmtop, mdin, mdout, mdcrd, mdtraj, restart, mdvels=None):
        if not os.path.isfile(prmtop):
            raise RuntimeError(f"Cannot find prmtop {prmtop} in {os.getcwd()}")
        if not os.path.isfile(mdin):
            raise RuntimeError(f"Cannot find mdin {mdin} in {os.getcwd()}")
        if not os.path.isfile(mdcrd):
            raise RuntimeError(f"Cannot find input crd {mdcrd} in {os.getcwd()}")
        if mdvels is None:
            os.system(
                f"{self.amberExe} -O -p {prmtop} -i {mdin} -o {mdout} -c {mdcrd} -x {mdtraj} -r {restart}"
            )
        else:
            os.system(
                f"{self.amberExe} -O -p {prmtop} -i {mdin} -o {mdout} -c {mdcrd} -x {mdtraj} -r {restart} -v {mdvels}"
            )
        if not os.path.isfile(restart):
            if self.amberExe == "pmemd.cuda":
                print("pmemd.cuda failed; trying MM sampling with pmemd")
                self.amberExe = "pmemd"
                runSander(pmrtop, mdin, mdout, mdcrd, mdtraj, restart, mdvels=mdvels)
            else:
                raise RuntimeError(f"MM dynamics with input {mdin} failed in {os.getcwd()}")

    def sample(self, frames, mdin):
        pdbIndex = 1
        for frame in frames:
            name = str(frame)
            copyfile(f"{name}.rst7", f"{name}_heat0.rst7")
            for j in range(1, self.heatCounter + 1):
                self.runSander(
                    os.path.join("..", self.prmtop),
                    os.path.join("..", f"heat{str(j)}.in"),
                    f"{name}_heat{str(j)}.out",
                    f"{name}_heat{str(j-1)}.rst7",
                    f"{name}_heat{str(j)}.nc",
                    f"{name}_heat{str(j)}.rst7",
                )
            self.runSander(
                os.path.join("..", self.prmtop),
                os.path.join("..", mdin),
                f"{name}_.out",
                f"{name}_heat{str(self.heatCounter)}.rst7",
                f"{name}.nc",
                f"{name}_md.rst7",
                mdvels=f"{name}_vel.nc",
            )
            with open("cpptraj.in", "w") as f:
                f.write(f"loadcrd {name}.nc coors\n")
                f.write(f"crdout coors {name}.pdb multi\n")
                f.write("exit\n")
            try:
                os.system(
                    f"cpptraj -p {os.path.join('..',self.prmtop)} -i cpptraj.in > cpptraj.out"
                )
            except Exception as e:
                print(e)
                raise RuntimeError(f"Error in trajectory postprocessing in {os.getcwd()}")
            for f in os.listdir():
                if ".pdb" in f and len(f.split(".")) > 2:
                    os.system(f"mv {f} {str(pdbIndex)}.pdb")
                    pdbIndex += 1


class ExternalOpenMMEngine(MMEngine):
    def __init__():
        super().__init__()
