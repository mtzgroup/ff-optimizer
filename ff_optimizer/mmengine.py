try:
    import sander
except:
    pass
from random import randint
from subprocess import check_output
from .utils import writeRst
import os
import GPUtil
from shutil import copyfile, rmtree


class MMEngine():

    def __init__(self, options):
        self.options = options
        self.coordPath = os.path.join(options['coordsDir'], options['coords'])
        self.startIndex, self.endIndex, self.splitIndex = self.getIndices()
        pass

    # 0-indexed
    # TODO: figure out what to do with this awful function
    # It really shouldn't be my responsibility to untangle an awful coors.xyz file
    def getIndices(self):
        start = self.options['start']
        end = self.options['end']
        split = self.options['split']
        self.coordIndex = []
        terachemFormat = True
        with open(self.coordPath, "r") as f:
            try:
                natoms = int(f.readline())
                if f.readline().split()[1] != "frame":
                    terachemFormat = False
                    print("Warning: coordinates file {self.options['coords']} does not have TeraChem formatted comment lines")
                    print("Assuming coordinates are contiguous and 0-indexed")
            except:
                raise RuntimeError("XYZ coordinates file {self.options['coords']} is not correctly formatted")
        if terachemFormat:
            splitIndex = None
            with open(self.coordPath, "r") as f:
                try:
                    for line in f.readlines():
                        if "frame" in line:
                            self.coordIndex.append(int(line.split()[2])+1)
                except:
                    raise RuntimeError("XYZ coordinates file {self.options['coords']} is not correctly formatted")
            nframes = len(self.coordIndex)
            startIndex = 0
            if start != None:
                while self.coordIndex[startIndex] < start and startIndex < nframes:
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
                while self.coordIndex[splitIndex] < split and splitIndex < nframes:
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
                    raise RuntimeError("XYZ coordinates file {self.options['coords']} is not correctly formatted")
            if end > nframes:
                end = nframes
            if split > nframes:
                raise ValueError("There must be frames after split")
            return start, end, split
            
    def getFrames(self, samplePath="."):
        frames = []
        if self.splitIndex is None:
            for i in range(self.options['nvalids']+1):
                frames.append(self.coordIndex[randint(self.startIndex, self.endIndex)])
        else:
            frames.append(self.coordIndex[randint(self.startIndex, self.splitIndex - 1)])
            for i in range(self.options['nvalids']):
                frames.append(self.coordIndex[randint(self.splitIndex, self.endIndex)])

        for frame in frames:
            self.getFrame(index, os.path.join(samplePath,f"{str(index)}.rst7"))
        return frames

    # 0-indexed 
    def getFrame(self, index, dest):
        frame = []
        atoms = []
        inFrame = False
        lineCounter = 0
        frameCounter = 1
        with open(self.coordPath, "r") as f:
            natoms = int(f.readline())
            for line in f.readlines():
                if int(lineCounter / (natoms + 2)) == index:
                    frame.append(line)
        writeRst(frame, natoms, dest)
            
    def getMMSamples(self):
        self.setup()
        frames = self.getFrames()
        name = frames[0].split('.')[0]
        if not os.path.isdir(name):
            os.mkdir(f"train_{name}")
        os.chdir(f"train_{name}")
        self.sample(frames[0], self.options["trainMdin"])
        with open("MMFinished.txt",'w') as f:
            f.write("MM sampling finished")
        os.chdir("..")
        for i in frames[1:]:
            name = frames[i].split('.')[0]
            if not os.path.isdir(name):
                os.mkdir(f"valid_{name}")
            os.chdir(f"valid_{name}")
            self.sample(frames[i], self.options["validMdin"])
            with open("MMFinished.txt",'w') as f:
                f.write("MM sampling finished")
            os.chdir("..")
            
    # This is the function that has to be implemented for every MMEngine
    def sample(self, rst, mdin):
        pass

    def setup(self):
        os.system(f"tleap -f {self.options['leap']} > leap.out")
        self.prmtop = None
        for f in os.listdir():
            if f.endswith(".prmtop"):
                self.prmtop = f
        if prmtop == None:
            raise RuntimeError(f"Tleap failed to create a new .prmtop file, check {os.path.join(os.getcwd(),'leap.out')} for more information")

    def restart(self):
        rsts = []
        for f in os.listdir(): 
            if f.endswith(".rst7"): 
                rsts.append(f)
        # Check to make sure we have the right number of ICs
        print(len(rsts))
        if len(rsts) != self.options['nvalids'] + 1:
            for rst in rsts:
                os.remove(rst)
            for f in os.listdir():
                if os.path.isdir(f):
                    rmtree(f)
            self.getMMSamples()
            return

        rsts = sorted(rsts)
        for i in range(len(rsts)):
            name = rsts[i].split('.')[0]
            folder = "."
            for f in os.listdir():
                if f.endswith(f"_{name}"):
                    folder = f
            # Skip finished jobs
            if os.path.isfile(os.path.join(folder,"MMFinished.txt")):
                continue
            if not os.path.isfile(os.path.join(folder,f"{name}.nc")):
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
                    self.sample(rsts[i], self.options['trainMdin'])
                else:
                    self.sample(rsts[i], self.options['validMdin'])
                os.chdir("..")
            # TODO: should distinguish between restarting from MD and from cpptraj
            else:
                os.chdir(folder)
                with open("cpptraj.in",'w') as f:
                    f.write(f"loadcrd {name}.nc coors\n")
                    f.write(f"crdout coors {name}.pdb multi\n")
                    f.write("exit\n")
                try:
                    os.system(f"cpptraj -p {self.prmtop} -i cpptraj.in > cpptraj.out")
                except Exception as e:
                    print(e)
                    raise RuntimeError(f"Error in trajectory postprocessing in {os.getcwd()}")
                os.chdir("..")
                        
class AmberEngine(MMEngine):
    
    def __init__(self, options):
        super().__init__(options)
        inputs = sander.gas_input()
        sander.setup(prmtop, options['coordinates'], None, inputs)
    
class ExternalAmberEngine(MMEngine):

    def __init__(self, options):
        self.options = options
        self.heatCounter = options['heatCounter']
        try:
            deviceIDs = GPUtil.getAvailable()
            if len(deviceIDs) > 0:
                print("Nvidia GPUs detected; defaulting to pmemd.cuda")
                self.amberExe = "pmemd.cuda"
            else:
                print("No Nvidia GPUs available; defaulting to sander")
                self.amberExe = "sander"
        except:
            print("No Nvidia GPUs available; defaulting to sander")
            self.amberExe = "pmemd"
        super().__init__(options)

    def runSander(self, prmtop, mdin, mdout, mdcrd, mdtraj, restart, mdvels=None):
        if mdvels is None:
            os.system(f"{self.amberExe} -O -p {prmtop} -i {mdin} -o {mdout} -c {mdcrd} -x {mdtraj} -r {restart}")
        else:
            os.system(f"{self.amberExe} -O -p {prmtop} -i {mdin} -o {mdout} -c {mdcrd} -x {mdtraj} -r {restart} -v {mdvels}")
            
    def sample(self, rst, mdin):
        name = rst.split(".")[0]
        if not os.path.isdir(name):
            os.mkdir(name)
        os.chdir(name)
        copyfile(os.path.join("..",rst),f"{name}_heat0.rst7")
        for j in range(1, self.heatCounter + 1):
            self.runSander(os.path.join("..",self.prmtop), os.path.join("..",f"heat{str(j)}.in"), f"{name}_heat{str(j)}.out", f"{name}_heat{str(j-1)}.rst7", f"{name}_heat{str(j)}.nc", f"{name}_heat{str(j)}.rst7")
        self.runSander(os.path.join("..",self.prmtop), os.path.join("..",mdin), f"{name}_.out", f"{name}_heat{str(self.heatCounter)}.rst7", f"{name}.nc", f"{name}_md.rst7", mdvels=f"{name}_vel.nc")
        with open("cpptraj.in",'w') as f:
            f.write(f"loadcrd {name}.nc coors\n")
            f.write(f"crdout coors {name}.pdb multi\n")
            f.write("exit\n")
        try:
            os.system(f"cpptraj -p {os.path.join('..',self.prmtop)} -i cpptraj.in > cpptraj.out")
        except Exception as e:
            print(e)
            raise RuntimeError(f"Error in trajectory postprocessing in {os.getcwd()}")
        for f in os.listdir():
            if ".pdb" in f and len(f.split(".")) > 2:
                label = f.split(".")[2]
                os.system(f"mv {f} {label}.pdb")
        os.chdir("../")

class ExternalOpenMMEngine(MMEngine):
    
    def __init__():
        uper().__init__()
