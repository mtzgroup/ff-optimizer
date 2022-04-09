try:
    import sander
except:
    pass
from random import randint
from subprocess import check_output
from utils import writeRst


class MMEngine():

    def __init__(self, homeDir, dynamicsDir, coords, start, end, split):
        self.dynamicsdir = dynamicsDir
        self.coords = coords
        self.homeDir = homeDir
        self.startIndex, self.endIndex, self.splitIndex = self.getIndices(start, end, split)
        pass

    def getIndices(self, start, end, split):
        self.coordIndex = []
        terachemFormat = True
        with open(os.path.join(self.dynamicsDir, self.coords), "r") as f:
            try:
                natoms = int(f.readline())
                if f.readline().split()[1] != "frame":
                    terachemFormat = False
                    print("Warning: coordinates file {self.coords} does not have TeraChem formatted comment lines")
                    print("Assuming coordinates are contiguous and 0-indexed")
            except:
                raise RuntimeError("XYZ coordinates file {self.coords} is not correctly formatted")
        if terachemFormat:
            splitIndex = None
            with open(os.path.join(self.dynamicsDir, self.coords), "r") as f:
                for line in f.readlines():
                    if "frame" in line:
                        self.coordIndex.append(int(line.split()[2])+1)
            startIndex = 0
            if start != None:
                while self.coordIndex[startIndex] < start:
                    startIndex += 1
            endIndex = len(coordIndex) - 1
            if end != None:
                while self.coordIndex[endIndex] > end:
                    endIndex -= 1
            splitIndex = 0
            if split != None:
                while self.coordIndex[splitIndex] < split:
                    splitIndex += 1
                if splitIndex == 0:
                    raise RuntimeError(f"There must be frames before {str(split)}")
                elif splitIndex >= len(self.coordIndex):
                    raise RuntimeError(f"There must be frames after {str(split)}")
            return startIndex, endIndex, splitIndex
        else:
            with open(os.path.join(self.dynamicsDir, self.coords), "r") as f:
            return start, end, split
            
    def getFrames(self, samplePath):
        frames = []
        if self.splitIndex is None:
            for i in range(self.options['nvalids']+1):
                frames.append(self.coordIndex[randint(self.startIndex, self.endIndex)])
        else:
            frames.append(self.coordIndex[randint(self.startIndex, self.splitIndex - 1)])
            for i in range(self.options['nvalids']):
                frames.append(self.coordIndex[randint(self.startIndex, self.splitIndex - 1)])

        for frame in frames:
            self.getFrame(index, os.path.join(self.options['dynamicsdir'], self.options['coors']), samplePath)
        return frames

    # 0-indexed 
    def getFrame(self, index, coordFile, dest):
        frame = []
        atoms = []
        inFrame = False
        lineCounter = 0
        frameCounter = 1
        with open(coordFile, "r") as f:
            natoms = int(f.readline())
            for line in f.readlines():
                if int(lineCounter / (natoms + 2)) == index:
                    frame.append(line)
        writeRst(frame, index, molSize, dest)
            
            

    # This is the function that has to be implemented for every MMEngine
    def getMMSamples(self):
        pass

class AmberEngine(MMEngine):
    
    def __init__(self, options):
        super().__init__(options)
        inputs = sander.gas_input()
        sander.setup(options['prmtop'], options['coordinates'], None, inputs)
    
class ExternalAmberEngine(MMEngine):

    def __init__(self, options):
        self.options = options
        ns = check_output(["nvidia-smi"])
        if "NVIDIA-SMI has failed" in ns:
            print("No Nvidia GPUs available; defaulting to sander")
            self.amberExe = "sander"
        else:
            print("Nvidia GPUs detected; defaulting to pmemd.cuda")
            self.amberExe = "pmemd.cuda"
            try:
                gpus = os.environ['CUDA_VISIBLE_DEVICES']
            except:
                print("Warning: CUDA_VISIBLE_DEVICES is not set!")
        super()__init__(options)

    def runSander(self, prmtop, mdin, mdout, mdcrd, mdtraj, restart, mdvels=None):
        if mdvels is None:
            os.system(f"{self.amberExe} -O -p {prmtop} -i {mdin} -o {mdout} -c {mdcrd} -x {mdtraj} -r {restart}")
        else:
            os.system(f"{self.amberExe} -O -p {prmtop} -i {mdin} -o {mdout} -c {mdcrd} -x {mdtraj} -r {restart} -v {mdvels}")
            
    def sample(self, rst, mdin):
        name = rst.split(".")[0]
        for j in range(1, heatCounter + 1):
            runSander(self.options['prmtop'], f"heat{str(j)}.in", f"{name}_heat{str(j)}.out", f"{name}_heat{str(j-1)}.rst7", f"{name}_heat{str(j)}.nc", f"{name}_heat{str(j)}.rst7")
        runSander(self.options['prmtop'], mdin, f"{name}_.out", f"{name}_heat{str(heatCounter)}.rst7", f"{name}.nc", f"{name}_md.rst7", mdvels=f"{name}_vel.nc")
        os.mkdir(name)
        os.system(f"mv {name}.nc {os.path.join(name,'.')}")
        os.chdir(name)
        with open("cpptraj.in",'w') as f:
            f.write(f"loadcrd {name}.nc coors\n")
            f.write(f"crdout coors {name}.pdb multi\n")
            f.write("exit\n")
        try:
            os.system(f"cpptraj -p {self.options['prmtop']} -i cpptraj.in > cpptraj.out")
        except:
            raise RuntimeError("Error in trajectory postprocessing in {os.getcwd()}")
        for f in os.listdir():
            if ".pdb" in f and len(f.split(".")) > 2:
                label = f.split(".")[2]
                os.system(f"mv {f} {label}.pdb")
        os.chdir("../")

    def getMMSamples(self):
        frames = self.getFrames(self.options['nvalids'])
        self.sample(frames[0], self.options["trainMdin"])
        name = frames[0].split('.')[0]
        os.system(f"mv {name} train_{name}")
        for i in frames[1:]:
            self.sample(frames[i], self.options["validMdin"])
            name = frames[i].split('.')[0]
            os.system(f"mv {name} valid_{str(i)_{name}")
            
class ExternalOpenMMEngine(MMEngine):
    
    def __init__():
        super()__init__()
