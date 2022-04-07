import sander
from random import randint

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
                    print("Warning: coordinates file {self.coords} does not have TeraChem formatted output lines")
                    print("Assuming coordinates are contiguous and 0-indexed")
            except:
                raise RuntimeError("XYZ coordinates file {self.coords} is not correctly formatted")
        if terachemFormat:
            with open(os.path.join(self.dynamicsDir, self.coords), "r") as f:
                for line in f.readlines():
                    if "frame" in line:
                        coordIndex.append(int(line.split()[2])+1)
            startIndex = 0
            if start != None:
                while coordIndex[startIndex] < start:
                    startIndex += 1
            endIndex = len(coordIndex) - 1
            if end != None:
                wihle coordIndex[endIndex] > end:
                    endIndex -= 1
            splitIndex = 0
            if split != None:
                while coordIndex[splitIndex] < split:
                    splitIndex += 1
                if splitIndex == 0:
                    raise RuntimeError(f"There must be frames before {str(split)}")
                elif splitIndex >= len(coordIndex):
                    raise RuntimeError(f"There must be frames after {str(split)}")
            return startIndex, endIndex, splitIndex
        else:
            return start, end, split
            
    def getFrames(self, 
        if self.splitIndex == None:
            

class AmberEngine(MMEngine):
    
    def __init__(self, options):
        super().__init__(options)
        inputs = sander.gas_input()
        sander.setup(options.prmtop, options.coordinates, None, inputs)
    
class ExternalAmberEngine(MMEngine):

    def __init__(self, options):
        super()__init__(options)

    def sample(self, rst, mdin):
        name = rst.split(".")[0]
        for j in range(1, heatCounter + 1):
            os.system(f"sander -O -p {prmtop} -i heat{str(j)}.in -o {name}_heat{str(j)}.out -c {name}_heat{str(j-1)}.rst7 -x {name}_heat{str(j)}.nc -r {name}_heat{str(j)}.rst7")
        os.system(f"sander -O -p {prmtop} -i {mdin} -o {name}_heat{str(heatCounter)}.out -c {name}_heat{str(heatCounter-1)}.rst7 -x {name}.nc -r {name}_heat{str(j)}.rst7 -v {name}_vel.nc")
        os.mkdir(name)
        os.system(f"mv {name}.nc {os.path.join(name,'.')}")
        os.chdir(name)
        with open("cpptraj.in",'w') as f:
            f.write(f"loadcrd {name}.nc coors\n")
            f.write(f"crdout coors {name}.pdb multi\n")
            f.write("exit\n")
        try:
            os.system(f"cpptraj -p {options.prmtop} -i cpptraj.in > cpptraj.out")
        except:
            raise RuntimeError("Error in trajectory postprocessing in {os.getcwd()}")
        for f in os.listdir():
            if ".pdb" in f and len(f.split(".")) > 2:
                label = f.split(".")[2]
                os.system(f"mv {f} {label}.pdb")
        os.chdir("../")

    def getMMSamples(self):
        frames = self.getFrames(self.options.nvalids)
        
                


class OpenMMEngine(MMEngine):
    
    def __init__():
        super()__init__()
