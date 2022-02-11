import numpy as np

def convertPDBtoXYZ(pdb):
    name = pdb.replace(".pdb",".xyz")
    xyzLines = []
    with open(pdb,'r') as f:
        for line in f.readlines():
            splitLine = line.split()
            if len(splitLine) == 0:
                continue
            if splitLine[0] == "ATOM" or splitLine[0] == "HETATM":
                xyzLines.append(f"{splitLine[10]}\t{splitLine[5]}\t{splitLine[6]}\t{splitLine[7]}\n")
    with open(name,'w') as f:
        f.write(f"{str(len(xyzLines))}\n")
        f.write("Converted from {pdb}\n")
        for line in xyzLines:
            f.write(line)
    return name

