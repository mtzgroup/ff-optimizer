from tccloud import TCClient
from tccloud.models import AtomicInput, Molecule
from os import listdir
from math import ceil
from time import sleep


def computeBatch(atomicInputs, batchSize):
    results = []
    for i in range(ceil(len(pdbs) / batchSize)):
        resultsBatch = []
        atomicInputsBatch = []
        for j in range(batchSize):
            if i * batchSize + j < len(atomicInputs):
                atomicInputsBatch.append(atomicInputs[i * batchSize + j])
        try:
            future_result_batch = client.compute(
                atomicInputsBatch, engine="terachem_pbs"
            )
            resultsBatch = future_result_batch.get()
        except:
            sleep(100)
            batchSizeResubmit = int(batchSize / 1.5)
            if batchSizeResubmit < 2:
                raise RuntimeError("Batch resubmission reached size 1")
            resultsBatch = computeBatch(atomicInputsBatch, batchSizeResubmit)
        for result in resultsBatch:
            results.append(result)
    return results


AU_TO_ANG = 0.5291772

files = listdir()
pdbs = []
for f in files:
    if ".pdb" in f:
        pdbs.append(f)

atomicInputs = []
client = TCClient()

for pdb in pdbs:
    name = pdb.replace(".pdb", ".xyz")
    coords = []
    atoms = []
    with open(pdb, "r", errors="ignore") as f:
        for line in f.readlines():
            if "ATOM" in line or "HETATM" in line:
                splitLine = line.split()
                coords.append(splitLine[5] + " " + splitLine[6] + " " + splitLine[7])
                atoms.append(splitLine[10])

    with open(name, "w") as f:
        f.write(str(len(coords)) + "\n")
        f.write("Converted from " + pdb + "\n")
        for i in range(len(coords)):
            f.write(atoms[i] + " " + coords[i] + "\n")

    mol = Molecule.from_file(name)
    atomicInput = AtomicInput(
        molecule=mol,
        model={"method": "B3LYP", "basis": "6-31gss"},
        driver="gradient",
        keywords={"closed": True, "restricted": True, "dftd": "d3"},
    )
    atomicInputs.append(atomicInput)

batchSize = 10
results = computeBatch(atomicInputs, batchSize)
retryInputs = []
failedIndices = []
for i in range(len(results)):
    if not results[i].success:
        retryInput = AtomicInput(
            molecule=atomicInputs[i].molecule,
            model=atomicInputs[i].model,
            driver=atomicInputs[i].driver,
            keywords={
                "closed": True,
                "restricted": True,
                "dftd": "d3",
                "threall": "1.0e-15",
                "diismaxvecs": 40,
            },
        )
        retryInputs.append(retryInput)
        failedIndices.append(i)
retryResults = computeBatch(retryInputs, batchSize)
print(failedIndices)
for i, result in zip(failedIndices, retryResults):
    print(atomicInputs[i].molecule.geometry)
    print(result)
    results[i] = result

i = 1
with open("qdata.txt", "w") as f:
    for result in results:
        f.write("JOB " + str(i) + "\n")
        coordLine = "COORDS "
        for atom in result.molecule.geometry:
            for coord in atom:
                coordLine = coordLine + str(round(coord * AU_TO_ANG, 3)) + " "
        gradLine = "FORCES "
        for atom in result.return_result:
            for coord in atom:
                gradLine = gradLine + str(coord) + " "
        f.write(coordLine + "\n")
        f.write("ENERGY " + str(result.properties.return_energy) + "\n")
        f.write(gradLine + "\n")
        f.write("\n")
        i += 1

with open("all.mdcrd", "w") as f:
    for result in results:
        tokenCounter = 1
        for atom in result.molecule.geometry:
            for coord in atom:
                f.write("%8.3f" % float(coord * AU_TO_ANG))
                if tokenCounter == 10:
                    f.write("\n")
                    tokenCounter = 1
                else:
                    tokenCounter += 1
        if tokenCounter != 1:
            f.write("\n")
