from .inputs import Input
from pathlib import Path
from shutil import copyfile
import os
from .utils import readXYZ, writePDB
from textwrap import dedent

def setup(inp, opt=True, mm=True, qm=True):
    home = Path(".").absolute()
    inp.optdir.mkdir(parents=True, exist_ok=True)
    inp.sampledir.mkdir(parents=True, exist_ok=True)
    # This setup code is meant for amber and chemcloud
    inp.initialtraining = False
    inp.mmengine = "amber"
    inp.qmengine = "chemcloud"
    inp.dynamicsdir = "."
    inp.coors = inp.easymode
    if opt:
        setupOpt(inp)
    if mm:
        setupMM(inp)
    if qm:
        setupQM(inp)
    # when we initialize the new input, it'll check for files and folder
    inp.easymode = None
    os.chdir(home)
    inp.toYaml("new_input.yaml")
    newInp = Input.fromYaml("new_input.yaml")
    return newInp

def setupOpt(inp):
    xyz = inp.easymode
    setupDir = inp.optdir / "setup"
    print(f"Making FF setup directory {str(setupDir)}")
    setupDir.mkdir(exist_ok=True)
    copyfile(xyz, setupDir / xyz)
    os.chdir(inp.optdir / "setup")
    frcmod, mol2, resname = setupFF(xyz, inp.optdir)
    os.chdir(inp.optdir)
    setupForceBalance(inp, frcmod, mol2, resname)
    editFrcmod(frcmod)

def setupMM(inp):
    print("Creating Amber MD input files")
    print("Head to https://ambermd.org/index.php for documentation")
    os.chdir(inp.sampledir)
    writeMDFiles()

def setupQM(inp):
    charge = getCharge(inp.optdir)
    os.chdir(inp.sampledir)
    writeTCFiles(inp, charge)

def getCharge(optdir):
    os.chdir(optdir)
    for f in os.listdir():
        if f.endswith("mol2"):
            mol2 = f
            break
    charge = 0.0
    inAtoms = False
    with open(mol2, "r") as f:
        for line in f.readlines():
            if inAtoms:
                splitLine = line.split()
                if len(splitLine) == 9:
                    charge += float(splitLine[-1])
                else:
                    inAtoms = False
            if "@<TRIPOS>ATOM" in line:
                inAtoms = True
    print(charge)
    return round(charge)

def writeTCFiles(inp, charge):
    print("Writing TeraChem input files")
    print("Head to http://www.petachem.com/doc/userguide.pdf for documentation")
    print("Defaulting to wB97X-D3/def2-svp and a spin multiplicity of 1")
    with open(inp.tctemplate, "w") as f:
        f.write(dedent(f"""\
            basis         def2-svp  
            method        wb97x-d3
            charge        {charge}        
            run           gradient 
            dftd          d3       
            spinmult      1        
        """))
    with open(inp.tctemplate_backup, "w") as f:
        f.write(dedent(f"""\
            basis         def2-svp  
            method        wb97x-d3
            charge        {charge}        
            run           gradient 
            dftd          d3       
            spinmult      1        
            threall       1.0e-14
            diismaxvecs   40
            maxit         200
        """))

def writeMDFiles():
    print("Writing inputs to equilibrate to 400K")
    for i in range(8):
        with open(f"heat{i+1}.in", "w") as f:
            f.write(dedent(f"""\
                heat molecule
                 &cntrl
                  imin=0,irest=0,ntx=1,
                  nstlim=10000,dt=0.001,
                  ntc=1,ntf=1,
                  cut=999, ntb=0,
                  ntpr=10000, ntwx=0,
                  ntwv=0,
                  ntt=3, gamma_ln=10.00,
                  tempi={i*50}, temp0={(i+1)*50},
                  ntr=0,
                  nmropt=0
                /
            """))
    print("Writing input for sampling trajectory at 400K")
    print("The number of samples is nstlim/ntwx = 200")
    print("(nstlim = number of MD steps, ntwx = frequency of writing out coordinates)")
    with open("md.in", "w") as f:
        f.write(dedent(f"""\
            sample conformations
            &cntrl
              imin=0,irest=1,ntx=5,
              nstlim=800000,dt=0.001,
              ntc=1,ntf=1,
              cut=999, ntb=0,
              ntpr=100, ntwx=4000,
              ntwv=1000,
              ntt=3, gamma_ln=10.0,
              tempi=400.0, temp0=400.0,
              ntr=0,
              nmropt=0
            /
        """))

def setupFF(xyz, optdir):
    name = xyz.replace(".xyz", "")
    resname = name[:3].upper()
    coords, atoms = readXYZ(xyz, True)
    pdb = name + ".pdb"
    mol2 = name + ".mol2"
    frcmod = name + ".frcmod"
    print(f"Converting xyz {xyz} to a pdb with residue name {resname}")
    writePDB(coords, pdb, atoms=atoms, resname=resname, template=None)
    
    print(f"Creating Amber forcefield files")
    print("Making mol2 file with antechamber")
    print("Using AM1-BCC to create charges; this should be done with HF/6-31g* RESP")
    os.system(f"antechamber -i {pdb} -fi pdb -o {mol2} -fo mol2 -c bcc")
    print("Making frcmod file with parmchk2")
    os.system(f"parmchk2 -i {mol2} -o {frcmod} -f mol2 -s 2 -a Y")
    copyfile(mol2, optdir / mol2)
    copyfile(frcmod, optdir / frcmod)
    copyfile(pdb, optdir / "conf.pdb")
    return frcmod, mol2, resname

def setupForceBalance(inp, frcmod, mol2, resname):
    print("Writing ForceBalance input files")
    print("Head to https://leeping.github.io/forcebalance/doc/html/index.html for documentation")
    with open(inp.opt0, "w") as f:
        f.write(dedent(f"""\
            $options
            jobtype             NewtonRaphson
            forcefield          {frcmod} {mol2}
            penalty_additive    0.01
            penalty_type        L2
            verbose_options     true
            backup              false
            constrain_charge    true
            priors
            BONDSK                               : 1.0e+05
            BONDSB                               : 1.0e-01
            ANGLESK                              : 2.0e+02
            ANGLESB                              : 1.0e+02
            PDIHS1K                              : 1.0e+02
            PDIHS2K                              : 4.0e+01
            PDIHS3K                              : 1.0e+01
            IDIHSK                               : 2.0e+01
            /priors
            $end 

            $target
            type                abinitio_amber
            name                train_1
            amber_leapcmd       setup.leap
            $end                              
            """))

    with open(inp.valid0, "w") as f:
        f.write(dedent(f"""\
            $options                                 
            jobtype             single               
            forcefield          {frcmod} {mol2}
            backup              false                
            $end                               
            """))

    with open("setup.leap", "w") as f:
        f.write(dedent(f"""\
            source leaprc.gaff
            loadamberparams {frcmod}
            {resname} = loadmol2 {mol2}
            x = loadpdb conf.pdb
            saveamberparm x ff.prmtop coord.inpcrd
            quit                                  
            """))

def editFrcmod(frcmod):
    print("Labeling frcmod file to optimize all bonds, angles, and torsions")
    inBond = False
    inAngle = False
    inDihedral = False
    inImproper = False
    with open(frcmod, "r") as inF:
        with open("temp.txt", "w") as outF:
            for line in inF.readlines():
                if len(line.split()) < 3:
                    inBond = False
                    inAngle = False
                    inDihedral = False
                    inImproper = False
                if inBond:
                    outF.write(line[:-1] + " # PRM 1 2\n")
                elif inAngle:
                    outF.write(line[:-1] + " # PRM 1 2\n")
                elif inDihedral:
                    outF.write(line[:-1] + " # PRM 2\n")
                elif inImproper:
                    outF.write(line[:-1] + " # PRM 1\n")
                else:
                    outF.write(line)
                    if "BOND" in line:
                        inBond = True
                    if "ANGLE" in line:
                        inAngle = True
                    if "DIHE" in line:
                        inDihedral = True
                    if "IMPROPER" in line:
                        inImproper = True
    os.rename("temp.txt", frcmod)
