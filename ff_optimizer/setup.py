import os
import subprocess
from pathlib import Path
from shutil import copyfile
from textwrap import dedent

from .inputs import Input
from .utils import *
from chemcloud import CCClient
from qcio import ProgramInput, Structure


class Setup:
    def __init__(self, xyz: str | Path, charge: int = 0, resp: str = "am1"):
        """
        Initialize Setup class and set some internal variables.

        Args:
            xyz (str | Path): xyz file for molecule to set up
            charge (int): molecular charge for that molecule
            resp (str): how to obtain charges
        """
        self.xyz = Path(xyz).absolute()
        readXYZ(
            self.xyz
        )  # easy way to check that xyz is formatted correctly and exists
        checkForAmber()
        self.home = Path(".").absolute()
        self.charge = charge
        self.spinMult = self.getSpinMult()
        self.resp = resp
        self.makeInput()

    def getSpinMult(self):
        """
        Obtains a guess for the spin multiplicity of the molecule based on 
        elements and charges. Defaults to the lowest reasonable option
        (singlet or doublet)

        Returns:
            spinmult (int): guess spin multiplicity for molecule
        """
        elements = loadElements()
        _, symbols = readXYZ(self.xyz, True)
        electrons = 0
        for symbol in symbols:
            atomicNumber = elements[symbol]["atomic_number"]
            electrons += atomicNumber
        electrons -= self.charge
        if electrons % 2 == 0:
            spinmult = 1
        else:
            spinmult = 2
        return spinmult

    def runResp(self):
        """
        Run RESP to obtain partial charges for the molecule. Then, change
        the charges in the mol2 file.
        """

        if self.resp == "am1":
            return
        elif self.resp == "chemcloud":
            charges = self.chemcloudResp()
        elif self.resp == "local":
            charges = self.localResp()
        else:
            raise ValueError("Invalid RESP method")
        self.changeCharges(charges, self.mol2)


    def chemcloudResp(self):
        """
        Run a ChemCloud job to obtain RESP results.
        """

        mod = {}
        mod["method"] = "hf"
        mod["basis"] = "6-31gs"
        keywords = {}
        keywords["resp"] = "yes"
        keywords["esp_grid_dens"] = 20.0
        tempMol = Structure.open(self.xyz)
        kwargs = tempMol.model_dump()
        kwargs["charge"] = self.charge
        #kwargs["spinmult"] = self.spinMult
        mol = Structure(**kwargs)
        inp = ProgramInput(structure=mol, model=mod, calctype="energy", keywords=keywords)
        client = CCClient()
        futureResult = client.compute("terachem", inp)
        result = futureResult.get()
        try:
            assert isinstance(result.stdout, str)
        except:
            raise RuntimeError("Chemcloud RESP calculation failed")
        charges = self.readCharges(result.stdout.split("\n"))
        return charges

    def localResp(self):
        """
        Run TeraChem locally to obtain RESP results.
        """

        with open("resp.in", "w") as f:
            f.write(
                dedent(
                    f"""\
                basis         6-31gs  
                method        hf
                charge        {self.charge}        
                run           energy 
                spinmult      {self.spinMult}        
                resp          yes
                esp_grid_dens 20.0
                coordinates   {self.xyz}
            """
                )
            )
        checkForTerachem(True)
        os.system("terachem resp.in > resp.out")
        try:
            with open("resp.out", "r") as f:
                lines = list(f.readlines())
        except:
            raise RuntimeError("Debug RESP calculation failed")
        charges = self.readCharges(lines)
        return charges
    
    def readCharges(self, lines):
        """
        Reads charges from lines of a TeraChem RESP output file.

        Args:
            lines (list[str]): lines of a TC RESP file
        Returns:
            charges (list[str]): RESP charges from TC output
        """

        charges = []
        i = 0
        for line in lines:                                                           
            if "ESP restrained charges:" in line:                                     
                break                                                                
            i = i + 1                                                                
        for line in lines[i+3:]:                                                     
            if "-----------------------------------------------------------" in line:
                break                                                                
            splitLine = line.split()                                                 
            charges.append(splitLine[4])
        return charges


    def changeCharges(self, charges, mol2):
        """
        Replaces charges in a mol2 file.

        Args:
            charges (list[str]): list of atomic partial charges
            mol2 (str | Path): mol2 file to have charges replaced
        """
        with open(mol2, 'r') as f:                                                                         
            lines = f.readlines()                                                                              
            i = 0                                                                                              
            start = 0                                                                                          
            end = 0                                                                                            
            for line in lines:                                                                                 
                if "@<TRIPOS>ATOM" in line:                                                                    
                    start = i + 1                                                                              
                if "@<TRIPOS>BOND" in line:                                                                    
                    end = i                                                                                    
                    break                                                                                      
                i = i + 1                                                                                      
            if end - start != len(charges):                                                                    
                raise RuntimeError("Size of molecule in mol2 file differs from size of molecule in tcout file")
        with open(mol2, 'w') as w:                                                                         
            for line in lines[:start]:                                                                         
                w.write(line)                                                                                  
            for line, charge in zip(lines[start:end],charges):                                                 
                oldCharge = line.split()[8]                                                                    
                w.write(line.replace(oldCharge,charge))                                                        
            for line in lines[end:]:                                                                           
                w.write(line)                                                                                  


    def makeInput(self):
        """ 
        Create Input object with defaults for setup. This object will be 
        returned later if the user wants to run an optimization with it.
        """
        self.inp = Input(**{"skipchecks": True})
        # This setup code is meant for amber and chemcloud
        self.inp.initialtraining = False
        self.inp.mmengine = "amber"
        self.inp.qmengine = "chemcloud"
        self.inp.dynamicsdir = "."
        self.inp.coors = self.xyz
        self.inp.conformers = self.inp.coors

    def setup(self, opt=True, mm=True, qm=True):
        """
        Set up the optimization, MM, and QM calculations.

        Args:
            opt (bool, optional): Whether to set up optimization. Defaults to True.
            mm (bool, optional): Whether to set up MM. Defaults to True.
            qm (bool, optional): Whether to set up QM. Defaults to True.

        Returns:
            Input: A new Input object with the setup configuration.
        """
        self.inp.optdir.mkdir(parents=True, exist_ok=True)
        self.inp.sampledir.mkdir(parents=True, exist_ok=True)
        if opt:
            self.setupOpt()
        if mm:
            self.setupMM()
        if qm:
            self.setupQM()
        # when we initialize the new input, it'll check for files and folder
        os.chdir(self.home)
        self.inp.toYaml("input.yaml")
        newInp = Input.fromYaml("input.yaml")
        return newInp

    def setupOpt(self):
        """
        Sets up the optimization directory for ff-opt. This includes running
        antechamber to make the mol2 and frcmod files and writing the 
        ForceBalance input files.
        """
        print("\n\n**** Setting up force field and ForceBalance files ****")
        setupDir = self.inp.optdir / "setup"
        print(f"Making FF setup directory {str(setupDir)}")
        setupDir.mkdir(exist_ok=True, parents=True)
        print(os.getcwd(), str(setupDir))
        print(str(self.xyz))
        copyfile(self.xyz, setupDir / self.xyz.name)
        os.chdir(self.inp.optdir / "setup")
        frcmod, mol2, resname = self.setupFF()
        os.chdir(self.inp.optdir)
        self.setupForceBalance(frcmod, mol2, resname)
        self.editFrcmod(frcmod)
        charge = self.getCharge(mol2)
        if charge != self.charge:
            raise ValueError(
                f"Specified charge {self.charge} does not match charge {charge} in mol2 file"
            )

    def setupMM(self):
        """
        Writes default sander input files for MD sampling.
        """
        print("\n\n**** Creating Amber MD input files ****")
        print("Head to https://ambermd.org/index.php for documentation")
        os.chdir(self.inp.sampledir)
        self.writeMDFiles()

    def setupQM(self):
        """
        Writes default TeraChem input files for QM calculations.
        """
        os.chdir(self.inp.sampledir)
        self.writeTCFiles()

    def getCharge(self, mol2) -> int:
        """
        Determines total molecular charge from a mol2 file.

        Args:
            mol2 (str | Path): mol2 file containing atomic partial charges

        Returns:
            int: the (rounded) molecular charge from the mol2 file.
        """
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
        diff = charge - round(charge)
        if abs(diff) > 1e-6:
            print(f"WARNING: residual charge of %13.9f in mol2 file" % diff)
            print("Consider adjusting a charge by %13.9f to cancel it" % -diff)
        return round(charge)

    def writeTCFiles(self):
        """
        Writes default TC input files.
        """
        print("\n\n**** Writing TeraChem input files ****")
        print("Head to http://www.petachem.com/doc/userguide.pdf for documentation")
        print(
            f"Defaulting to wB97X-D3/def2-svp and a spin multiplicity of {self.spinMult}"
        )
        with open(self.inp.tctemplate, "w") as f:
            f.write(
                dedent(
                    f"""\
                basis         def2-svp  
                method        wb97xd3
                charge        {self.charge}        
                run           gradient 
                dftd          d3       
                spinmult      {self.spinMult}        
            """
                )
            )
        with open(self.inp.tctemplate_backup, "w") as f:
            f.write(
                dedent(
                    f"""\
                basis         def2-svp  
                method        wb97xd3
                charge        {self.charge}        
                run           gradient 
                dftd          d3       
                spinmult      {self.spinMult}        
                threall       1.0e-14
                diismaxvecs   40
                maxit         200
            """
                )
            )

    def writeMDFiles(self):
        """
        Writes default sander input files.
        """
        print("Writing inputs to equilibrate to 400K")
        for i in range(8):
            with open(f"heat{i+1}.in", "w") as f:
                f.write(
                    dedent(
                        f"""\
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
                """
                    )
                )
        print("Writing input for sampling trajectory at 400K")
        print("The number of samples is nstlim/ntwx = 200")
        print(
            "(nstlim = number of MD steps, ntwx = frequency of writing out coordinates)"
        )
        with open("md.in", "w") as f:
            f.write(
                dedent(
                    f"""\
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
            """
                )
            )

    def setupFF(self) -> tuple:
        """
        Runs antechamber to produce frcmod and mol2 files from the xyz file.

        Returns:
            tuple: a tuple containing the names of the frcmod file, the mol2
            file, and the residue id for the molecule.
        """
        name = self.xyz.name.replace(".xyz", "")
        resname = name[:3].upper()
        coords, atoms = readXYZ(self.xyz, True)
        pdb = name + ".pdb"
        mol2 = name + ".mol2"
        self.mol2 = mol2
        frcmod = name + ".frcmod"
        print(f"Converting xyz {self.xyz} to a pdb with residue name {resname}")
        writePDB(coords, pdb, atoms=atoms, resname=resname, template=None)

        print(f"Creating Amber forcefield files")
        print("Making mol2 file with antechamber")
        if self.resp == "am1":
            method = "bcc"
            print(
                "Using AM1-BCC to create charges; this should be done with HF/6-31g* RESP"
            )
        else:
            method = "dc"
        subprocess.check_output(
            [
                "antechamber",
                "-i",
                pdb,
                "-fi",
                "pdb",
                "-o",
                mol2,
                "-fo",
                "mol2",
                "-c",
                method,
                "-nc",
                str(self.charge),
            ]
        )
        self.runResp()
        print("Making frcmod file with parmchk2")
        os.system(f"parmchk2 -i {mol2} -o {frcmod} -f mol2 -s 2 -a Y")
        copyfile(mol2, self.inp.optdir / mol2)
        copyfile(frcmod, self.inp.optdir / frcmod)
        copyfile(pdb, self.inp.optdir / "conf.pdb")
        return frcmod, mol2, resname

    def setupForceBalance(self, frcmod: str | Path, mol2: str | Path, resname: str):
        """
        Writes ForceBalance and tleap input files.

        Args:
            frcmod (str | Path): frcmod file for molecule to be optimized
            mol2 (str | Path): mol2 file for molecule to be optimized
            resname (str): residue ID used for molecule in pdb files. 
        """
        print("Writing ForceBalance input files")
        print(
            "Head to https://leeping.github.io/forcebalance/doc/html/index.html for documentation"
        )
        with open(self.inp.opt0, "w") as f:
            f.write(
                dedent(
                    f"""\
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
                """
                )
            )

        with open(self.inp.valid0, "w") as f:
            f.write(
                dedent(
                    f"""\
                $options                                 
                jobtype             single               
                forcefield          {frcmod} {mol2}
                backup              false                
                $end                               
                """
                )
            )

        with open("setup.leap", "w") as f:
            f.write(
                dedent(
                    f"""\
                source leaprc.gaff
                loadamberparams {frcmod}
                {resname} = loadmol2 {mol2}
                x = loadpdb conf.pdb
                saveamberparm x ff.prmtop coord.inpcrd
                quit                                  
                """
                )
            )

    def editFrcmod(self, frcmod):
        """
        Edit the frcmod file to optimize all bonds, angles, and torsions.

        Args:
            frcmod (str): Name of the frcmod file.
        """
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
