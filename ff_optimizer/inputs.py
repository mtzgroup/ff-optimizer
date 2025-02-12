import errno
import os
from dataclasses import dataclass, field
from os import strerror
from pathlib import Path

import yaml

from . import utils


def checkForFile(f: Path, isFile: bool = True):
    """
    Check if a file or directory exists.

    Args:
        f (Path): Path to the file or directory to check.
        isFile (bool): If True, check for a file; if False, check for a directory.

    Raises:
        FileNotFoundError: If the file or directory does not exist.
    """
    if isFile:
        if not f.is_file():
            raise FileNotFoundError(errno.ENOENT, strerror(errno.ENOENT), f.absolute())
    else:
        if not f.is_dir():
            raise FileNotFoundError(errno.ENOENT, strerror(errno.ENOENT), f.absolute())


def checkDirectory(directory: Path, fs: list):
    """
    Check if a directory exists and contains specified files.

    Args:
        directory (Path): Path to the directory to check.
        fs (list): List of filenames to check for in the directory.

    Raises:
        FileNotFoundError: If the directory or any of the specified files do not exist.
    """
    if isinstance(directory, str):
        directory = Path(directory)
    checkForFile(directory, False)
    for f in fs:
        checkForFile(directory / f, True)


@dataclass
class Input:
    """
    Input class containing configuration parameters for the force field optimization process.
    """

    # General parameters
    dynamicsdir: Path = field(
        default=None,
        metadata={
            "comment": """
        Specifies the directory containing QM energies,
        gradients, and coordinates from an MD trajectory. If not
        provided, ff-opt will default to using conf.pdb in optdir
        to provide starting coordinates for MM sampling.
        """
        },
    )
    optdir: Path = field(
        default=Path("1_opt"),
        metadata={
            "comment": """
        Specifies the directory where ForceBalance optimization will
        be run; contains FB inputs and initial parameters.
        """
        },
    )
    sampledir: Path = field(
        default=Path("2_sampling"),
        metadata={
            "comment": """
        Specifies the directory where MM sampling and QM 
        calculations will be run; contains MM and QM inputs.
        """
        },
    )
    coors: str = field(
        default="coors.xyz",
        metadata={
            "comment": """
        Specifies the name of the coordinates file in dynamicsdir.
        """
        },
    )
    start: int = field(
        default=None,
        metadata={
            "comment": """
        Specifies the first frame of the coors file to sample from
        to create the initial dataset. 0-indexed. If not provided, 
        default is the first frame.
        """
        },
    )
    end: int = field(
        default=None,
        metadata={
            "comment": """
        Specifies the last frame of the coors file to sample from. 
        0-indexed. If not provided, defaults to last frame in coors file.
        """
        },
    )
    split: int = field(
        default=None,
        metadata={
            "comment": """
        Divides the coors file in two at the provided frame number.
        If provided, initial coordinates for MM MD for the training
        set will come from before split, and ICs for the validation
        set will come from after. 0-indexed.
        """
        },
    )
    resp: float = field(
        default=0,
        metadata={
            "comment": """
        Specifies the weight of RESP fits for the charges in the
        objective function. If nonzero, QM calculations will
        include RESP fits, and parameter optimization will attempt
        to match the ESP in those QM calculations. Can severely
        slow down force field fitting.
        """
        },
    )
    nvalids: int = field(
        default=1,
        metadata={
            "comment": """
        Specifies how many validation sets to produce and evaluate.
        """
        },
    )
    restart: bool = field(
        default=False,
        metadata={
            "comment": """
        Specifies whether or not to restart the force field
        optimization. If True, it will automatically determine
        where to restart.
        """
        },
    )
    activelearning: int = field(
        default=1,
        metadata={
            "comment": """
        Specifies how many different models to train. If
        activelearning > 1, then an active learning procedure is
        used to generate the new data points. In the MM MD step,
        each model generates a large set of candidate points. For
        each set of candidates, every model evaluates forces.
        Uncertainty (standard deviation) in the forces is computed
        for each set. Half the geometries chosen for QM
        calculations have the highest uncertainty of the set, and
        the other half are sampled randomly. These sets become the
        training set for the model which did the MM MD, and a
        validation set for the other models.
        """
        },
    )
    maxcycles: int = field(
        default=30,
        metadata={
            "comment": """
        Specifies how many cycles to run before stopping.
        """
        },
    )
    # MM-specific parameters
    mmengine: str = field(
        default="amber",
        metadata={
            "comment": """
         specifies which software will run MM MD. Options are
        "amber" and "openmm".
        If "amber", provided MD input files must be valid inputs
        for sander/pmemd/pmemd.cuda. The number of frames sampled
        from the dynamics is determined by nstlim / ntwx.
        If "openmm", provided MD input files must be valid python
        scripts which run OpenMM and satisfy the following:
        1. They accept as arguments (in order) the .prmtop and
            .rst7 files
        2. They perform all necessary equilibrations internally
            (the heat*in files for equilibration are skipped)
        3. They produce a {name}.nc file containing the sampled
            geometries, where {name}.rst7 is the rst7 file passed
            to the script. The coordinates in this netcdf file will
            be passed to the QM engine to augment the dataset.
        An example can be found in examples/openmm/2_sampling.
        """
        },
    )
    trainmdin: str = field(
        default="md.in",
        metadata={
            "comment": """
        Specifies the input file for MM MD. This input will be
        used to generate the training set.
        """
        },
    )
    validmdin: str = field(
        default="md.in",
        metadata={
            "comment": """
        Specifies the input file for MM MD. This input will be
        used to generate validation sets.
        """
        },
    )
    conformers: str = field(
        default=None,
        metadata={
            "comment": """
        Specifies the name of a conformers coordinates file in 
        dynamicsdir. If provided, initial coordinates for MM sampling 
        will be drawn from this file instead of the coordinates file.
        """
        },
    )
    conformersperset: int = field(
        default=1,
        metadata={
            "comment": """
        Specifies how many unique MM MD trajectories to run per dataset.
        """
        },
    )
    # QM-specific parameters
    qmengine: str = field(
        default="chemcloud",
        metadata={
            "comment": """
        Specifies which qmengine to use to run QM calculations.
        Options are "chemcloud", "slurm", and "debug".
        If "chemcloud", options are read from a valid TeraChem input
        file (tc_template.in) and sent to chemcloud.
        If "queue", a sbatch_template.sh file must be provided. This
        file allows the QM TeraChem calculations to be run by a
        SLURM scheduler.                                            
        If "debug", TeraChem will be run locally, one job at a time.
        """
        },
    )
    # opt-specific parameters
    resppriors: int = field(
        default=0,
        metadata={
            "comment": """
        Specifies an alternative to using RESP calculations
        without incorporating the RESP fit into the objective
        function. Instead, a prior is calculated for each atomic
        charge, and a penalty function is added to the objective.
        If 0, don't use RESP priors for optimizing charges.
        If 1, the width for each atom is determined by the
        distribution of RESP charges for that atom in the training
        set.
        If 2, the width will be determined based on the RESP and ESP
        charges for that atom in the training set.
        """
        },
    )
    stride: int = field(
        default=50,
        metadata={
            "comment": """
        Specifies stride in sampling initial training data from the
        QM trajectory in coors.
        """
        },
    )
    tcout: str = field(
        default="tc.out",
        metadata={
            "comment": """
        Specifies TeraChem output file with energies and gradients
        from the QM trajectory in coors.
        """
        },
    )
    initialtraining: bool = field(
        default=False,
        metadata={
            "comment": """
        Specifies whether or not to sample initial training data from 
        a QM trajectory and optimize parameters before beginning MM 
        sampling in the first iteration. If False, an xyz file
        must still be provided in dynamicsdir as starting coordinates
        for MM sampling, but tcout is not used.
        """
        },
    )
    validinitial: bool = field(
        default=False,
        metadata={
            "comment": """
        If true, compute validation set performance on the initial
        parameters at every iteration.
        """
        },
    )
    dryrun: bool = field(
        default=False,
        metadata={
            "comment": """
        If true, check all folders and files, initialize everything,
        then stop before running any calculations.
        """
        },
    )
    skipchecks: bool = field(
        default=False,
        metadata={
            "comment": """
        If true, skip checking for files and folders. Useful for
        debugging or if you just want to initialize an Input full
        of defaults for something. If you do this, you'll need to 
        call Input.__post_init__() later.
        """
        },
    )
    # there's no reason to change these, but it's useful to have them accessible
    # from this class)
    tctemplate: str = field(
        default="tc_template.in",
        metadata={
            "comment": """
        Name of template TeraChem input file in sampledir. Should contain
        all relevant keywords, but ff-opt will take care of the 
        coordinates keyword for you.
        """
        },
    )
    tctemplate_backup: str = field(
        default="tc_template_backup.in",
        metadata={
            "comment": """
        Name of backup template TeraChem input file in sampledir. Only
        used if the calculation with tctemplate fails first.
        """
        },
    )
    sbatchtemplate: str = field(
        default="sbatch_template.sh",
        metadata={
            "comment": """
        Name of slurm sbatch template submission script. Should contain
        all relevant SBATCH options, as well as a script which runs
        the TeraChem input files. See examples/slurm/2_sampling for an
        example.
        """
        },
    )
    opt0: str = field(
        default="opt_0.in",
        metadata={
            "comment": """
        Name of the template ForceBalance parameter optimization 
        input file. Should contain all relevant keywords and an 
        example target section which ff-opt will use as targets are added.
        """
        },
    )
    valid0: str = field(
        default="valid_0.in",
        metadata={
            "comment": """
        Name of the template ForceBalance validation set evaluation
        input file. An example target section is unnecessary here.
        """
        },
    )
    batchsize: int = field(
        default=10,
        metadata={
            "comment": """
        Number of calculations to bundle together to send to a 
        ChemCloud server. You shouldn't mess with this number, except
        to decrease it if you're running a very large system.
        """
        },
    )
    retries: int = field(
        default=3,
        metadata={
            "comment": """
        Number of times failed calculations on ChemCloud are retried.
        """
        },
    )
    patience: int = field(
        default=5,
        metadata={
            "comment": """
        Number of consecutive iterations where the convergence criterion 
        must be satisfied before ending the calculation early. 
        """
        },
    )
    cutoff: float = field(
        default=-1.0,
        metadata={
            "comment": """
        Cutoff for percentage change in validation criterion:
        If (\Chi_V(r_j, \Theta_j) - \Chi_V(r_j, \Theta_{j-1}) / 
            \Chi_V(r_j, \Theta_j) * 100 > cutoff, then we enter the
        patience period. This criterion should always be negative
        if the force field is improving, so cutoff should always be 
        a small, negative number.
        """
        },
    )

    @classmethod
    def fromYaml(cls, inputs: str | Path) -> "Input":
        """
        Create an Input object from a YAML file.

        Args:
            inputs (str | Path): Path to the YAML file.

        Returns:
            Input: An Input object with parameters from the YAML file.
        """
        inputs = Path(inputs)  # Convert to Path if string
        if inputs.is_file():
            with open(inputs, "r") as f:
                di = yaml.safe_load(f)
        else:
            di = {}
        # If using all defaults
        if di is None:
            di = {}
        return cls(**di)

    def toYaml(self, dest: str | Path):
        """
        Write the Input object parameters to a YAML file.

        Args:
            dest (str | Path): Path to the destination YAML file.
        """
        data = vars(self)
        # convert pathlib Paths to strs for YAML serialization
        for key, value in data.items():
            if isinstance(value, Path):
                data[key] = str(value)
        with open(dest, "w") as f:
            yaml.dump(data, f)

    def __post_init__(self):
        """
        Perform post-initialization checks and setup.
        """
        self.checkParams()
        self.pathify()
        if not self.skipchecks:
            self.checkFiles()
        if self.dryrun:
            self.maxcycles = -1

    def pathify(self):
        """
        Convert directory strings to absolute Path objects.
        """
        if self.dynamicsdir:
            self.dynamicsdir = Path(self.dynamicsdir).absolute()
        self.optdir = Path(self.optdir).absolute()
        self.sampledir = Path(self.sampledir).absolute()

    def setupDynamicsFolder(self):
        """
        If dynamicsdir is not provided, set it to optdir and set up
        coordinates within it.
        """
        if self.dynamicsdir:
            return
        self.dynamicsdir = self.optdir
        home = os.getcwd()
        os.chdir(self.dynamicsdir)
        utils.convertPDBtoXYZ("conf.pdb")
        os.chdir(home)
        self.coors = "conf.xyz"

    def checkFiles(self):
        """
        Check for the existence of required files and directories.

        Raises:
            FileNotFoundError: If required files or directories are missing.
        """
        optFiles = ["conf.pdb", "setup.leap", self.opt0, self.valid0]
        checkDirectory(self.optdir, optFiles)
        dynamicsFiles = []
        if self.conformers:
            dynamicsFiles.append(self.conformers)
        if self.initialtraining:
            dynamicsFiles.append(self.tcout)
        else:
            self.setupDynamicsFolder()
        dynamicsFiles.append(self.coors)
        checkDirectory(self.dynamicsdir, dynamicsFiles)
        sampleFiles = [
            self.tctemplate,
            self.tctemplate_backup,
            self.trainmdin,
            self.validmdin,
        ]
        if self.qmengine == "slurm":
            sampleFiles.append(self.sbatchtemplate)
        checkDirectory(self.sampledir, sampleFiles)

    def checkParams(self):
        """
        Validate input parameters.

        Raises:
            ValueError: If any input parameters are invalid.
        """
        if self.initialtraining and not self.dynamicsdir:
            raise ValueError(
                "Initial training requires providing a dynamics directory (via keyword dynamicsdir) with a coordinates and terachem output file"
            )
        if self.resppriors != 1 and self.resppriors != 2 and self.resppriors != 0:
            raise ValueError("RESP prior mode must be either 0, 1, or 2")
        if self.stride < 1:
            raise ValueError(
                "Stride for creating initial data from dynamics must be at least 1"
            )
        qmengines = ["chemcloud", "slurm", "debug"]
        if self.qmengine not in qmengines:
            raise ValueError(f"QMEngine {self.qmengine} has not been implemented")
        mmengines = ["amber", "openmm"]
        if self.mmengine not in mmengines:
            raise ValueError(f"MMEngine {self.mmengine} is not implemented")
        if self.conformersperset < 1:
            raise ValueError("Must be at least one conformer per set")
        if self.start is not None:
            if self.start < 0:
                raise ValueError(f"Start frame {self.start} must be > 0")
            if self.end is not None:
                if self.end <= self.start:
                    raise ValueError(
                        f"Start frame {self.start} must be less than end frame {self.end}"
                    )
            if self.split is not None:
                if self.split <= self.start:
                    raise ValueError(
                        f"Start frame {self.start} must be less than split frame {self.split}"
                    )

        if self.end is not None:
            if self.end < 0:
                raise ValueError(f"End frame {self.end} must be > 0")
            if self.split is not None:
                if self.split > self.end:
                    raise ValueError(
                        f"End frame {self.end} must be >= split frame {self.split}"
                    )
        if self.resp > 1 or self.resp < 0:
            raise ValueError("RESP weight must be in [0, 1]")
        if self.nvalids < 1:
            raise ValueError("Must use at least one validation set (for now)")
        if self.activelearning < 1:
            raise ValueError("Must have at least one model")
