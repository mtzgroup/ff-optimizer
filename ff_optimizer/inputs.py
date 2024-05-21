import errno
from dataclasses import dataclass
from os import strerror
from pathlib import Path

import yaml


def checkForFile(f, isFile=True):
    if isFile:
        if not f.is_file():
            raise FileNotFoundError(errno.ENOENT, strerror(errno.ENOENT), f.absolute())
    else:
        if not f.is_dir():
            raise FileNotFoundError(errno.ENOENT, strerror(errno.ENOENT), f.absolute())


def checkDirectory(directory, fs):
    if type(directory) is str:
        directory = Path(directory)
    checkForFile(directory, False)
    for f in fs:
        checkForFile(directory / f, True)


# Flat list of options for now
@dataclass
class Input:
    # General parameters
    dynamicsdir: str = "0_dynamics"
    optdir: str = "1_opt"
    sampledir: str = "2_sampling"
    coors: str = "coors.xyz"
    start: int = None
    end: int = None
    split: int = None
    resp: float = 0
    nvalids: int = 1
    restart: bool = False
    activelearning: int = 1
    maxcycles: int = 30
    # MM-specific parameters
    mmengine: str = "amber"
    trainmdin: str = "md.in"
    validmdin: str = "md.in"
    conformers: str = None
    conformersperset: int = 1
    # QM-specific parameters
    qmengine: str = "chemcloud"
    # opt-specific parameters
    resppriors: int = 0
    stride: int = 50
    tcout: str = "tc.out"
    initialtraining: bool = True
    validinitial: bool = False
    easymode: str = None
    setuponly: bool = False

    # there's no reason to change these, but they're here just in case
    # (so they're not documented, but it's useful to have them accessible
    # from this class)
    tctemplate: str = "tc_template.in"
    tctemplate_backup: str = "tc_template_backup.in"
    sbatchtemplate: str = "sbatch_template.sh"
    opt0: str = "opt_0.in"
    valid0: str = "valid_0.in"
    batchsize: int = 10
    retries: int = 3
    patience: int = 5
    cutoff: float = -1.0

    @classmethod
    def fromYaml(cls, inputs):
        if Path(inputs).is_file():
            with open(inputs, "r") as f:
                di = yaml.safe_load(f)
        else:
            di = {}
        # If using all defaults
        if di is None:
            di = {}
        return cls(**di)

    def toYaml(self, dest):
        data = vars(self)
        # convert pathlib Paths back to strs
        for key in data.keys():
            if isinstance(data[key], Path):
                data[key] = str(data[key])
        with open(dest, "w") as f:
            yaml.dump(vars(self), f)

    def __post_init__(self):
        # in easymode, files/folders aren't set up yet
        self.pathify()
        if self.easymode is None:
            self.checkParams()
            self.checkFiles()
        else:
            checkDirectory(".", [self.easymode])
        if self.setuponly:
            self.maxcycles = -1

    def pathify(self):
        self.dynamicsdir = Path(self.dynamicsdir).absolute()
        self.optdir = Path(self.optdir).absolute()
        self.sampledir = Path(self.sampledir).absolute()

    def checkFiles(self):
        if self.initialtraining:
            dynamicsFiles = [self.tcout, self.coors, self.conformers]
            checkDirectory(self.dynamicsdir, dynamicsFiles)
        optFiles = ["conf.pdb", "setup.leap", "opt_0.in", "valid_0.in"]
        checkDirectory(self.optdir, optFiles)
        sampleFiles = ["tc_template.in", "tc_template_backup.in"]
        if self.qmengine == "slurm":
            sampleFiles.append(self.sbatchtemplate)
        checkDirectory(self.sampledir, sampleFiles)

    def checkParams(self):
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
        if self.conformers is None:
            self.conformers = self.coors
        else:
            checkForFile(self.sampledir, self.conformers)
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
        if self.maxcycles < 0:
            raise ValueError(
                "So you just don't want me to do anything? Check maxcycles."
            )
