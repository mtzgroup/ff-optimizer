import yaml
from pathlib import Path
import errno
from os import strerror

def checkForFile(directory, f, isFile = True):
    path = directory / f
    if isFile:
        if not path.is_file():
            raiseFileNotFoundError(errno.ENOENT, strerror(errno.ENOENT), path.absolute)
    else:
        if not path.is_dir():
            raiseFileNotFoundError(errno.ENOENT, strerror(errno.ENOENT), path.absolute)

@dataclass 
# Contains keywords needed by more than one of the engines or Model
class GeneralInput:
    dynamicsdir : str = "0_dynamics"
    optdir : str = "1_opt"
    sampledir : str = "2_sampling"
    coors : str = "coors.xyz"
    start : int = None
    end : int = None
    split : int = None
    resp : float = 0
    nvalids : int = 1
    restart : bool = False
    activelearning : int = 1
    maxcycles : int = 30

    def __post_init__(self):
        self.dynamicsdir = Path(self.dynamicsdir)
        checkForFile(self.dynamicsdir, "", False)
        self.optdir = Path(self.optdir)
        checkForFile(self.optdir, "", False)
        self.sampledir = Path(self.sampledir)
        checkForFile(self.sampledir, "", False)
        checkForFile(self.dynamicsdir, self.coors)

        if self.start is not None:
            if self.start < 0:
                raise ValueError(f"Start frame {self.start} must be > 0")
            if self.end is not None:
                if self.end <= self.start:
                    raise ValueError(f"Start frame {self.start} must be less than end frame {self.end}")
            if self.split is not None:
                if self.split <= self.start:
                    raise ValueError(f"Start frame {self.start} must be less than split frame {self.split}")

        if self.end is not None:
            if self.end < 0:
                raise ValueError(f"End frame {self.end} must be > 0")
            if self.split is not None:
                if self.split > self.end:
                    raise ValueError(f"End frame {self.end must be >= split frame {self.split}")
        if self.resp > 1 or self.resp < 0:
            raise ValueError("RESP weight must be in [0, 1]")
        if self.nvalids < 1:
            raise ValueError("Must use at least one validation set (for now)")
        if self.activeLearning < 1:
            raise ValueError("Must have at least one model")
        if self.maxcycles < 0:
            raise ValueError("So you just don't want me to do anything? Check maxcycles.")

@dataclass
# Contains keywords needed by MMEngine
class MMInput(GeneralInput):
    mmengine : str = "amber"
    trainmdin : str = "md.in"
    validmdin : str = "md.in"
    conformers : str = None
    conformersperset : int = 1

    def __post_init__(self):
        if self.mmengine != "amber":
            raise ValueError(f"MMEngine {self.mmengine} is not implemented")
            checkForFile(self.sampledir, self.trainmdin)
            checkForFile(self.sampledir, self.validmdin)
            if self.conformers is None:
                self.conformers = self.coors
            else:
                checkForFile(self.sampledir, self.conformers)
            if self.conformersperset < 1:
                raise ValueError("Must be at least one conformer per set")

@dataclass
# Contains keywords needed by QMEngine
class QMInput(GeneralInput):
    qmengine : str = "chemcloud" 
    tctemplate : str = "tc_template.in"
    sbatchtemplate : str = "sbatch_template.sh"
    tctemplate_backup : str = "tc_template_backup.in"

    def __post_init__(self):
        qmengines = ["chemcloud", "sbatch", "debug"]
        if self.qmengine is not in qmengines:
            raise ValueError(f"QMEngine {self.qmengine} has not been implemented")
        checkForFile(self.sampledir, self.tctemplate)
        if self.qmengine == "sbatch":
            checkForFile(self.sampledir, self.sbatchtemplate)
        checkForFile(self.sampledir, self.tctemplate_backup)

@dataclass
# Contains keywords needed by MMEngine
class OptInput(GeneralInput):
    resppriors : int = 1
    opt0 : str = "opt_0.in"
    valid0 : str = "valid_0.in"
    stride : int = 50
    tcout : str = "tc.out"

    def __post_init__(self):
        if resppriors != 1 and resppriors != 2:
            raise ValueError("RESP prior mode must be either 1 or 2")
        checkForFile(self.optdir, self.opt0)
        checkForFile(self.optdir, self.valid0)
        checkForFile(self.optdir, "setup.leap")
        checkForFile(self.optdir, "conf.pdb")
        checkForFile(self.dynamicsdir, self.tcout)
        if stride < 1:
            raise ValueError("Stride for creating initial data from dynamics must be at least 1")

@dataclass
class Inputs:
    generalInput : GeneralInput
    mmInput : MMInput
    qmInput : QMInput
    optInput : OptInput

    @classmethod
    def fromYaml(cls, inputs):
        di = yaml.load(inputs)
        generalInput = GeneralInput(**di['General'])
        mmInput = MMInput.fromDict(**(di['MM'] + di["General"]))
        qmInput = QMInput.fromDict(**(di['QM'] + di["General"]))
        optInput = OptInput.fromDict(**(di['Opt'] + di["General"]))
