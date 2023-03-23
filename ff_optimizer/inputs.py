import yaml
from pathlib import Path
import errno
from os import strerror
from dataclasses import dataclass

def checkForFile(directory, f, isFile = True):
    if type(directory) is str:
        directory = Path(directory)
    path = directory / f
    if isFile:
        if not path.is_file():
            raise FileNotFoundError(errno.ENOENT, strerror(errno.ENOENT), path.absolute())
    else:
        if not path.is_dir():
            raise FileNotFoundError(errno.ENOENT, strerror(errno.ENOENT), path.absolute())

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

    def pathify(self):
        self.dynamicsdir = Path(self.dynamicsdir)
        self.optdir = Path(self.optdir)
        self.sampledir = Path(self.sampledir)

    def __post_init__(self):
        self.pathify()
        checkForFile(self.dynamicsdir, "", False)
        checkForFile(self.optdir, "", False)
        checkForFile(self.sampledir, "", False)
        checkForFile(self.dynamicsdir, self.coors)
        self.checkParams()

    def checkParams(self):
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
                    raise ValueError(f"End frame {self.end} must be >= split frame {self.split}")
        if self.resp > 1 or self.resp < 0:
            raise ValueError("RESP weight must be in [0, 1]")
        if self.nvalids < 1:
            raise ValueError("Must use at least one validation set (for now)")
        if self.activelearning < 1:
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
        super().pathify()
        checkForFile(self.sampledir, self.trainmdin)
        checkForFile(self.sampledir, self.validmdin)
        self.checkParams()

    def checkParams(self):
        if self.mmengine != "amber":
            raise ValueError(f"MMEngine {self.mmengine} is not implemented")
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
        super().pathify()
        checkForFile(self.sampledir, self.tctemplate)
        if self.qmengine == "sbatch":
            checkForFile(self.sampledir, self.sbatchtemplate)
        checkForFile(self.sampledir, self.tctemplate_backup)
        self.checkParams()

    def checkParams(self):
        qmengines = ["chemcloud", "sbatch", "debug"]
        if self.qmengine not in qmengines:
            raise ValueError(f"QMEngine {self.qmengine} has not been implemented")

@dataclass
# Contains keywords needed by MMEngine
class OptInput(GeneralInput):
    resppriors : int = 1
    opt0 : str = "opt_0.in"
    valid0 : str = "valid_0.in"
    stride : int = 50
    tcout : str = "tc.out"

    def __post_init__(self):
        super().pathify()
        checkForFile(self.optdir, self.opt0)
        checkForFile(self.optdir, self.valid0)
        checkForFile(self.optdir, "setup.leap")
        checkForFile(self.optdir, "conf.pdb")
        checkForFile(self.dynamicsdir, self.tcout)
        self.checkParams()

    def checkParams(self):
        if self.resppriors != 1 and self.resppriors != 2:
            raise ValueError("RESP prior mode must be either 1 or 2")
        if self.stride < 1:
            raise ValueError("Stride for creating initial data from dynamics must be at least 1")
        

@dataclass
class Inputs:
    generalInput : GeneralInput
    mmInput : MMInput
    qmInput : QMInput
    optInput : OptInput

    @classmethod
    def fromYaml(cls, inputs):
        if Path(inputs).is_file():
            with open(inputs, 'r') as f:
                di = yaml.safe_load(f)
        else:
            di = {}
        # If using all defaults
        if di is None:
            di = {}
        # If using all defaults for one set of inputs
        keys = ["general", "mm", "qm", "opt"]
        for key in keys:
            if key not in di.keys():
                di[key] = {}
        generalInput = GeneralInput(**di['general'])
        mmInput = MMInput(**{**di['mm'], **di["general"]})
        qmInput = QMInput(**{**di['qm'], **di["general"]})
        optInput = OptInput(**{**di['opt'], **di["general"]})
        return cls(generalInput, mmInput, qmInput, optInput)
