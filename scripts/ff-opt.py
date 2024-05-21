#!/usr/bin/env python

import sys
from textwrap import dedent
from time import perf_counter

from ff_optimizer import active_learning, inputs, model


def printHelp():
    summary = dedent(
        """\
    Summary of necessary directories and files
    Dynamics directory ("dynamicsdir", 0_dynamics) contains:
        -- QM dynamics coordinates ("coors", coors.xyz)
        -- TeraChem output file ("tcout", tc.out)
        -- Conformers file ("conformers", coors.xyz)
    FB Optimization directory ("optdir", 1_opt) contains:
        -- FB optimization input file (opt_0.in)
        -- FB validation input file (valid_0.in)
        -- PDB conformation file (conf.pdb)
        -- Tleap setup file (setup.leap)
        -- Parameter frcmod file (*.frcmod)
        -- Parameter mol2 file (*.mol2)
    MM sampling directory ("sampledir", 2_sampling) contains:
        -- Equilibration inputs (heat*in)
        -- MM sampling input(s) ("trainmdin", "validmdin", md.in)
        -- TeraChem input file for fire (tc_template.in)
        -- Backup TeraChem input file if job fails (tc_template_long.in)
        -- sbatch template file for queue (sbatch_template.sh)
    """
    )
    keywords = dedent(
        """\
    Keywords are to be provided in a .yaml file. Example:
    {
    "restart": False,
    "mmengine": "amber",
    "maxcycles": 30
    }
    For all keywords, the syntax is simply "keyword": value. Take care that 
    keywords have the the right type: strings should be in quotes, ints should
    not.

    List of all keywords:
    "dynamicsdir": str, specifies the directory containing QM energies, 
                    gradients, and coordinates from an MD trajectory.
                    Default = "0_dynamics"

    "optdir": str, specifies the directory where ForceBalance optimization will
                    be run; contains FB inputs and initial parameters.
                    Default = "1_opt"

    "sampledir": str, specifies the directory where MM sampling and QM 
                    calculations will be run; contains MM and QM inputs.
                    Default = "2_sampling"

    "easymode": str, specifies an xyz file for parameter optimization. The 
                    setup populates all files/folders with reasonable defaults,
                    and prepares everything for an amber/chemcloud calculation 
                    without initial training.
                    Default = None

    "initialtraining": bool, specifies whether or not to perform a parameter optimization
                    prior to beginning the first round of sampling. If so, the program
                    uses the MD trajectory in dynamicsdir to create an initial target
                    for optimization.
                    Default = True

    "validinitial": bool, specifies whether or not to evaluate the initial 
                    parameters on the validation set at each iteration.
                    Default = False

    "coors": str, specifies the name of the coordinates file in dynamicsdir.
                    Default = "coors.xyz"

    "start": int, specifies the first frame of the coors file to sample from 
                    to create the initial dataset. 0-indexed.
                    Default = 0

    "end": int, specifies the last frame of the coors file to sample from.
                    0-indexed.
                    If not provided, defaults to last frame in coors file.

    "split": int, divides the coors file in two at the provided frame number.
                    If provided, initial coordinates for MM MD for the training
                    set will come from before split, and ICs for the validation
                    set will come from after. 0-indexed.
                    Default = None

    "resp": float, specifies the weight of RESP fits for the charges in the 
                    objective function. If nonzero, QM calculations will
                    include RESP fits, and parameter optimization will attempt
                    to match the ESP in those QM calculations. Can severely
                    slow down force field fitting.
                    Default = 0

    "nvalids": int, specifies how many validation sets to produce and evaluate.
                    Default = 1

    "restart": bool, specifies whether or not to restart the force field 
                    optimization. If True, it will automatically determine
                    where to restart.
                    Default = False

    "activelearning": int, specifies how many different models to train. If
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
                    Default = 1

    "maxcycles": int, specifies how many cycles to run before stopping. 
                    Default = 30

    "mmengine": str, specifies which software will run MM MD. Options are 
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
                    An example can be found in tests/mmengine/md.py.

    "trainmdin": str, specifies the input file for MM MD. This input will be 
                    used to generate the training set.
                    Default = "md.in"

    "validmdin": str, specifies the input file for MM MD. This input will be
                    used to generate validation sets.
                    Default = "md.in"

    "conformers": str, specifies the name of a conformers coordinates file. If
                    provided, initial coordinates for MM sampling will be drawn
                    from this file instead of the coordinates file.
                    Default = None

    "conformersperset": int, specifies how many unique MM MD trajectories to 
                    run per dataset.
                    Default = 1

    "qmengine": str, specifies which qmengine to use to run QM calculations. 
                    Options are "chemcloud", "slurm", and "debug".
                    If "chemcloud", options are read from a valid TeraChem input
                    file (tc_template.in) and sent to chemcloud.
                    If "queue", a sbatch_template.sh file must be provided. This
                    file allows the QM TeraChem calculations to be run by a 
                    SLURM scheduler. 

    "resppriors": int, specifies an alternative to using RESP calculations
                    without incorporating the RESP fit into the objective
                    function. Instead, a prior is calculated for each atomic 
                    charge, and a penalty function is added to the objective.
                    If 0, don't use RESP priors for optimizing charges.
                    If 1, the width for each atom is determined by the 
                    distribution of RESP charges for that atom in the training 
                    set.
                    If 2, the width will be determined based on the RESP and ESP
                    charges for that atom in the training set.
                    Default = 0

    "stride": int, specifies stride in sampling initial training data from the
                    QM trajectory in coors.
                    Default = 50

    "tcout": str, specifies TeraChem output file with energies and gradients
                    from the QM trajectory in coors. 
                    Default = "tc.out"
    """
    )
    print(summary)
    print(keywords)


def createModel(inp):
    if inp.activelearning > 1:
        ffModel = active_learning.ActiveLearningModel(inp)
    else:
        ffModel = model.Model(inp)
    return ffModel


def getRestartCycle(inp):
    if inp.restart:
        restartCycle = ffModel.restartCycle
        print(f"Restarting optimization at cycle {restartCycle}")
    else:
        restartCycle = -1
    return restartCycle


if __name__ == "__main__":
    if len(sys.argv) == 1:
        print("No input file provided, printing summary message instead")
        printHelp()
        sys.exit(0)
    if sys.argv[1] == "-h" or sys.argv[1] == "--help":
        printHelp()
        sys.exit(0)

    inputFile = sys.argv[1]
    inp = inputs.Input.fromYaml(inputFile)
    ffModel = createModel(inp)
    restartCycle = getRestartCycle(inp)

    # First optimization cycle is not necessary if restarting from somewhere later
    if restartCycle < 0:
        ffModel.initialCycle()

    # Begin sampling/optimization cycling
    print(
        "%7s%15s%15s%22s%25s%8s%8s%8s"
        % (
            "Epoch",
            "Training",
            "Validation",
            "Validation change / %",
            "Validation, initial params",
            "MM time",
            "QM time",
            "FB time",
        )
    )
    for i in range(1, inp.maxcycles + 1):

        if i < restartCycle:
            continue

        mmStart = perf_counter()
        ffModel.doMMSampling(i)
        mmEnd = perf_counter()
        mmTime = mmEnd - mmStart

        ffModel.doQMCalculations(i)
        qmEnd = perf_counter()
        qmTime = qmEnd - mmEnd

        ffModel.doParameterOptimization(i)
        optResults = ffModel.optResults
        fbEnd = perf_counter()
        fbTime = fbEnd - qmEnd

        if inp.validinitial:
            print(
                "%7d%15.8f%15.8f%22.8f%25.8f%8.1f%8.1f%8.1f"
                % (
                    i,
                    optResults[0],
                    optResults[1],
                    optResults[2],
                    optResults[3],
                    mmTime,
                    qmTime,
                    fbTime,
                )
            )
        else:
            print(
                "%7d%15.8f%15.8f%22.8f%25.8f%8.1f%8.1f%8.1f"
                % (
                    i,
                    optResults[0],
                    optResults[1],
                    optResults[2],
                    "",
                    mmTime,
                    qmTime,
                    fbTime,
                )
            )

        if ffModel.converged:
            break
