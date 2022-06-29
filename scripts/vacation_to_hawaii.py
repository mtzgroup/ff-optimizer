#!/usr/bin/env python

import argparse
import errno
import os
from shutil import copyfile, rmtree
from time import perf_counter

from ff_optimizer import mmengine, optengine, qmengine

from textwrap import dedent

# Some helper functions
def die():
    raise RuntimeError("die here")


def convertTCtoFB(
    tcout, coors, stride, start=None, end=None, qdata="qdata.txt", mdcrd="all.mdcrd"
):

    molSize = 0
    coords = []
    frame = []
    energies = []
    grads = []
    coordIndex = []
    gradIndex = []
    lineCounter = 0
    frameStart = -1
    coorEnergies = []

    with open(coors, "r") as f:
        maxFrame = -1
        for line in f.readlines():
            splitLine = line.split()
            if lineCounter == 0:
                molSize = int(splitLine[0])
            if lineCounter == 1:
                index = int(splitLine[2]) + 1  # coor files are 0-indexed
                energy = splitLine[0]
                if frameStart == -1:
                    frameStart = index
            elif lineCounter > 1 and lineCounter < molSize + 2:
                for token in splitLine[1:]:
                    frame.append(token)
            if lineCounter == molSize + 1:
                if index > maxFrame:
                    maxFrame = index
                    coords.append(frame)
                    coordIndex.append(index)
                    coorEnergies.append(energy)
                else:
                    j = len(coords) - 1
                    while coordIndex[j] > index:
                        j -= 1
                    coords[j] = frame
                    coordIndex[j] = index
                    coorEnergies[j] = energy
                frame = []
                lineCounter = -1
            lineCounter = lineCounter + 1

    gradCounter = molSize + 37
    index = 0
    frameStart = -1

    with open(tcout, "r") as f:
        maxFrame = -1
        for line in f.readlines():
            if "MD STEP" in line:
                index = int(line.split()[4])
                if frameStart == -1:
                    if len(grads) > 0:
                        frameStart = index - 1  # first step is unnumbered "step 0"
                        gradIndex[0] = index - 1
                    else:
                        frameStart = index
            if "FINAL ENERGY" in line:
                energy = float(line.split()[2])
            if "Gradient units" in line:
                gradCounter = 0
                frame = []
            if gradCounter < molSize + 2 and gradCounter > 2:
                for token in line.split():
                    frame.append(token)
            if gradCounter == molSize + 2:
                for token in line.split():
                    frame.append(token)
                if index > maxFrame:
                    maxFrame = index
                    grads.append(frame)
                    gradIndex.append(index)
                    energies.append(energy)
                else:
                    j = len(grads) - 1
                    while gradIndex[j] > index:
                        j -= 1
                    grads[j] = frame
                    gradIndex[j] = index
                    energies[j] = energy
            gradCounter = gradCounter + 1

    if start == None:
        start = gradIndex[0]
    if end == None:
        end = gradIndex[-1]

    jobCounter = 0
    indices = []
    precision = len(coorEnergies[0].split(".")[1])
    for i in range(len(gradIndex)):
        for j in range(len(coordIndex)):
            if gradIndex[i] == coordIndex[j]:
                eFormat = "%." + str(precision) + "f"
                gradEnergy = eFormat % energies[i]
                if coorEnergies[j] != gradEnergy:
                    # print("Mismatched energies in step " + str(gradIndex[i]))
                    # raise RuntimeError("Mismatched energies from " + tcout + " and " + coors)
                    break
                indices.append([i, j])

    usedIndices = []
    lastFrame = -args.stride - 37
    for i in range(len(indices)):
        if gradIndex[indices[i][1]] >= start and gradIndex[indices[i][1]] <= end:
            if gradIndex[indices[i][1]] - lastFrame >= stride:
                lastFrame = gradIndex[indices[i][1]]
                usedIndices.append(indices[i])

    with open(qdata, "w") as f:
        for j in usedIndices:
            f.write("\n")
            f.write("JOB " + str(jobCounter) + "\n")
            coordLine = "COORDS "
            for coord in coords[j[1]]:
                coordLine = coordLine + str(coord) + " "
            gradLine = "FORCES "
            for grad in grads[j[0]]:
                gradLine = gradLine + str(grad) + " "
            f.write(coordLine + "\n")
            f.write("ENERGY " + str(energies[j[0]]) + "\n")
            f.write(gradLine + "\n")
            jobCounter = jobCounter + 1

    with open(mdcrd, "w") as f:
        f.write("converted from TC with convertTCMDtoFB.py\n")
        for j in usedIndices:
            tokenCounter = 1
            for coord in coords[j[1]]:
                f.write("%8.3f" % float(coord))
                if tokenCounter == 10:
                    f.write("\n")
                    tokenCounter = 1
                else:
                    tokenCounter = tokenCounter + 1
            if tokenCounter != 1:
                f.write("\n")
    return jobCounter

def checkArgs(args):
    # Check for necessary folders and files
    # TODO: Make these the proper errors, and just put it in a function
    if not os.path.isdir(args.dynamicsdir):
        raise RuntimeError("Dynamics directory " + args.dynamicsdir + " does not exist")
    if not os.path.isfile(os.path.join(args.dynamicsdir, args.coors)):
        raise RuntimeError(
            "XYZ coordinates (from QM dynamics) "
            + args.coors
            + " does not exist in "
            + args.dynamicsdir
        )
    if not os.path.isfile(os.path.join(args.dynamicsdir, args.tcout)):
        raise RuntimeError(
            "TC output file (from QM dynamics "
            + args.tcout
            + " does not exist in "
            + args.dynamicsdir
        )
    if not os.path.isdir(args.optdir):
        raise RuntimeError(
            "ForceBalance optimization directory " + args.optdir + " does not exist"
        )
    if not os.path.isfile(os.path.join(args.optdir, "conf.pdb")):
        raise RuntimeError(
            "Prototype PDB coordinates conf.pdb does not exist in " + args.optdir
        )
    if not os.path.isfile(os.path.join(args.optdir, "setup.leap")):
        raise RuntimeError("Tleap input file setup.leap does not exist in " + args.optdir)
    if not os.path.isfile(os.path.join(args.optdir, args.opt0)):
        raise RuntimeError(
            "Initial ForceBalance optimization input file "
            + args.opt0
            + " does not exist in "
            + args.optdir
        )
    if not os.path.isfile(os.path.join(args.optdir, args.valid0)):
        raise RuntimeError(
            "Initial ForceBalance validation input file "
            + args.valid0
            + " does not exist in "
            + args.optdir
        )
    if not os.path.isdir(args.sampledir):
        raise RuntimeError("MM sampling directory " + args.sampledir + " does not exist")
    if not os.path.isfile(os.path.join(args.sampledir, "heat1.in")):
        raise RuntimeError(
            "No sander input file for equilibration named heat1.in provided in "
            + args.sampledir
        )
    if not os.path.isfile(os.path.join(args.sampledir, "md.in")):
        raise RuntimeError(
            "No sander input file for sampling named md.in provided in " + args.sampledir
        )
    if not os.path.isfile(os.path.join(args.sampledir, "cpptraj.in")):
        raise RuntimeError("No cpptraj input file provided in " + args.sampledir)
    
    if args.qmengine != "queue" and args.qmengine != "debug" and args.qmengine != "tccloud":
        raise RuntimeError("Engine " + args.qmengine + " is not implemented")
    if args.qmengine == "queue":
        if not os.path.isfile(os.path.join(args.sampledir, args.sbatch)):
            raise RuntimeError(
                "Sbatch template " + args.sbatch + " does not exist in " + args.sampledir
            )
        if not os.path.isfile(os.path.join(args.sampledir, args.tctemplate)):
            raise RuntimeError(
                "TC input template "
                + args.tctemplate
                + " does not exist in "
                + args.sampledir
            )
        if not os.path.isfile(os.path.join(args.sampledir, args.tctemplate_long)):
            raise RuntimeError(
                "TC input template (for resubmitting jobs) "
                + args.tctemplate_long
                + " does not exist in "
                + args.sampledir
            )
    
        if (
            args.qmengine != "queue"
            and args.qmengine != "debug"
            and args.qmengine != "tccloud"
        ):
            raise RuntimeError("Engine " + args.qmengine + " is not implemented")
        if args.qmengine == "queue":
            if not os.path.isfile(os.path.join(args.sampledir, args.sbatch)):
                raise RuntimeError(
                    "Sbatch template "
                    + args.sbatch
                    + " does not exist in "
                    + args.sampledir
                )
            if not os.path.isfile(os.path.join(args.sampledir, args.tctemplate)):
                raise RuntimeError(
                    "TC input template "
                    + args.tctemplate
                    + " does not exist in "
                    + args.sampledir
                )
            if not os.path.isfile(os.path.join(args.sampledir, args.tctemplate_long)):
                raise RuntimeError(
                    "TC input template (for resubmitting jobs) "
                    + args.tctemplate_long
                    + " does not exist in "
                    + args.sampledir
                )
    
        if args.mmengine != "amber":
            raise ValueError(f"MM Engine {args.mmengine} is unsupported!")
        if args.nvalids < 1:
            raise ValueError(f"Must use at least one validation set for now")
        if not os.path.isfile(os.path.join(args.sampledir, args.trainMdin)):
            raise FileNotFoundError(
                errno.ENOENT,
                os.srterror(errno.ENOENT),
                os.path.join(args.sampledir, args.trainMdin),
            )
        if not os.path.isfile(os.path.join(args.sampledir, args.validMdin)):
            raise FileNotFoundError(
                errno.ENOENT,
                os.srterror(errno.ENOENT),
                os.path.join(args.sampledir, args.validMdin),
            )
        if args.conformers is not None:
            if not os.path.isfile(os.path.join(args.dynamicsdir,args.conformers)):
                raise RuntimeError(f"Conformer XYZ file {args.conformers} does not exist in {args.dynamicsdir}")
        if args.conformersPerSet < 1:
            raise ValueError("conformersPerSet must be a positive integer")

def initializeOptEngine(args):
    optOptions = {}
    optOptions["optdir"] = args.optdir
    optOptions["sampledir"] = args.sampledir
    optOptions["respPriors"] = args.respPriors
    optOptions["resp"] = args.resp
    optOptions["maxCycles"] = args.maxcycles
    optOptions["restart"] = args.restart
    optEngine = optengine.OptEngine(optOptions)
    restartCycle = optEngine.restartCycle
    return optEngine
    
def initializeQMEngine(args):
    if args.qmengine == "debug":
        qmEngine = qmengine.DebugEngine(
            os.path.join(args.sampledir, args.tctemplate),
            os.path.join(args.sampledir, args.tctemplate_long),
            doResp=doResp,
        )
    elif args.qmengine == "queue":
        qmEngine = qmengine.SbatchEngine(
            os.path.join(args.sampledir, args.tctemplate),
            os.path.join(args.sampledir, args.tctemplate_long),
            os.path.join(args.sampledir, args.sbatch),
            os.getenv("USER"),
            doResp=doResp,
        )
    elif args.qmengine == "tccloud":
        qmEngine = qmengine.TCCloudEngine(
            os.path.join(args.sampledir, args.tctemplate),
            os.path.join(args.sampledir, args.tctemplate_long),
            doResp=doResp,
        )
    return qmEngine

def initializeMMEngine(args):
    mmOptions = {}
    if args.conformers is None:
        mmOptions["start"] = args.start
        mmOptions["end"] = args.end
        mmOptions["split"] = args.split
        mmOptions["coordPath"] = os.path.join(args.dynamicsdir, args.coors)
    else:
        mmOptions["start"] = None
        mmOptions["end"] = None
        mmOptions["split"] = None
        mmOptions["coordPath"] = os.path.join(args.dynamicsdir, args.conformers)
    mmOptions["conformers"] = args.conformersPerSet
    mmOptions["nvalids"] = args.nvalids
    mmOptions["trainMdin"] = args.trainMdin
    mmOptions["validMdin"] = args.validMdin
    mmOptions["leap"] = "setup.leap"
    mmOptions["heatCounter"] = heatCounter
    if args.mmengine == "amber":
        mmEngine = mmengine.ExternalAmberEngine(mmOptions)
    return mmEngine

if __name__ == "__main__":
    # Summary stuff
    summary = dedent(
        """\
    Summary of necessary directories and files
    Dynamics directory (--dynamicsdir, 0_dynamics) contains:
        -- QM dynamics coordinates (--coors, coors.xyz)
        -- TeraChem output file (--tcout, tc.out)
        -- Conformers file (--conformers, coors.xyz)
    FB Optimization directory (--optdir, 1_opt) contains:
        -- FB optimization input file (--opt0, opt_0.in)
        -- FB validation input file (--valid0, valid_0.in)
        -- PDB conformation file (conf.pdb)
        -- Tleap setup file (setup.leap)
        -- Parameter frcmod file (*.frcmod)
        -- Parameter mol2 file (*.mol2)
    MM sampling directory (--sampledir, 2_sampling) contains:
        -- Amber equilibration inputs (heat*in)
        -- Amber MM sampling input (md.in)
        -- Cpptraj input file (cpptraj.in)
        -- TeraChem input file for fire (tc_template.in)
        -- Backup TeraChem input file if job fails (tc_template_long.in)
        -- sbatch template file for queue (sbatch_template.sh)
    """

    )
    parser = argparse.ArgumentParser(
        epilog=summary, formatter_class=argparse.RawDescriptionHelpFormatter
    )
    parser.add_argument(
        "--dynamicsdir",
        help="Folder containing TC output file from MD trajectory and XYZ coord file from MD trajectory, default is 0_dynamics",
        type=str,
        default="0_dynamics",
    )
    parser.add_argument(
        "--tcout",
        help="TeraChem output file from MD trajectory, default is tc.out",
        type=str,
        default="tc.out",
    )
    parser.add_argument(
        "--coors",
        help="XYZ file containing coords from TeraChem MD trajectory, default is coors.xyz",
        type=str,
        default="coors.xyz",
    )
    parser.add_argument(
        "--start",
        help="Lower bound for sampling of frames from MD trajectory, default is the first frame",
        type=int,
        default=None,
    )
    parser.add_argument(
        "--end",
        help="Upper bound for sampling of frames from MD trajectory, default is the last frame",
        type=int,
        default=None,
    )
    parser.add_argument(
        "--stride",
        help="Stride between successive samples of frames from MD trajectory, default is 50",
        type=int,
        default=50,
    )
    parser.add_argument(
        "--optdir",
        help="Directory where ForceBalance optimization is performed, default is 1_opt",
        type=str,
        default="1_opt",
    )
    parser.add_argument(
        "--maxcycles",
        help="Maximum number of optimization/sampling cycles to be performed, default is 10",
        type=int,
        default=10,
    )
    parser.add_argument(
        "--sampledir",
        help="Directory where MM sampling is performed, default is 2_sampling",
        type=str,
        default="2_sampling",
    )
    parser.add_argument(
        "--split",
        help="If set, sample one geometry from before this frame and one from after",
        type=int,
        default=None,
    )
    parser.add_argument(
        "--opt0",
        help="ForceBalance input file for first optimization, with a single target named dynamics, default is opt_0.in",
        type=str,
        default="opt_0.in",
    )
    parser.add_argument(
        "--valid0",
        help="Template ($options section only) for ForceBalance input files for validation single points, default is valid_0.in",
        type=str,
        default="valid_0.in",
    )
    parser.add_argument(
        "--qmengine",
        help="Engine for performing QM calculations, either queue, debug, or tccloud",
        type=str.lower,
        default="tccloud",
    )
    parser.add_argument(
        "--sbatch",
        help="Sbatch template file for submitting jobs to fire, default is sbatch_template.sh",
        type=str,
        default="sbatch_template.sh",
    )
    parser.add_argument(
        "--tctemplate",
        help="TC input template file for submitting jobs to fire, default is tc_template.in",
        type=str,
        default="tc_template.in",
    )
    parser.add_argument(
        "--tctemplate_long",
        help="TC input template file for resubmitting failed jobs to fire, default is tc_template_long.in",
        type=str,
        default="tc_template_long.in",
    )
    parser.add_argument(
        "--restart",
        help="Restart partially complete forcefield optimization",
        action="store_true",
    )
    parser.add_argument(
        "--resp",
        help="Weight of RESP in the objective function. If 0, don't do RESP. Default is 0.",
        type=float,
        default=0,
    )
    
    parser.add_argument(
        "--mmengine",
        help="Package for running MM sampling. Default is amber.",
        type=str.lower,
        default="amber",
    )
    parser.add_argument(
        "--nvalids",
        help="Number of validation sets to create. Default is 1.",
        type=int,
        default=1,
    )
    parser.add_argument(
        "--trainMdin",
        help="MD input file for MM sampling for training set. Default is md.in",
        type=str,
        default="md.in",
    )
    parser.add_argument(
        "--validMdin",
        help="MD input file for MM sampling for validation set. Default is md.in",
        type=str,
        default="md.in",
    )
    parser.add_argument(
        "--respPriors",
        help="Use RESP to compute priors for the charges. Mode 1: RESP st. dev. Mode 2: Uses ESP-RESP difference. Default is not.",
        type=int,
        default=0,
    )
    parser.add_argument("--conformers",help="XYZ file containing conformers to use for initial MM sampling conditions. If not provided, will use the initial dynamics XYZ.",type=str,default=None)
    parser.add_argument("--conformersPerSet",help="Number of conformers to do MM sampling on per training/validation dataset, default is 1",type=int, default=1)
    
    args = parser.parse_args()
    checkArgs(args)
    
    # Set some miscellaneous variables
    home = os.getcwd()
    os.rename(
        os.path.join(args.optdir, args.valid0), os.path.join(args.optdir, "valid_0.in")
    )
    os.rename(os.path.join(args.optdir, args.opt0), os.path.join(args.optdir, "opt_0.in"))
    mdFiles = []
    heatCounter = 0
    for f in os.listdir(args.sampledir):
        if os.path.isfile(os.path.join(args.sampledir, f)):
            mdFiles.append(f)
            if f.startswith("heat"):
                heatCounter += 1
    if args.resp != 0 or args.respPriors != 0:
        doResp = True
    else:
        doResp = False

    optEngine = initializeOptEngine(args)
    qmEngine = initializeQMEngine(args)
    mmEngine = initializeMMEngine(args)
    
    # First optimization cycle is not necessary if restarting from somewhere later
    if restartCycle < 0:
    
        # Create initial target data from dynamics
        with open(os.path.join(args.optdir, args.opt0)) as f:
            for line in f.readlines():
                splitLine = line.split()
                if len(splitLine) > 1:
                    if splitLine[0] == "name":
                        initialTarget = splitLine[1]
        if not os.path.isdir(os.path.join(args.optdir, "targets")):
            os.mkdir(os.path.join(args.optdir, "targets"))
        path = os.path.join(args.optdir, "targets", initialTarget)
        if not os.path.isdir(path):
            os.mkdir(path)
        l = convertTCtoFB(
            os.path.join(args.dynamicsdir, args.tcout),
            os.path.join(args.dynamicsdir, args.coors),
            args.stride,
            args.start,
            args.end,
            os.path.join(path, "qdata.txt"),
            os.path.join(path, "all.mdcrd"),
        )
        for f in ["setup.leap", "conf.pdb", "setup_valid_initial.leap"]:
            copyfile(os.path.join(args.optdir, f), os.path.join(path, f))
    
        os.chdir(args.optdir)
        optEngine.optimizeForcefield(0)
        os.chdir(home)
    
    # Begin sampling/optimization cycling
    print(
        "%7s%15s%15s%20s%23s%8s%8s%8s"
        % (
            "Epoch",
            "Validation",
            "Valid ratio",
            "Current-Previous",
            "Current-last Current",
            "MM time",
            "QM time",
            "FB time",
        )
    )
    for i in range(1, args.maxcycles + 1):
    
        if i <= restartCycle:
            continue
    
        mmStart = perf_counter()
        # Make sampling directory and copy files into it
        sampleName = str(i) + "_cycle_" + str(i)
        samplePath = os.path.join(args.sampledir, sampleName)
        if not os.path.isdir(samplePath):
            os.mkdir(samplePath)
        elif restartCycle == -1:
            rmtree(samplePath)
            os.mkdir(samplePath)
        src = os.path.join(args.optdir, "result", "opt_" + str(i - 1), "*")
        dest = os.path.join(samplePath, ".")
        os.system(f"cp {src} {dest}")
        src = os.path.join(args.optdir, "conf.pdb")
        os.system(f"cp {src} {dest}")
        src = os.path.join(args.optdir, "setup.leap")
        os.system(f"cp {src} {dest}")
        for f in mdFiles:
            src = os.path.join(args.sampledir, f)
            os.system(f"cp {src} {dest}")
        os.chdir(samplePath)
        # Do MM sampling
        if i == restartCycle + 1:
            mmEngine.restart()
        else:
            mmEngine.getMMSamples()
        mmEnd = perf_counter()
        mmTime = mmEnd - mmStart
    
        # Run QM calculations for each sampling trajectory
        for f in os.listdir():
            if (f.startswith("train") or f.startswith("valid")) and os.path.isdir(f):
                os.chdir(f)
                if i == restartCycle + 1:
                    qmEngine.restart(".")
                else:
                    pdbs = []
                    for g in os.listdir():
                        if g.endswith(".pdb"):
                            pdbs.append(g)
                    qmEngine.getQMRefData(pdbs, ".")
                os.chdir("..")
        os.chdir(home)
        qmEnd = perf_counter()
        qmTime = qmEnd - mmEnd
    
        # Set up new ForceBalance optimization
    
        # Copy new QM data into appropriate folders
        src = os.path.join(args.optdir, "targets", initialTarget, "*")
        trainFolder = os.path.join(args.optdir, "targets", "train_" + str(i))
        validFolder = os.path.join(args.optdir, "targets", "valid_" + str(i))
        if not os.path.isdir(trainFolder):
            os.mkdir(trainFolder)
        if not os.path.isdir(validFolder):
            os.mkdir(validFolder)
        dest = os.path.join(args.optdir, "targets", "train_" + str(i), ".")
        os.system(f"cp {src} {dest}")
        dest = os.path.join(args.optdir, "targets", "valid_" + str(i), ".")
        os.system(f"cp {src} {dest}")
    
    
        # NOTE: currently only supports a single validation set
        valids = []
        for f in os.listdir(os.path.join(args.sampledir, f"{str(i)}_cycle_{str(i)}")):
            if f.startswith("train_") and os.path.isdir(
                os.path.join(args.sampledir, f"{str(i)}_cycle_{str(i)}", f)
            ):
                trainFolder = f
            elif f.startswith("valid_") and os.path.isdir(
                os.path.join(args.sampledir, f"{str(i)}_cycle_{str(i)}", f)
            ):
                valids.append(f)
    
        src = os.path.join(
            args.sampledir, str(i) + "_cycle_" + str(i), trainFolder, "all.mdcrd"
        )
        dest = os.path.join(args.optdir, "targets", "train_" + str(i), ".")
        os.system(f"cp {src} {dest}")
        src = os.path.join(
            args.sampledir, str(i) + "_cycle_" + str(i), trainFolder, "qdata.txt"
        )
        os.system(f"cp {src} {dest}")
        src = os.path.join(
            args.sampledir,
            str(i) + "_cycle_" + str(i),
            valids[0],
            "all.mdcrd",
        )
        dest = os.path.join(args.optdir, "targets", "valid_" + str(i), ".")
        os.system(f"cp {src} {dest}")
        src = os.path.join(
            args.sampledir,
            str(i) + "_cycle_" + str(i),
            valids[0],
            "qdata.txt",
        )
        os.system(f"cp {src} {dest}")
    
        # Run ForceBalance on each input
        os.chdir(args.optdir)
        optEngine.optimizeForcefield(i)
        os.chdir(home)
        fbEnd = perf_counter()
        fbTime = fbEnd - qmEnd
    
        if i == 1:
            print(
                "%7d%15.8f%15.8f%20.8f%23s%8.1f%8.1f%8.1f"
                % (
                    i,
                    optEngine.valid[-1],
                    optEngine.valid[-1] / optEngine.validInitial[-1],
                    optEngine.valid[-1] - optEngine.validPrevious[-1],
                    "",
                    mmTime,
                    qmTime,
                    fbTime,
                )
            )
        else:
            print(
                "%7d%15.8f%15.8f%20.8f%23.8f%8.1f%8.1f%8.1f"
                % (
                    i,
                    optEngine.valid[-1],
                    optEngine.valid[-1] / optEngine.validInitial[-1],
                    optEngine.valid[-1] - optEngine.validPrevious[-1],
                    optEngine.valid[-1] - optEngine.valid[-2],
                    mmTime,
                    qmTime,
                    fbTime,
                )
            )
