#!/usr/bin/env python

import argparse
import errno
from time import perf_counter
from ff_optimizer import model
from textwrap import dedent

# Some helper functions

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

if __name__ == "main":
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
    
    ffModel = model.Model(args)
    
    # First optimization cycle is not necessary if restarting from somewhere later
    if restartCycle < 0:
        ffModel.initialCycle() 
    
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
        ffModel.doMMSampling(i)
        mmEnd = perf_counter()
        mmTime = mmEnd - mmStart
    
        ffModel.doQMCalculations(i)
        qmEnd = perf_counter()
        qmTime = qmEnd - mmEnd
    
        optResults = ffModel.doParameterOptimization(i)
        fbEnd = perf_counter()
        fbTime = fbEnd - qmEnd
    
        if i == 1:
            print(
                "%7d%15.8f%15.8f%20.8f%23s%8.1f%8.1f%8.1f"
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
        else:
            print(
                "%7d%15.8f%15.8f%20.8f%23.8f%8.1f%8.1f%8.1f"
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
