#!/usr/bin/env python

from dataclasses import fields
from textwrap import dedent
from time import perf_counter
from typing import get_type_hints

import typer

from ff_optimizer import active_learning, inputs, model
from ff_optimizer import setup as st

app = typer.Typer()


def print_info(cls=inputs.Input):
    type_hints = get_type_hints(cls)
    for field in fields(cls):
        fname = field.name
        ftype = type_hints.get(fname, None)
        ftypestr = ftype.__name__ if ftype else "None"
        fdefault = field.default if field.default != field.default_factory else None
        fcomment = field.metadata.get("comment", "")
        print(f"{fname}: {ftypestr}, Default = {str(fdefault)}")
        print(fcomment.strip("\n"))


@app.command()
def print_manual():
    summary = dedent(
        """\
    ### Summary of necessary directories and files ###

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
    keywordIntro = dedent(
        """\
    ### Input files ###

    Keywords are to be provided in a .yaml (or .json) file. Example:
    restart: False
    mmengine: amber
    maxcycles: 30

    For all keywords, the syntax is simply keyword: value. If using a json, 
    take care that keywords have the the right type: strings should be in 
    quotes, ints should not. If using a yaml file, this doesn't matter.

     ### List of all keywords ###
     """
    )
    print(summary)
    print(keywordIntro)
    print_info()


def createModel(inp):
    if inp.activelearning > 1:
        ffModel = active_learning.ActiveLearningModel(inp)
    else:
        ffModel = model.Model(inp)
    return ffModel


def getRestartCycle(inp, ffModel):
    if inp.restart:
        restartCycle = ffModel.restartCycle
        print(f"Restarting optimization at cycle {restartCycle}")
    else:
        restartCycle = -1
    return restartCycle


@app.command()
def setup(xyz: str):
    st.setup(xyz)


@app.command()
def optimize(input_file: str):
    inp = inputs.Input.fromYaml(input_file)
    ffModel = createModel(inp)
    restartCycle = getRestartCycle(inp, ffModel)

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
                "%7d%15.8f%15.8f%22.8f%25.8s%8.1f%8.1f%8.1f"
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


def main():
    app()


if __name__ == "__main__":
    app()
