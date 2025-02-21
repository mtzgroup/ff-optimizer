from pathlib import Path

from qcio import *
from qcparse import parse


def grep(filename, string):
    with open(filename, "r") as f:
        for line in f.readlines():
            if string in line:
                return line


# Makes a ProgramInput with all the keywords in a TeraChem input file
def inputFromTCin(inp, extras={}):
    keywords = {}
    charge = None
    spinmult = None
    calc = "energy"
    with open(inp, "r") as f:
        for line in f.readlines():
            split = line.split()
            if len(split) == 0:
                continue
            keyword = split[0].lower()
            if keyword.startswith("!") or keyword.startswith("#"):
                continue
            value = " ".join(split[1:])  # for the rare multitoken values
            if keyword == "method":
                method = value
            elif keyword == "basis":
                basis = value
            elif keyword == "charge":
                charge = value
            elif keyword == "spinmult":
                spinmult = value
            elif keyword == "coordinates":
                coords = value
            elif keyword == "run":
                calc = value
            else:
                keywords[keyword] = value
    if calc == "frequencies":
        calc = "hessian"
    if calc == "ts":
        calc = "transition_state"
    if calc == "minimize":
        calc = "optimization"
    model = {"method": method, "basis": basis}
    mol = Structure.open(coords, charge=charge, multiplicity=spinmult)
    inp = ProgramInput(
        structure=mol, calctype=calc, model=model, keywords=keywords, extras=extras
    )
    return inp


# Reads the BARE MINIMUM from a TeraChem output to make a valid ProgramInput
# Requires coordinates file to be where the output file says it is
def inputFromTCout(output, extras={}):
    coordinates = grep(output, "XYZ coordinates").split()[2]
    method = grep(output, "Method:").split()[1]
    basis = grep(output, "Using basis set:").split()[3]
    mol = Structure.open(coordinates)
    model = {"method": method, "basis": basis}
    if grep(output, "SINGLE POINT GRADIENT CALCULATIONS"):
        calc = "gradient"
    # man I hope they don't fix the typo
    elif grep(output, "RUNNING GEOMETRY OPTMIZATION"):
        calc = "optimization"
    elif grep(output, "FREQUENCY ANALYSIS"):
        calc = "hessian"
    # Neglect MD, etc. for now
    else:
        calc = "energy"
    inp = ProgramInput(structure=mol, calctype=calc, model=model, extras=extras)
    return inp


def tclinesFromInput(inp):
    name = inp.calctype.name
    if name == "hessian":
        calctype = "frequencies"
    elif name == "optimization":
        calctype = "minimize"
    elif name == "transition_state":
        calctype = "ts"
    else:
        calctype = name
    keys = ["run", "coordinates", "charge", "spinmult", "method", "basis"]
    vals = [
        calctype,
        "geometry.xyz",
        inp.structure.charge,
        inp.structure.multiplicity,
        inp.model.method,
        inp.model.basis,
    ]
    lines = ""
    for key, val in zip(keys, vals):
        lines += "%-20s %1s\n" % (key, val)
    for key in inp.keywords.keys():
        lines += ("%-20s %1s\n" % (key, inp.keywords[key])).rstrip(" ")
    return lines


def readLines(file):
    try:
        with open(file, "r") as f:
            lines = list(f.readlines())
        lines = "".join(lines)
    except:
        with open(file, "rb") as f:
            lines = list(f.readlines())
        lines = b"".join(lines)
    return lines


def collectFiles(scratch, inp):
    scratch = Path(scratch)
    files = {}
    if scratch.is_dir():
        name = scratch.name.split(".")[1]
        for scratchfile in sorted(scratch.iterdir()):
            lines = readLines(scratchfile)
            files["scr.geometry/" + scratchfile.name.replace(name, "geometry")] = lines
    tclines = tclinesFromInput(inp)
    coords = inp.structure.to_xyz()
    files["tc.in"] = tclines
    files["geometry.xyz"] = coords
    return files


def programOutputFromTCout(tcout, inp=None, getFiles=False):
    if not inp:
        inp = inputFromTCout(tcout)
    if grep(tcout, "Job finished"):
        success = True
        traceback = ""
    else:
        success = False
        traceback = "oh no!"
    try:
        result = parse(tcout, "terachem")
        resultDict = result.model_dump(mode="json", exclude_unset=True)
        if getFiles:
            scr = Path(grep(tcout, "Scratch directory:").split()[2])
            files = collectFiles(scr, inp)
        else:
            files = {}
        resultDict["files"] = files
        result = SinglePointResults(**resultDict)
    except:
        result = Files(files={})
    with open(tcout, "r") as f:
        lines = list(f.readlines())
    stdout = "".join(lines)
    extras = {}
    provenance = Provenance(program="terachem")
    output = ProgramOutput(
        input_data=inp,
        success=success,
        results=result,
        provenance=provenance,
        extras=extras,
        stdout=stdout,
        traceback=traceback,
    )
    return output


class MockFutureResult:
    def __init__(self, outputs):
        self.outputs = outputs

    def get(self):
        return self.outputs

    @classmethod
    def fromOutputFiles(cls, files):
        if isinstance(files, list):
            outputs = [programOutputFromTCout(f) for f in files]
        else:
            outputs = programOutputFromTCout(f)
        return cls(outputs)


def patcher(monkeypatch, target, outputs):

    def mockCompute(
        program, inputs, return_future=False, collect_files=False, outputs=outputs
    ):
        if not isinstance(outputs, list):
            outputs = [outputs]
        if not isinstance(inputs, list):
            inputs = [inputs]
        if isinstance(outputs[0], ProgramOutput):
            programOutputs = outputs
        else:
            if collect_files:
                programOutputs = [
                    programOutputFromTCout(out, inp=inp, getFiles=True)
                    for inp, out in zip(inputs, outputs)
                ]
            else:
                programOutputs = [
                    programOutputFromTCout(out, inp=inp)
                    for inp, out in zip(inputs, outputs)
                ]
        if len(programOutputs) == 1:
            programOutputs = programOutputs[0]
        if return_future:
            result = MockFutureResult(programOutputs)
        else:
            result = programOutputs
        return result

    monkeypatch.setattr(target, "compute", mockCompute)
