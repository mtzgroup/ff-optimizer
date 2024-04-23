from chemcloud import CCClient
from qcio import Molecule, ProgramInput, SinglePointOutput

water = Molecule(
    symbols=["O", "H", "H"],
    geometry=[
        [0.0000, 0.00000, 0.0000],
        [0.2774, 0.89290, 0.2544],
        [0.6067, -0.23830, -0.7169],
    ],
)

client = CCClient()

prog_inp = ProgramInput(
    molecule=water,
    model={"method": "b3lyp", "basis": "6-31g"},
    calctype="energy",  # Or "gradient" or "hessian"
    keywords={},
)
future_result = client.compute("terachem", prog_inp, collect_files=True)
output: SinglePointOutput = future_result.get()
# SinglePointOutput object containing all returned data
print(output.stdout)
# The energy value requested
print(output.return_result)
