import sys

import simtk.openmm as mm
from parmed.openmm import NetCDFReporter
from simtk import unit
from simtk.openmm import app

prmtopFile = sys.argv[1]
rst = sys.argv[2]
name = rst.split(".")[0]

coordinates = app.AmberInpcrdFile(rst)
prmtop = app.AmberPrmtopFile(prmtopFile)
system = prmtop.createSystem(
    nonbondedMethod=app.NoCutoff,
    nonbondedCutoff=2 * unit.nanometer,
    removeCMMotion=True,
    constraints=None,
    rigidWater=True,
)
integrator = mm.VerletIntegrator(0.0002 * unit.picoseconds)
integrator.setConstraintTolerance(1e-05)

force = mm.CustomExternalForce("k*max(0, r-0.563)^2; r=sqrt(x*x+y*y+z*z)")
force.addGlobalParameter("k", 10.0 * unit.kilocalories_per_mole / unit.angstroms**2)
system.addForce(force)
for i in range(system.getNumParticles()):
    force.addParticle(i, [])
system.setParticleMass(0, 0)

platform = mm.Platform.getPlatformByName("CPU")
simulation = app.Simulation(prmtop.topology, system, integrator, platform)
simulation.context.setPositions(coordinates.positions)
simulation.minimizeEnergy(maxIterations=10)

simulation.reporters.append(NetCDFReporter(f"{name}.nc", 40))
simulation.step(320)
