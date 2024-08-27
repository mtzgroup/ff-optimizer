import sys
from simtk.openmm import app
from simtk import unit
import simtk.openmm as mm
from parmed.openmm import NetCDFReporter

prmtopFile = sys.argv[1]
rst = sys.argv[2]
coordinates = app.AmberInpcrdFile(rst)
prmtop = app.AmberPrmtopFile(prmtopFile)
system = prmtop.createSystem(nonbondedMethod=app.NoCutoff, nonbondedCutoff=1*unit.nanometer, removeCMMotion=True, constraints=None, rigidWater=True)

integrator = mm.LangevinIntegrator(30*unit.kelvin, 1.0/unit.picosecond, 0.001*unit.picoseconds)

integrator.setConstraintTolerance(1e-05)
force = mm.CustomExternalForce('k*max(0, r-0.563)^2; r=sqrt(x*x+y*y+z*z)')
force.addGlobalParameter("k", 10.0*unit.kilocalories_per_mole/unit.angstroms**2)
system.addForce(force)
for i in range(system.getNumParticles()):
    force.addParticle(i,[])

system.setParticleMass(0, 1e10)
try:
    platform = mm.Platform.getPlatformByName("CUDA")
    simulation = app.Simulation(prmtop.topology, system, integrator, platform)
except:
    platform = mm.Platform.getPlatformByName("CPU")
    simulation = app.Simulation(prmtop.topology, system, integrator, platform)


simulation.context.setPositions(coordinates.positions)
simulation.minimizeEnergy()
simulation.context.setVelocitiesToTemperature(30 * unit.kelvin)
for i in range(10):
    integrator.setTemperature(30 * (i+1) * unit.kelvin)
    simulation.step(10000)
    
name = rst.split(".")[0]
simulation.reporters.append(NetCDFReporter(f"{name}.nc",10000))
simulation.step(900000)
