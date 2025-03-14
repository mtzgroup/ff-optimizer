import sys
from simtk.openmm import app
from simtk import unit
import simtk.openmm as mm
from parmed.openmm import NetCDFReporter

# The OpenMM python file must satisfy the following conditions:    
# 1. It accepts as arguments (in order) the .prmtop and .rst7 files
# 2. It performs all necessary equilibrations internally           
# 3. It produces a {name}.nc file containing the sampled geometries

# read prmtop and rst7 file names from command line
prmtopFile = sys.argv[1]
rst = sys.argv[2]

# Load Amber files
coordinates = app.AmberInpcrdFile(rst)
prmtop = app.AmberPrmtopFile(prmtopFile)

# Create OpenMM system
system = prmtop.createSystem(nonbondedMethod=app.NoCutoff, nonbondedCutoff=1*unit.nanometer, removeCMMotion=True, constraints=None, rigidWater=True)
integrator = mm.LangevinIntegrator(30*unit.kelvin, 1.0/unit.picosecond, 0.001*unit.picoseconds)
integrator.setConstraintTolerance(1e-05)

# Add spherical, harmonic restraining potential
force = mm.CustomExternalForce('k*max(0, r-0.563)^2; r=sqrt(x*x+y*y+z*z)')
force.addGlobalParameter("k", 10.0*unit.kilocalories_per_mole/unit.angstroms**2)
system.addForce(force)
for i in range(system.getNumParticles()):
    force.addParticle(i,[])

# Fix one molecule in the center of the system
system.setParticleMass(0, 1e10)

# Figure out whether we're running on CPU or GPU
try:
    platform = mm.Platform.getPlatformByName("CUDA")
    simulation = app.Simulation(prmtop.topology, system, integrator, platform)
except:
    platform = mm.Platform.getPlatformByName("CPU")
    simulation = app.Simulation(prmtop.topology, system, integrator, platform)

# Minimize system
simulation.context.setPositions(coordinates.positions)
simulation.minimizeEnergy()

# Equilibrate system
simulation.context.setVelocitiesToTemperature(30 * unit.kelvin)
for i in range(10):
    integrator.setTemperature(30 * (i+1) * unit.kelvin)
    simulation.step(10000)
    
# Run sampling dynamics    
name = rst.split(".")[0]
# Sample geometries every 10 ps
simulation.reporters.append(NetCDFReporter(f"{name}.nc",10000))
simulation.step(900000)
