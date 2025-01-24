#  ff optimizer

Uses ForceBalance, Amber, and TeraChem to optimize forcefield parameters.

### Installation ###

1. Clone the ff-opt repo: `git clone git@github.com:mtzgroup/ff-optimizer.git`
2. With your python environment manager of choice, create an environment with a python version compatible with ff-opt (currently ^3.11)
3. Install poetry into this environment: `pip install poetry`
4. Ensure that the poetry executable is on your PATH: `which poetry`
5. Move into the ff-opt folder: `cd ff-optimizer`
6. Run `poetry install`

   ##### Installing Amber #####
   ff-opt depends on ForceBalance, which uses Amber's pysander interface. For ff-opt to work, Amber's pysander interface
   must be compiled with the same version of python used by ff-opt. As of right now, the python environment bundled with Amber
   is python3.8, which is deprecated. Follow the instructions below to install Amber with an up-to-date version of python.

   1. Create a python environment with the desired python version (currently ^3.11) and the following packages
       - numpy
       - scipy
       - cython
       - setuptools
   2. Activate this environment before installing Amber.
   3. In Amber's `run_cmake` file, set `DOWNLOAD_MINICONDA` to `FALSE`.
   4. Compile Amber.

   This will compile Amber's pysander interface with a newer version of python.

### Running ff-opt ###

1. Head to the directory where ff-opt was installed.
2. Create a new shell with ff-opt on PATH: `poetry shell`
3. Add Amber to the PATH with, e.g., `module load Amber`
4. Move to the desired folder and run ff-opt!

ff-opt currently has three commands:
- `ff-opt optimize [input file]`:
Optimizes force field parameters according to [input file]. [Here's](#summary-of-necessary-directories-and-files) a summary of the
directory structure and files required by the parameter optimizer. The input file specification can be found [here](#input-files), and
a full list of keywords can be found [here](#list-of-all-keywords). Examples are in `/ff-optimizer/examples`.
- `ff-opt setup [xyz file]`:
Sets up a force field optimization for [xyz file]. Creates all necessary files and folders required by `ff-opt optimize` and populates
them with reasonable defaults. Currently, it obtains atom types and starting parameters from GAFF, charges from AM1-BCC, and sets up an
optimization of all bonded parameters using wB97X-D3/def2-svp.
- `ff-opt print-manual`:
Prints this file, then quits.


### Summary of necessary directories and files ###

- Dynamics directory ("dynamicsdir", 0_dynamics) contains:
    - QM dynamics coordinates ("coors", coors.xyz)
    - TeraChem output file ("tcout", tc.out)
    - Conformers file ("conformers", coors.xyz)
- FB Optimization directory ("optdir", 1_opt) contains:
    - FB optimization input file (opt_0.in)
    - FB validation input file (valid_0.in)
    - PDB conformation file (conf.pdb)
    - Tleap setup file (setup.leap)
    - Parameter frcmod file (*.frcmod)
    - Parameter mol2 file (*.mol2)
- MM sampling directory ("sampledir", 2_sampling) contains:
    - Equilibration inputs (heat*in)
    - MM sampling input(s) ("trainmdin", "validmdin", md.in)
    - TeraChem input file for fire (tc_template.in)
    - Backup TeraChem input file if job fails (tc_template_long.in)
    - sbatch template file for queue (sbatch_template.sh)

### Input files ###

Keywords are to be provided in a .yaml (or .json) file. Example yaml format:
```
restart: False
mmengine: amber
maxcycles: 30
```

For all keywords, the syntax is simply keyword: value. If using a json, 
take care that keywords have the the right type: strings should be in 
quotes, ints should not. If using a yaml file, this doesn't matter.

 ### List of all keywords ###

dynamicsdir: str, Default = 0_dynamics
        Specifies the directory containing QM energies,
        gradients, and coordinates from an MD trajectory.
        
optdir: str, Default = 1_opt
        Specifies the directory where ForceBalance optimization will
        be run; contains FB inputs and initial parameters.
        
sampledir: str, Default = 2_sampling
        Specifies the directory where MM sampling and QM 
        calculations will be run; contains MM and QM inputs.
        
coors: str, Default = coors.xyz
        Specifies the name of the coordinates file in dynamicsdir.
        
start: int, Default = None
        Specifies the first frame of the coors file to sample from
        to create the initial dataset. 0-indexed. If not provided, 
        default is the first frame.
        
end: int, Default = None
        Specifies the last frame of the coors file to sample from. 
        0-indexed. If not provided, defaults to last frame in coors file.
        
split: int, Default = None
        Divides the coors file in two at the provided frame number.
        If provided, initial coordinates for MM MD for the training
        set will come from before split, and ICs for the validation
        set will come from after. 0-indexed.
        
resp: float, Default = 0
        Specifies the weight of RESP fits for the charges in the
        objective function. If nonzero, QM calculations will
        include RESP fits, and parameter optimization will attempt
        to match the ESP in those QM calculations. Can severely
        slow down force field fitting.
        
nvalids: int, Default = 1
        Specifies how many validation sets to produce and evaluate.
        
restart: bool, Default = False
        Specifies whether or not to restart the force field
        optimization. If True, it will automatically determine
        where to restart.
        
activelearning: int, Default = 1
        Specifies how many different models to train. If
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
        
maxcycles: int, Default = 30
        Specifies how many cycles to run before stopping.
        
mmengine: str, Default = amber
         specifies which software will run MM MD. Options are
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
        An example can be found in examples/openmm/2_sampling.
        
trainmdin: str, Default = md.in
        Specifies the input file for MM MD. This input will be
        used to generate the training set.
        
validmdin: str, Default = md.in
        Specifies the input file for MM MD. This input will be
        used to generate validation sets.
        
conformers: str, Default = None
        Specifies the name of a conformers coordinates file in 
        dynamicsdir. If provided, initial coordinates for MM sampling 
        will be drawn from this file instead of the coordinates file.
        
conformersperset: int, Default = 1
        Specifies how many unique MM MD trajectories to run per dataset.
        
qmengine: str, Default = chemcloud
        Specifies which qmengine to use to run QM calculations.
        Options are "chemcloud", "slurm", and "debug".
        If "chemcloud", options are read from a valid TeraChem input
        file (tc_template.in) and sent to chemcloud.
        If "queue", a sbatch_template.sh file must be provided. This
        file allows the QM TeraChem calculations to be run by a
        SLURM scheduler.                                            
        If "debug", TeraChem will be run locally, one job at a time.
        
resppriors: int, Default = 0
        Specifies an alternative to using RESP calculations
        without incorporating the RESP fit into the objective
        function. Instead, a prior is calculated for each atomic
        charge, and a penalty function is added to the objective.
        If 0, don't use RESP priors for optimizing charges.
        If 1, the width for each atom is determined by the
        distribution of RESP charges for that atom in the training
        set.
        If 2, the width will be determined based on the RESP and ESP
        charges for that atom in the training set.
        
stride: int, Default = 50
        Specifies stride in sampling initial training data from the
        QM trajectory in coors.
        
tcout: str, Default = tc.out
        Specifies TeraChem output file with energies and gradients
        from the QM trajectory in coors.
        
initialtraining: bool, Default = True
        Specifies whether or not to sample initial training data from 
        a QM trajectory and optimize parameters before beginning MM 
        sampling in the first iteration. If False, an xyz file
        must still be provided in dynamicsdir as starting coordinates
        for MM sampling, but tcout is not used.
        
validinitial: bool, Default = False
        If true, compute validation set performance on the initial
        parameters at every iteration.
        
dryrun: bool, Default = False
        If true, check all folders and files, initialize everything,
        then stop before running any calculations.
        
skipchecks: bool, Default = False
        If true, skip checking for files and folders. Useful for
        debugging or if you just want to initialize an Input full
        of defaults for something. If you do this, you'll need to 
        call Input.__post_init__() later.
        
tctemplate: str, Default = tc_template.in
        Name of template TeraChem input file in sampledir. Should contain
        all relevant keywords, but ff-opt will take care of the 
        coordinates keyword for you.
        
tctemplate_backup: str, Default = tc_template_backup.in
        Name of backup template TeraChem input file in sampledir. Only
        used if the calculation with tctemplate fails first.
        
sbatchtemplate: str, Default = sbatch_template.sh
        Name of slurm sbatch template submission script. Should contain
        all relevant SBATCH options, as well as a script which runs
        the TeraChem input files. See examples/slurm/2_sampling for an
        example.
        
opt0: str, Default = opt_0.in
        Name of the template ForceBalance parameter optimization 
        input file. Should contain all relevant keywords and an 
        example target section which ff-opt will use as targets are added.
        
valid0: str, Default = valid_0.in
        Name of the template ForceBalance validation set evaluation
        input file. An example target section is unnecessary here.
        
batchsize: int, Default = 10
        Number of calculations to bundle together to send to a 
        ChemCloud server. You shouldn't mess with this number, except
        to decrease it if you're running a very large system.
        
retries: int, Default = 3
        Number of times failed calculations on ChemCloud are retried.
        
patience: int, Default = 5
        Number of consecutive iterations where the convergence criterion 
        must be satisfied before ending the calculation early. 
        
cutoff: float, Default = -1.0
        Cutoff for percentage change in validation criterion:
        If (\Chi_V(r_j, \Theta_j) - \Chi_V(r_j, \Theta_{j-1}) / 
            \Chi_V(r_j, \Theta_j) * 100 > cutoff, then we enter the
        patience period. This criterion should always be negative
        if the force field is improving, so cutoff should always be 
        a small, negative number.
        
