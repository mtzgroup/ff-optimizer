#!/bin/bash

# Slurm parameters, you may need to adapt these for your cluster

#SBATCH -t 1:00:00
#SBATCH -J FB_ref_gradient_JOBID
#SBATCH --gres gpu:1
#SBATCH -p gpu_shortq
#SBATCH --qos=gpu_short
#SBATCH --mem=8GB

# To write the sbatch files used by ff-opt, it replaces certain
# all-caps variables with the appropriate value: 
# JOBID - the ID of the xyz file used for a gradient calculation
# TCTEMPLATE - name of the terachem template input
# TCTEMPLATEBACKUP - name of the backup terachem template input

# The Martinez lab cluster uses a scratch folder system where
# $SCRATCH is set in the environment. You can remove this if
# your cluster doesn't do something similar. 
cp JOBID.xyz TCTEMPLATE TCTEMPLATEBACKUP $SCRATCH/.
cd $SCRATCH

# Load the appropriate terachem module
module load TeraChem
terachem TCTEMPLATE > tc_JOBID.out
echo JOBID

# Re-run using backup tc input if the job failed
if [ $(grep -c "Job finished" tc_JOBID.out) -ne 1 ]
then
  terachem TCTEMPLATEBACKUP > tc_JOBID.out
fi
