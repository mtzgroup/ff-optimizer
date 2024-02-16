#!/bin/bash

#SBATCH -t 1:00:00
#SBATCH -J FB_ref_gradient_JOBID
#SBATCH --gres gpu:1
#SBATCH -p gpuq
#SBATCH --mem=16GB
cp JOBID.xyz TCTEMPLATE TCTEMPLATEBACKUP $SCRATCH/.
cd $SCRATCH

ml TeraChem/2023.11-intel-2023.2.0-CUDA-12.2.1
terachem TCTEMPLATE > tc_JOBID.out
if [ $(grep -c "Job finished" tc_JOBID.out) -ne 1 ]
then
  rm -rf scr*
  terachem TCTEMPLATEBACKUP > tc_JOBID.out
fi
