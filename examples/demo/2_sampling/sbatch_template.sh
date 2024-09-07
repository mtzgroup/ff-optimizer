#!/bin/bash

#SBATCH -t 1:00:00
#SBATCH -J FB_ref_gradient_JOBID
#SBATCH --gres gpu:1
#SBATCH -p gpu_shortq
#SBATCH --qos=gpu_short
#SBATCH --mem=8GB
cp JOBID.xyz TCTEMPLATE TCTEMPLATEBACKUP $SCRATCH/.
cd $SCRATCH

ml TeraChem
terachem TCTEMPLATE > tc_JOBID.out
echo JOBID
if [ $(grep -c "Job finished" tc_JOBID.out) -ne 1 ]
then
  terachem TCTEMPLATEBACKUP > tc_JOBID.out
fi
