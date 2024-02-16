#!/bin/bash

#SBATCH -t 2:00:00
#SBATCH -J FB_ref_gradient_13
#SBATCH --gres gpu:2
#SBATCH -p gpu_shortq
#SBATCH --qos=gpu_short
#SBATCH --mem=16GB
cp 13.xyz tc_13.in tc_backup_13.in $SCRATCH/.
cd $SCRATCH

ml TeraChem/2023.11-intel-2023.2.0-CUDA-12.2.1
terachem tc_13.in > tc_13.out
if [ $(grep -c "Job finished" tc_13.out) -ne 1 ]
then
  terachem tc_backup_13.in > tc_13.out
fi
