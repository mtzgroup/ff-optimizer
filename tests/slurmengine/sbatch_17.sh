#!/bin/bash

#SBATCH -t 2:00:00
#SBATCH -J FB_ref_gradient_17
#SBATCH --gres gpu:2
#SBATCH -p gpu_shortq
#SBATCH --qos=gpu_short
#SBATCH --mem=16GB
cp 17.xyz tc_17.in tc_backup_17.in $SCRATCH/.
cd $SCRATCH

ml TeraChem/2023.11-intel-2023.2.0-CUDA-12.2.1
terachem tc_17.in > tc_17.out
if [ $(grep -c "Job finished" tc_17.out) -ne 1 ]
then
  terachem tc_backup_17.in > tc_17.out
fi
