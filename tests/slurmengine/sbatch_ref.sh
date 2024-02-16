#!/bin/bash

#SBATCH -t 2:00:00
#SBATCH -J FB_ref_gradient_23
#SBATCH --gres gpu:2
#SBATCH -p gpu_shortq
#SBATCH --qos=gpu_short
#SBATCH --mem=16GB
cp 23.xyz tc_23.in tc_backup_23.in $SCRATCH/.
cd $SCRATCH

ml TeraChem/2023.11-intel-2023.2.0-CUDA-12.2.1
terachem tc_23.in > tc_23.out
if [ $(grep -c "Job finished" tc_23.out) -ne 1 ]
then
  terachem tc_backup_23.in > tc_23.out
fi
