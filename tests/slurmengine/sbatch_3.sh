#!/bin/bash

#SBATCH -t 2:00:00
#SBATCH -J FB_ref_gradient_3
#SBATCH --gres gpu:2
#SBATCH -p gpu_shortq
#SBATCH --qos=gpu_short
#SBATCH --mem=16GB
cp 3.xyz tc_3.in tc_backup_3.in $SCRATCH/.
cd $SCRATCH

ml TeraChem/2023.11-intel-2023.2.0-CUDA-12.2.1
terachem tc_3.in > tc_3.out
if [ $(grep -c "Job finished" tc_3.out) -ne 1 ]
then
  terachem tc_backup_3.in > tc_3.out
fi
