#!/bin/bash

#SBATCH -t 2:00:00
#SBATCH -J FB_ref_gradient_1
#SBATCH --gres gpu:2
#SBATCH -p gpu_shortq
#SBATCH --qos=gpu_short
#SBATCH --mem=16GB
cp 1.xyz tc_1.in tc_backup_1.in $SCRATCH/.
cd $SCRATCH

ml TeraChem/2023.11-intel-2023.2.0-CUDA-12.2.1
terachem tc_1.in > tc_1.out
if [ $(grep -c "Job finished" tc_1.out) -ne 1 ]
then
  terachem tc_backup_1.in > tc_1.out
fi
