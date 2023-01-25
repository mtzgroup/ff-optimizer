#!/bin/bash

#SBATCH -t 2:00:00
#SBATCH --gres gpu:2
#SBATCH --mem=16GB
#SBATCH --exclude=fire-11-01
#SBATCH -J FB_ref_gradient_23
#SBATCH --fin=tc_23.in,23.xyz,tc_23_backup.in
#SBATCH --fout=tc_23.out,esp_23.xyz
cd $SCRATCH

module unload MPICH2
#ml TeraChem/2021.02-intel-2017.8.262-CUDA-9.0.176
ml TeraChem/2019.12-intel-2017.0.098-MPICH2-1.4.1p1-CUDA-9.0.176
terachem tc_23.in > tc_23.out
if [ $(grep -c "Job finished" tc_23.out) -ne 1 ]
then
  terachem tc_23_backup.in > tc_23.out
fi
mv scr.23/esp.xyz esp_23.xyz
