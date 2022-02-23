#!/bin/bash

#SBATCH -t 2:00:00
#SBATCH -J FB_ref_gradient_XXX
#SBATCH --gres gpu:2
#SBATCH --mem=16GB
#SBATCH --fin=tc_XXX.in,XXX.pdb,tc_XXX_long.in
#SBATCH --fout=tc_XXX.out
#SBATCH --exclude=fire-11-01
cd $SCRATCH

ml TeraChem/2021.02-intel-2017.8.262-CUDA-9.0.176
#ml TeraChem/2019.12-intel-2017.0.098-MPICH2-1.4.1p1-CUDA-9.0.176
terachem tc_XXX.in > tc_XXX.out
if [ $(grep -c "Job finished" tc_XXX.out) -ne 1 ]
then
  terachem tc_XXX_long.in > tc_XXX.out
fi
