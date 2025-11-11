#$ -S /bin/bash                                
#$ -N FB_ref_gradient_23                               
#$ -cwd                                        
#$ -l h_rt=8:00:00                             
#$ -l h_rss=8G                                 
#$ -q small                                    
#$ -l gpus=1                                   
#$ -pe smp 1                                   
# -fin tc_23.in
# -fin tc_backup_23.in                                  
# -fin 23.xyz                                                          
# -fout scr/                                   
                                               
module purge                                   
module load terachem/terachem-1.9H-2025.03-dev
export OMP_NUM_THREADS=1                              
terachem tc_23.in > $SGE_O_WORKDIR/tc_23.out
if [ $(grep -c "Job finished" $SGE_O_WORKDIR/tc_23.out) -ne 1 ]
then
  terachem tc_backup_23.in > $SGE_O_WORKDIR/tc_23.out
fi
sleep 60                                       