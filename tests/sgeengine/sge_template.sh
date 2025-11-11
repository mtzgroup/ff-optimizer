#$ -S /bin/bash                                
#$ -N FB_ref_gradient_JOBID                               
#$ -cwd                                        
#$ -l h_rt=8:00:00                             
#$ -l h_rss=8G                                 
#$ -q small                                    
#$ -l gpus=1                                   
#$ -pe smp 1                                   
# -fin TCTEMPLATE
# -fin TCTEMPLATEBACKUP                                  
# -fin JOBID.xyz                                                          
# -fout scr/                                   
                                               
module purge                                   
module load terachem/terachem-1.9H-2025.03-dev
export OMP_NUM_THREADS=1                              
terachem TCTEMPLATE > $SGE_O_WORKDIR/tc_JOBID.out
if [ $(grep -c "Job finished" $SGE_O_WORKDIR/tc_JOBID.out) -ne 1 ]
then
  terachem TCTEMPLATEBACKUP > $SGE_O_WORKDIR/tc_JOBID.out
fi
sleep 60                                       