#!/bin/bash                                                                                                                                
                                                                                                                                           
#SBATCH --job-name=02F_Downsamp                                                                                                            
                                                                                                                                           
#SBATCH --time=16:00:00                                                                                                                   
                                                                                                                                       
#SBATCH --mem=400GB                                                                                                                       
                                                                                                                                        
#SBATCH --cpus-per-task=2                                                                                                                  
                                                                                                                                        
#SBATCH --partition=bigmem                                                                                                             
                                                                                                                                       
#SBATCH --out=out.02F_Downsamp                                                                                                            
ml biology
ml samtools






cd ~/dpetrov/MarkB/Orchard2021Data/18_HAFPipe_Validation/Bams.08x.A/
for i in ./*.bam;do
    samtools view -bs 0.25 $i > ${i/.08x.A.bam/.04.bam}
done