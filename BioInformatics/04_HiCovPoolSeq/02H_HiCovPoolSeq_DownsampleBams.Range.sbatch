#!/bin/bash                                                                                                                                                                                                        
#SBATCH --job-name=02F_Downsamp                                                                                                                                                                                    
#SBATCH --time=16:00:00                                                                                                                                                                                            
#SBATCH --mem=400GB                                                                                                                                                                                                
#SBATCH --cpus-per-task=2                                                                                                                                                                                          
#SBATCH --partition=bigmem                                                                                                                                                                                         
#SBATCH --out=out.02F_Downsamp                                                                                                                                                                                     

ml biology
ml samtools


cd ~/dpetrov/MarkB/Orchard2021Data/18_HAFPipe_Validation/Bams.100x/
for i in ./*.bam;do
    samtools view -bs 0.25 $i > ${i/.100x.bam/.25x.bam}
done
mv *.25x.bam ../Bams.25x/

cd ~/dpetrov/MarkB/Orchard2021Data/18_HAFPipe_Validation/Bams.100x/
for i in ./*.bam;do
    samtools view -bs 0.50 $i > ${i/.100x.bam/.50x.bam}
done

mv *.50x.bam ../Bams.5x/

cd ~/dpetrov/MarkB/Orchard2021Data/18_HAFPipe_Validation/Bams.100x/
for i in ./*.bam;do
    samtools view -bs 0.75 $i > ${i/.100x.bam/.75x.bam}
done

mv *.75x.bam ../Bams.75x/

