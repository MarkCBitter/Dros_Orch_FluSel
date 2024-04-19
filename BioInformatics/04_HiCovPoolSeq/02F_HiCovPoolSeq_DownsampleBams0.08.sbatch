#!/bin/bash                                                                                                                                                                                  
#SBATCH --job-name=02F_Downsamp                                                                                                                                                              
#SBATCH --time=16:00:00                                                                                                                                                                      
#SBATCH --mem=400GB                                                                                                                                                                           
#SBATCH --cpus-per-task=2                                                                                                                                                                   
#SBATCH --partition=bigmem                                                                                                                                                                      
#SBATCH --out=out.02F_Downsamp 


ml biology
ml samtools

cd ~/dpetrov/MarkB/Orchard2021Data/18_HAFPipe_Validation/Bams/
samtools view -bs 0.8094217 tp11_F1_E1.HiCov.bam > tp11_F1_E1.100x.bam
samtools view -bs 0.8076435 tp11_F1_E3.HiCov.bam > tp11_F1_E3.100x.bam
samtools view -bs 0.6059908 tp11_F1_E4.HiCov.bam > tp11_F1_E4.100x.bam
samtools view -bs 0.6468640 tp11_F1_E5.HiCov.bam > tp11_F1_E5.100x.bam
samtools view -bs 0.8696938 tp11_F1_E7.HiCov.bam > tp11_F1_E7.100x.bam
samtools view -bs 0.9698192 tp11_F1_E8.HiCov.bam > tp11_F1_E8.100x.bam
samtools view -bs 0.6447827 tp3_F1_E1.HiCov.bam > tp3_F1_E1.100x.bam
samtools view -bs 0.9761050 tp3_F1_E10.HiCov.bam > tp3_F1_E10.100x.bam
samtools view -bs 0.9724316 tp3_F1_E11.HiCov.bam > tp3_F1_E11.100x.bam
#samtools view -bs 0.76183 tp3_F1_E12.HiCov.bam > tp3_F1_E12.100x.bam
samtools view -bs 0.7032299 tp3_F1_E2.HiCov.bam > tp3_F1_E2.100x.bam
samtools view -bs 0.7240081 tp3_F1_E3.HiCov.bam > tp3_F1_E3.100x.bam
samtools view -bs 0.6964516 tp3_F1_E4.HiCov.bam > tp3_F1_E4.100x.bam
samtools view -bs 0.5960079 tp3_F1_E6.HiCov.bam > tp3_F1_E6.100x.bam
samtools view -bs 0.6987144 tp3_F1_E7.HiCov.bam > tp3_F1_E7.100x.bam
samtools view -bs 0.9735866 tp3_F1_E8.HiCov.bam > tp3_F1_E8.100x.bam

mv *.100x.bam ../Bams.100x/


cd ~/dpetrov/MarkB/Orchard2021Data/18_HAFPipe_Validation/Bams.100x/
for i in ./*.bam;do
    samtools view -bs 0.08 $i > ${i/.100x.bam/.08x.A.bam} 
done


mv *.A.bam ../Bams.08x.A/

cd ~/dpetrov/MarkB/Orchard2021Data/18_HAFPipe_Validation/Bams.100x//
for i in ./*.bam;do
    samtools view -bs 43.08 $i > ${i/.100x.bam/.08x.B.bam} 
done

mv *.B.bam ../Bams.08x.B/
