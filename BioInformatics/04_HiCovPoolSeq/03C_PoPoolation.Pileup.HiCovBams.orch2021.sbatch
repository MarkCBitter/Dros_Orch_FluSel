#!/bin/bash                                                                                                                                                                             
#SBATCH --job-name=03C_GeneratePileups                                                                                                                                                  
#SBATCH --time=05:00:00                                                                                                                                                                 
#SBATCH --mem=300GB                                                                                                                                                                     
#SBATCH --cpus-per-task=2                                                                                                                                                                
#SBATCH --partition=bigmem                                                                                                                                                               
#SBATCH --out=out.03C_GeneratePileups                                                                                                                                                    
module load biology
module load samtools perl java
Dir=${1}
Cov=${2}

cd ${Dir}

ls *.bam > samps.txt

samtools mpileup -B $(<samps.txt) -f ../all_dmel.fasta -l ../orch2021.sites.txt > hi.cov.orch21sites.${Cov}.mpileup
