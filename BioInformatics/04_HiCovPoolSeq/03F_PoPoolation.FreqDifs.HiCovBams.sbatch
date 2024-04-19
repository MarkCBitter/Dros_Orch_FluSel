#!/bin/bash                                                                                                                                                                                              
#SBATCH --job-name=03D_GenerateFreqs                                                                                                                                                                      
#SBATCH --time=05:00:00                                                                                                                                                                                  
#SBATCH --mem=700GB                                                                                                                                                                                      
#SBATCH --cpus-per-task=2                                                                                                                                                                                
#SBATCH --partition=bigmem                                                                                                                                                                               
#SBATCH --out=out.03E_GenerateFreqs                                                                                                                                                                      
module load biology
module load samtools perl java
Dir=${1}
cd ${Dir}

perl ~/dpetrov/SOFTWARE/popoolation2_1201/snp-frequency-diff.pl --input hi.cov.sync --output-prefix hi.cov --min-count 6 --min-coverage 50 --max-coverage 200
