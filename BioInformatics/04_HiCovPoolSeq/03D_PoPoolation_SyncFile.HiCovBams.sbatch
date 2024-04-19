#!/bin/bash                                                                                                                                                                                              
#SBATCH --job-name=03D_GenerateSync                                                                                                                                                                   
#SBATCH --time=01:00:00                                                                                                                                                                                  
#SBATCH --mem=300GB                                                                                                                                                                                     
#SBATCH --cpus-per-task=2                                                                                                                                                                               
#SBATCH --partition=bigmem                                                                                                                                                                              
#SBATCH --out=out.03D_GenerateSync                                                                                                                                                                   
module load biology
module load samtools perl java
Dir=${1}
Cov=${2}

cd ${Dir}


java -ea -Xmx7g -jar ~/dpetrov/SOFTWARE/popoolation2_1201/mpileup2sync.jar --input hi.cov.orch21sites.${Cov}.mpileup --output hi.cov.orch21sites.${Cov}.sync --fastq-type   sanger --min-qual 20 --threads 2
