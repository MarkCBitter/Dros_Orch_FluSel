#!/bin/bash                                                                                                                                                                                              
#SBATCH --job-name=03E_GenerateFreqs                                                                                                                                                                      
#SBATCH --time=00:25:00                                                                                                                                                                                  
#SBATCH --mem=70GB                                                                                                                                                                                      
#SBATCH --cpus-per-task=2                                                                                                                                                                                
#SBATCH --partition=hns                                                                                                                                                                               
#SBATCH --out=out.03E_GenerateFreqs                                                                                                                                                                      
module load biology
module load samtools perl java
Dir=${1}
cd ${Dir}

./calfreq -i hi.cov.orch21sites.sync -f 0 -d 0,999 -b 0 -m 0.05
