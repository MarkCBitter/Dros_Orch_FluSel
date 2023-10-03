#!/bin/bash
#SBATCH --job-name=SubsetRData
#SBATCH --time=06:00:00
#SBATCH --mem=300GB
#SBATCH --cpus-per-task=2
#SBATCH --partition=bigmem
#SBATCH --out=out.SubsetRData

ml R

Script=${1}

Rscript ${Script}

