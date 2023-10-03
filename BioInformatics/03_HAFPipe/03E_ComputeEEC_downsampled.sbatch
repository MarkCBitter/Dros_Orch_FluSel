#!/bin/bash
#SBATCH --job-name=ComputeEEC_downsampled
#SBATCH --time=03:00:00
#SBATCH --mem=300GB
#SBATCH --cpus-per-task=2
#SBATCH --partition=bigmem
#SBATCH --out=out.ComputeEEC_downsampled 

ml R

Rscript ComputeEEC_Downsampled.R

