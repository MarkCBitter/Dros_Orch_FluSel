# Dros_Orch_FluSel
Code and allele frequency data repository for analysis of 2021 Drosophila melanogaster mesocosm data
Each subdirectory contains a seaparate README file containing specific information regarding the data and code within the directory.

The bioinformatics directory contains code for processing raw reads and ultimately producing haplotype-informed allele frequencies.
The software used for these analyses include:
1) bwa: https://bio-bwa.sourceforge.net/
2) samtools: https://www.htslib.org/
3) picard: https://broadinstitute.github.io/picard/
4) bamtools: https://github.com/pezmaster31/bamtools/wiki
5) HAFPipe: https://github.com/petrov-lab/HAFpipe-line (which necessitates dependencies: HARP (download from https://bitbucket.org/dkessner/harp); Python 3 (including numpy for the 'npute' method of Task 2); and tabix and bgzip (http://www.htslib.org/download/)

The remaining directories are all sets of scripts (each containing their own README) for analysis of allele frequency data. 
These analyses were conducited in the R v. 3.5.6: https://www.r-project.org/
The orch2021_functions.R script in the main directory is source code that contains functions that are called upon in many of the scripts and notebooks used in the various subdirectories.
The Baseline_NaturalPop_comparison directory contains a notebook for assessing differentiation between the founder/baseline population and natural popualtions of D. melanogaster.
The Fst directory contains notebooks and scripts for assessing patterns of Fst throughout the experiment.
The GLM directory contains all analyses associated with the use of generlized linear models to infer systematic allele frequency change across replicates. These include the leave-out-out cross validation of GLM results (located in the LOO subdirectory) as well as the scripts for running GLM across all replicates and using these results to identify clusters of unlinked loci underpinning the GLM signal (NonLOO subdirectory)
      -Notably, the clustering of unlinked loci draws upon code/repository found here: https://github.com/greensii/dros-adaptive-tracking
      -The analysis of correlations in allele frequency movement between SNPs was conducted using code deposited here:    
       https://github.com/EgorLappo/temporal-covariances/tree/main/for_sharing

The allele frequency files used in these analyses are too large to upload to github. 
But you can access the raw reads for the samples collected for this experiment on the NCBI SRA: PRJNA1031645
And can access the allele frequency data here: https://doi.org/10.5061/dryad.xd2547dpv
You can also run the code on a small subset of the data for the evolved allele frequency samples and baseline samples provided in the Data directory







