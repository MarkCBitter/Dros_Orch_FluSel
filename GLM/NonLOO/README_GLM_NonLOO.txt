###README GLM's that *ARE NOT* run using the leave-one-out approach

##This directory contains a series of R scripts to generate glm results for all possible interval lengths/comparison.
##Each script runs all possible glm's for a distinct interval length - e.g. 'ExpandingPhases.BiWeekly.glm.ecage.orch21.R' runs glm's for all two timepoint intervals (e.g. T1->3...T10->12).


##The subdirectory GLM_Scripts_NonLOA_DiscreteTime has code for GLM's run using just the first and final time point of the interval in the regression. So, for regression fo t1->3, just using allele frequency data from t1 and t3. 

##The notebook  GetSigSItes_ExpPhas_NonLOA_Final.ipynb contains code for taking in glm results (generated via ExpandingPHases.XWeekly.glm* result files) and identifes significant sites for each comparison 

##The notebook GLMSigSiteAssessment.orch2021.Final.ipynb contains code analyzing the number of significant sites for all GLM's

##The notebook orch2021_ManhattanPlots_Final.ipynb generates manhattan plots for all interval length/comparison glm's

##The script for clustering of unlinked loci from glm results is Clustering.R 
###IMPORTANT - this clustering code is modeled off of code written and provided by S. Greenblum here - https://github.com/greensii/dros-adaptive-tracking
##This script specifically generates results for 10 timepoint intervals, which includes the Spring to Fall interval (T1->11)
##This code was similarly used to generate results for the T1->7, 3->9, 5->12, and 7->9 intervals for the comparison to Rudman et al. (2022) results

