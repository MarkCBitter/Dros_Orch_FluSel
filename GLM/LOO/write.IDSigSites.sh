#!/bin/bash                                                                                                                                                                                                            
LOC=${1}
PhaseLength=${2}

echo \
"  
source('/home/users/mcbitter/OrchardProject/Code/plotting_functions.R')
source('/home/users/mcbitter/OrchardProject/Code/workflow_functions.R')
source('/home/users/mcbitter/OrchardProject/Code/general_cage_functions.R')
source('/home/users/mcbitter/OrchardProject/Code/Orchard2021_Scripts/LOA/loa.orch2021.functions.R')
setwd('/home/users/mcbitter/dpetrov/MarkB/Orchard2021Data/09_ExpandingPhases/LOA/GLMRes/')
fdrThreshs=c(.2,.05,.01)
esThreshs=c(0,0.02,0.02)
#Load leave-one-out RDATA
load(paste0('~/dpetrov/MarkB/Orchard2021Data/RData/Downsampled/LOA_RData/orch2021_ECage_LOA_', '${LOC}','.RData'))###ALTER

#Load leave-one-out GLM results file
load(paste0('./ExpandingPhases.LOC', '${LOC}.', '${PhaseLength}' , 'Week.glm.ecage.orch2021.RData'))

sites = df.glm %>% dplyr::select(chrom, pos)

comparisons = grep('coef.', names(df.glm), value = TRUE)
for (i in (1:length(comparisons))) {
    comparisons[i] = strsplit(comparisons[i], '[.]')[[1]][2]
} 

af.shifts <- get_af_shifts(afmat, samps, cage_set = NULL, comparisons) 

sites = df.glm %>% dplyr::select(chrom, pos) 
FDR = get_glm_FDR.V2(df.glm) 

df.sig <- get_sig_sites(df.glm, comparisons, FDR, afShifts = af.shifts, fdrThreshs, esThreshs)

write.csv(df.sig, paste0('df.sig.ExpandingPhases.LOC', '${LOC}', '.' ,'${PhaseLength}', 'Week.orch21.csv'), row.names = FALSE)
                                                                                                                                                                                                                                                                                                                                                                                 
"  > ExpandingPhases.IDSigSites.LOC${LOC}.${PhaseLength}Week.R
