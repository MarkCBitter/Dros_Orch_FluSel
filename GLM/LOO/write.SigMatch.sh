#!/bin/bash                                                                                                                                                                              \
                                                                                                                                                                                          
LOC=${1}
PhaseLength=${2}

echo \
"                                                                                                                                                                                         
source('/home/users/mcbitter/OrchardProject/Code/plotting_functions.R')                                                                                                                   
source('/home/users/mcbitter/OrchardProject/Code/workflow_functions.R')                                                                                                                   
source('/home/users/mcbitter/OrchardProject/Code/general_cage_functions.R')                                                                                                               
source('/home/users/mcbitter/OrchardProject/Code/Orchard2021_Scripts/LOA/loa.orch2021.functions.R')                                                                                       
setwd('/home/users/mcbitter/dpetrov/MarkB/Orchard2021Data/09_ExpandingPhases/LOA/SigSites/')                                                                                                

matched.sites = fread('../../../RData/Downsampled/ShuffledSites20X/Orch2021.ShuffledSites20x.csv')
df.sig = fread(paste0('df.sig.ExpandingPhases.LOC', '${LOC}', '.' ,'${PhaseLength}', 'Week.orch21.csv'))
comps = (df.sig %>% distinct(comparison)) %>% pull(comparison)


df.sig.matched = data.frame()
for (comp in comps){
    df.c = df.sig %>% filter(comparison == comp)
    df.c = df.c %>% mutate(snp = paste0(chrom, pos))
    df.c = df.c %>% filter(FDR < 0.2 & FDR >= 0.125 & (abs(afShift) > 0.005))
    mc.cores <- 14  # Set the number of cores you want to use here
    result_list <- mclapply(1:nrow(df.c), function(i) select.matched(df.c[i, ], matched.sites, df.c), mc.cores = mc.cores)
    # Combine the individual result DataFrames into a single DataFrame
    result_df <- do.call(rbind, result_list)
    df.sig.matched = rbind(df.sig.matched, result_df)
    }
                                                                                                                                                                                          
fwrite(df.sig.matched, paste0('df.sig.Matched.ExpandingPhases.LOC', '${LOC}', '.','${PhaseLength}', 'Week.FDR2.ES01.csv'))                                                  
                                                                                                                                                                                                                                                                                                                                                                                  \
                                                                                                                                                                                          
"  > ExpandingPhases.SigMatch.LOC${LOC}.${PhaseLength}Week.FDR2.ES01.R
