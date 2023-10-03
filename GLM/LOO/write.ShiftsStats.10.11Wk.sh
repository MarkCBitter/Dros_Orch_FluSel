#!/bin/bash
LOC=${1}
PhaseLength=${2}
echo \
"
source('/home/users/mcbitter/OrchardProject/Code/plotting_functions.R')
source('/home/users/mcbitter/OrchardProject/Code/workflow_functions.R')
source('/home/users/mcbitter/OrchardProject/Code/general_cage_functions.R')
source('/home/users/mcbitter/OrchardProject/Code/Orchard2021_Scripts/LOA/loa.orch2021.functions.R')
setwd('/home/users/mcbitter/dpetrov/MarkB/Orchard2021Data/09_ExpandingPhases/LOA/')
HAFsFile = paste0('~/dpetrov/MarkB/Orchard2021Data/RData/Downsampled/RData_ByCage/orch2021_Downsampled_', '${LOC}', '_Filtered.RData')
BaseData = '~/dpetrov/MarkB/Orchard2021Data/RData/Downsampled/orch2021_Baseline_Downsampled_Filtered_Mean.RData'
WeeklyShifts = paste0('../WeeklyShifts_ByCage/af.shifts.Weekly.ecage.orch21.', '${LOC}', '.csv') 

#Load this cage/interval lengths significant sites file 
df = read.csv(paste0('./SigSites/SigMatched/df.sig.Matched.ExpandingPhases.LOC', '${LOC}', '.', '${PhaseLength}', 'Week.csv'))
#Loop through each comparison and get weekly shifts
comparisons = as.character((df %>% distinct(comparison))[,1])
df.shifts = data.frame()
for (comp in comparisons) {
    df.sig.matched.comp = df %>% filter(comparison == comp)
    df.shifts.comp = get.WeeklyShifts(df.sig.matched.comp,
         HAFsFile, WeeklyShifts) 
    df.shifts = rbind(df.shifts, df.shifts.comp)
}
#Subset shifts file into only those snps with FDR < 0.1 and removing spring to fall snps
df.shifts = df.shifts %>% mutate(snp = paste0(chrom, pos))
#Get statistics for all different shift dataframes
df.stats = get.shift.stats(df.shifts)
#Write output files                                                                                                                                                                                                 
write.csv(df.shifts, paste0('./FDR2/df.shifts.ExpandingPhases.LOC', '${LOC}', '.', '${PhaseLength}', 'Week.FDR.2.csv'), row.names = FALSE)                                                     write.csv(df.shifts, paste0('./FDR2/df.shifts.ExpandingPhases.LOC', '${LOC}', '.', '${PhaseLength}', 'Week.BERem.FDR.2.csv'), row.names = FALSE)
write.csv(df.stats, paste0('./FDR2/df.stats.ExpandingPhases.LOC', '${LOC}', '.', '${PhaseLength}', 'Week.FDR.2.csv'), row.names = FALSE)
write.csv(df.stats, paste0('./FDR2/df.stats.ExpandingPhases.LOC', '${LOC}', '.', '${PhaseLength}', 'Week.BERem.FDR.2.csv'), row.names = FALSE)
"  > ExpandingPhases.ShiftsStats.PhasesSeparate.LOC${LOC}.${PhaseLength}Week.R
