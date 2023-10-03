#Getting LOA RData
source('/home/users/mcbitter/OrchardProject/Code/plotting_functions.R')
source('/home/users/mcbitter/OrchardProject/Code/workflow_functions.R')
source('/home/users/mcbitter/OrchardProject/Code/general_cage_functions.R')
source('/home/users/mcbitter/OrchardProject/Code/Orchard2021_Scripts/LOA/loa.orch2021.functions.R')
setwd('/home/users/mcbitter/dpetrov/MarkB/Orchard2021Data/RData/Downsampled/')

load('./orch2021_Downsampled_ECage_Filtered.RData', verbose = TRUE)
cages = as.character(unique(samps$cage))
ncol.samps = ncol(samps)
df = cbind(samps, t(afmat))
df.eec = cbind(samps, t(eec))

for (c in cages){
    df.lo = df %>% filter(cage != c)
    df.eec.lo = df.eec %>% filter(cage != c)
    samps = df.lo[,1:ncol.samps]
    afmat = df.lo[,-c(1:ncol.samps)]
    afmat = t(afmat)
    eec = df.eec.lo[,-c(1:ncol.samps)]
    eec = t(eec)
    save(samps, sites, afmat, eec, file = paste0('./LOA_RData/orch2021_ECage_LOA_', c, '.RData'))
}
