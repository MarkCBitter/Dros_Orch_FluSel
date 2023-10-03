source("/home/users/mcbitter/OrchardProject/Code/config.R")
source("/home/users/mcbitter/OrchardProject/Code/helper_functions.R")
source("/home/users/mcbitter/OrchardProject/Code/load_packages.R")
source("/home/users/mcbitter/OrchardProject/Code/plotting_functions.R")
source("/home/users/mcbitter/OrchardProject/Code/workflow_functions.R")
source("/home/users/mcbitter/OrchardProject/Code/general_cage_functions.R")
setwd('~/dpetrov/MarkB/Orchard2021Data/')

load('./RData/Downsampled/orch2021_Downsampled_ECage_Filtered.RData')
df.glm = sites


tps = 1:9

for (tp in tps){
    i = tp 
    i.2 = i+3
    vec = c(i:i.2)
    t1 = vec[1]
    tf = vec[length(vec)]
    vec = c(t1, tf)
    load('./RData/Downsampled/orch2021_Downsampled_ECage_Filtered.RData')
    afmat = t(afmat)
    eec = t(eec)
    afmat = cbind(samps, afmat)
    eec = cbind(samps, eec)
    afmat = afmat %>%
        filter(tpt %in% vec)
    eec = eec %>%
        filter(tpt %in% vec)
    samps = afmat[,1:ncol(samps)]
    samps$tpt = as.factor(samps$tpt)
    afmat = afmat[,-c(1:ncol(samps))]
    afmat = t(afmat)
    afmat = as.data.frame(afmat)
    eec = eec[,-c(1:ncol(samps))]
    eec = t(eec)
    eec = as.data.frame(eec)
    res = fit_GLM_ContinuousTime(afMatrix = afmat ,rdMatrix = eec, vec = vec, sampleData = samps, model.vars = 'tpt', poolCt=100,ncores = 8)
    df.glm = cbind(df.glm, res)
    save(df.glm, file = "./05_GLM/ExpandingPhases_TriWeekly_glm_Discrete.ecage.orch2021.RData")
    }
