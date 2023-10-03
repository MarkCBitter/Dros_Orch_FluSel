#!/bin/bash

LOC=${1}
PhaseLength=${2}
LastTP=${3}


echo \
"source('/home/users/mcbitter/OrchardProject/Code/config.R')
source('/home/users/mcbitter/OrchardProject/Code/helper_functions.R')
source('/home/users/mcbitter/OrchardProject/Code/load_packages.R')
source('/home/users/mcbitter/OrchardProject/Code/plotting_functions.R')
source('/home/users/mcbitter/OrchardProject/Code/workflow_functions.R')
source('/home/users/mcbitter/OrchardProject/Code/general_cage_functions.R')
setwd('~/dpetrov/MarkB/Orchard2021Data/')

load(paste0('./RData/Downsampled/LOA_RData/orch2021_ECage_LOA_', '${LOC}', '.RData'))

df.glm = sites


tps = 1:${LastTP}

for (tp in tps){
    i = tp
    i.2 = i + ${PhaseLength}
    vec = c(i:i.2)
    load(paste0('./RData/Downsampled/LOA_RData/orch2021_ECage_LOA_', '${LOC}', '.RData'))
    afmat = t(afmat)
    eec = t(eec)
    afmat = cbind(samps, afmat)
    eec = cbind(samps, eec)
    afmat = afmat %>%
        filter(tpt %in% vec)
    eec = eec %>%
        filter(tpt %in% vec)
    samps = afmat[,1:ncol(samps)]
    samps = samps %>% mutate(tpt = as.numeric(tpt))
    afmat = afmat[,-c(1:ncol(samps))]
    afmat = t(afmat)
    afmat = as.data.frame(afmat)
    eec = eec[,-c(1:ncol(samps))]
    eec = t(eec)
    eec = as.data.frame(eec)
    res = fit_GLM_ContinuousTime(afMatrix = afmat ,rdMatrix = eec, vec = vec, sampleData = samps, model.vars = 'tpt', poolCt=100,ncores = 16)
    df.glm = cbind(df.glm, res)
    save(df.glm, file = paste0('./09_ExpandingPhases/LOA/GLMRes/ExpandingPhases.LOC', '${LOC}' , '.', '${PhaseLength}', 'Week.glm.ecage.orch2021.RData'))
    }



" > ExpandingPhases.LOA.${LOC}.${PhaseLength}Week.R
