###Fst between technical replicates - run on cluster TechReps_Fst.R
setwd('/scratch/groups/dpetrov/MarkB/Orchard2021Data/RData/Downsampled/')
source("/home/users/mcbitter/OrchardProject/Code/config.R")
source("/home/users/mcbitter/OrchardProject/Code/general_cage_functions.R")

load('./orch2021_Downsampled_ECage_TechRep.RData')
rep.samps = unique(samps$sample)

data = data.frame()
for (samp in rep.samps){
    load('./orch2021_Downsampled_ECage_TechRep.RData')
    df = cbind(samps, t(afmat))
    df = df %>%
        filter(sample == samp)
    df = df[,-(1:ncol(samps))]
    df = t(df)
    rep = samp
    fst.mat = Fst.mat(df)
    fst.val = fst.mat[1, 2]
    fst.val = cbind(rep, fst.val)
    data = rbind(data, fst.val)
    
}


write.csv(data,"../../FST_Analysis/TechReps_Fst_orch2021.csv", row.names = FALSE)