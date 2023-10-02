###Fst between biological replicates - run on cluster BiolReps_ECage_Fst.R
setwd('/scratch/groups/dpetrov/MarkB/Orchard2021Data/RData/Downsampled/')
source("/home/users/mcbitter/OrchardProject/Code/config.R")
source("/home/users/mcbitter/OrchardProject/Code/general_cage_functions.R")

load('./orch2021_Downsampled_BiolReps_ECages.RData')
rep.samps = unique(samps$sample)


data = data.frame()
for (samp in rep.samps){
    load('./orch2021_Downsampled_BiolReps_ECages.RData')
    df = cbind(samps, t(afmat))
    df = df %>%
        filter(sample == samp)
    df = df[,-(1:ncol(samps))]
    df = t(df)
    rep = samp
    fst.mat = Fst.mat(df)
    fst.val = fst.mat[1, 2]
    fst.df = cbind(rep, fst.val)
    fst.df = fst.df %>% as.data.frame(fst.df) %>%
                rename(fst = fst.val, cage = rep) %>%
                mutate(type = "biol.rep") %>%
                mutate(time.seg = "biol.rep") %>%
                dplyr::select(cage, fst, time.seg, type)
    data = rbind(data, fst.df)    
}


write.csv(data, "../../FST/BiolReps_ECage_Fst_orch2021.csv", row.names = FALSE)