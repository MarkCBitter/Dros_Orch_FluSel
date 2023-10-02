##Getting baseline vs. Ecage

setwd('/scratch/groups/dpetrov/MarkB/Orchard2021Data/RData/Downsampled/')
source("/home/users/mcbitter/OrchardProject/Code/config.R")
source("/home/users/mcbitter/OrchardProject/Code/general_cage_functions.R")

TP = c(1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12)


data = data.frame()

for (tp in TP){
    load('./orch2021_Downsampled_ECage_Filtered.RData')
    df = cbind(samps, t(afmat))
    df = df %>%
    filter(tpt == tp)
    df = df[,-(1:ncol(samps))]
    df = t(df)
    load('./orch2021_Baseline_Downsampled_Filtered.RData')
    df.fst = cbind(afmat.base, df)
    fst.mat = Fst.mat(df.fst)
    fst.mat = (fst.mat[-c(5:nrow(fst.mat)),-c(1:4)])
    fst.vec = as.vector(fst.mat)
    cages = colnames(df)
    time.seg = paste0("0_", tp)
    type = 'Base.TP'
    fst.df = as.data.frame(fst.vec)
    fst.df$cage = rep(cages, 4)
    fst.df$time.seg = time.seg
    fst.df$type = type
    fst.df = fst.df %>%
        rename(fst = fst.vec) %>%
        dplyr::select(cage, fst, time.seg, type)
    data = rbind(data, fst.df)

}
write.csv(data, '../../FST/BaselineTP_Fst_ECage_Orch2021.csv', row.names = FALSE)