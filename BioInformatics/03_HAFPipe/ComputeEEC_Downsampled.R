##Eec for time point samples   ****Real code run on cluster: 02E_ComputeEEC_Downsampled.sbatch

library('data.table')
library('tidyverse')

setwd('/scratch/groups/dpetrov/MarkB/Orchard2021Data/03_bams/DP_AlignStats/')

load("../../RData_V2/orch2021_Downsampled_META_RAW_V2.RData")

pct.missing = fread('../../founders/inbredv2.filtered.Orch21.PctMissing.csv')
gens = fread('/home/users/mcbitter/OrchardProject/Code/Orchard2021_Scripts/BioInformatics/03_HAFPipe_V2/HAF.params.csv')
colnames(gens) = c('tpt', 'gens')
gens = gens %>%
    mutate(tpt = as.character(tpt))
chroms = c("2L", "2R", "3L", "3R", "X")

n.snps = c(nrow(filter(sites, chrom == "2L")), nrow(filter(sites, chrom == "2R")), nrow(filter(sites, chrom == "3L")), nrow(filter(sites, chrom == "3R")), nrow(filter(sites, chrom == "X")))

chrom_length <- 23000000
recomb_rate <- 0.0000000239


#EEC function from S. Greenblum
calc_expected_ec=function(rd,gen,pct_missing,nof_snps,chrom_length,recomb_rate){
  
#example input: rd=5,gen=20,pct_missing=2,nof_snps=283438
# rd: actual read depth
# gen: nof generations since population founding
# pct_missing: percent of founder genotype calls that are missing (%, not a fraction)
# nof_snps: average number of snps per chromosome
# chrom_length: average length of a chromosome
# recomb_rate: average recombination rate

q=18

mycoeffs=data.frame(a=0.5199118,b=-0.6909052,c=0.3553630)

    
winSize=round(qexp(q/100,1/((chrom_length)/((recomb_rate)*(chrom_length)*(gen)+1)))/1000)

    
nReadsPerWin=(rd)*(nof_snps)*(winSize)*1000/(chrom_length)

    
ec = 10^(mycoeffs$a * log10( nReadsPerWin ) + mycoeffs$b * log10(1+pct_missing) + mycoeffs$c )

ec

}


##EEC matrix generation


library(data.table)
library(tidyverse)



files <- list.files(pattern="*MeanChromDP*", full.names = TRUE, recursive=FALSE)


data <- data.frame()
for (file in files){
    clean.file <- strsplit(file, "./")[[1]][2]
    full.sample.name <- strsplit(clean.file, ".MeanChromDP.txt")[[1]][1]
    chrom.dp = read.csv(file)
    chrom.dp$sample = full.sample.name
    chrom.dp = chrom.dp %>%
        dplyr::select(sample, chrom, mean_dp)
    data = rbind(data, chrom.dp)
            
}

samps.ext = samps
samps.ext$n.times = 5
samps.ext = as.data.frame(lapply(samps.ext, rep, samps.ext$n.times))
samps.ext = as.data.frame(samps.ext[,-2])
colnames(samps.ext) = "sample"
samps.ext = samps.ext %>%
    mutate(sample = as.character(sample))
data = data[order(match(data[,1],samps.ext[,1])),]


data = data %>%
    separate(sample, into = c("TP" , NA), sep = "_F1", remove = FALSE) %>%
    separate(TP, into = c(NA, 'tpt'), sep = "tp")
data = left_join(data, gens)
data = data %>% 
    mutate(num.snps = rep(as.numeric(as.character(n.snps)), length(unique(data$sample))),
           pct.missing = rep(pct.missing$pct, length(unique(data$sample))),
            chrom.length = chrom_length, 
           recomb.rate = recomb_rate)

data = data %>% rowwise() %>%
    mutate(eec = calc_expected_ec(mean_dp, gens, pct.missing, num.snps, chrom.length, recomb.rate))


sp.data = data %>%
    select(sample, chrom, eec) %>%
    spread(sample, eec)
sp.data = sp.data %>%
    mutate(n.times = as.numeric(as.character(n.snps)))
sp.data = as.data.frame(lapply(sp.data, rep, sp.data$n.times))

eec.sites = sp.data %>%
    select(chrom)
sp.data = sp.data %>%
    dplyr::select(-c(chrom, n.times))

eec.samps = as.data.frame(names(sp.data))


eec = sp.data


##adding eec matrix to orch 2021 baseline .RDATA
setwd("/home/users/mcbitter/dpetrov/MarkB/Orchard2021Data/RData_V2/")
load('orch2021_Downsampled_META_RAW_V2.RData')


df = samps
df = df %>%
    separate(sample.name, into = c("tpt",NA), sep = "_F1", remove = FALSE) %>%
    separate(tpt, into = c(NA, 'tpt'), sep = 'tp')
df = df %>%
    separate(sample.name, into = c(NA, 'cage'), sep = "F1_", remove = FALSE) %>%
    separate(cage, into = c('cage', 'seq.rep'), sep = "_", extra = "merge")
df = df %>%
    mutate(treatment = substr(cage, 1, 1))
df = df %>%
    separate(seq.rep, into = c('seq.rep', 'down'), sep = "down")
biol.reps = df %>%
    filter(grepl('rep', seq.rep))
tech.reps = df %>%
    filter(grepl('Rd', seq.rep))
df = df %>% mutate(biol.rep = ifelse(sample.name %in% biol.reps$sample.name, yes = "Yes", no = "No"))
df = df %>% mutate(tech.rep = ifelse(sample.name %in% tech.reps$sample.name, yes = "Yes", no = "No"))

df = df %>%
    rename(full.sample.name = sample.name)

df$tp = "tp"

df = df %>% rowwise() %>%
    mutate(sample = paste0(tp, tpt, "_", cage)) %>%
    dplyr::select(-c(down, tp, )) %>%
    dplyr::select(full.sample.name, sample, tpt, treatment, biol.rep, tech.rep)


samps = df

save(sites, samps, afmat, eec, file = './orch2021_Downsampled_META_RAW_V2.RData')