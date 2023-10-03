##Code for clustering on intervals 10 time point intervals - including the Spring to Fall (T1->11) interval 
###This script uses code and functions provided by S. Greemblum - https://github.com/greensii/dros-adaptive-tracking 
source('/home/users/mcbitter/OrchardProject/Code/config.R')
source('/home/users/mcbitter/OrchardProject/Code/helper_functions.R')
source('/home/users/mcbitter/OrchardProject/Code/load_packages.R')
source('/home/users/mcbitter/OrchardProject/Code/plotting_functions.R')
source('/home/users/mcbitter/OrchardProject/Code/workflow_functions.R')
source('/home/users/mcbitter/OrchardProject/Code/general_cage_functions.R')
setwd('~/dpetrov/MarkB/Orchard2021Data/')

## globals
chroms=c('2L','2R','3L','3R','X')
cages=1:12

## parameters
fdrThreshs=c(.2,.05,.01) ## maximum fdr-corrected pvalue for difference in allele frequency between treatments for a site to be considered significantly diverged
esThreshs=c(0,0.02,0.02) ## minimum mean allele frequency difference between treatments for a site to be considered significantly diverged (combined with p-value)
windowSize=500
windowShift=100
maxClusterGap=100
linkageThresh=0.03
maxSNPPairDist=1000000
poolSize=100

##Allele frequency data
HAFsFile='./RData/Downsampled/orch2021_Downsampled_ECage_Filtered.RData'
load(HAFsFile)

#Interval length glm to generate clustering results
lengths = c('DecWeekly')

for (length in lengths)
{
    glmFile=paste0('./09_ExpandingPhases/NonLOA/ExpandingPhases_', length , '_glm.ecage.orch2021.RData')
    load(glmFile)
    df.sig = read.csv(paste0('./09_ExpandingPhases/NonLOA/df.sig.ExpandingPhases.', length, '.orch21.csv'))
    comparisons = as.character(unique(df.sig$comparison))
    

    #######score windows in data in which site labels are shuffled and real data
    df.wins=score_wins(df.sig ,sites,windowSize,windowShift)
    df.sig.shuff=df.sig %>% group_by(comparison) %>% mutate(ix=sample(1:nrow(sites),n())) %>%
          ungroup() %>% mutate(chrom=sites$chrom[ix],pos=sites$pos[ix])
    df.wins.shuff=score_wins(df.sig.shuff,sites,windowSize,windowShift)
    df.winfdr=get_win_fdr(df.wins,df.wins.shuff)


    ##########Generating Clusters with funciton that difers slightly from S. Greenblum workflow 
    df.wins=df.wins %>%  arrange(chrom,comparison,winstartSNP,winSize)
    chroms=unique(df.wins$chrom)
    comparisons=levels(df.wins$comparison)
    cc_pairs=expand.grid(chroms,comparisons)

    df.clust = clust_wins.v2(df.wins, df.winfdr, maxBreak)


    write.csv(df.clust, paste0('./12_Clustering/AllPhases/df.clust.', length ,'.orch2021.UnMerged.csv'))


    df.snps = df.sig %>% filter(sigLevel>1)

    df.snps = do.call(rbind,lapply(1:nrow(df.clust),function(cl){
            df.snps %>% merge(df.clust[cl,] %>% dplyr::select(cl,comparison,chrom,startPos,endPos),by=c("comparison","chrom")) %>%
            filter(pos>=startPos,pos<=endPos)
            }))
    
    ###Randomly sampling snps for each comparison to make total number < 200K (proportionally across clusters/chromsomal arms) 
    ##this allows merging linked clusters to be computationally feasible
     d.snps.meta = data.frame()
        for (comp in as.character(unique(df.snps$comparison))){
            d.snps.comp = df.snps %>% filter(comparison == comp)
            if(nrow(d.snps.comp) > 200000){
                prct.down = 200000/nrow(d.snps.comp)
                d.snps.comp.sub = data.frame()
                for (i in as.character(unique(d.snps.comp$cl))){
                    d.comp.cl = d.snps.comp %>% filter(cl == i)
                    d.comp.cl.sub = data.frame()
                    for (chr in as.character(unique(d.comp.cl$chrom))) {
                        d.comp.cl.chr = d.comp.cl %>% filter(chrom == chr)
                        Nrow = nrow(d.comp.cl.chr)
                        d.comp.cl.chr.sub = d.comp.cl.chr %>% sample_n(ceiling((Nrow*prct.down)))
                        d.comp.cl.sub = rbind(d.comp.cl.sub, d.comp.cl.chr.sub)
                        }
                    d.snps.comp.sub = rbind(d.snps.comp.sub, d.comp.cl.sub)
                }
                }else{d.snps.comp.sub = d.snps.comp}
                d.snps.meta = rbind(d.snps.meta, d.snps.comp.sub)
            }
    
        
    df.snps = d.snps.meta
        
  
    
    
    ####Identify SNP pairs and calculate linkage between them
    df.pairs = df.snps %>% find_snp_pairs(maxDist=maxSNPPairDist) %>% filter(pairType=="inter")

    if (nrow(df.pairs) > 0){
        data = data.frame()
            for (chr in as.character(unique(df.pairs$chrom))){
                snpFile= paste0("~/dpetrov/MarkB/Orchard2021Data/founders/snptables/Orchard2021/inbredv2_withHets.orch2021.", chr, ".snpTable.numeric")
                df.pairs.chrom = df.pairs %>% filter(chrom == chr)
                snppairs = calc_Rsq_for_snp_pairs(df.pairs.chrom, ncores = 16, snpFile)
                snppairs = snppairs %>%
                    dplyr::select(-snp1.cl,-snp2.cl)
                data = rbind(data, snppairs)
            }
        snppairs = data
        df.final = merge_linked_clusters.v2(snppairs = snppairs, df.clust = df.clust, df.sig = df.sig, Rsq.thresh = linkageThresh)
        write.csv(df.final, paste0('./12_Clustering/AllPhases/df.clust.', length ,'.orch2021.csv'), row.names = FALSE)
      } else{write.csv(df.clust,paste0('./12_Clustering/AllPhases/df.clust.', length,'.orch2021.csv'), row.names = FALSE)}    
    }