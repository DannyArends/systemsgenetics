#
# test.mappings.R
#
# Copyright (c) 2013-2013 GBIC: Danny Arends, Harm-Jan Westra, Lude Franke and Ritsert C. Jansen
# last modified Jul, 2013
# first written Jun, 2013
#

setwd("~/Github/systemsgenetics/eqtl-mapping-pipeline/src/main/scripts")
source("helper.functions.R")

setwd("~/Github/Juha/")
# Load the Illumina celltype data
Illu <- read.illumina.celltypes()
Illu <- annotate.illumina.celltypes(Illu)

WBillu <- read.illumina.wb()
WBillu <- add.illumina.probes.information(WBillu)
        
translation <- read.affy.translation()
WBaffy <- read.csv("GPL570_WholeBlood.txt", sep='\t', row.names=1)
WBaffy <- annotate.affy.by.rownames(WBaffy, translation)
        
genes <- unique(as.character(WBillu[,1]))
ilGene <- as.character(WBillu[,1])
afGene <- as.character(WBaffy[,1])
MeanMatrix <- NULL
        cat("",file="WB.Affy.txt")
        cat("",file="WB.Illu.txt")
cnt <- 1
for(gene in genes){
    illuID <- which(ilGene == gene)
    affyID <- which(afGene == gene)
    if(length(illuID) > 0 && length(affyID) > 0){
        MeanMatrix <- rbind(MeanMatrix, c(gene, mean(unlist(WBaffy[affyID,-1])), mean(unlist(WBillu[illuID, -1]))))
        for(id in affyID){ cat(unlist(c(gene, as.character(WBaffy[id,-1]),"\n")), sep='\t',file="WB.Affy.txt", append=TRUE) }
        for(id in illuID){ cat(unlist(c(gene, as.character(WBillu[id,-1]),"\n")), sep='\t',file="WB.Illu.txt", append=TRUE) }
    }
    if(cnt %% 100 == 0) cat(cnt,length(genes),"\n")
    cnt <- cnt+1
}