#
# create.tstat.matrix.illumina.R
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

WB <- read.illumina.wb()
WBAnnot <- add.illumina.probes.information(WB)

IlluWBHigh <- cor(WBAnnot[,-1])
IlluWBCor <- names(which(apply(IlluWBHigh, 1, mean, na.rm=T) > 0.85))

WBAnnot <- WBAnnot[,IlluWBCor]

WBAnnot <- WBAnnot[which(rownames(WBAnnot) %in% rownames(Illu)),]
Illu <- Illu[which(rownames(Illu) %in% rownames(WBAnnot)),]

IlluSort <- match(rownames(WBAnnot), rownames(Illu))
Illu <- Illu[IlluSort,]

cellTypes <- unique(colnames(Illu))
probes <- rownames(Illu)
if(!file.exists("tstat.matrix.illumina.txt")){
  cat("",cellTypes, file="tstat.matrix.illumina.txt",sep='\t')
  cat("\n", file="tstat.matrix.illumina.txt", append=TRUE)
  signLVL <- 0.05/nrow(Illu)
  for(p in 1:nrow(Illu)){
    data <- NULL  
    for(celltype in cellTypes){
      cols <-  which(colnames(Illu)==celltype)
      tC  <- t.test(unlist(Illu[p,cols]), log2(unlist(WBAnnot[p,-1])))
      sC  <- tC$statistic
#      if(tC$p.value > signLVL) sC  <- NA
      data <- cbind(data, sC)
    }
    if(p %% 100 == 0) cat("Done",p,"Out of",nrow(Illu),"\n")
    cat(probes[p], data, file="tstat.matrix.illumina.txt", append=TRUE,sep='\t')
    cat("\n", file="tstat.matrix.illumina.txt", append=TRUE)
  }
}else{
  full <- read.csv("tstat.matrix.illumina.txt", sep='\t', row.names = 1)
}

onlyonce <- which(apply(full,1,function(x){ sum(!is.na(x)) })==1)
onceCell <- full[onlyonce,]

