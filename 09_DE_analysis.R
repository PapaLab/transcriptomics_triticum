### Find the differential expression gene (DEG) 
### Input: the raw count of each transcripts
### Three DEG detection approach applied: linear model based implemented in the R package limma, 
### and two Poisson-based models implemented in the R packages edgeR and DESeq2. 
### Output: DEG in three approaches and the intersection. 
### Contact: tong@mpimp-golm.mpg.de

#######################################################################
### set path and packages
#######################################################################
### set path
dir <- ("~/Nextcloud/10_co-work/wheat/1_Diff/")
setwd(dir)
rm(list=ls())

### load library
library(limma)
library(edgeR)
library(DESeq2)

#######################################################################
### read data set 
#######################################################################
### design information 
### species, genotype, condition
infodata <- read.table("../0_data/exp_design.csv",head=T,sep=",")
speciesinfo <- infodata[,3]
genotypeinfo <- infodata[,4]
conditioninfo <- infodata[,5]

genotypeinfo <- gsub("-", ".", genotypeinfo)
genotypeinfo <- gsub("/", ".", genotypeinfo)

### raw count
rawdata <- read.table("../0_data/all_samples_counts_matrix",head=T,sep="\t")

#sum(colnames(rawdata) != infodata[,1])

#######################################################################
### filter genes based on CPM
#######################################################################
### only the expressed genes remained 
### keep the genes if at least 10 samples of CPM > 1

### This is calculated on raw counts of each species separately!!!
speciesname <- c("T.dicoccoides","T.dicoccum","T.durum")

for (i in c(1:3)) {
  ss <- speciesname[i]
  datai <- rawdata[,which(speciesinfo==ss)]
  
  dataicpm <- cpm(datai)
  dataif <- datai[rowSums(dataicpm>1) >= 10,]
  idi <- as.numeric(rowSums(dataicpm>1) >= 10)
  idi <- which(idi==1)
  
  write.table(dataif,paste0("rawcount_filter_",ss,".csv"),sep=',',row.names=T,col.names=T,quote=F)
  write.table(idi,paste0("id_filter_",ss,".csv"),sep=',',row.names=F,col.names=F,quote=F)

}

######################
### combine the three by intersection
id1 <- read.table("id_filter_T.dicoccoides.csv",head=F,sep=",")
id2 <- read.table("id_filter_T.dicoccum.csv",head=F,sep=",")
id3 <- read.table("id_filter_T.durum.csv",head=F,sep=",")
idall <- Reduce(intersect, list(as.matrix(id1),as.matrix(id2),as.matrix(id3)))
rawdata_filter <- rawdata[idall,]

write.table(rawdata_filter,"rawcount_filter_all.csv",sep=',',row.names=T,col.names=T,quote=F)
write.table(cbind(idall,rownames(rawdata_filter)),"id_filter_all.csv",sep=',',row.names=F,col.names=F,quote=F)

#######################################################################
### compare between conditions of each species
#######################################################################

rawdata_filter <- read.table("rawcount_filter_all.csv",head=T,sep=",")

###############################################
### limma
###############################################
### This is calculated of each species separately
speciesname <- c("T.dicoccoides","T.dicoccum","T.durum")

for (i in c(1:3)) {
  ss <- speciesname[i]
  accid <- which(speciesinfo==ss)
  
  ### make the list
  counts <- rawdata_filter[,accid]
  dge <- DGEList(counts=counts,genes=rownames(counts))
  group <- as.factor(conditioninfo[accid])
  dge$samples$group <- group
  species <- as.factor(speciesinfo[accid])
  dge$samples$species <- species
  genotype <- as.factor(genotypeinfo[accid])
  dge$samples$genotype <- genotype
  
  ### normalization by trimmed mean of M-values (TMM) 
  dge <- calcNormFactors(dge)
  
  ### set design matrix 
  design <- model.matrix(~0+group+genotype)
  colnames(design) <- gsub("group", "", colnames(design))
  rownames(design) <- colnames(counts)
  
  ### set contrasts
  mc=makeContrasts('High-Low',levels=design)
  
  ### dge <- estimateDisp(dge, design, robust=TRUE)
  
  ### voom normalization
  v <- voom(dge,design,plot=F) 
  
  ### find DEG
  fit <- lmFit(v,design)
  c.fit=contrasts.fit(fit,mc)
  # calculate P-values with moderate t-statistic with eBayes
  eb=eBayes(c.fit)
  
  ### save ouput
  # adjust p-value
  pres <- apply(eb$p.value,2,p.adjust,method="fdr")
  #adjpres <- apply(pres,2,function(col){names(col)[which(col<0.05)]})
  
  degout <- data.frame(LFC = eb$coefficients, 
                     p.value= eb$p.value, adj.pvalue = pres)
  colnames(degout) <- c("LFC","p.value","adj.p.value")
  rownames(degout) <- rownames(counts)
  
  write.table(degout,paste0("DEG_limma_between.condition_in.",ss,".csv"),sep=',',row.names=T,col.names=T,quote=F)
  
}

###############################################
### edgeR
###############################################
### This is calculated of each species separately
speciesname <- c("T.dicoccoides","T.dicoccum","T.durum")

for (i in c(1:3)) {
  ss <- speciesname[i]
  accid <- which(speciesinfo==ss)
  
  ### make the list
  counts <- rawdata_filter[,accid]
  dge <- DGEList(counts=counts,genes=rownames(counts))
  group <- as.factor(conditioninfo[accid])
  dge$samples$group <- group
  species <- as.factor(speciesinfo[accid])
  dge$samples$species <- species
  genotype <- as.factor(genotypeinfo[accid])
  dge$samples$genotype <- genotype
  
  ### normalization by trimmed mean of M-values (TMM) 
  dge <- calcNormFactors(dge)
  
  ### set design matrix 
  design <- model.matrix(~0+group+genotype)
  colnames(design) <- gsub("group", "", colnames(design))
  rownames(design) <- colnames(counts)
  
  ### set contrasts
  mc=makeContrasts('High-Low',levels=design)
  
  ### estimate dispersions
  dge <- estimateGLMCommonDisp(dge, design=design)
  dge <- estimateGLMTrendedDisp(dge, design=design)
  dge <- estimateGLMTagwiseDisp(dge, design=design)
  
  ### find DEG
  fit <- glmFit(dge, design)
  lrt <- glmLRT(fit, contrast=mc)
  
  ### save output
  # adjust p-value
  pres <- apply(matrix(lrt$table$PValue,ncol=1),2,p.adjust,method="fdr")
  #adjpres <- apply(pres,2,function(col){names(col)[which(col<0.05)]})
  
  degout <- data.frame(LFC = lrt$table$logFC, 
                       p.value= lrt$table$PValue, adj.pvalue = pres)
  colnames(degout) <- c("LFC","p.value","adj.p.value")
  rownames(degout) <- rownames(counts)
  
  write.table(degout,paste0("DEG_edgeR_between.condition_in.",ss,".csv"),sep=',',row.names=T,col.names=T,quote=F)
  
}

###############################################
### DESeq2
###############################################
### This is calculated of each species separately
speciesname <- c("T.dicoccoides","T.dicoccum","T.durum")

for (i in c(1:3)) {
  ss <- speciesname[i]
  accid <- which(speciesinfo==ss)
  
  ### set up data
  counts <- rawdata_filter[,accid]
  group <- as.factor(conditioninfo[accid])
  species <- as.factor(speciesinfo[accid])
  genotype <- as.factor(genotypeinfo[accid])
  coldata <- data.frame(group,species,genotype)
  
  ### set design matrix 
  #design <- model.matrix(~0+group+genotype)
  #colnames(design) <- gsub("group", "", colnames(design))
  #colnames(design) <- gsub("genotype", "", colnames(design))
  #rownames(design) <- colnames(counts)
  #dds <- DESeqDataSetFromMatrix(countData=counts, colData=coldata, design=design)
  
  dds <- DESeqDataSetFromMatrix(countData=counts, colData=coldata, design=~group+genotype)
  # The factor of interest should go at the end of the formula
  
  ### set control condition as reference
  #dds$group <- relevel(dds$group, ref = "High")
  
  ####cds <- estimateSizeFactors(dds)
  #cdsB <-estimateDispersions(cds)
  
  ### find DEG
  resdds <- DESeq(dds)
  res <- results(resdds, contrast=c("group","High","Low"))
  #res <- results(resdds, contrast=list("High","Low"))
  
  ### save output
  # adjust p-value
  #pres <- apply(matrix(res$pvalue,ncol=1),2,p.adjust,method="fdr")
  #adjpres <- apply(pres,2,function(col){names(col)[which(col<0.05)]})
  
  degout <- data.frame(LFC=res$log2FoldChange, 
                       p.value=res$pvalue, adj.p.value=res$padj)
  colnames(degout) <- c("LFC","p.value","adj.p.value")
  rownames(degout) <- rownames(counts)
  
  write.table(degout,paste0("DEG_DESeq2_between.condition_in.",ss,".csv"),sep=',',row.names=T,col.names=T,quote=F)
  
}

###############################################
### overlap between the three
###############################################
### set the significant level as adj.p.vale < alpha
### report the intersection of significant DEG in three approaches

pp <- 0.001
speciesname <- c("T.dicoccoides","T.dicoccum","T.durum")

for (i in c(1:3)) {
  ss <- speciesname[i]
  results_limma <- read.table(paste0("DEG_limma_between.condition_in.",ss,".csv"),head=T,sep=",")
  results_edgeR <- read.table(paste0("DEG_edgeR_between.condition_in.",ss,".csv"),head=T,sep=",")
  results_DEseq <- read.table(paste0("DEG_DESeq2_between.condition_in.",ss,".csv"),head=T,sep=",")
  
  ### significant genes
  sig_limma <- results_limma[which(results_limma[,"adj.p.value"]<pp),]
  sig_edgeR <- results_edgeR[which(results_edgeR[,"adj.p.value"]<pp),]
  sig_DESeq <- results_DEseq[which(results_DEseq[,"adj.p.value"]<pp),]
  
  write.table(sig_limma,paste0("DEG_limma_between.condition_in.",ss,"_adjP.",pp,".csv"),sep=',',row.names=T,col.names=T,quote=F)
  write.table(sig_edgeR,paste0("DEG_edgeR_between.condition_in.",ss,"_adjP.",pp,".csv"),sep=',',row.names=T,col.names=T,quote=F)
  write.table(sig_DESeq,paste0("DEG_DESeq2_between.condition_in.",ss,"_adjP.",pp,".csv"),sep=',',row.names=T,col.names=T,quote=F)
  
  ### overlap
  #sig_all <- Reduce(intersect, list(rownames(sig_limma),rownames(sig_edgeR)))
  sig_all <- Reduce(intersect, list(rownames(sig_limma),rownames(sig_edgeR),rownames(sig_DESeq)))
  write.table(sig_all,paste0("DEG_overlap_between.condition_in.",ss,"_adjP.",pp,".csv"),sep=',',row.names=F,col.names=F,quote=F)
  
}

#######################################################################
### compare between species of each condition
#######################################################################

rawdata_filter <- read.table("rawcount_filter_all.csv",head=T,sep=",")

###############################################
### limma
###############################################
### This is calculated of each condition separately
conditionname <- c("High","Low")

for (i in c(1:2)) {
  ss <- conditionname[i]
  accid <- which(conditioninfo==ss)
  
  ### make the list
  counts <- rawdata_filter[,accid]
  dge <- DGEList(counts=counts,genes=rownames(counts))
  group <- as.factor(conditioninfo[accid])
  dge$samples$group <- group
  species <- as.factor(speciesinfo[accid])
  dge$samples$species <- species
  genotype <- as.factor(genotypeinfo[accid])
  dge$samples$genotype <- genotype
  
  ### normalization by trimmed mean of M-values (TMM) 
  dge <- calcNormFactors(dge)
  
  ### set design matrix 
  # correct the nested issue
  dge$samples$species.n <- factor(as.numeric(dge$samples$species))
  genotype.n <- NULL
  for (s in 1:3){
    idd <- which(dge$samples$species.n==s)
    g.n <- as.character(dge$samples$genotype[idd])
    g.n1 <- as.numeric(as.factor(g.n))
    genotype.n <- c(genotype.n, g.n1)
  }
  dge$samples$genotype.n <- factor(genotype.n)
  
  design <- model.matrix(~0+species+species:genotype.n)
  colnames(design) <- gsub("species", "", colnames(design))
  colnames(design) <- gsub(":", ".", colnames(design))
  rownames(design) <- colnames(counts)
  
  ### set contrasts
  mc=makeContrasts('T.dicoccoides-T.dicoccum',
                   'T.dicoccoides-T.durum',
                   'T.dicoccum-T.durum',
                   levels=design)
  
  ####dge <- estimateDisp(dge, design, robust=TRUE)
  
  ### voom normalization
  v <- voom(dge,design,plot=F) 
  
  ### find DEG
  fit <- lmFit(v,design)
  c.fit=contrasts.fit(fit,mc)
  # calculate P-values with moderate t-statistic with eBayes
  eb=eBayes(c.fit)
  
  ### save ouput
  # adjust p-value
  pres <- apply(eb$p.value,2,p.adjust,method="fdr")
  #adjpres <- apply(pres,2,function(col){names(col)[which(col<0.05)]})
  
  degout <- data.frame(LFC = eb$coefficients, 
                       p.value= eb$p.value, adj.p.value = pres)
  #colnames(degout) <- c("LFC","p.value","adj.p.value")
  rownames(degout) <- rownames(counts)
  
  write.table(degout,paste0("DEG_limma_between.species_in.",ss,".csv"),sep=',',row.names=T,col.names=T,quote=F)
  
}

###############################################
### edgeR
###############################################
### This is calculated of each condition separately
conditionname <- c("High","Low")

for (i in c(1:2)) {
  ss <- conditionname[i]
  accid <- which(conditioninfo==ss)
  
  ### make the list
  counts <- rawdata_filter[,accid]
  dge <- DGEList(counts=counts,genes=rownames(counts))
  group <- as.factor(conditioninfo[accid])
  dge$samples$group <- group
  species <- as.factor(speciesinfo[accid])
  dge$samples$species <- species
  genotype <- as.factor(genotypeinfo[accid])
  dge$samples$genotype <- genotype
  
  ### normalization by trimmed mean of M-values (TMM) 
  dge <- calcNormFactors(dge)
  
  ### set design matrix 
  # correct the nested issue
  dge$samples$species.n <- factor(as.numeric(dge$samples$species))
  genotype.n <- NULL
  for (s in 1:3){
    idd <- which(dge$samples$species.n==s)
    g.n <- as.character(dge$samples$genotype[idd])
    g.n1 <- as.numeric(as.factor(g.n))
    genotype.n <- c(genotype.n, g.n1)
  }
  dge$samples$genotype.n <- factor(genotype.n)
  
  design <- model.matrix(~0+species+species:genotype.n)
  colnames(design) <- gsub("species", "", colnames(design))
  colnames(design) <- gsub(":", ".", colnames(design))
  rownames(design) <- colnames(counts)
  
  ### set contrasts
  mc1=makeContrasts('T.dicoccoides-T.dicoccum',levels=design)
  mc2=makeContrasts('T.dicoccoides-T.durum',levels=design)
  mc3=makeContrasts('T.dicoccum-T.durum',levels=design)
  
  ### estimate dispersions
  dge <- estimateGLMCommonDisp(dge, design=design)
  dge <- estimateGLMTrendedDisp(dge, design=design)
  dge <- estimateGLMTagwiseDisp(dge, design=design)
  
  ### find DEG
  fit <- glmFit(dge, design)
  lrt1 <- glmLRT(fit, contrast=mc1)
  lrt2 <- glmLRT(fit, contrast=mc2)
  lrt3 <- glmLRT(fit, contrast=mc3)

  ### save ouput
  # adjust p-value
  pres1 <- apply(matrix(lrt1$table$PValue,ncol=1),2,p.adjust,method="fdr")
  pres2 <- apply(matrix(lrt2$table$PValue,ncol=1),2,p.adjust,method="fdr")
  pres3 <- apply(matrix(lrt3$table$PValue,ncol=1),2,p.adjust,method="fdr")
  comp <- c('T.dicoccoides-T.dicoccum','T.dicoccoides-T.durum','T.dicoccum-T.durum')
  
  ### output
  LFC <- cbind(lrt1$table$logFC,lrt2$table$logFC,lrt3$table$logFC)
  p.value <- cbind(lrt1$table$PValue,lrt2$table$PValue,lrt3$table$PValue)
  adj.p.value <- cbind(pres1,pres2,pres3)
  colnames(LFC) <- comp
  colnames(p.value) <- comp
  colnames(adj.p.value) <- comp
  
  degout <- data.frame(LFC=LFC,p.value=p.value,adj.p.value=adj.p.value)
  rownames(degout) <- rownames(counts)
  
  write.table(degout,paste0("DEG_edgeR_between.species_in.",ss,".csv"),sep=',',row.names=T,col.names=T,quote=F)
  
}

###############################################
### DESeq2
###############################################
### This is calculated of each condition separately
conditionname <- c("High","Low")

for (i in c(1:2)) {
  ss <- conditionname[i]
  accid <- which(conditioninfo==ss)
  
  ### set up data
  counts <- rawdata_filter[,accid]
  group <- as.factor(conditioninfo[accid])
  species <- as.factor(speciesinfo[accid])
  genotype <- as.factor(genotypeinfo[accid])
  coldata <- data.frame(group,species,genotype)
  
  ### set design matrix 
  # correct the nested issue
  coldata$species.n <- factor(as.numeric(coldata$species))
  genotype.n <- NULL
  for (s in 1:3){
    idd <- which(coldata$species.n==s)
    g.n <- as.character(coldata$genotype[idd])
    g.n1 <- as.numeric(as.factor(g.n))
    genotype.n <- c(genotype.n, g.n1)
  }
  coldata$genotype.n <- factor(genotype.n)
  
  design <- model.matrix(~0+species+species:genotype.n)
  colnames(design) <- gsub("species", "", colnames(design))
  colnames(design) <- gsub(":", ".", colnames(design))
  rownames(design) <- colnames(counts)
  
  dds <- DESeqDataSetFromMatrix(countData=counts, colData=coldata, design=design)
  # The factor of interest should go at the end of the formula
  
  ### set control condition as reference
  #dds$group <- relevel(dds$group, ref = "High")
  
  ####cds <- estimateSizeFactors(dds)
  #cdsB <-estimateDispersions(cds)
  
  ### find DEG
  resdds <- DESeq(dds)
  comp <- c('T.dicoccoides-T.dicoccum','T.dicoccoides-T.durum','T.dicoccum-T.durum')
  res1 <- results(resdds, contrast=list("T.dicoccoides","T.dicoccum"))
  res2 <- results(resdds, contrast=list("T.dicoccoides","T.durum"))
  res3 <- results(resdds, contrast=list("T.dicoccum","T.durum"))
  
  ### save output
  # adjust p-value
  #pres <- apply(matrix(res$pvalue,ncol=1),2,p.adjust,method="fdr")
  #adjpres <- apply(pres,2,function(col){names(col)[which(col<0.05)]})
  
  LFC <- cbind(res1$log2FoldChange,res2$log2FoldChange,res3$log2FoldChange)
  p.value <- cbind(res1$pvalue,res2$pvalue,res3$pvalue)
  adj.p.value <- cbind(res1$padj,res2$padj,res3$padj)
  colnames(LFC) <- comp
  colnames(p.value) <- comp
  colnames(adj.p.value) <- comp
  
  degout <- data.frame(LFC=LFC,p.value=p.value,adj.p.value=adj.p.value)
  rownames(degout) <- rownames(counts)
  
  write.table(degout,paste0("DEG_DESeq2_between.species_in.",ss,".csv"),sep=',',row.names=T,col.names=T,quote=F)
  
}

###############################################
### overlap between the three
###############################################
### set the significant level as adj.p.vale < alpha
### report the intersection of significant DEG in three approaches

pp <- 0.001
conditionname <- c("High","Low")
bets <- c('T.dicoccoides-T.dicoccum','T.dicoccoides-T.durum', 'T.dicoccum-T.durum')

for (i in c(1:2)) {
  ss <- conditionname[i]
  results_limma <- read.table(paste0("DEG_limma_between.species_in.",ss,".csv"),head=T,sep=",")
  results_edgeR <- read.table(paste0("DEG_edgeR_between.species_in.",ss,".csv"),head=T,sep=",")
  results_DEseq <- read.table(paste0("DEG_DESeq2_between.species_in.",ss,".csv"),head=T,sep=",")
  
  for (j in c(1:3)){
    ### significant genes
    sig_limma <- results_limma[which(results_limma[,j+6]<pp),]
    sig_edgeR <- results_edgeR[which(results_edgeR[,j+6]<pp),]
    sig_DESeq <- results_DEseq[which(results_DEseq[,j+6]<pp),]
    
    write.table(sig_limma,paste0("DEG_limma_between.",bets[j],"_in.",ss,"_adjP.",pp,".csv"),sep=',',row.names=T,col.names=T,quote=F)
    write.table(sig_edgeR,paste0("DEG_edgeR_between.",bets[j],"_in.",ss,"_adjP.",pp,".csv"),sep=',',row.names=T,col.names=T,quote=F)
    write.table(sig_DESeq,paste0("DEG_DESeq2_between.",bets[j],"_in.",ss,"_adjP.",pp,".csv"),sep=',',row.names=T,col.names=T,quote=F)
    
    ### overlap
    #sig_all <- Reduce(intersect, list(rownames(sig_limma),rownames(sig_edgeR)))
    sig_all <- Reduce(intersect, list(rownames(sig_limma),rownames(sig_edgeR),rownames(sig_DESeq)))
    write.table(sig_all,paste0("DEG_overlap_between.",bets[j],"_in.",ss,"_adjP.",pp,".csv"),sep=',',row.names=F,col.names=F,quote=F)
  }

}

#######################################################################
### the end
#######################################################################
