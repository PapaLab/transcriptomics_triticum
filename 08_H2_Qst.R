### BLUP estimation for phenotypic traits across multiple environments
### Variance components of G, E, and GxE
### Broad-sense heritability estimation across multiple environments
### Contact: tong@mpimp-golm.mpg.de

#######################################################################
### set path and packages
#######################################################################
### set path
dir <- ("~/Nextcloud/10_co-work/wheat/2_H2_Qst/")
setwd(dir)
rm(list=ls())

### load library
library(limma)
library(edgeR)
library(DESeq2)
library(lme4)

#######################################################################
### read data set 
#######################################################################
### raw count
rawdata_filter <- read.table("../1_Diff/rawcount_filter_all.csv",head=T,sep=",")
### other information 
### species, genotype, condition
infodata <- read.table("../0_data/exp_design.csv",head=T,sep=",")
speciesinfo <- infodata[,3]
genotypeinfo <- infodata[,4]
conditioninfo <- infodata[,5]

genotypeinfo <- gsub("-", ".", genotypeinfo)
genotypeinfo <- gsub("/", ".", genotypeinfo)

#sum(colnames(rawdata) != infodata[,1])

#######################################################################
### normalization
#######################################################################

###############################################
### normalization of each species and condition

speciesname <- c("T.dicoccoides","T.dicoccum","T.durum")
conditionname <- c("High","Low")

for (i in c(1:3)) {
  for (j in c(1:2)){
    ssi <- speciesname[i]
    ssj <- conditionname[j]
    accid <- intersect(which(speciesinfo==ssi),which(conditioninfo==ssj))
    
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
    design <- model.matrix(~0+genotype)
    colnames(design) <- gsub("genotype", "", colnames(design))
    rownames(design) <- colnames(counts)
    
    ### 
    #dge <- estimateDisp(dge, design, robust=TRUE)
    
    ### voom normalization
    v <- voom(dge,design,plot=F) 
    
    #E	is numeric matrix of normalized expression values on the log2 scale
    ### save output
    write.table(v$E,paste0("TMM_voom_logE_",ssi,"_",ssj,".csv"),sep=',',row.names=T,col.names=T,quote=F)
    
  }
}

#######################################################################
### linar mixed model
#######################################################################

###############################################
### input data format: row for each line, column for each trait

speciesname <- c("T.dicoccoides","T.dicoccum","T.durum")
conditionname <- c("High","Low")

dataall <- NULL
for (i in c(1:3)) {
  for (j in c(1:2)){
    ssi <- speciesname[i]
    ssj <- conditionname[j]
    dataij <- read.table(paste0("TMM_voom_logE_",ssi,"_",ssj,".csv"),head=T,sep=",")
    dataall <- cbind(dataall,as.matrix(dataij))
  }   
}

dataall_f <- dataall[,infodata[,1]]
dataall_f <- t(dataall_f)

p <- ncol(dataall_f) #phenotype number
pname <- colnames(dataall_f)

xa <- cbind(speciesinfo,genotypeinfo,conditioninfo,dataall_f)
colnames(xa) <- c("species","genotype","condition",pname)
write.table(xa,"expression-data-all.csv",sep=',',quote=F,row.names=F)

###############################################
### estimate variance components and H2
###############################################
### run on Calculon !!!!

### read data
xa <- read.table("expression-data-all.csv",sep=",",header=T)

#line.blup <- c(1:length(unique(genotypeinfo)))
heritability <- NULL

for(i in 4:ncol(xa)){
  
  #print(i-3)
  
  control=lmerControl(check.conv.singular = .makeCC(action = "ignore",  tol = 1e-4),
                      check.nobs.vs.nlev = "ignore",
                      check.nobs.vs.rankZ = "ignore",
                      check.nobs.vs.nRE="ignore")
  # varcomp <- lmer(x[,i]~(1|LINE)+(1|LOC)+(1|LINE:LOC),data=x,control=control) 
  # #isSingular(varcomp, tol = 1e-4)
  
  varcomp <- lmer(xa[,i]~(1|condition)+(1|species:condition)+(1|species:genotype:condition)+(1|species/genotype),
                  #+(1|species)+(1|species/genotype),
                  data=xa,  control=control) 

  #summary(varcomp)
  var.trans <- lme4::VarCorr(varcomp)
  species.var <- as.numeric(var.trans$species)
  genotype.var <- as.numeric(var.trans$`genotype:species`)
  condition.var <- as.numeric(var.trans$condition)
  species.condition.var <- as.numeric(var.trans$`species:condition`)
  genotype.condition.var <- as.numeric(var.trans$`species:genotype:condition`)
  # residual standard deviation is stored as attribute "sc"
  e.var <- attr(var.trans,'sc')^2

  h2 <- c(species.var,genotype.var,condition.var,species.condition.var,genotype.condition.var,e.var,
          (species.var+genotype.var)/(species.var+genotype.var+condition.var+(species.condition.var/2)+(genotype.condition.var/2)+(e.var/2)))
  heritability <- rbind(heritability,h2)
  
  #f <- fixef(varcomp)
  #r <- ranef(varcomp)
  #blup <- rep(f,128)+c(rep(r$species[[1]][1],36),rep(r$species[[1]][2],43),rep(r$species[[1]][3],49))+r$`genotype:species`
  #line.blup <- cbind(line.blup,blup)
}

#colnames(line.blup) <- c('line',colnames(x)[idd+2])
heritability <- cbind(colnames(xa)[-c(1:3)],heritability)
colnames(heritability) <- c("gene","species.var","genotype.var","condition.var","species.condition.var","genotype.condition.var","residual","h2")

#write.table(line.blup,"all-gene-BLUP.csv",row.names=F,sep=",")
write.table(heritability,"all-gene-h2.csv",row.names=F,col.names=T,sep=",")

#######################################################################
### selective sweep - Qst
#######################################################################

### read data
xa <- read.table("expression-data-all.csv",sep=",",header=T)

### This is calculated of each condition separately
### For each condition, the three species are compared pairwise

###############################################
### for high / low N condition
###############################################

conditionname <- c("High","Low")

for (n in 1:2){
  
  nn <- conditionname[n]
  xacon <- xa[which(xa[,3]==nn),]
  
  print(nn)
  
  speciesname1 <- c("T.dicoccoides","T.dicoccum","T.durum")
  speciesname2 <- c("T.dicoccum","T.durum","T.dicoccoides")
  
  for (s in c(1:3)) {
      ss1 <- speciesname1[s]
      ss2 <- speciesname2[s]
      accid1 <- which(xacon[,1]==ss1)
      accid2 <- which(xacon[,1]==ss2)
      xaconi <- xacon[c(accid1,accid2),]
      
      print(s)
      
      #line.blup <- c(1:length(unique(genotypeinfo)))
      Qst <- NULL
      
      for(i in 4:ncol(xaconi)){
        
        #print(i-3)
        
        control=lmerControl(check.conv.singular = .makeCC(action = "ignore",  tol = 1e-4),
                            check.nobs.vs.nlev = "ignore",
                            check.nobs.vs.rankZ = "ignore",
                            check.nobs.vs.nRE="ignore")
        # varcomp <- lmer(x[,i]~(1|LINE)+(1|LOC)+(1|LINE:LOC),data=x,control=control) 
        # #isSingular(varcomp, tol = 1e-4)
        
        varcomp <- lmer(xaconi[,i]~1+(1|species/genotype),
                        data=xaconi, control=control) 
        
        #summary(varcomp)
        var.trans <- lme4::VarCorr(varcomp)
        species.var <- as.numeric(var.trans$species)
        genotype.var <- as.numeric(var.trans$`genotype:species`)
        # residual standard deviation is stored as attribute "sc"
        e.var <- attr(var.trans,'sc')^2
        
        qst <- c(species.var,genotype.var,e.var,
                (species.var)/(species.var+genotype.var))
        Qst <- rbind(Qst,qst)
        
        #f <- fixef(varcomp)
        #r <- ranef(varcomp)
        #blup <- rep(f,128)+c(rep(r$species[[1]][1],36),rep(r$species[[1]][2],43),rep(r$species[[1]][3],49))+r$`genotype:species`
        #line.blup <- cbind(line.blup,blup)
      }
      
      #colnames(line.blup) <- c('line',colnames(x)[idd+2])
      Qst <- cbind(colnames(xa)[-c(1:3)],Qst)
      colnames(Qst) <- c("gene","species.var","genotype.var","residual","Qst")
      
      #write.table(line.blup,"all-gene-BLUP.csv",row.names=F,sep=",")
      write.table(Qst,paste0("all-gene-Qst_",nn,"_",ss1,"_",ss2,".csv"),row.names=F,col.names=T,sep=",")
  
  }

}

#######################################################################
### the end
#######################################################################
