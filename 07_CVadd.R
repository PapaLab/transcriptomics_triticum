### load library
library(DESeq2)

### design information 
### species, genotype, condition
infodata <- read.table("./0_data/exp_design.csv",head=T,sep=",")
speciesinfo <- infodata[,3]
genotypeinfo <- infodata[,4]
conditioninfo <- infodata[,5]

genotypeinfo <- gsub("-", ".", genotypeinfo)
genotypeinfo <- gsub("/", ".", genotypeinfo)

### filter genes based on CPM
### only the expressed genes remained 
### keep the genes if at least 10 samples of CPM > 1
### This is calculated on raw counts of each species separately!!!

### read data set 
rawdata_filter <- read.table("./1_Diff/rawcount_filter_all.csv",head=T,sep=",")

###create dds
### This is calculated of each condition separately
conditionname <- c("High","Low")

###########
#in HIGH N#
for (i in c(1)) {
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
  
  dds_high <- DESeqDataSetFromMatrix(countData=counts, colData=coldata, design=design)
  # The factor of interest should go at the end of the formula
  
}

###transform
vsd_high <- vst(dds_high, blind=FALSE)


###Function to calculate the COEFFICIENT OF VARIATION
co.var <- function(x) ( sd(x)/mean(x) )

df_vsd_high <- as.data.frame(assay(vsd_high))


dicoccoides_vsd_high <- df_vsd_high[,c(1:20)]
dicoccum_vsd_high <- df_vsd_high[,c(21:41)]
durum_vsd_high <- df_vsd_high[,c(42:66)]

####DICOCCOIDES####

dicoccoides_vsd_high$BAZ382 <- apply(dicoccoides_vsd_high[,1:2] , 1, mean)               
dicoccoides_vsd_high$PI233288 <- apply(dicoccoides_vsd_high[,3:4] , 1, mean)               
dicoccoides_vsd_high$PI428014 <- apply(dicoccoides_vsd_high[,5:6] , 1, mean)               
dicoccoides_vsd_high$PI466991 <- apply(dicoccoides_vsd_high[,7:8] , 1, mean)               
dicoccoides_vsd_high$IG46504 <- apply(dicoccoides_vsd_high[,9:10] , 1, mean)               
dicoccoides_vsd_high$PI538656 <- apply(dicoccoides_vsd_high[,11:12] , 1, mean)               
dicoccoides_vsd_high$PI554584 <- apply(dicoccoides_vsd_high[,13:14] , 1, mean)               
dicoccoides_vsd_high$PI343446 <- apply(dicoccoides_vsd_high[,15:16] , 1, mean)               
dicoccoides_vsd_high$PI352324 <- apply(dicoccoides_vsd_high[,17:18] , 1, mean)               
dicoccoides_vsd_high$PI355459 <- apply(dicoccoides_vsd_high[,19:20] , 1, mean)               

dicoccoides_vsd_high$CV_add <- apply(dicoccoides_vsd_high[,22:31] , 1, co.var)               
mean(dicoccoides_vsd_high$CV_add)


CV_dicoccoides <- as.data.frame(dicoccoides_vsd_high[,c(21,32)])
CV_dicoccoides$gene <- row.names(dicoccoides_vsd_high)

####DICOCCUM####

dicoccum_vsd_high$MG5473 <- apply(dicoccum_vsd_high[,1:2] , 1, mean)               
dicoccum_vsd_high$Lucanica <- apply(dicoccum_vsd_high[,3:4] , 1, mean)               
dicoccum_vsd_high$MG5350 <- apply(dicoccum_vsd_high[,5:6] , 1, mean)               
dicoccum_vsd_high$MG5293.1 <- dicoccum_vsd_high[,7]
dicoccum_vsd_high$PI74106 <- apply(dicoccum_vsd_high[,8:9] , 1, mean)               
dicoccum_vsd_high$CItr17675 <- apply(dicoccum_vsd_high[,10:11] , 1, mean)               
dicoccum_vsd_high$PI254169 <- apply(dicoccum_vsd_high[,12:13] , 1, mean)               
dicoccum_vsd_high$PI352347 <- apply(dicoccum_vsd_high[,14:15] , 1, mean)               
dicoccum_vsd_high$PI470739 <- apply(dicoccum_vsd_high[,16:17] , 1, mean)               
dicoccum_vsd_high$MoliseSel.Colli <- apply(dicoccum_vsd_high[,18:21] , 1, mean)               

dicoccum_vsd_high$CV_add <- apply(dicoccum_vsd_high[,23:32] , 1, co.var)               
mean(dicoccum_vsd_high$CV_add)

CV_dicoccum <- as.data.frame(dicoccum_vsd_high[,c(22,33)])
CV_dicoccum$gene <- row.names(dicoccum_vsd_high)



####DURUM####

durum_vsd_high$Appulo <- apply(durum_vsd_high[,1:2] , 1, mean)               
durum_vsd_high$Capeiti.8 <- apply(durum_vsd_high[,3:4] , 1, mean)               
durum_vsd_high$Cappelli <- apply(durum_vsd_high[,5:6] , 1, mean)               
durum_vsd_high$Cirillo <- apply(durum_vsd_high[,7:8] , 1, mean)
durum_vsd_high$Creso <- apply(durum_vsd_high[,9:10] , 1, mean)               
durum_vsd_high$Neodur <- apply(durum_vsd_high[,11:12] , 1, mean)               
durum_vsd_high$Ofanto <- apply(durum_vsd_high[,13:14] , 1, mean)               
durum_vsd_high$Pedroso <- durum_vsd_high[,15]            
durum_vsd_high$PR22D89 <- apply(durum_vsd_high[,16:17] , 1, mean)               
durum_vsd_high$Timilia <- apply(durum_vsd_high[,18:19] , 1, mean)  
durum_vsd_high$Trinakria <- apply(durum_vsd_high[,20:21] , 1, mean)
durum_vsd_high$Simeto <- apply(durum_vsd_high[,22:25] , 1, mean)

durum_vsd_high$CV_add <- apply(durum_vsd_high[,27:38] , 1, co.var)               
mean(durum_vsd_high$CV_add)

CV_durum <- as.data.frame(durum_vsd_high[,c(26,39)])
CV_durum$gene <- row.names(durum_vsd_high)



###########
#in low N#
for (i in c(2)) {
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
  
  dds_low <- DESeqDataSetFromMatrix(countData=counts, colData=coldata, design=design)
  # The factor of interest should go at the end of the formula
  
}

###transform
vsd_low <- vst(dds_low, blind=FALSE)

###Function to calculate the COEFFICIENT OF VARIATION
co.var <- function(x) ( sd(x)/mean(x) )

df_vsd_low <- as.data.frame(assay(vsd_low))

dicoccoides_vsd_low <- df_vsd_low[,1:16]
dicoccum_vsd_low <- df_vsd_low[,17:38]
durum_vsd_low <- df_vsd_low[,39:62]

####DICOCCOIDES####


dicoccoides_vsd_low$BAZ382 <- apply(dicoccoides_vsd_low[,c("Ti55","Ti90")] , 1, mean)               
dicoccoides_vsd_low$PI233288 <- dicoccoides_vsd_low[,c("Ti13")]              
dicoccoides_vsd_low$PI428014 <- apply(dicoccoides_vsd_low[,c("Ti38","Ti59")] , 1, mean)               
dicoccoides_vsd_low$PI466991 <- dicoccoides_vsd_low[,c("Ti5")]              
dicoccoides_vsd_low$IG46504 <- dicoccoides_vsd_low[,c("Ti89")]              
dicoccoides_vsd_low$PI538656 <- dicoccoides_vsd_low[,c("Ti95")]              
dicoccoides_vsd_low$PI554584 <- apply(dicoccoides_vsd_low[,c("Ti42","Ti94")] , 1, mean)               
dicoccoides_vsd_low$PI343446 <- apply(dicoccoides_vsd_low[,c("Ti126","Ti143")] , 1, mean)               
dicoccoides_vsd_low$PI352324 <- apply(dicoccoides_vsd_low[,c("Ti35","Ti88")] , 1, mean)               
dicoccoides_vsd_low$PI355459 <- apply(dicoccoides_vsd_low[,c("Ti49","Ti96")] , 1, mean)

dicoccoides_vsd_low$CV_add <- apply(dicoccoides_vsd_low[,18:27] , 1, co.var)               
mean(dicoccoides_vsd_low$CV_add)

CV_dicoccoides <- as.data.frame(dicoccoides_vsd_low[,c(17,28)])
CV_dicoccoides$gene <- row.names(dicoccoides_vsd_low)


####DICOCCUM####

dicoccum_vsd_low$MG5473 <- apply(dicoccum_vsd_low[,c("Ti12","Ti123")] , 1, mean)               
dicoccum_vsd_low$Lucanica <- apply(dicoccum_vsd_low[,c("Ti7","Ti132")] , 1, mean)               
dicoccum_vsd_low$MG5350 <- apply(dicoccum_vsd_low[,c("Ti56","Ti129")] , 1, mean)               
dicoccum_vsd_low$MG5293.1 <- apply(dicoccum_vsd_low[,c("Ti3","Ti140")] , 1, mean)
dicoccum_vsd_low$PI74106 <- apply(dicoccum_vsd_low[,c("Ti66","Ti139")] , 1, mean)               
dicoccum_vsd_low$CItr17675 <- apply(dicoccum_vsd_low[,c("Ti20","Ti92")] , 1, mean)               
dicoccum_vsd_low$PI254169 <- apply(dicoccum_vsd_low[,c("Ti63","Ti79")] , 1, mean)               
dicoccum_vsd_low$PI352347 <- apply(dicoccum_vsd_low[,c("Ti26","Ti76")] , 1, mean)               
dicoccum_vsd_low$PI470739 <- apply(dicoccum_vsd_low[,c("Ti64","Ti111")] , 1, mean)               
dicoccum_vsd_low$MoliseSel.Colli <- apply(dicoccum_vsd_low[,c("Ti91","Ti121","Ti278","Ti239")] , 1, mean) 

dicoccum_vsd_low$CV_add <- apply(dicoccum_vsd_low[,24:33] , 1, co.var)               
mean(dicoccum_vsd_low$CV_add)

CV_dicoccum <- as.data.frame(dicoccum_vsd_low[,c(23,34)])
CV_dicoccum$gene <- row.names(dicoccum_vsd_low)

####DURUM####

durum_vsd_low$Appulo <- apply(durum_vsd_low[,1:2] , 1, mean)               
durum_vsd_low$Capeiti.8 <- apply(durum_vsd_low[,3:4] , 1, mean)               
durum_vsd_low$Cappelli <- apply(durum_vsd_low[,5:6] , 1, mean)               
durum_vsd_low$Cirillo <- apply(durum_vsd_low[,7:8] , 1, mean)
durum_vsd_low$Creso <- durum_vsd_low[,9]     
durum_vsd_low$Neodur <- apply(durum_vsd_low[,10:11] , 1, mean)               
durum_vsd_low$Ofanto <- apply(durum_vsd_low[,12:13] , 1, mean)               
durum_vsd_low$Pedroso <- apply(durum_vsd_low[,14:15] , 1, mean)               
durum_vsd_low$PR22D89 <- apply(durum_vsd_low[,16:17] , 1, mean)               
durum_vsd_low$Timilia <- apply(durum_vsd_low[,18:19] , 1, mean)  
durum_vsd_low$Trinakria <- durum_vsd_low[,20]
durum_vsd_low$Simeto <- apply(durum_vsd_low[,21:24] , 1, mean)

durum_vsd_low$CV_add <- apply(durum_vsd_low[,26:37] , 1, co.var)               
mean(durum_vsd_low$CV_add)

CV_durum <- as.data.frame(durum_vsd_low[,c(25,38)])
CV_durum$gene <- row.names(durum_vsd_low)