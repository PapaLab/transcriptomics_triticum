#this script was written to extract Fst estimates for each gene. The same script was run for each comparison (dicoccoides vs. dicoccum, dicoccum vs. durum, dicoccoides vs. durum)

fst = read.table("dicoccoides_dicoccum_FOLDED_fst.txt", header=F)
genes = read.table("HC_genes_coord.txt")

genes = read.table("HC_genes_coord.txt")

colnames(fst) <- c("chr", "pos", "A1", "B1")
colnames(genes) <- c ("chr", "start", "end", "name")

genes$Fst_dicoccoides_dicoccum = NA

#loop over genes and calculate fst

for (i in 1:nrow(genes))
{
  sites = (fst$chr == genes[i,1] & fst$pos >= genes[i,2] & fst$pos <= genes[i,3])
  
  genes$Fst_dicoccoides_dicoccum[i]  = sum(fst$A1[sites])/sum(fst$B1[sites])
  
}


write.table(genes, "Genes_Fst_dicoccoides_dicoccum.txt", sep = "\t", quote = F)