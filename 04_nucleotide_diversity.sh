#this script runs ANGSD for nucleotide diveristy estimates. The same script was run for the three subspecies (dicoccoides, dicoccum and durum)

#Step 1: Finding a 'global estimate' of the SFS
    #First estimate the site allele frequency likelihood
#dicoccoides
/home/ANGSD/angsd/angsd -bam dicoccoides_ANGSD_bam.filelist -doSaf 1 -anc ./../../iwgsc_refseqv2.1_assembly_ONLY_AB.fa -GL 1 -out dicoccoides -P 5 -minMapQ 30 -minQ 20

#Step 2: Calculate the thetas for each site
/home/ANGSD/angsd/misc/realSFS saf2theta dicoccoides.saf.idx -sfs dicoccoides.sfs -outname dicoccoides_theta

#Step 3a: Estimate thetas
/home/ANGSD/angsd/misc/thetaStat do_stat dicoccoides_theta.thetas.idx
