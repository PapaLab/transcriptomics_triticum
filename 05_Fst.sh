#this script runs ANGSD for Fst estimate. The same script was run for each comparison (dicoccoides vs. dicoccum, dicoccum vs. durum, dicoccoides vs. durum)

#prepare the fst
/home/ANGSD/angsd/misc/realSFS fst index ./../dicoccoides.saf.idx ./../dicoccum.saf.idx -sfs ./../dicoccoides.dicoccum.ml -fold 1 -fstout dicoccoides_dicoccum_FOLDED_fst

#get the global estimate
/home/ANGSD/angsd/misc/realSFS fst stats dicoccoides_dicoccum_FOLDED_fst.fst.idx > dicoccoides_dicoccum_FOLDED_fst_stats.txt

/home/ANGSD/angsd/misc/realSFS fst print dicoccoides_dicoccum_FOLDED_fst.fst.idx > dicoccoides_dicoccum_FOLDED_fst.txt
