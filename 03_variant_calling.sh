# this script runs bcftools for variant calling

/m100_work/IscrB_GenoCB/apieri00/sw/bcftools-1.15/bcftools mpileup \
					--fasta-ref /m100_work/IscrB_GenoCB/apieri00/triticum_genomes/iwgsc_refseqv2.1_assembly_ONLY_AB.fa \
					--bam-list /m100_work/IscrB_GenoCB/apieri00/SNP_calling/scripts/all_bam.txt -Ou -o all_var.raw.bcf
					

/m100_work/IscrB_GenoCB/apieri00/sw/bcftools-1.15/bcftools call \
					-m -v -Oz -o all_var.raw.vcf all_var.raw.bcf
