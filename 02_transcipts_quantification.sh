# this script runs salmon for transcripts quantification

for i in $(ls /gatk/data/Alice_RNAseq/alignments_toTranscr_CS_Singleend/*.Aligned.toTranscriptome.out.bam); do

FULLSAMPLE=${i/.Aligned.toTranscriptome.out.bam/}; SAMPLE=$(basename $FULLSAMPLE)

${Salmon} quant -t ./../iwgsc_refseqv2.1_gene_annotation_200916/iwgsc_refseqv2.1_annotation_200916_HC_mrna_ONLY_AB.fasta \
			-l A \
			-a ${i} \
			-p 8 \
			-o ${SAMPLE}
			
done
