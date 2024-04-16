#this script runs STAR for reads alignment

for i in $(ls /m100_work/IscrB_GenoCB/apieri00/reads_fastq/*_R1.cleaned.paired.fastq.gz); do

    THR=18
    FULLSAMPLE=${i/_R1.cleaned.paired.fastq.gz/}; SAMPLE=$(basename $FULLSAMPLE)
    

    mkdir -p ${WORKDIR}/mappings_Singleend/${SAMPLE} ${WORKDIR}/mappings_Singleend/${SAMPLE}/logs
    cd ${WORKDIR}/mappings_Singleend/${SAMPLE}
    
   
    STAR --runThreadN ${THR} --genomeDir /m100_work/IscrB_GenoCB/apieri00/triticum_genomes/ChineseSpring_index \
        --readFilesIn ${FULLSAMPLE}_R1.cleaned.paired.fastq.gz ${FULLSAMPLE}_R2.cleaned.paired.fastq.gz \
        --outFileNamePrefix ${WORKDIR}/mappings_Singleend/${SAMPLE}/${SAMPLE}. --outSAMtype BAM SortedByCoordinate --outSAMunmapped Within --quantMode TranscriptomeSAM --quantTranscriptomeBan Singleend \
        --outFilterMultimapNmax 20 --outBAMsortingThreadN ${THR} \
        --readFilesCommand zcat --twopassMode Basic >${WORKDIR}/mappings_Singleend/${SAMPLE}/logs/mapping_${SAMPLE}.out 2>${WORKDIR}/mappings_Singleend/${SAMPLE}/logs/mapping_${SAMPLE}.err

done

