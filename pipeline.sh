#!/bin/bash

/mnt/ds2413p/daniyar/Pipeline/bwa/bwa-0.7.17/bwa mem -M -t 16 /mnt/ds2413p/daniyar/hg38/Homo_sapiens_assembly38.fasta /mnt/ds2413p/daniyar/Botai/fastq_gz/ERR3148631.fastq.gz > /mnt/ds2413p/daniyar/Botai/botai_6/ERR3148631/ERR3148631.mem.sam &&
/mnt/ds2413p/daniyar/Pipeline/JavaScript/k8-0.2.5/k8-Linux /mnt/ds2413p/daniyar/Pipeline/bwa/bwa-0.7.17/bwakit/bwa-postalt.js /mnt/ds2413p/daniyar/hg38/Homo_sapiens_assembly38.fasta.alt /mnt/ds2413p/daniyar/Botai/botai_6/ERR3148631/ERR3148631.mem.sam  > /mnt/ds2413p/daniyar/Botai/botai_6/ERR3148631/ERR3148631.mem_postalt.sam &&

/mnt/ds2413p/daniyar/Pipeline/java/jdk1.8.0_101/bin/java -jar /mnt/ds2413p/daniyar/Pipeline/picard/picard_2_18_7/picard.jar \
RevertSam \
    I=/mnt/ds2413p/daniyar/Botai/botai_6/ERR3148631/ERR3148631.mem_postalt.sam O=/mnt/ds2413p/daniyar/Botai/botai_6/ERR3148631/ERR3148631_u.bam \
    ATTRIBUTE_TO_CLEAR=XS ATTRIBUTE_TO_CLEAR=XA &&
	
/mnt/ds2413p/daniyar/Pipeline/java/jdk1.8.0_101/bin/java -jar /mnt/ds2413p/daniyar/Pipeline/picard/picard_2_18_7/picard.jar AddOrReplaceReadGroups \
    I=/mnt/ds2413p/daniyar/Botai/botai_6/ERR3148631/ERR3148631_u.bam O=/mnt/ds2413p/daniyar/Botai/botai_6/ERR3148631/ERR3148631_rg.bam \
    RGID=altalt RGSM=altalt RGLB=wgsim RGPU=shlee RGPL=illumina &&
	
/mnt/ds2413p/daniyar/Pipeline/java/jdk1.8.0_101/bin/java -jar /mnt/ds2413p/daniyar/Pipeline/picard/picard_2_18_7/picard.jar MergeBamAlignment \
    ALIGNED=/mnt/ds2413p/daniyar/Botai/botai_6/ERR3148631/ERR3148631.mem_postalt.sam UNMAPPED=/mnt/ds2413p/daniyar/Botai/botai_6/ERR3148631/ERR3148631_rg.bam O=/mnt/ds2413p/daniyar/Botai/botai_6/ERR3148631/ERR3148631_m.bam \
    R=/mnt/ds2413p/daniyar/hg38/Homo_sapiens_assembly38.fasta \
    SORT_ORDER=unsorted CLIP_ADAPTERS=false \
    ADD_MATE_CIGAR=true MAX_INSERTIONS_OR_DELETIONS=-1 \
    PRIMARY_ALIGNMENT_STRATEGY=MostDistant \
    UNMAP_CONTAMINANT_READS=false \
    ATTRIBUTES_TO_RETAIN=XS ATTRIBUTES_TO_RETAIN=XA &&
	
/mnt/ds2413p/daniyar/Pipeline/java/jdk1.8.0_101/bin/java -jar /mnt/ds2413p/daniyar/Pipeline/picard/picard_2_18_7/picard.jar MarkDuplicates \
    INPUT=/mnt/ds2413p/daniyar/Botai/botai_6/ERR3148631/ERR3148631_m.bam OUTPUT=/mnt/ds2413p/daniyar/Botai/botai_6/ERR3148631/ERR3148631_md.bam METRICS_FILE=/mnt/ds2413p/daniyar/Botai/botai_6/ERR3148631/ERR3148631_md.bam.txt \
    OPTICAL_DUPLICATE_PIXEL_DISTANCE=2500 ASSUME_SORT_ORDER=queryname && 
	
set -o pipefail
/mnt/ds2413p/daniyar/Pipeline/java/jdk1.8.0_101/bin/java -jar /mnt/ds2413p/daniyar/Pipeline/picard/picard_2_18_7/picard.jar SortSam \
    INPUT=/mnt/ds2413p/daniyar/Botai/botai_6/ERR3148631/ERR3148631_md.bam OUTPUT=/dev/stdout SORT_ORDER=coordinate | \
    /mnt/ds2413p/daniyar/Pipeline/java/jdk1.8.0_101/bin/java -jar /mnt/ds2413p/daniyar/Pipeline/picard/picard_2_18_7/picard.jar SetNmMdAndUqTags  \
    INPUT=/dev/stdin OUTPUT=/mnt/ds2413p/daniyar/Botai/botai_6/ERR3148631/ERR3148631_snaut.bam \
    CREATE_INDEX=true R=/mnt/ds2413p/daniyar/hg38/Homo_sapiens_assembly38.fasta &&
	
/mnt/ds2413p/daniyar/Pipeline/samtools/samtools-1.4.1/samtools view -h /mnt/ds2413p/daniyar/Botai/botai_6/ERR3148631/ERR3148631_snaut.bam | gawk '{printf "%s\t", $1; if(and($2,0x1)) {t=$2-0x1}else{t=$2}; printf "%s\t" , t; for (i=3; i<NF; i++){printf "%s\t", $i} ; printf "%s\n",$NF}'| /mnt/ds2413p/daniyar/Pipeline/samtools/samtools-1.4.1/samtools view -Sb - > /mnt/ds2413p/daniyar/Botai/botai_6/ERR3148631/ERR3148631_se.bam &&
/mnt/ds2413p/daniyar/Pipeline/samtools/samtools-1.4.1/samtools index /mnt/ds2413p/daniyar/Botai/botai_6/ERR3148631/ERR3148631_se.bam &&

/mnt/ds2413p/daniyar/Pipeline/java/jdk1.8.0_101/bin/java -jar /mnt/ds2413p/daniyar/Pipeline/gatk/GATK_3_8_1/GenomeAnalysisTK.jar -T HaplotypeCaller \
    -R /mnt/ds2413p/daniyar/hg38/Homo_sapiens_assembly38.fasta \
    -I /mnt/ds2413p/daniyar/Botai/botai_6/ERR3148631/ERR3148631_se.bam -o /mnt/ds2413p/daniyar/Botai/botai_6/ERR3148631/ERR3148631_hc.g.vcf \
    -ERC GVCF \
    --emitDroppedReads -bamout /mnt/ds2413p/daniyar/Botai/botai_6/ERR3148631/ERR3148631_hc.bam &&

##EXTRACT MULTISAMPLE 

#/mnt/ds2413p/daniyar/Pipeline/java/jdk1.8.0_101/bin/java -jar /mnt/ds2413p/daniyar/Pipeline/gatk/GATK_3_8_1/GenomeAnalysisTK.jar -T GenotypeGVCFs \
    -R /mnt/ds2413p/daniyar/hg38/Homo_sapiens_assembly38.fasta -o /mnt/ds2413p/daniyar/Multisample/multisample.vcf \
    --variant /mnt/ds2413p/daniyar/Botai/botai_6/ERR3148630/ERR3148630_hc.g.vcf --variant /mnt/ds2413p/daniyar/Botai/botai_6/ERR3148631/ERR3148631_hc.g.vcf --variant /mnt/ds2413p/daniyar/Botai/botai_6/ERR3148632/ERR3148632_hc.g.vcf --variant /mnt/ds2413p/daniyar/Botai/botai_6/ERR3148633/ERR3148633_hc.g.vcf  --variant /mnt/ds2413p/daniyar/Botai/botai_6/ERR3380195/ERR3380195_hc.g.vcf --variant /mnt/ds2413p/daniyar/Botai/botai_6/ERR3380196/ERR3380196_hc.g.vcf --variant /mnt/ds2413p/daniyar/KAZ_WG/KAZ_WG_5_hg38/KAZ_WG2/KAZ_WG2_hc.g.vcf --variant /mnt/ds2413p/daniyar/KAZ_WG/KAZ_WG_5_hg38/KAZ_WG4/KAZ_WG4_hc.g.vcf --variant /mnt/ds2413p/daniyar/KAZ_WG/KAZ_WG_5_hg38/KAZ_WG5/KAZ_WG5_hc.g.vcf --variant /mnt/ds2413p/daniyar/KAZ_WG/KAZ_WG_5_hg38/KAZ_WG6/KAZ_WG6_hc.g.vcf --variant /mnt/ds2413p/daniyar/KAZ_WG/KAZ_WG_5_hg38/KAZ_WG7/KAZ_WG7_hc.g.vcf 
	