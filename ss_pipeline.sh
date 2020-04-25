#!/bin/bash



mkdir \
-p /mnt/ds2413p/projects/ancient_dna_on_hg38/alignment/ERR3148631



bwa mem \
-M \
-t 4 /mnt/ds2413p/PublicData/ftp.broadinstitute.org/bundle/hg38/Homo_sapiens_assembly38.fasta /mnt/ds2413p/daniyar/Botai/fastq_gz/ERR3148631.fastq.gz > /mnt/ds2413p/projects/ancient_dna_on_hg38/alignment/ERR3148631/ERR3148631.bwamem.sam



picard RevertSam I=/mnt/ds2413p/projects/ancient_dna_on_hg38/alignment/ERR3148631/ERR3148631.bwamem.sam O=/mnt/ds2413p/projects/ancient_dna_on_hg38/alignment/ERR3148631/ERR3148631.picard_RS.bam ATTRIBUTE_TO_CLEAR=XS ATTRIBUTE_TO_CLEAR=XA



picard AddOrReplaceReadGroups INPUT=/mnt/ds2413p/projects/ancient_dna_on_hg38/alignment/ERR3148631/ERR3148631.picard_RS.bam OUTPUT=/mnt/ds2413p/projects/ancient_dna_on_hg38/alignment/ERR3148631/ERR3148631.picard_ARRG.bam RGID=ERR3148630 RGLB=ERR3148630 RGPL=ILLUMINA RGPU=SureSelectV4 RGSM=ERR3148630 RGCN=NLA



picard MergeBamAlignment ALIGNED=/mnt/ds2413p/projects/ancient_dna_on_hg38/alignment/ERR3148631/ERR3148631.bwamem.sam UNMAPPED=/mnt/ds2413p/projects/ancient_dna_on_hg38/alignment/ERR3148631/ERR3148631.picard_ARRG.bam O=/mnt/ds2413p/projects/ancient_dna_on_hg38/alignment/ERR3148631/ERR3148631.picard_MBA.bam R=/mnt/ds2413p/PublicData/ftp.broadinstitute.org/bundle/hg38/Homo_sapiens_assembly38.fasta SORT_ORDER=unsorted CLIP_ADAPTERS=false ADD_MATE_CIGAR=true MAX_INSERTIONS_OR_DELETIONS=-1 PRIMARY_ALIGNMENT_STRATEGY=MostDistant UNMAP_CONTAMINANT_READS=false ATTRIBUTES_TO_RETAIN=XS ATTRIBUTES_TO_RETAIN=XA



picard MarkDuplicates INPUT=/mnt/ds2413p/projects/ancient_dna_on_hg38/alignment/ERR3148631/ERR3148631.picard_MBA.bam OUTPUT=/mnt/ds2413p/projects/ancient_dna_on_hg38/alignment/ERR3148631/ERR3148631.picard_MD.bam METRICS_FILE=/mnt/ds2413p/projects/ancient_dna_on_hg38/alignment/ERR3148631/ERR3148631.picard_MD.metrix.txt OPTICAL_DUPLICATE_PIXEL_DISTANCE=2500 ASSUME_SORT_ORDER=queryname



set \
-o pipefail && picard SortSam INPUT=/mnt/ds2413p/projects/ancient_dna_on_hg38/alignment/ERR3148631/ERR3148631.picard_MD.bam OUTPUT=/dev/stdout SORT_ORDER=coordinate | picard SetNmAndUqTags INPUT=/dev/stdin OUTPUT=/mnt/ds2413p/projects/ancient_dna_on_hg38/alignment/ERR3148631/ERR3148631.picard_SS_SNAUT.bam CREATE_INDEX=true R=/mnt/ds2413p/PublicData/ftp.broadinstitute.org/bundle/hg38/Homo_sapiens_assembly38.fasta



gatk \
-T BaseRecalibrator \
-R /mnt/ds2413p/PublicData/ftp.broadinstitute.org/bundle/hg38/Homo_sapiens_assembly38.fasta \
-I /mnt/ds2413p/projects/ancient_dna_on_hg38/alignment/ERR3148631/ERR3148631.picard_SS_SNAUT.bam \
-knownSites /mnt/ds2413p/PublicData/ftp.broadinstitute.org/bundle/hg38/dbsnp_138.hg38.vcf \
-knownSites /mnt/ds2413p/PublicData/ftp.broadinstitute.org/bundle/hg38/Mills_and_1000G_gold_standard.indels.hg38.vcf \
-knownSites /mnt/ds2413p/PublicData/ftp.broadinstitute.org/bundle/hg38/1000G_phase1.snps.high_confidence.hg38.vcf \
-o /mnt/ds2413p/projects/ancient_dna_on_hg38/alignment/ERR3148631/ERR3148631.gatk_BR.table



gatk \
-T BaseRecalibrator \
-R /mnt/ds2413p/PublicData/ftp.broadinstitute.org/bundle/hg38/Homo_sapiens_assembly38.fasta \
-I /mnt/ds2413p/projects/ancient_dna_on_hg38/alignment/ERR3148631/ERR3148631.picard_SS_SNAUT.bam \
-knownSites /mnt/ds2413p/PublicData/ftp.broadinstitute.org/bundle/hg38/dbsnp_138.hg38.vcf \
-knownSites /mnt/ds2413p/PublicData/ftp.broadinstitute.org/bundle/hg38/Mills_and_1000G_gold_standard.indels.hg38.vcf \
-knownSites /mnt/ds2413p/PublicData/ftp.broadinstitute.org/bundle/hg38/1000G_phase1.snps.high_confidence.hg38.vcf \
-BQSR /mnt/ds2413p/projects/ancient_dna_on_hg38/alignment/ERR3148631/ERR3148631.gatk_BR.table \
-o /mnt/ds2413p/projects/ancient_dna_on_hg38/alignment/ERR3148631/ERR3148631.gatk_BR_BQSR.table



gatk \
-T AnalyzeCovariates \
-R /mnt/ds2413p/PublicData/ftp.broadinstitute.org/bundle/hg38/Homo_sapiens_assembly38.fasta \
-before /mnt/ds2413p/projects/ancient_dna_on_hg38/alignment/ERR3148631/ERR3148631.gatk_BR.table \
-after /mnt/ds2413p/projects/ancient_dna_on_hg38/alignment/ERR3148631/ERR3148631.gatk_BR_BQSR.table \
-plots /mnt/ds2413p/projects/ancient_dna_on_hg38/alignment/ERR3148631/ERR3148631.gatk_AC.pdf



gatk \
-T PrintReads \
-R /mnt/ds2413p/PublicData/ftp.broadinstitute.org/bundle/hg38/Homo_sapiens_assembly38.fasta \
-I /mnt/ds2413p/projects/ancient_dna_on_hg38/alignment/ERR3148631/ERR3148631.picard_SS_SNAUT.bam \
-BQSR /mnt/ds2413p/projects/ancient_dna_on_hg38/alignment/ERR3148631/ERR3148631.gatk_BR_BQSR.table \
-o /mnt/ds2413p/projects/ancient_dna_on_hg38/alignment/ERR3148631/ERR3148631.gatk_PR_BR_BQSR.bam



gatk \
-T HaplotypeCaller \
-R /mnt/ds2413p/PublicData/ftp.broadinstitute.org/bundle/hg38/Homo_sapiens_assembly38.fasta \
-I /mnt/ds2413p/projects/ancient_dna_on_hg38/alignment/ERR3148631/ERR3148631.gatk_PR_BR_BQSR.bam \
--dbsnp /mnt/ds2413p/PublicData/ftp.broadinstitute.org/bundle/hg38/dbsnp_138.hg38.vcf \
-o /mnt/ds2413p/projects/ancient_dna_on_hg38/alignment/ERR3148631/ERR3148631.gatk_HC.vcf \
-ERC NONE \
--max_alternate_alleles 3 \
--read_filter OverclippedRead \
--emitDroppedReads \
-bamout /mnt/ds2413p/projects/ancient_dna_on_hg38/alignment/ERR3148631/ERR3148631.gatk_HC.bam \
-variant_index_type LINEAR \
-variant_index_parameter 128000


# added echo to skip
echo gatk \
-T VariantRecalibrator \
-R /mnt/ds2413p/PublicData/ftp.broadinstitute.org/bundle/hg38/Homo_sapiens_assembly38.fasta \
-input /mnt/ds2413p/projects/ancient_dna_on_hg38/alignment/ERR3148631/ERR3148631.gatk_HC.vcf \
-resource:hapmap,known=false,training=true,truth=true,prior=15.0 /mnt/ds2413p/PublicData/ftp.broadinstitute.org/bundle/hg38/hapmap_3.3.hg38.vcf \
-resource:omni,known=false,training=true,truth=true,prior=12.0 /mnt/ds2413p/PublicData/ftp.broadinstitute.org/bundle/hg38/1000G_omni2.5.hg38.vcf \
-resource:1000G,known=false,training=true,truth=false,prior=10.0 /mnt/ds2413p/PublicData/ftp.broadinstitute.org/bundle/hg38/1000G_phase1.snps.high_confidence.hg38.vcf \
-resource:dbsnp,known=true,training=false,truth=false,prior=2.0 /mnt/ds2413p/PublicData/ftp.broadinstitute.org/bundle/hg38/dbsnp_138.hg38.vcf \
-an DP \
-an QD \
-an FS \
-an MQRankSum \
-an ReadPosRankSum \
-mode SNP \
-tranche 100.0 \
-tranche 99.9 \
-tranche 99.0 \
-tranche 90.0 \
--maxGaussians 4 \
-recalFile /mnt/ds2413p/projects/ancient_dna_on_hg38/alignment/ERR3148631/ERR3148631.gatk_VR_SNP.recal \
-tranchesFile /mnt/ds2413p/projects/ancient_dna_on_hg38/alignment/ERR3148631/ERR3148631.gatk_VR_SNP.tranches \
-rscriptFile /mnt/ds2413p/projects/ancient_dna_on_hg38/alignment/ERR3148631/ERR3148631.gatk_VR_SNP_plots.Rscript


# added echo to skip
echo gatk \
-T ApplyRecalibration \
-R /mnt/ds2413p/PublicData/ftp.broadinstitute.org/bundle/hg38/Homo_sapiens_assembly38.fasta \
-input /mnt/ds2413p/projects/ancient_dna_on_hg38/alignment/ERR3148631/ERR3148631.gatk_HC.vcf \
-mode SNP \
--ts_filter_level 99.0 \
-recalFile /mnt/ds2413p/projects/ancient_dna_on_hg38/alignment/ERR3148631/ERR3148631.gatk_VR_SNP.recal \
-tranchesFile /mnt/ds2413p/projects/ancient_dna_on_hg38/alignment/ERR3148631/ERR3148631.gatk_VR_SNP.tranches \
-o /mnt/ds2413p/projects/ancient_dna_on_hg38/alignment/ERR3148631/ERR3148631.gatk_AR_SNP.vcf


# added echo to skip this step
echo gatk \
-T VariantRecalibrator \
-R /mnt/ds2413p/PublicData/ftp.broadinstitute.org/bundle/hg38/Homo_sapiens_assembly38.fasta \
-input /mnt/ds2413p/projects/ancient_dna_on_hg38/alignment/ERR3148631/ERR3148631.gatk_HC.vcf \
-resource:mills,known=true,training=true,truth=true,prior=12.0 /mnt/ds2413p/PublicData/ftp.broadinstitute.org/bundle/hg38/Mills_and_1000G_gold_standard.indels.hg38.vcf \
-an DP \
-an FS \
-an MQRankSum \
-an ReadPosRankSum \
-mode INDEL \
-tranche 100.0 \
-tranche 99.9 \
-tranche 99.0 \
-tranche 90.0 \
--maxGaussians 4 \
-recalFile /mnt/ds2413p/projects/ancient_dna_on_hg38/alignment/ERR3148631/ERR3148631.gatk_VR_INDEL.recal \
-tranchesFile /mnt/ds2413p/projects/ancient_dna_on_hg38/alignment/ERR3148631/ERR3148631.gatk_VR_INDEL.tranches \
-rscriptFile /mnt/ds2413p/projects/ancient_dna_on_hg38/alignment/ERR3148631/ERR3148631.gatk_VR_INDEL_plots.Rscript


# added echo to skip this step
echo gatk \
-T ApplyRecalibration \
-R /mnt/ds2413p/PublicData/ftp.broadinstitute.org/bundle/hg38/Homo_sapiens_assembly38.fasta \
-input /mnt/ds2413p/projects/ancient_dna_on_hg38/alignment/ERR3148631/ERR3148631.gatk_AR_SNP.vcf \
-mode INDEL \
--ts_filter_level 99.0 \
-recalFile /mnt/ds2413p/projects/ancient_dna_on_hg38/alignment/ERR3148631/ERR3148631.gatk_VR_INDEL.recal \
-tranchesFile /mnt/ds2413p/projects/ancient_dna_on_hg38/alignment/ERR3148631/ERR3148631.gatk_VR_INDEL.tranches \
-o /mnt/ds2413p/projects/ancient_dna_on_hg38/alignment/ERR3148631/ERR3148631.gatk_AR_INDEL.vcf



gatk \
-T SelectVariants \
-R /mnt/ds2413p/PublicData/ftp.broadinstitute.org/bundle/hg38/Homo_sapiens_assembly38.fasta \
-V /mnt/ds2413p/projects/ancient_dna_on_hg38/alignment/ERR3148631/ERR3148631.gatk_HC.vcf \
-selectType SNP \
-o /mnt/ds2413p/projects/ancient_dna_on_hg38/alignment/ERR3148631/ERR3148631.gatk_SV_SNP_raw.vcf



gatk \
-T VariantFiltration \
-R /mnt/ds2413p/PublicData/ftp.broadinstitute.org/bundle/hg38/Homo_sapiens_assembly38.fasta \
-V /mnt/ds2413p/projects/ancient_dna_on_hg38/alignment/ERR3148631/ERR3148631.gatk_SV_SNP_raw.vcf \
--filterExpression "QD < 2.0 || FS > 60.0 || MQ < 40.0 || MQRankSum < -12.5 || ReadPosRankSum < -8.0" \
--filterName "SNP_FAIL" \
-o /mnt/ds2413p/projects/ancient_dna_on_hg38/alignment/ERR3148631/ERR3148631.gatk_SV_SNP_fil.vcf



gatk \
-T SelectVariants \
-R /mnt/ds2413p/PublicData/ftp.broadinstitute.org/bundle/hg38/Homo_sapiens_assembly38.fasta \
-V /mnt/ds2413p/projects/ancient_dna_on_hg38/alignment/ERR3148631/ERR3148631.gatk_HC.vcf \
-selectType INDEL \
-o /mnt/ds2413p/projects/ancient_dna_on_hg38/alignment/ERR3148631/ERR3148631.gatk_SV_INDEL_raw.vcf



gatk \
-T VariantFiltration \
-R /mnt/ds2413p/PublicData/ftp.broadinstitute.org/bundle/hg38/Homo_sapiens_assembly38.fasta \
-V /mnt/ds2413p/projects/ancient_dna_on_hg38/alignment/ERR3148631/ERR3148631.gatk_SV_INDEL_raw.vcf \
--filterExpression "QD < 2.0 || FS > 200.0 || ReadPosRankSum < -20.0" \
--filterName "INDEL_FAIL" \
-o /mnt/ds2413p/projects/ancient_dna_on_hg38/alignment/ERR3148631/ERR3148631.gatk_SV_INDEL_fil.vcf



vcf-concat /mnt/ds2413p/projects/ancient_dna_on_hg38/alignment/ERR3148631/ERR3148631.gatk_SV_SNP_fil.vcf /mnt/ds2413p/projects/ancient_dna_on_hg38/alignment/ERR3148631/ERR3148631.gatk_SV_INDEL_fil.vcf > /mnt/ds2413p/projects/ancient_dna_on_hg38/alignment/ERR3148631/ERR3148631.vcftools_concat.FINAL.vcf



cat /mnt/ds2413p/projects/ancient_dna_on_hg38/alignment/ERR3148631/ERR3148631.vcftools_concat.FINAL.vcf | vcf-sort > /mnt/ds2413p/projects/ancient_dna_on_hg38/alignment/ERR3148631/ERR3148631.vcftools_sorted.FINAL.vcf

