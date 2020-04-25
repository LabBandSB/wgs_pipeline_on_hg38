# wgs_pipeline_on_hg38

whole human genome sequence alignment/mapping pipeline on hg38 reference genome based on https://gatk.broadinstitute.org/hc/en-us/articles/360037498992--How-to-Map-reads-to-a-reference-with-alternate-contigs-like-GRCH38

reference genome and public databases  were downloaded from broadinstitute

wget -m -N ftp://gsapubftp-anonymous@ftp.broadinstitute.org/bundle/hg38

gatk version used in this pipeline is 3.8-1-0-gf15c1c3ef

### installation

conda env create --file environment.yml

### list of files:

- .gitignore
- CHANGELOG.md - list of all changes 
- environment.yml - instructions for new conda environment
- pipeline.sh - old version of pipeline 
- README.md - this file
- ss_pipeline.sh - sample of script generated by wgs_pipeline_on_hg38.py   
- wgs_pipeline_on_hg38.py - full version of pipeline
- wgs_pipeline_on_hg38_WO_VQSR.py - pipeline version without VQSR

 