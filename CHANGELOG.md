
# Changelog
All notable changes to this project will be documented in this file.

## [Unreleased]

## 0.5.13 - 2020-04-011

RGID fix using copy.deepcopy 


## 0.5.12 - 2020-04-011

added tmp_dir option for picard 

## 0.5.11 - 2020-04-05

added CHANGELOG.md

added envorionment.yml

updated README.md 

## 0.5.10

added new pipeline "wgs_pipeline_on_hg38_WO_VQSR.py" without VQSR step

## 0.5.9

-D changed to --dbsnp

## 0.5.8

add_tokens default value is False

## 0.5.7

added -D {dbsnp} to gatk_HC for adding rs to vcf file

added VQSR steps

## 0.5.6

added to pipeline get_cmd_gatk_AC_pdf

-ERC GVCF changed to -ERC NONE to show only 0/1 and 1/1 in vcf
