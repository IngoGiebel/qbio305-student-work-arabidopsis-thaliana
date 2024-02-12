#!/bin/bash

cd data

vcftools --gzvcf 1001genomes_snp-short-indel_only_ACGTN.vcf.gz --keep accessions.tsv --bed at_drought_cold_heat_stress_only.bed --minDP 10 --minGQ 20 --minQ 30 --max-missing 0.80 --remove-indels --max-alleles 2 --recode --recode-INFO-all --stdout | bgzip -c > final_1001genomes_snp-short-indel_only_ACGTN_Dp10GQ20Q30_NoIndel_Bialleleic_80PcMissing.vcf.gz

tabix -p vcf final_1001genomes_snp-short-indel_only_ACGTN_Dp10GQ20Q30_NoIndel_Bialleleic_80PcMissing.vcf.gz
