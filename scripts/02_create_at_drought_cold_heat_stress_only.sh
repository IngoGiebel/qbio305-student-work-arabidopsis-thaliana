#!/bin/bash

cd data

grep -e "drought-responsive" -e "RD29A" -e "Drought-responsive family protein" -e "DREB2A" -e "Calmodulin-independent protein kinase" -e "EXS" -e "drought-induced" -e "ABA overly sensitive mutant 3" -e "Drought-repressed 4" -e "CDSP32" -e "WD40 repeat family" -e "Drought" -e "cold" -e "heat" TAIR10_GFF3_genes_annot.gff > at_drought_cold_heat_stress_only.gff

awk '{OFS="\t"; print $1, $4, $5}' at_drought_cold_heat_stress_only.gff > at_drought_cold_heat_stress_only.bed
