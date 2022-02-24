#!/bin/bash

# this assumes you have the following tools on your system
# picard 
# gatk (we used version 3.8.1)
# bcftools 
# vcftools 
# samtools

# CREATE FASTA INDEX
printf "*\n*\n*\nNow creating FASTA index with samtools faidx\n"
samtools faidx genomeAssembly.fa \
	 -o genomeAssembly.fa.fai


# CREATE FASTA SEQUENCE DICTIONARY FILE
printf "*\n*\n*\nNow creating FASTA sequence dictionary using Picard\n"
java -Xmx64g \
	 -jar picard.jar CreateSequenceDictionary \
	 R=genomeAssembly.fa \
	 O=genomeAssembly.dict

# EXTRACT SNPS and INDELS INTO SEPARATE FILES
#
printf "*\n*\n*\nNow extracting SNPS\n"
java -Xmx64g \
     -jar GenomeAnalysisTK.jar \
     -T SelectVariants \
     -R genomeAssembly.fa \
     -V pythonfiltered.vcf \
     -selectType SNP \
     -o raw_snps.vcf

printf "*\n*\n*\nNow extracting indels\n"
java -Xmx64g \
     -jar GenomeAnalysisTK.jar \
     -T SelectVariants \
     -R genomeAssembly.fa \
     -V pythonfiltered.vcf \
     -selectType INDEL \
     -o raw_indels.vcf


# HARD FILTERING SNPS and INDELS ACCORDING TO GATK BEST PRACTICES
#
printf "*\n*\n*\nNow HARD filtering SNPs\n"
java -Xmx64g \
    -jar GenomeAnalysisTK.jar \
    -T VariantFiltration \
    -R genomeAssembly.fa \
    -V raw_snps.vcf \
    --filterExpression "QD < 2.0 || FS > 60.0 || MQ < 40.0 || MQRankSum < -12.5 || ReadPosRankSum < -8.0" \
    --filterName "GATKBP-snp" \
    -o snps.GATKBP-flagged.vcf

printf "*\n*\n*\nNow HARD filtering indels\n"
java -Xmx64g \
    -jar GenomeAnalysisTK.jar \
    -T VariantFiltration \
    -R genomeAssembly.fa \
    -V raw_indels.vcf \
    --filterExpression "QD < 2.0 || FS > 200.0 || ReadPosRankSum < -20.0" \
    --filterName "GATKBP-indel" \
    -o indels.GATKBP-flagged.vcf


# EXTRACT VARIANTS THAT PASSED
printf "*\n*\n*\nNow extracting variants that passed the filters\n"
bcftools view -f "PASS" snps.GATKBP-flagged.vcf > snps.GATKBP-passed.vcf
rm snps.GATKBP-flagged.vcf
bcftools view -f "PASS" indels.GATKBP-flagged.vcf > indels.GATKBP-passed.vcf
rm indels.GATKBP-flagged.vcf


# CONCATENATE THE TWO FILES, THEN SORT THEM
vcf-concat snps.GATKBP-passed.vcf indels.GATKBP-passed.vcf > GATKBP-passed.unsorted.vcf
vcf-sort GATKBP-passed.unsorted.vcf > GATKBP-passed.vcf

# CLEAN UP
rm GATKBP-passed.unsorted.vcf
