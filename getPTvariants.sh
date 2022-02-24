#!/bin/bash



# this will likely be specific to your system
module load picard gatk/3.8.1 gcc/8.2.0 bcftools vcftools samtools

# CREATE FASTA INDEX
printf "*\n*\n*\nNow creating FASTA index with samtools faidx\n"
samtools faidx startingfiles/hsym.wildtype.canu.dovetail.pbjelly.arrow.pilon.primary.polished.renamed.fa \
	 -o startingfiles/hsym.wildtype.canu.dovetail.pbjelly.arrow.pilon.primary.polished.renamed.fa.fai


# CREATE FASTA SEQUENCE DICTIONARY FILE
printf "*\n*\n*\nNow creating FASTA sequence dictionary using GATK\n"
java -Xmx64g \
	 -jar /ihome/crc/install/picard/2.18.12/picard.jar CreateSequenceDictionary \
	 R=startingfiles/hsym.wildtype.canu.dovetail.pbjelly.arrow.pilon.primary.polished.renamed.fa \
	 O=startingfiles/hsym.wildtype.canu.dovetail.pbjelly.arrow.pilon.primary.polished.renamed.dict

# EXTRACT SNPS and INDELS INTO SEPARATE FILES
#
# The path to the .jar file will likely be specific to your system, and may not need to be specified.
printf "*\n*\n*\nNow extracting SNPS\n"
java -Xmx64g \
     -jar /ihome/crc/install/gatk/GenomeAnalysisTK-3.8-1/GenomeAnalysisTK.jar \
     -T SelectVariants \
     -R startingfiles/hsym.wildtype.canu.dovetail.pbjelly.arrow.pilon.primary.polished.renamed.fa \
	 -V pythonfiltered.vcf \
     -selectType SNP \
     -o raw_snps.vcf

printf "*\n*\n*\nNow extracting indels\n"
java -Xmx64g \
     -jar /ihome/crc/install/gatk/GenomeAnalysisTK-3.8-1/GenomeAnalysisTK.jar \
     -T SelectVariants \
     -R startingfiles/hsym.wildtype.canu.dovetail.pbjelly.arrow.pilon.primary.polished.renamed.fa \
     -V pythonfiltered.vcf \
     -selectType INDEL \
     -o raw_indels.vcf


# HARD FILTERING SNPS and INDELS ACCORDING TO GATK BEST PRACTICES
#
# The path the to the .jar file will likely be specific to your system, and may not need to be specified.
printf "*\n*\n*\nNow HARD filtering SNPs\n"
java -Xmx64g \
    -jar /ihome/crc/install/gatk/GenomeAnalysisTK-3.8-1/GenomeAnalysisTK.jar \
    -T VariantFiltration \
    -R startingfiles/hsym.wildtype.canu.dovetail.pbjelly.arrow.pilon.primary.polished.renamed.fa \
    -V raw_snps.vcf \
    --filterExpression "QD < 2.0 || FS > 60.0 || MQ < 40.0 || MQRankSum < -12.5 || ReadPosRankSum < -8.0" \
    --filterName "GATKBP-snp" \
    -o snps.GATKBP-flagged.vcf

printf "*\n*\n*\nNow HARD filtering indels\n"
java -Xmx64g \
    -jar /ihome/crc/install/gatk/GenomeAnalysisTK-3.8-1/GenomeAnalysisTK.jar \
    -T VariantFiltration \
    -R startingfiles/hsym.wildtype.canu.dovetail.pbjelly.arrow.pilon.primary.polished.renamed.fa \
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
rm GATKBP-passed.unsorted.vcf
