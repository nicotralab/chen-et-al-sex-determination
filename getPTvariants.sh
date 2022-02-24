#!/bin/bash

# In our dataset, the female parent is in column 100 of the VCF file 
# and the male parent is in column 101 of the VCF file. 

# CREATE FILES OF "aa x ab", "ab x aa", or "ab x ab" VARIANTS
awk '$100~/0\/1/ && $101~/0\/0/ {print}' GATKBP-passed.vcf > GATKBP-passed.femaleHetMale00.txt
awk '$100~/0\/0/ && $101~/0\/1/ {print}' GATKBP-passed.vcf > GATKBP-passed.maleHetFemale00.txt
awk '$100~/1\/1/ && $101~/0\/1/ {print}' GATKBP-passed.vcf > GATKBP-passed.maleHetFemale11.txt
#awk '$100~/0\/1/ && $101~/1\/1/ {print}' GATKBP-passed.vcf > GATKBP-passed.femaleHetMale11.txt
#awk '$100~/0\/1/ && $101~/0\/1/ {print}' GATKBP-passed.vcf > GATKBP-passed.bothHet.txt


# reattach the header
grep '^#' GATKBP-passed.vcf > GATKBP-passed.vcf.header.txt
cat GATKBP-passed.vcf.header.txt GATKBP-passed.femaleHetMale00.txt > GATKBP-passed.femaleHetMale00.vcf
cat GATKBP-passed.vcf.header.txt GATKBP-passed.maleHetFemale00.txt > GATKBP-passed.maleHetFemale00.vcf
cat GATKBP-passed.vcf.header.txt GATKBP-passed.maleHetFemale11.txt > GATKBP-passed.maleHetFemale11.vcf
#cat GATKBP-passed.vcf.header.txt GATKBP-passed.femaleHetMale11.txt > GATKBP-passed.femaleHetMale11.vcf
#cat GATKBP-passed.vcf.header.txt GATKBP-passed.bothHet.txt > GATKBP-passed.bothHet.vcf
rm GATKBP-passed.femaleHet*.txt GATKBP-passed.maleHet*.txt # GATKBP-passed.bothHet.txt

# calculate statistics for each dataset
bcftools stats -s - GATKBP-passed.femaleHetMale00.vcf > GATKBP-passed.femaleHetMale00.vcf.STATS
bcftools stats -s - GATKBP-passed.maleHetFemale00.vcf > GATKBP-passed.maleHetFemale00.vcf.STATS
bcftools stats -s - GATKBP-passed.maleHetFemale11.vcf > GATKBP-passed.maleHetFemale11.vcf.STATS
#bcftools stats -s - GATKBP-passed.femaleHetMale11.vcf > GATKBP-passed.femaleHetMale11.vcf.STATS
#bcftools stats -s - GATKBP-passed.bothHet.vcf > GATKBP-passed.bothHet.vcf.STATS

# now assess how many impossible genotypes are present:
#
# If only one parent is heterozygous (e.g. 0/0 and 0/1) no offspring should be homozygous
# for the minor allele at that locus (e.g. no 1/1).
#
# This code will count the number of "impossible" genotypes in the male and female datasets
# and then save a "*.hist" file that can be used to create a histogram of counts of impossible genotype per variant

# for the femaleHetmale00:
echo "Calculating histogram for femaleHetMale00"
while read i ; do echo $i | grep -o '1/1' | wc -l ; done < GATKBP-passed.femaleHetMale00.vcf > counts.txt
for i in {0..92} ; do grep -c "^$i$" counts.txt ; done > GATKBP-passed.femaleHetMale00.11counts.hist.txt
rm counts.txt

# for the maleHetfemale00:
echo "Calculating histogram for maleHetFemale00"
while read i ; do echo $i | grep -o '1/1' | wc -l ; done < GATKBP-passed.maleHetFemale00.vcf > counts.txt
for i in {0..92} ; do grep -c "^$i$" counts.txt ; done > GATKBP-passed.maleHetFemale00.11counts.hist.txt
rm counts.txt

# for the maleHetfemale11:
echo "Calculating histogram for maleHetFemale11"
while read i ; do echo $i | grep -o '0/0' | wc -l ; done < GATKBP-passed.maleHetFemale11.vcf > counts.txt
for i in {0..92} ; do grep -c "^$i$" counts.txt ; done > GATKBP-passed.maleHetFemale11.00counts.hist.txt
rm counts.txt

# looking at these histograms identifies variants that have a few impossible genotypes (likely genotyping errors) or many impossible genotypes.
# These create a bell curve around 0.25, which is where you'd expect a marker to appear if it was actually an "ab x ab" type.
# Therefore, all impossible genotypes were changed to missing data using sed, and variants with more than 10% missing data were removed.

# for the femaleHetMale00:
#   change 1/1 to ./.
#   use vcf-tools to remove sites with >10% missing data.
sed 's/1\/1/.\/./g' GATKBP-passed.femaleHetMale00.vcf > GATKBP-passed.femaleHetMale00.11toMissing.vcf
vcftools --vcf GATKBP-passed.femaleHetMale00.11toMissing.vcf --max-missing 0.90 --recode --recode-INFO-all --out GATKBP-passed.femaleHetMale00.abxabRemoved
mv GATKBP-passed.femaleHetMale00.abxabRemoved.recode.vcf GATKBP-passed.femaleHet.abxabRemoved.vcf

# for the maleHet datasets
#   change 1/1 or 0/0 to ./., as approprite
sed 's/1\/1/.\/./g' GATKBP-passed.maleHetFemale00.vcf > GATKBP-passed.maleHetFemale00.11toMissing.vcf
sed 's/0\/0/.\/./g' GATKBP-passed.maleHetFemale11.vcf > GATKBP-passed.maleHetFemale11.00toMissing.vcf

# now concatenate, sort and remove the SNPs with missing data
vcf-concat GATKBP-passed.maleHetFemale00.11toMissing.vcf GATKBP-passed.maleHetFemale11.00toMissing.vcf > GATKBP-passed.maleHet.00toMissing.unsorted.vcf
vcf-sort GATKBP-passed.maleHet.00toMissing.unsorted.vcf > GATKBP-passed.maleHet.00toMissing.vcf
vcftools --vcf GATKBP-passed.maleHet.00toMissing.vcf --max-missing 0.90 --recode --recode-INFO-all --out GATKBP-passed.maleHet.abxabRemoved
mv GATKBP-passed.maleHet.abxabRemoved.recode.vcf GATKBP-passed.maleHet.abxabRemoved.vcf

# generate stats files
bcftools stats -s - GATKBP-passed.maleHet.abxabRemoved.vcf > GATKBP-passed.maleHet.abxabRemoved.vcf.stats
bcftools stats -s - GATKBP-passed.femaleHet.abxabRemoved.vcf > GATKBP-passed.femaleHet.abxabRemoved.vcf.stats
[mln29@login0b HY526_filteringVariantsForPaper]$
