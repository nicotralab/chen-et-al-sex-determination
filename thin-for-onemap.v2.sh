#!/bin/bash


############
#  FEMALE  #
############

# thin to 5kb
vcftools --gzvcf femalePT.vcf.gz --thin 5000 --recode --recode-INFO-all --out femalePT.5kb

mv femalePT.5kb.recode.vcf femalePT.5kb.vcf

# Create onemap file
bcftools query -f '%CHROM\.%POS\t[%GT\t]\n' femalePT.5kb.vcf > femalePT.5kb.tab


# For the female dataset, changed genotypes as follows:
# renamed them to *Fcontignumber.position
# ./. to -
# 0/0 to a
# 0/1 to ab
sed 's/^HyS/\*HyS/g' femalePT.5kb.tab | \
	sed 's/\.\/\./-/g' | \
	sed 's/0\/0/a/g' | \
	sed 's/0\/1/ab/g' \
	> femalePT.5kb.genotypes4onemap.tab

# Now used awk to remove last two columns, then add D1.10 as second column
awk 'NF{NF-=2};1' < femalePT.5kb.genotypes4onemap.tab | \
	awk '$1 = $1 FS "D1.10"' > femalePT.5kb.genotypes4onemap.raw.noheader


# add header:
# ...which can be done by grabbing the CHR line from the vcf file, lopping off the last two lines, and then prepending that to the file, along with the data type info, etc.
IND="$(grep 'CHR' femalePT.5kb.vcf | sed 's/#CHR.*FORMAT\t//g' | sed 's/\tPARENT.*WT//g' | sed 's/\t/ /g')"

# now determine the number of markers in the file
NUMMARKER="$(wc -l femalePT.5kb.genotypes4onemap.raw.noheader | awk '{print $1}')"

printf "data type outcross \n90 $NUMMARKER 0 0 0\n$IND\n" > header.txt

cat header.txt femalePT.5kb.genotypes4onemap.raw.noheader > femalePT.5kb.raw


##########
#  MALE  #
##########

# thin to 5kb
vcftools --gzvcf malePT.vcf.gz --thin 5000 --recode --recode-INFO-all --out malePT.5kb

mv malePT.5kb.recode.vcf malePT.5kb.vcf

# Create onemap file
bcftools query -f '%CHROM\.%POS\t[%GT\t]\n' malePT.5kb.vcf > malePT.5kb.tab


# For the male dataset, changed genotypes as follows:
# renamed them to *Fcontignumber.position
# ./. to -
# 0/0 to a
# 0/1 to ab
sed 's/^HyS/\*HyS/g' malePT.5kb.tab | \
	sed 's/\.\/\./-/g' | \
	sed 's/0\/0/a/g' | \
	sed 's/0\/1/ab/g' | \
        sed 's/1\/1/a/g' \
	> malePT.5kb.genotypes4onemap.tab

# Now used awk to remove last two columns, then add D1.10 as second column
awk 'NF{NF-=2};1' < malePT.5kb.genotypes4onemap.tab | \
	awk '$1 = $1 FS "D2.15"' > malePT.5kb.genotypes4onemap.raw.noheader


# add header:
# ...which can be done by grabbing the CHR line from the vcf file, lopping off the last two lines, and then prepending that to the file, along with the data type info, etc.
IND="$(grep 'CHR' malePT.5kb.vcf | sed 's/#CHR.*FORMAT\t//g' | sed 's/\tPARENT.*WT//g' | sed 's/\t/ /g')"

# now determine the number of markers in the file
NUMMARKER="$(wc -l malePT.5kb.genotypes4onemap.raw.noheader | awk '{print $1}')"

printf "data type outcross \n90 $NUMMARKER 0 0 0\n$IND\n" > header.txt

cat header.txt malePT.5kb.genotypes4onemap.raw.noheader > malePT.5kb.raw


#clean up
rm *noheader *tab *log header.txt
