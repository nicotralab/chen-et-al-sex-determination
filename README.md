## Overview

This document describes how we generated linkage maps of the _Hydractinia symbiolongicarpus_ genome and used them to identify sex-linked markers. It provides large files that we could not include as supplemental data in the original publication, along with scripts and some other details.

In theory, you should be able to replicate our analysis perfectly. However, the details of running each analysis tool may differ depending on your platform and version. 


## Sequencing

1. The paternal genome assembly can be downloaded at **TODO: add link to paternal genome**. 

2. Male parent 291-10 was previously sequenced to high coverage. To generate a dataset of coverage equivalent to that obtained for the female parent and F1 progeny, we downsampled the original files (**TODO: add link to SRA for paternal genome files**) with seqtk as follows:

   ```bash
   $ seqtk sample -s100 160429_CRATOS_HKM5MBCXX.1.12268798.1.fq.gz 0.2 \
       | gzip > 160429_CRATOS_HKM5MBCXX.lane1.12268798.read1.DOWNSAMPLE_0.2.fq.gz
   
   $ seqtk sample -s100 160429_CRATOS_HKM5MBCXX.1.12268798.2.fq.gz 0.2 \
       | gzip > 160429_CRATOS_HKM5MBCXX.lane1.12268798.read2.DOWNSAMPLE_0.2.fq.gz
   
   ```

    The downsampled files can be downloaded from **TODO: upload data and add link**

3. Female parent 295-8 and F1 progeny were sequenced as described in the paper. Sequences can be downloaded from the SRA at **TODO: upload data and add link**



## Alignment of sequencing reads to the paternal genome

1. For each animal, including the parents, the Raw reads were mapped to an assembly of the paternal genome ([accession]()) with BWA-MEM with mapping parameters '-M -t 8'. The resulting .sam file was converted to .bam format  and then sorted with sam tools. An generic command is below:

   ```bash
   $ bwa mem -M -t 8 genomeAssembly.fa.gz \
       colonyXXX-fwdReads_R1_001.fastq.gz \
       colonyXXX-revReads_R2_001.fastq.gz \
       | samtools view -h -b - \
       | samtools sort -m 6G -@ 4 \
       > colonyXXX.sorted.bam
   ```

2. Use Picard to fix read groups

   ```bash
   $ java -jar $PICARDJARPATH/picard.jar AddOrReplaceReadGroups \
       I=colonyXXX.sorted.bam \
       O= colonyXXX.sorted.rg.bam \
       RGID=1 \
       RGLB=1 \
       RGPL=illumina \
       RGPU=1 \
       RGSM=colonyXXX.sorted.rg.bam
   ```

3. Use Picard to mark duplicates

   ```bash
   $ java -jar $PICARDJARPATH/picard.jar MarkDuplicates \
       I=colonyXXX.sorted.rg.bam \
       O=colonyXXX.sorted.rg.md.bam
   ```



## Calling Variants

1. Create indexes for GATK on reference:

   ```bash
   $ java -jar CreateSequenceDictionary.jar genomeAssembly.fa O=genomeAssembly.dict
   $ sametools faidx genomeAssembly.fa
   ```

2. Pre-call variants in each sample with GATK HaplotypeCaller

   ```bash
   $ java -jar $GATK_JAR \
     -T HaplotypeCaller \
     -R genomeAssembly.fa \
     -I colonyXXX.sorted.rg.md.bam \
     --emitRefConfidence GVCF \
     -o colonyXXX.out.rawsnps.indels.g.vcf
   ```

   

3. Create a textfile that lists filenames for pre-called variants for each sample. For example a file called ```pre-called.list```:

   ```
   colony001.out.rawsnps.indels.g.vcf
   colony002.out.rawsnps.indels.g.vcf
   colony003.out.rawsnps.indels.g.vcf
   ...
   colonyNNN.out.rawsnps.indels.g.vcf
   ```

   

4. Perform joint genotyping on the pre-called variants with GenotypeGVCFs:

   ```bash
   $ java -Xmx24g -jar $GATK_JAR \
     -T GenotypeGVCFs \
     -R genomeAssembly.fa \
     --variant pre-called.list \  
     -o rawvariants.90f1.vcf
   	
   ```



5. File ```rawvariants.90f1.vcf``` contains genotypes for all samples, including both parents in .vcf format. This file can be downloaded from **TODO: upload file and provide link**. 



## Quality Filtering Variants

We performed two rounds of variant filtering to ensure that only the most well-supported genotypes were used for linkage mapping.

#### Custom Filtering

We first used a custom python script [qualityfilter.py](https://github.com/nicotralab/chen-et-al-sex-determination/blob/main/qualityfilter.py) to ensure that:

1. No samples were missing data
2. All samples genotyped as homozygousall samples genotyped as homozygous reference (0/0) had no more than two mapped reads corresponding to the alternative allele and also had more than ten mapped reads corresponding to the reference
3. All samples genotyped as homozygous alternative (1/1) had no more than two mapped reads corresponding to the reference allele and also had more than ten mapped reads corresponding to the alternate allele
4. All samples genotyped as heterozygous (0/1) had an alternate allele read count percentage of greater than 0.3 or less than 0.7.
5. The result of this was a temporary file, referred to as ```pythonfiltered.vcf``` in the next section



#### Filtering according to GATK best practices

After custom filtering, we filtered variants according to GATK best practices. The process is outlined in the bash script [filterbyGATKBP.sh](https://github.com/nicotralab/chen-et-al-sex-determination/blob/main/filterbyGATKBP.sh). 

The output of this step is a file called ```GATKBP-passed.vcf```. It can be downloaded from **TODO: upload and add link**





## Identifying Pseudotestcross Variants

Next, we extracted variants (SNPs and indels) from ```GATKBP-passed.vcf``` that would meet the criteria for a pseudotestcross. The approach was to:

1. Extract variants with the correct genotype (i.e. heterozygous in one parent and homozygous in the other)
2. Identify and remove variants at which some animals have "unexpected" genotypes. These are
   - Animals that are homozygous for the wrong allele (e.g. 1/1 when the parents are 0/0 and 0/1)
   - Variants with 0/0:0/1:1/1 genotypes in a 1:2:1 ratio, indicating a likely misgenotype of one parent. 
3. Create two datafiles of variants, one to map recombination events in the male parent, and one to map recombination events in the female parent. 

This was achieved with the bash script [getPTvariants.sh](https://github.com/nicotralab/chen-et-al-sex-determination/blob/main/getPTvariants.sh).

The resulting datafiles, `GATKBP-passed.maleHet.abxabRemoved.vcf` and `GATKBP-passed.femaleHet.abxabRemoved.vcf` can be downloaded from **TODO upload file and provide link**



## Generating a linkage maps of the maternal and paternal genomes


#### Preparing data for OneMap

Prior to linkage mapping, the variants in each file were filtered to remove those that were distorted using the script [removedistorted.pl](https://github.com/nicotralab/chen-et-al-sex-determination/blob/main/removeDistorted.pl). The command used was:

```bash
$ removeDistorted.pl --vcf GATKBP-passed.maleHet.abxabRemoved.vcf \
    --parent1 PARENT_295_8 \
    --parent2 WT \
    --out-file GATKBP-passed.maleHet.abxabRemoved.nodist-p-00001.vcf \
    --siglevel 0.00001 \
    --bonf true

$ removeDistorted.pl --vcf GATKBP-passed.femaleHet.abxabRemoved.vcf \
    --parent1 PARENT_295_8 \
    --parent2 WT \
    --out-file GATKBP-passed.femaleHet.abxabRemoved.nodist-p-00001.vcf \
    --siglevel 0.00001 \
    --bonf true
```



These were then thinned to 1 variant per 5 kb and converted from .vcf to OneMap's with vcftools and command line tools via [thin-for-onemap.v2.sh](https://github.com/nicotralab/chen-et-al-sex-determination/blob/main/thin-for-onemap.v2.sh).

The resulting files were [femalePT.5kb.raw](https://github.com/nicotralab/chen-et-al-sex-determination/blob/main/femalePT.5kb.raw) and [malePT.5kb.raw](https://github.com/nicotralab/chen-et-al-sex-determination/blob/main/malePT.5kb.raw) and are also provided as supplemental datasets in Chen et al. 



#### Linkage mapping in OneMap

R was run in R studio. Each linkage map (paternal and maternal) was created in a separate Rstudio project:
- [Rstudio-femaleMap](https://github.com/nicotralab/chen-et-al-sex-determination/blob/main/Rstudio-femaleMap.zip)
- [Rstudio-maleMap](https://github.com/nicotralab/chen-et-al-sex-determination/blob/main/Rstudio-maleMap.zip)

A full description of the analysis is documented in the R script `create_maps.r`, which is provided in the project folder.

Note that part of the analysis must be run with the perl scripts [identify_redundant_markers.pl](https://github.com/nicotralab/chen-et-al-sex-determination/blob/main/identify_redundant_markers.pl), [unbin_markers_in_map.pl](https://github.com/nicotralab/chen-et-al-sex-determination/blob/main/unbin_markers_in_map.pl), and [calculate_unbinned_map_statistics.pl](https://github.com/nicotralab/chen-et-al-sex-determination/blob/main/calculate_unbinned_map_statistics.pl). 


## QTL analysis

Loci linked to sex were identified in R using R/qtl. The entire analysis can be recreated by downloading the Rstudio project and following the R script `QTLanalysis_binaryTraits.R`





