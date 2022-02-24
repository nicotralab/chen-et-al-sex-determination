# This script was written by Justin Paschall of the NHGRI at NIH.
# This script was modified by Matt Nicotra to deal with an error caused by values of AD for which there was not data (i.e. a "." instead of "number,number"
# Specifically, we added an if than statement to reject variants for which this was true.
# Takes as input a file in VCF format and outputs each line from that file after removing genotypes that do not meet filtering criteria.


import sys

for line in sys.stdin:
  if( not line.startswith("#")):
    line = line.rstrip()
    linelist = line.split("\t")
    #print "Processing variant "
    reject = 0


    if("," in linelist[4]):
      reject = 1  

    for geno in linelist[9:]:
      genolist = geno.split(":")
      GT = genolist[0]
      AD = genolist[1]
      DP = genolist[2]
      GQ = genolist[3]
      #print AD
      alleles = AD.split(",")
      if( len(alleles) > 1 ): 
         allele1 = int(alleles[0])
         allele2 = int(alleles[1])

         if( (GT == "0/0" or GT == "0|0") and (allele2 > 2 or allele1 < 10)):
            #print "geno error 0/0: " + str(genolist)
            reject = 1

         if( (GT == "1/1" or GT == "1|1") and (allele1 > 2 or allele2 < 10)):
            #print "geno error 0/0: " + str(genolist)
            reject = 1

         total_alleles = allele1 + allele2
      
         if((allele1+allele2) > 0):
            allele2perc = float(allele2) / float(allele1+allele2)
 
         if( ( (allele1+allele2) > 0 ) and ( (GT == "0/1" or GT == "0|1") and ((allele2 < 10 or allele1 < 10) or (  allele2perc < 0.3 or allele2perc > 0.7 )) )) :
            #print "geno error 0/1: " + str(genolist) + str(allele2perc)
            reject = 1
    
      else: 
         reject = 1

    if not reject:
      print line
       
       

    #print line
