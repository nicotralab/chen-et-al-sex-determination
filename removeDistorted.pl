#!usr/bin/perl
#
# Matt Nicotra
# removeDistorted.pl
#
# Version: 03 for linux
#
#
#

use warnings;
use strict;
use Cwd;
use Statistics::Distributions;

# change into current working directory
my $cwd = getcwd();
chdir "$cwd";

# check whether the appropriate user options are there
unless( @ARGV == 12 ){
	print "USAGE: removeDistorted.pl --vcf <vcf file name> --parent1 <string> --parent2 <string> --out-file <string> --siglevel <float> --bonf <true,false>\n";
	print "All command line options must be used in this order.\n";	
	exit;
}


# initialize some global variables
my $numsamples = 0;
my $parent1position = 0;
my $parent2position = 0;
my $nondistortedcount = 0;

# determine number of variants
my $numvariants = `grep -cv '^#' $ARGV[1]`;
chomp( $numvariants );

# read user supplied commands into variables
my $filename = $ARGV[1];
print "File name is $filename.\n";
my $parent1 = $ARGV[3];
print "Parent 1 is $parent1.\n";
my $parent2 = $ARGV[5];
print "Parent 1 is $parent2.\n";
my $outfilename = $ARGV[7];
print "Output file name is $outfilename.\n";
my $siglevel = $ARGV[9];
print "P-value cutoff is $siglevel.\n";
my $bonf = $ARGV[11];
print "bonf = $bonf\n";
if( $bonf eq "true" ){
	print "Bonferroni correction will be used.\n";
	$siglevel = $siglevel / $numvariants;
	print "P-value cutoff now $siglevel.\n";
} elsif( $bonf eq "false" ){
	print "Bonferroni correction will not be used.\n";
} else {
	print "ERROR: Bad value for --bonf.\n";
	exit;
}


# open the VCF file for reading
open(my $fh, '<:encoding(UTF-8)', $filename)
  or die "Could not open file '$filename' $!";

  
#open the output file for writing
open(my $outfh, '>', $outfilename)
  or die "Could not create file '$outfilename' $!";
  
while (my $row = <$fh>) {

	# if it's a comment, write it to the output file and go to next line
	if( $row =~ /^##/ ){
		print $outfh $row;
		next;
	}

	# determine number of samples and position of the parents
	if( $row =~ /#CHROM/ ){
		print $outfh $row;
	
		#determine the number of samples
		my @fields = split( "\t", $row );
		chomp( @fields );
		
		my $numfields = @fields;
		$numsamples = $numfields - 9;
		
		print "There are $numsamples samples.\n";
	
		# determine position of each parent.
		for( my $i=0; $i< $numfields; $i++ ){
			if( $fields[$i] =~ $parent1 ){ 
				$parent1position = $i;
			}
		}

		for( my $i=0; $i< $numfields; $i++ ){
			if( $fields[$i] =~ $parent2 ){ 
				$parent2position = $i;
			}
		}

		print "Parent 1 is $parent1 and is in position $parent1position.\n";
		print "Parent 2 is $parent2 and is in position $parent2position.\n";
		
		next;
	} 
	
	#get the line
	#print "The next line is\n: $row\n"; # This is only for debugging
	
	# Turn line into array, then get genotypes of parents
	my @fields = split( "\t", $row);
	chomp( @fields );
	
	# store genotype of parent1 and parent2
	my $parent1GT = substr( $fields[$parent1position], 0, 3);
	$fields[$parent1position] = "parent1";
	my $parent2GT = substr( $fields[$parent2position], 0, 3);
	$fields[$parent2position] = "parent2";

	
	
	my $markertype = "none";
	# determine marker type:
	if( $parent1GT eq "0/1" && $parent2GT eq "0/0" ){
		$markertype = "oneHetOneRef";
	} elsif( $parent1GT eq "0/0" && $parent2GT eq "0/1") {
		$markertype = "oneHetOneRef";
	} elsif( $parent1GT eq "0/1" && $parent2GT eq "1/1") {
		$markertype = "oneHetOneAlt";
	} elsif( $parent1GT eq "1/1" && $parent2GT eq "0/1") {
		$markertype = "oneHetOneAlt";
	} elsif( $parent1GT eq "0/1" && $parent2GT eq "0/1") {
		$markertype = "bothHet";
	}
		
	#print "\n\nGenotypes:", $parent1GT, " ", $parent2GT, "\n";
	#print "Marker type: $markertype\n";

	if( $markertype eq "oneHetOneRef") {
	
		#count genotypes
		my $aa = "0/0";
		my $ab = "0/1";
		my $aa_count = 0;
		my $ab_count = 0;
		foreach my $genotype ( @fields ){
			if( substr( $genotype, 0, 3) eq $aa ){ $aa_count++; }
			if( substr( $genotype, 0, 3) eq $ab ){ $ab_count++; }
		}
		#print "There are $aa_count aa genotypes and $ab_count ab genotypes.\n";
	
		#determine expected numbers
		my $expected_aa_count = ($aa_count + $ab_count)/2;
		my $expected_ab_count = $expected_aa_count;
	
		# calculate chisq value for one-way goodness of fit test
		my $chisqvalue = (($aa_count - $expected_aa_count )**2 / $expected_aa_count ) + (($ab_count - $expected_ab_count )**2 / $expected_ab_count );
		#print "Chi sq value is $chisqvalue\n";
		
		# calculate significance of this chi sq value
		my $probofchisqvalue = Statistics::Distributions::chisqrprob (1, $chisqvalue);
		#print "prob = $prob\n";
		
		# determine whether marker is distorted or not
		if( $probofchisqvalue > $siglevel ){
			print $outfh $row;
			$nondistortedcount++;
		}
		
	} elsif ( $markertype eq "oneHetOneAlt") {
	
		#count genotypes
		my $aa = "1/1";
		my $ab = "0/1";
		my $aa_count = 0;
		my $ab_count = 0;
		foreach my $genotype ( @fields ){
			if( substr( $genotype, 0, 3) eq $aa ){ $aa_count++; }
			if( substr( $genotype, 0, 3) eq $ab ){ $ab_count++; }
		}
		#print "There are $aa_count aa genotypes and $ab_count ab genotypes.\n";
	
		#determine expected numbers
		my $expected_aa_count = ($aa_count + $ab_count)/2;
		my $expected_ab_count = $expected_aa_count;
	
		# calculate chisq value for one-way goodness of fit test
		my $chisqvalue = (($aa_count - $expected_aa_count )**2 / $expected_aa_count ) + (($ab_count - $expected_ab_count )**2 / $expected_ab_count );
		#print "Chi sq value is $chisqvalue\n";
		
		# calculate significance of this chi sq value
		my $probofchisqvalue = Statistics::Distributions::chisqrprob (1, $chisqvalue);
		#print "prob = $prob\n";
		
		# determine whether marker is distorted or not
		if( $probofchisqvalue > $siglevel ){
			print $outfh $row;
			$nondistortedcount++;
		}
		
	} elsif ( $markertype eq "bothHet") {
	
		#count genotypes
		my $aa = "0/0";
		my $ab = "0/1";
		my $bb = "1/1";
		my $aa_count = 0;
		my $ab_count = 0;
		my $bb_count = 0;
		foreach my $genotype ( @fields ){
			if( substr( $genotype, 0, 3) eq $aa ){ $aa_count++; }
			if( substr( $genotype, 0, 3) eq $ab ){ $ab_count++; }
			if( substr( $genotype, 0, 3) eq $bb ){ $bb_count++; }
		}
	
		#print "There are $aa_count to $ab_count to $bb_count.\n";

	
		#determine expected numbers
		my $expected_aa_count = ($aa_count + $ab_count + $bb_count)/4;
		my $expected_bb_count = $expected_aa_count;
		my $expected_ab_count = 2 * $expected_aa_count;
	
		# calculate chisq value for one-way goodness of fit test
		my $chisqvalue = (($aa_count - $expected_aa_count )**2 / $expected_aa_count ) + (($ab_count - $expected_ab_count )**2 / $expected_ab_count ) + (($bb_count - $expected_bb_count )**2 / $expected_bb_count );
		#print "Chi sq value is $chisqvalue\n";
		
		# calculate significance of this chi sq value
		my $probofchisqvalue = Statistics::Distributions::chisqrprob (1, $chisqvalue);
		#print "prob = $prob\n";
		
		# determine whether marker is distorted or not
		if( $probofchisqvalue > $siglevel ){
			print $outfh $row;
			$nondistortedcount++;
		}
	}	
}

print "$nondistortedcount of $numvariants were not distorted in this file.\n";
