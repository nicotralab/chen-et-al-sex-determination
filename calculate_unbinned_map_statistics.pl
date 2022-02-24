#!/usr/bin/perl
# Matt Nicotra
# calculate_unbinned_map_statistics.pl
#
#Version: 01 - 2020-01-07 for linux
#

use warnings;
use strict;
use Cwd;
use Data::Dumper;

# change into current working directory
my $cwd = getcwd();
chdir "$cwd";


# check whether the appropriate user options are there
unless( @ARGV == 1 ){
	print "USAGE: calculate_unbinned_map_statistics.pl mapfile\n\n";
	exit;
}

#make sure there are the correct EOL characters.
my $answer = "";
print "Your text files must have Unix-style EOL characters. Is it correctly formatted? (y/n)\n";

do {
	$answer = <STDIN>;
	chomp($answer);

	if( $answer eq "N" ){
		print "Please reformat the data.\n";
		exit;
	} elsif( $answer ne "Y" ){
		print "Please answer \"y\" or \"n\"\n";
	}
} until( $answer eq "y" | $answer eq "n");


#create hash for the bin and contig info
my %contighash = ();

#the hash will be an hash of arrays with a key for each contig and the value an array with the same number of elements as the number of linkage groups.
# format of this hash will be
# %binhash = (
	# contig1  => [value1, value2... valuei] (where i = number of linkage groups and value is the count of markers in that linkage group)
	# contig2  => [value1, value2... valuei] (where i = number of linkage groups and value is the count of markers in that linkage group)
# )


# open the map file
open( my $mapfh, '<', $ARGV[0])
	or die "Could not open mapfile: $ARGV[1]\n";

#determine number of the largest linkage group
print "What is the number of the LARGEST linkage group in the map?(answer must be an integery\n";
do{
	$answer = <STDIN>;
	chomp($answer);
} until( $answer =~ /\d+/);

while( my $line = <$mapfh> ){

	chomp($line);
	#print $line, "\n";

	if( $line =~ /^\d/ ){
	

		$line =~ /(^\d+?) (.+?) (.*)/;
		my $currentlg = $1;
		my $currentmarker = $2;
		my $currentpos = $3;
		
		#print $currentlg, "\n";
		#print $currentmarker, "\n";
		
		$currentmarker =~ /(.*?)\.(.*?)/;
		my $currentcontig = $1;
		#print $currentcontig, "\n";
		#<STDIN>;
		
		
		# populate the hash with the data from this marker.
		if( exists( $contighash{$currentcontig} )){
			
			#print "EXISTS";
			#add one to the value for that linkage group
			$contighash{$currentcontig}[$currentlg] = $contighash{$currentcontig}[$currentlg] + 1;
			#print "value is now", $contighash{$currentcontig}[$currentlg], "\n";
		
		} else {
		
			#print "DOES NOT EXIST";
			#populate the array for this one.
			for( my $i=0; $i <= $answer; $i++ ){
				push( @{$contighash{$currentcontig}}, 0 );
			}
			$contighash{$currentcontig}[$currentlg] = $contighash{$currentcontig}[$currentlg] + 1;			
			#<STDIN>;
		}
	}
	
}


# summarize by linkage groups
@LGs
@LGs = ();















############################
# OUTPUT
############################

#print matrix of all contigs by all linkage groups

open( my $outfh1, '>', "$ARGV[0]\.matrix.tab");

my @sortedcontignames = sort( keys( %contighash ));

# print column headers
print $outfh1 "Contig Name";
for( my $i=1; $i <= $answer; $i++ ){
	print $outfh1 "\t", $i;
}
print $outfh1 "\n";

foreach my $contigname ( @sortedcontignames ){

	print $outfh1 "$contigname";
		for( my $i=1; $i <= $answer; $i++ ){
			print $outfh1 "\t", $contighash{$contigname}[$i];
		}
	print $outfh1 "\n";

}
close( $outfh1 );


# print summary 