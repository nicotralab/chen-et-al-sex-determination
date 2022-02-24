#!/usr/bin/perl
# Matt Nicotra
# identify_redundant_markers.pl
#
#Version: 01 - 2019-11-26 for linux
#

use warnings;
use strict;
use Cwd;
use Data::Dumper;

# change into current working directory
my $cwd = getcwd();
chdir "$cwd";

# check whether the appropriate user options are there
unless( @ARGV == 2 ){
	print "USAGE: identify_redundant_markers.pl binsfile mapfile\n";
    print "      This script creates three files:\n";
    print "        - [binsfile]nonredundant.tab = summary of which contigs are in each bin\n";
    print "        - markerswithdifferentcontigs.txt = list of markers that are at the same cM position, but are from different contigs\n";
    print "        - redundantmarkers.txt = markers that are at the same cM position and from same contig\n";
	exit;
}

# read bins file and create a hash with each marker and its constituent contigs



#create hash for the bin and contig info
my %binhash = ();


 



open( my $binsfilehandle, '<', $ARGV[0])
	or die "Could not open binsfile: $ARGV[0]\n";

while (my $line = <$binsfilehandle>){

	chomp $line;
	
	#ignore lines that do not have marker info
	if( $line =~ /^\$/ or $line =~/^\:/ or $line =~/^\s/ or $line=~/^\-/ or $line =~ /^\[/ or $line =~ /^at/ or $line =~/^$/){
		#print "skipped line: $line \n";
		next;
	}

	#extract contig name and bin number from the line
	#print "now doing $line\n";
	
	$line =~ /^(\S+?)\s+?(\d+?)\s+?\d+?/;
	my $currentmarker = $1;
	my $currentbin = $2;

	$currentmarker =~ /(^\S+?)\./;  # this gets rid of the bp position on the marker name
	my $currentcontig = $1;

	#print "bin = $currentbin, contig = $currentcontig, marker = $currentmarker\n";
	#<STDIN>;
	
	# we can take advantage of the fact that the markers within the bins file are sorted by marker name.
	# Therefore we can just add contig names onto the string in the hash.
	# First, test whether this bin exists
	
	if( exists $binhash{$currentbin} ){
			
		#print "this is not a new bin\n";
		
		# if the bin exists, test whether this contig exists in this bin
		my $isnewcontig = 1;
		foreach my $item ( @{$binhash{$currentbin}} ){
			if( $item eq $currentcontig ){
				$isnewcontig = 0;
			}
		}
		
		# if the contig is new to the bin, push it to the end of the bin
		if( $isnewcontig ){
			push( @{$binhash{$currentbin}}, $currentcontig );
		}
	
	} else {
	
		#print "This is a new bin\n";
		${$binhash{$currentbin}}[0] = $currentcontig;
	
	}
		

}

#close( $binsfilehandle );
#<STDIN>;
#print Dumper(\%binhash);
#<STDIN>;




# now print a summary file of the bins
my $outfilename = $ARGV[0]."nonredundant.tab";
open( my $outfh, '>', $outfilename);

foreach my $key ( keys %binhash ){

	#print $key;
	print $outfh "$key";
	foreach my $element ( @{$binhash{$key}} ){
		print $outfh "\t", $element;
	}
	print $outfh "\n";
}

close( $outfh );


# now read map file into an array
open( my $mapfh, '<', $ARGV[1])
	or die "Could not open mapfile: $ARGV[1]\n";

my @lines = <$mapfh>;

close( $mapfh );

chomp( @lines );


my $counter = 0;

open( my $redundantmarkerfh, '>', "redundantmarkers.txt");
open( my $markerswithdiffcontigsfh, '>', "markerswithdifferentcontigs.txt");

while( $counter <= @lines ){

	if( $lines[$counter] =~ /^\[/ ){
		print $redundantmarkerfh "Chromosome $lines[$counter]\n";
		print $markerswithdiffcontigsfh "Chromosome $lines[$counter]\n";
		$counter = $counter + 6;
	}
	
	
	
	# get this line and the next line
	my $currentline = $lines[$counter];
	my $nextindex = $counter + 1;
	my $nextline = $lines[$nextindex];
	
	print "Currentline: $currentline\n";
	print "nextline:    $nextline\n";
	
	if( $currentline =~ /(\d+)\s+\S+\s+(\S+)/ and $nextline =~ /(\d+)\s+\S+\s+(\S+)/ ){

		$currentline =~ /(\d+)\s+\S+\s+(\S+)/;
		my $currentbin = $1;
		my $currentpos = $2;
		
		$nextline =~ /(\d+)\s+\S+\s+(\S+)/;
		my $nextbin = $1;
		my $nextpos = $2;
		
		if( $currentpos eq $nextpos ){
		
			print "These markers are redundant\n";

			# now test if they have the same contigs
			my @contigsincurrentbin = sort( @{$binhash{$currentbin}} );
			my @contigsinnextbin = sort ( @{$binhash{$nextbin}} );
			
			print @contigsincurrentbin, "\n";
			print @contigsinnextbin, "\n";
			
			if( @contigsincurrentbin ~~ @contigsinnextbin ){
			
				print "These markers have same contigs\n";
				print $redundantmarkerfh $currentbin, "\n";
			
			} else {
			
				print "These markers do not have same contigs\n";
				print $markerswithdiffcontigsfh $currentbin, "\t", $nextbin, "\n";
				
			}
		

	
	}
	
	
	}
	
	$counter++;
}


close( $redundantmarkerfh, $markerswithdiffcontigsfh );













