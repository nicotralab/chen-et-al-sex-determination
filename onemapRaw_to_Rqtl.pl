#!/usr/bin/perl
# Matt Nicotra
# onemapRaw_to_Rqtl.pl
#
#Version: 01 - 2020-01-03 for linux
#

use warnings;
use strict;
use Cwd;
use Data::Dumper;

# change into current working directory
my $cwd = getcwd();
chdir "$cwd";

# check whether the appropriate user options are there
unless( @ARGV == 3 ){
	print "USAGE: onemapRaw_to_Rqtl.pl rawfile mapfile qtlfilename\n";
	exit;
}

my $answer = "";
print "Your text files must have Unix-style EOL characters. Is it correctly formatted? (y/n)\n";

do {
	$answer = <STDIN>;
	chomp($answer);

	if( $answer eq "n" ){
		print "Please reformat the data.\n";
		exit;
	} elsif( $answer ne "y" ){
		print "Please answer \"y\" or \"n\"\n";
	}
} until( $answer eq "y" | $answer eq "n");






#create array for the marker names
my @markersinmap = ();


############################################################
# create an array with all of the raw data
############################################################

open( my $rawfh, '<', $ARGV[0])
	or die "Could not open rawfile: $ARGV[0]\n";

my @rawfilelines = <$rawfh>;
chomp( @rawfilelines );

close( $rawfh );



#######################################################################
# open the map file and create an array with lines for the marker data
#######################################################################

open( my $mapfh, '<', $ARGV[1])
	or die "Could not open mapfile: $ARGV[1]\n";

my @linesforRqtl = ();

while (my $line = <$mapfh>){

	if( $line =~/^.* .* .*/ ) {
		
		chomp( $line );
		
		#print $line."\n";
		
		$line =~ /(.*) (.* )(.*)/;
		my $chromosome = $1;
		my $markerid = $2;
		my $position = $3;
	
		
		my @matches = grep( /$markerid/, @rawfilelines );
		
		if( @matches > 1 ){
			print "ERROR: more than one marker with name $markerid in rawfile\n";
			exit;
		} elsif ( @matches < 1 ){
			print "ERROR: no marker with name $markerid in rawfile\n";
			exit;
		}
			
		my $linefromrawfile = $matches[0];

		#print $linefromrawfile, "\n";

		$linefromrawfile =~ /^.*? .*? (.*?)$/;
		
		my $genotypes = $1;

		# translate genotypes into Rqtl symbols
		$genotypes =~ s/ab/H/g;
		$genotypes =~ s/a/A/g;
		
		#remove that extra space on the marker id
		$markerid =~ s/ //g;
		
		my $lineforRqtl = "$markerid $chromosome $position $genotypes\n"; 
		
		$lineforRqtl =~ s/ /,/g;
		
		push( @linesforRqtl, $lineforRqtl );
		
	}

}




#######################################################################
# create file for Rqtl
#######################################################################
open( my $qtlfh, '>', $ARGV[2])
	or die "Could not create Rqtl datafile: $ARGV[2]\n";

# print header line
my $ids = $rawfilelines[2];
$ids =~ s/ /,/g;
print $qtlfh "IDS,,,", $ids, "\n";

foreach my $line ( @linesforRqtl ){ print $qtlfh $line};

close( $qtlfh );
