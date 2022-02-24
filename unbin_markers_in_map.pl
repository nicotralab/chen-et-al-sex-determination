#!/usr/bin/perl
# Matt Nicotra
# unbin_markers_in_map.pl
#
#Version: 01 - 2019-12-13 for linux
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
	print "USAGE: unbin_markers_in_map.pl binsfile mapfile\n";
	exit;
}


#create hash for the bin and contig info
my %binhash = ();

#the hash will be an hash of arrays, with each element of the array being a marker ID from the bin. 
# format of this hash will be
# %binhash = (
	# 1  => [F0001.1, F0001.2, F0001.3]
	# 2  => [F0001.4, F0001.5, F0001.6]
# )


open( my $binsfilehandle, '<', $ARGV[0])
	or die "Could not open binsfile: $ARGV[0]\n";

while (my $line = <$binsfilehandle>){

	chomp $line;
	
	#print "$line\n";
	#<STDIN>;
	
	# if this is the beginning of a bin, create an entry in the hash that has all of the constituent markers
	if( $line =~/^\:/ ){
		
		$line = <$binsfilehandle>;
		$line = <$binsfilehandle>;
		
		do {
		
			$line =~ /^(\S+?)\s+?(\d+?)\s+?\d+?/;
			my $currentmarker = $1;
			my $currentbin = $2;

			push( @{$binhash{$currentbin}}, $currentmarker );
		
			$line = <$binsfilehandle>;
		
		} until ( $line =~ /^\-\-\-/ or $line =~ /\$info/ );

	}
	
}



close( $binsfilehandle );


# because this is not an elegant script, we end up with an undefined value in the hash. This is removed here:
delete($binhash{''});


#for error checking
#print Dumper(\%binhash);
#<STDIN>;




# now print a new map file with all of the unbinned markers
my $outfilename = $ARGV[1].".unbinned.txt";
open( my $outfh, '>', $outfilename);


# now open the map file and use it to print the new map file. 
open( my $mapfh, '<', $ARGV[1])
	or die "Could not open mapfile: $ARGV[1]\n";

while( my $line = <$mapfh> ){

	if( $line =~ /^\d/ ){
	
		#print $line;

		$line =~ /(^\d+?) (.+?) (.*)/;
		my $currentlg = $1;
		my $currentmarker = $2;
		my $currentpos = $3;
		
		#print $currentlg, "\n";
		#print $currentmarker, "\n";
		#print $currentpos, "\n";
		#<STDIN>;

		
		#determine which bin this marker is from
		my @bins = keys %binhash;
		
		
		my $currentbin = '';
		
		foreach my $bin ( @bins ){
		
			#print "Now looking at bin $bin\n";
			#<STDIN>;
			foreach my $markerinbin ( @{$binhash{$bin}} ) {
			
				if( $markerinbin eq $currentmarker ){
				
					#print "It Matched bin $bin\n";
					$currentbin = $bin;
				}
			}
		}
		#print "Current lg =  $currentlg\n";
		#print "Current marker =  $currentmarker\n";
		#print "Current pos =  $currentpos\n";
		#print "Current Bin = $currentbin \n";
		#<STDIN>;
		
		foreach my $binnedmarker ( @{$binhash{$currentbin}} ){
		
			print $outfh "$currentlg $binnedmarker $currentpos\n";

		}
	
	}
}

close( $outfh );
close( $mapfh );
