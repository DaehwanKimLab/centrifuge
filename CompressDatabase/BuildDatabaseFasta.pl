#!/bin/perl

# Read the genomes organized according to a nice directory structure.
# Then merge the sequence for the chosen level

use strict ;
use threads ;

my $usage = "perl a.pl path_to_the_superkingdom [-level species -o compressed.fa -bss ./ -t 1]\n" ;

my $level = "species" ;
my $path = $ARGV[0] ;
my $i ;
my $output = "compressed.fa" ;
my $bssPath = "./" ;
my $numOfThreads = 1 ;
die $usage if ( scalar( @ARGV ) == 0 ) ;

for ( $i = 1 ; $i < scalar( @ARGV ) ; ++$i )
{
	if ( $ARGV[$i] eq "-level" )
	{
		$level = $ARGV[$i + 1] ;
		++$i ;
	}
	elsif ( $ARGV[$i] eq "-o" )
	{
		$output = $ARGV[$i + 1] ;
		++$i ;
	}
	elsif ( $ARGV[$i] eq "-bss" )
	{
		$bssPath = $ARGV[$i + 1] ;
		++$i ;
	}
	elsif ( $ARGV[$i] eq "-t" )
	{
		$numOfThreads = $ARGV[$i + 1] ;
		++$i ;
	}
	else
	{
		die $usage ;
	}
}


# recursively go through the directories.
sub solve
{
	my $pwd = $_[0] ;
	#print "$pwd\n" ;
	if ( $pwd =~ /$level/ ) # reached the desired level
	{
		# Extracts the files
		#print "find $pwd -name *.fna > tmp.list\n" ;
		my $tid = threads->tid() - 1 ;
		`find $pwd -name *.fna > tmp_$tid.list` ;
		#`find $pwd -name *.fa >> tmp.list` ; # Just in case
		# Build the header
		my $genusId ;
		my $speciesId ;
		my $speciesName ;
		my $paths = `head -1 tmp_$tid.list` ;
		#print $paths ;	
		if ( $paths =~ /genus_(.*?)\// )
		{
			$genusId = $1 ;
		}
		else
		{
			return ;
		}

		if ( $paths =~ /species_(.*?)_(.*?)\// )
		{
			$speciesId = $2 ;
		}
		else
		{
			return ;
		}

		if ( $paths =~ /species_.*\/(.*?)\.fna$/ )
		{
			$speciesName = ( split /_/, $1 )[0]."_". ( split /_/, $1 )[1] ;
		}
		else
		{
			return ;
		}

		return if ( ( $genusId * 17 + $speciesId ) % $numOfThreads != $tid ) ; 

		my $id = ( $speciesId << 32 ) | $genusId ;
		my $header = ">$id $speciesName" ;
		print $header, "\n" ;
		#return ;
		# Build the sequence
		`perl $bssPath/BuildSharedSequence.pl tmp_$tid.list -prefix tmp_${tid}_$id` ;
		# Merge all the fragmented sequence into one big chunk.
		`cat tmp_${tid}_${id}_*.fa > tmp_${tid}_$id.fa` ;

		open FP1, "tmp_${tid}_$id.fa" ;
		my $seq = "" ;
		while ( <FP1> )
		{
			chomp ;
			next if ( /^>/ ) ;
			$seq .= $_ ;
		}
		close FP1 ;

		open fpOut, ">>${output}_${tid}" ;
		print fpOut "$header\n$seq\n" ;
		close fpOut ;

		`rm tmp_${tid}_*` ;	
		return ;
	}

	my $ls = `ls -d $pwd/*` ;

	return if ( $ls eq "" ) ; 
	my @dirs = split /\s+/, $ls ;
	my $name ;
	foreach $name (@dirs)
	{
		solve( "$name" ) ;
	}
}


sub threadWrapper
{
	my $tid = threads->tid() - 1 ;
	`rm ${output}_${tid}` ;
	my $superKingdom = `ls  -d $path/superkingdom_*` ;
	my @array = split /\s+/, $superKingdom ;
	my $name ;
	foreach $name (@array)
	{
		solve( "$name" ) ;
	}
}

my @threads ;
for ( $i = 0 ; $i < $numOfThreads ; ++$i )
{
	push @threads, $i ;
}

foreach (@threads)
{
	$_ = threads->create( \&threadWrapper ) ;
}

foreach (@threads)
{
	$_->join() ;
}

# merge the files generated from each threads
`cat ${output}_* > tmp_output.fa ; rm ${output}_*` ;

# Remove the Ns from the file
`perl $bssPath/RemoveN.pl tmp_output.fa > $output` ;

`rm tmp_*` ;
