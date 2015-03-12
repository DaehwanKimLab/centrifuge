#!/bin/perl

# Read .
# Then merge the sequence for the chosen level

use strict ;
use threads ;
use threads::shared ;

my $usage = "perl a.pl path_to_download_files path_to_taxnonomy [-o compressed.fa -bss ./ -t 1]\n" ;

my $level = "species" ;
my $path = $ARGV[0] ;
my $taxPath = $ARGV[1] ;
my $i ;
my $output = "compressed.fa" ;
my $bssPath = "./" ;
my $numOfThreads = 1 ;

my %gidToFile ;
my %gidUsed ;
my %tidToGid ; # each tid can corresponds several gid
my %gidToTid ;
my %speciesList ; # hold the tid in the species
my %species ; # hold the species tid
my %genus ; # hold the genus tid
my %speciesToGenus ;
my %fileUsed : shared ;
my $fileUsedLock : shared ;
my %taxTree ;

my @speciesListKey ;
my $speciesUsed : shared ;
my $debug: shared ;

die $usage if ( scalar( @ARGV ) == 0 ) ;

for ( $i = 2 ; $i < scalar( @ARGV ) ; ++$i )
{
	if ( $ARGV[$i] eq "-level" )
	{
		# No use now.
		#$level = $ARGV[$i + 1] ; 
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

#Extract the gid we met in the file
`find $path -name *.fna > tmp_fna.out` ;
#print "find $path -name *.fna > tmp_fna.out\n" ;
open FP1, "tmp_fna.out" ;
while ( <FP1> )
{
	chomp ;
	my $file = $_ ;
	open FP2, $file ;
	my $head = <FP2> ;
	close FP2 ;
	#print $file, "\n" ;
	my $gid ;
	if ( $head =~ /^>gi\|([0-9]+)?\|/ )
	{
		$gid = $1 ;
		print "$gid $file\n" ;
		$gidUsed{ $gid } = 1 ;
		$gidToFile{ $gid } = $file ;
		$fileUsed{ $file } = 0 ;
	}
	else
	{
		die "$file has wrong header.\n"
	}
}
close FP1 ;
#print "1 over\n" ;
# Extract the tid that are associated with the gids
open FP1, "$taxPath/gi_taxid_nucl.dmp" ;
while ( <FP1> )
{
	chomp ;
	my @cols = split ;		
	#print $cols[0], "\n" ;	
	if ( defined( $gidUsed{ $cols[0] } ) )
	{
		push @{ $tidToGid{ $cols[1] } }, $cols[0] ;	
		print $cols[0], " ", $cols[1], " ", $gidToFile{ $cols[0] }, "\n" ;
	}
}
close FP1 ;
#print "2 over\n" ;
# Organize the tree
open FP1, "$taxPath/nodes.dmp" ;
while ( <FP1> )
{
	chomp ;

	my @cols = split ;
	#next if ( !defined $tidToGid{ $cols[0] } ) ;
	
	my $tid = $cols[0] ;
	my $parentTid = $cols[2] ;
	my $rank = $cols[4] ;
	#print "subspecies: $tid $parentTid\n" ;
	#push @{ $species{ $parentTid } }, $tid ;
	#$tidToSpecies{ $tid } = $parentTid ;

	$taxTree{ $cols[0] } = $cols[2] ;
	#print $cols[0], "=>", $cols[2], "\n" ;
	if ( $rank eq "species" )	
	{
		$species{ $cols[0] } = 1 ;
	}
	elsif ( $rank eq "genus" )
	{
		$genus{ $cols[0] } = 1 ;
	}
}
close FP1 ;

# Organize put the sub-species tid into the corresponding species.
for $i ( keys %tidToGid )
{
	my $p = $i ;
	my $found = 0 ;
	while ( 1 )
	{
		last if ( $p <= 1 ) ;
		if ( defined $species{ $p } ) 
		{
			$found = 1 ;
			last ;
		}
		if ( defined $taxTree{ $p } )
		{
			$p = $taxTree{ $p } ;
		}
		else
		{
			last ;
		}
	}

	if ( $found == 1 )
	{
		push @{ $speciesList{ $p } }, $i ;
	}
}

#exit ;

# Compress one species.
sub solve
{
# Extracts the files
#print "find $pwd -name *.fna > tmp.list\n" ;
	my $tid = threads->tid() - 1 ;
	#`find $pwd -name *.fna > tmp_$tid.list` ;

#`find $pwd -name *.fa >> tmp.list` ; # Just in case
# Build the header
	my $genusId ;
	my $speciesId = $_[0] ;
	my $speciesName ;
	my $i ;
	my $file ;
	my @subspeciesList = @{ $speciesList{ $speciesId } } ;
	$genusId = $taxTree{ $speciesId } ;
	while ( 1 )
	{
		if ( !defined $genusId || $genusId <= 1 )
		{
			$genusId = 0 ;
			last ;
		}

		if ( defined $genus{ $genusId } )
		{
			last ;
		}
		else
		{
			$genusId = $taxTree{ $genusId } ;
		}
	}

	my $FP1 ;
	open FP1, ">tmp_$tid.list" ;
	foreach $i ( @subspeciesList )
	{
		foreach my $j ( @{$tidToGid{ $i } } )
		{
			{
				lock( $debug ) ;
				print "Merge ", $gidToFile{ $j }, "\n" ;
			}
			$file =  $gidToFile{ $j } ;
			{
				lock( $fileUsedLock ) ;
				$fileUsed{ $file } = 1 ;
			}
			print FP1 $file, "\n" ;
		}
	}
	close FP1 ;
		
	#print $file, "\n" ;	
	if ( $file =~ /\/(\w*?)uid/ )
	{
		$speciesName = ( split /_/, $1 )[0]."_". ( split /_/, $1 )[1] ;
	}
	else
	{
		$speciesName = "" ;
	}
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
}


sub threadWrapper
{
	my $tid = threads->tid() - 1 ;
	`rm ${output}_${tid}` ;

	while ( 1 )
	{
		my $u ;
		{
		lock $speciesUsed ;
		$u = $speciesUsed ;
		++$speciesUsed ;
		}
		last if ( $u >= scalar( @speciesListKey ) ) ;
		solve( $speciesListKey[ $u ] ) ;
	}
}


my @threads ;
@speciesListKey = keys %speciesList ; 
$speciesUsed = 0 ;
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

#Merge the unused files
open FP2, ">>", "tmp_output.fa" ;
my $seq = "" ;
foreach $i ( keys %fileUsed )
{
	if ( $fileUsed{ $i } == 0 )
	{
#print $i, "\n" ;
#`cat $i >> tmp_output.fa` ;
		my $speciesName ;
		open FP1, $i ;

		print "Unused file: $i\n" ;
		while ( <FP1> )
		{
			chomp ;
			next if ( /^>/ ) ;
			$seq .= $_ ;
		}
		
		close FP1 ;
	}
}
print FP2 ">0 unknown_taxno_genomes\n$seq\n" ;
close FP2 ;

# Remove the Ns from the file
`perl $bssPath/RemoveN.pl tmp_output.fa > $output` ;

`rm tmp_*` ;
