#!/usr/bin/perl

# Read and merge the sequence for the chosen level

use strict ;
use warnings ;

use threads ;
use threads::shared ;
use FindBin qw($Bin);
use File::Basename;
use File::Find;
use Getopt::Long;
use Cwd;

my $usage = "USAGE: perl ".basename($0)." path_to_download_files path_to_taxnonomy [-o compressed.fa -bss . -noCompress -t 1 -maxG -1]\n" ;

my $level = "species" ;
my $output = "compressed.fa" ;
my $bssPath = $Bin; # take path of binary as script directory
my $numOfThreads = 1 ;
my $noCompress = 0 ;
my $verbose = 0;
my $maxGenomeSizeForCompression = -1 ;

GetOptions ("level|l=s" => \$level,
			"output|o=s" => \$output,
			"bss=s" => \$bssPath,
			"threads|t=i" => \$numOfThreads,
			"maxG=i" => \$maxGenomeSizeForCompression, 
            "verbose|v" => \$verbose,
			"noCompress" => \$noCompress)
or die("Error in command line arguments. \n\n$usage");

die $usage unless @ARGV == 2;

my $path = $ARGV[0] ;
my $taxPath = $ARGV[1] ;

my $i ;

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
my %speciesIdToName : shared ;

#Extract the gid we met in the file


print STDERR "Step 1: Collecting all .fna files in $path and getting gids\n";
find ( sub {
    return unless -f;        # Must be a file
    return unless -s;        # Must be non-zero

    if ( !( /\.f[nf]?a$/ || /\.ffn$/ || /\.fasta$/ ) )
    {
    	return ;
    }

    my $fullfile = $File::Find::name; ## file with full path, but the CWD is actually the file's path
    my $file = $_; ## file name
	open FP2, $file or die "Error opening $file: $!";
	my $head = <FP2> ;
	close FP2 ;
	if ( $head =~ /^>gi\|([0-9]+)?\|/ ) {
		my $gid = $1 ;
		print "gid=$gid $file\n" if $verbose;
		if ( defined $gidUsed{ $gid } )
		{
			print "Repeated gid $gid\n" if $verbose ;
			$fileUsed{ $fullfile } = 1 ;
		}
		else
		{
			$fileUsed{ $fullfile } = 0 ;
			$gidToFile{ $gid } = $fullfile ;
		}

		$gidUsed{ $gid } = 1 ;
	} elsif ( $head =~ /^>kraken:taxid\|([0-9]+)?\|/ ) {
		my $tid = $1 ;
		my $dummyGid = "centrifuge_gid_".$fullfile."_$1" ;
		$gidUsed{ $dummyGid } = 1 ;
		$gidToFile{ $dummyGid } = $fullfile ;
		$fileUsed{ $fullfile } = 0 ;
		push @{ $tidToGid{ $tid } }, $dummyGid ;	
		
		print "tid=$tid $file\n" if $verbose;
	} else {
		print STDERR "Excluding $fullfile: Wrong header.\n";
	}

}, $path );

# Extract the tid that are associated with the gids
print STDERR "Step 2: Extract the taxonomy ids that are associated with the gids\n";
open FP1, "$taxPath/gi_taxid_nucl.dmp" ;
while ( <FP1> )
{
	chomp ;
	my @cols = split ;		
	#print $cols[0], "\n" ;	
	if ( defined( $gidUsed{ $cols[0] } ) )
	{
		push @{ $tidToGid{ $cols[1] } }, $cols[0] ;	
		print "gid=", $cols[0], " tid=", $cols[1], " ", $gidToFile{ $cols[0] }, "\n" if $verbose;
	}
}
close FP1 ;


# Organize the tree
print STDERR "Step 3: Organizing the taxonomical tree\n";
open FP1, "$taxPath/nodes.dmp" ;
while ( <FP1> ) {
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

# Put the sub-species taxonomy id into the corresponding species.
print STDERR "Step 4: Putting the sub-species taxonomy id into the corresponding species\n";
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

print STDERR "Step 5: Reading the name of the species\n";
open FP1, "$taxPath/names.dmp" ;
while ( <FP1> )
{
	next if (!/scientific name/ ) ;
	my @cols = split /\t/ ;
	
	if ( defined $species{ $cols[0] } )
	{
		$speciesIdToName{ $cols[0] } = $cols[2] ;
	}
}
close FP1 ;
#exit ;

sub GetGenomeSize
{
	open FPfa, $_[0] ;
	my $size = 0 ;
	while ( <FPfa> )
	{
		next if ( /^>/ ) ;
		$size += length( $_ ) - 1 ;
	}
	close FPfa ;
	return $size ;
}

# Compress one species.
sub solve
{
# Extracts the files
#print "find $pwd -name *.fna > tmp.list\n" ;
	my $tid = threads->tid() - 1 ;
	#system_call("find $pwd -name *.fna > tmp_$tid.list") ;

#system_call("find $pwd -name *.fa >> tmp.list") ; # Just in case
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
	my $genomeSize = 0 ;
	foreach $i ( @subspeciesList )
	{
		foreach my $j ( @{$tidToGid{ $i } } )
		{
			#{
			#	lock( $debug ) ;
			#	print "Merge ", $gidToFile{ $j }, "\n" ;
			#}
			$file =  $gidToFile{ $j } ;
			{
				lock( $fileUsedLock ) ;
				$fileUsed{ $file } = 1 ;
			}
			print FP1 $file, "\n" ;

			$genomeSize += GetGenomeSize( $file ) ;
		}
	}
	close FP1 ;

	$genomeSize = int( $genomeSize / scalar( @subspeciesList ) ) ;
		
	#print $file, "\n" ;	
	#if ( $file =~ /\/(\w*?)uid/ )
	#{
	#	$speciesName = ( split /_/, $1 )[0]."_". ( split /_/, $1 )[1] ;
	#}
	if ( defined $speciesIdToName{ $speciesId } )
	{
		$speciesName = $speciesIdToName{ $speciesId } ;
		$speciesName =~ s/ /_/g ;
	}
	else
	{
		$speciesName = "Unknown_species_name" ;
	}
	my $id = ( $speciesId << 32 ) | $genusId ;
	my $header = ">$id $speciesName $genomeSize ".scalar( @subspeciesList ) ;
	print $header, "\n" ;

#return ;
# Build the sequence
	my $seq = "" ;
	if ( $noCompress == 0 &&  ( $maxGenomeSizeForCompression < 0 || $genomeSize <= $maxGenomeSizeForCompression ) ) #$genomeSize < 50000000 )
	{
		system_call("perl $bssPath/BuildSharedSequence.pl tmp_$tid.list -prefix tmp_${tid}_$id") ;

# Merge all the fragmented sequence into one big chunk.
		system_call("cat tmp_${tid}_${id}_*.fa > tmp_${tid}_$id.fa");

		open FP1, "tmp_${tid}_$id.fa" ;
		while ( <FP1> )
		{
			chomp ;
			next if ( /^>/ ) ;
			$seq .= $_ ;
		}
		close FP1 ;
	}
	else
	{
		foreach $i ( @subspeciesList )
		{
			foreach my $j ( @{$tidToGid{ $i } } )
			{
				$file =  $gidToFile{ $j } ;
				open FP1, $file ;
				while ( <FP1> )
				{
					chomp ;
					next if ( /^>/ ) ;
					$seq .= $_ ;
				}
				close FP1 ;
			}
		}
	}
	open fpOut, ">>${output}_${tid}" ;
	print fpOut "$header\n$seq\n" ;
	close fpOut ;

	unlink glob("tmp_${tid}_*");	
}

sub system_call {
    print STDERR "SYSTEM CALL: ".join(" ",@_);
	system(@_) == 0
	  or die "system @_ failed: $?";
    print STDERR " finished\n";
}

sub threadWrapper
{
	my $tid = threads->tid() - 1 ;
	unlink("${output}_${tid}");

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


print STDERR "Step 6: Merging sub-species\n";
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
system_call("cat ${output}_* > tmp_output.fa");
unlink glob("${output}_*");

#print unused files
foreach $i ( keys %fileUsed )
{
	if ( $fileUsed{ $i } == 0 )
	{
#print $i, "\n" ;
#`cat $i >> tmp_output.fa` ;
		print "Unused file: $i\n" ;
	}
}

# Remove the Ns from the file
system_call("$bssPath/RemoveN tmp_output.fa > $output") ;

unlink glob("tmp_*") ;
