#!/bin/perl

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

my $usage = "USAGE: perl ".basename($0)." --taxonomy=<taxnonomy directory>  <fasta file directories>* [-o compressed.fa -bss . -noCompress -t 1]\n" ;

my $level = "species" ;
my $output = "compressed.fa" ;
my $bssPath = $Bin; # take path of binary as script directory
my $numOfThreads = 1 ;
my $noCompress = 0 ;
my $verbose = 0;
my $useEutils = 0;

my $taxPath;
GetOptions ("level|l=s" => \$level,
			"output|o=s" => \$output,
			"bss=s" => \$bssPath,
			"threads|t=i" => \$numOfThreads,
            "verbose|v" => \$verbose,
            "taxonomy=s" => \$taxPath,
            "use-eutils|e" => \$useEutils,
			"noCompress" => \$noCompress)
or die("Error in command line arguments. \n\n$usage");

die $usage unless @ARGV > 1;
die $usage unless defined $taxPath;

my @paths = @ARGV;

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

###############################################################################
print STDERR "Step 1: Collecting all .fna files in @paths and getting gids\n";
find ( sub {
    return unless -f;        # Must be a file
    return unless /\.fna$/;  # Must end with `.fna` suffix
    return unless -s;        # Must be non-zero

    my $fullfile = $File::Find::name; ## file with full path, but the CWD is actually the file's path
    my $file = $_; ## file name
	open FP2, $file or die "Error opening $file: $!";
	my $head = <FP2> ;
	close FP2 ;
	if ( $head =~ /^>gi\|([0-9]+)?\|/ ) 
    {
		my $gid = $1 ;
		print "$gid $file\n" if $verbose;
		$gidUsed{ $gid } = 1 ;
		$gidToFile{ $gid } = $fullfile ;
		$fileUsed{ $fullfile } = 0 ;
	} elsif ( $head =~ /^>[a-zA-Z]*:taxid\|([0-9]+)?\|/ ) 
    {
		my $tid = $1 ;
		$gidUsed{ $tid } = 1 ;
		$gidToFile{ $tid } = $fullfile ;
		$fileUsed{ $fullfile } = 0 ;
		push @{ $tidToGid{ $tid } }, $tid ;	

	} else 
    {
		print STDERR "Excluding $fullfile: Wrong header.\n";
	}

}, @paths );


###############################################################################
print STDERR "Step 2: Extract the mapping taxonomy id -> gids\n";

if ($useEutils) {
    # First, build the Eutils query
    my $utils = 'http://eutils.ncbi.nlm.nih.gov/entrez/eutils'; # Base URL for searches
    my $db = 'nuccore'; # Default to PubMed database; this may be changed.

    my $esearch = "$utils/esummary.fcgi?db=$db";

    foreach my $gid (keys %gidUsed) {
        my $esearch_result = get("$esearch&id=$gid");
        my @matches = ($esearch_result =~ m(<Item Name="TaxId" Type="Integer">([0-9]*)</Item>)g);
        die "so many or few matches for $gid: @matches" if (scalar @matches != 1);
	    push @{ $tidToGid{ $matches[0] } }, $gid;
    }
} else {
    open my $FP1,"<","$taxPath/gi_taxid_nucl.dmp" or die "Error opening $taxPath/gi_taxid_nucl.dmp: $!";
    while ( <$FP1> )
    {
    	chomp ;
    	my ($gid,$tid) = split ;
	    push @{ $tidToGid{ $tid } }, $gid;
    }
    close ($FP1);
}


###############################################################################
print STDERR "Step 3: Organizing the taxonomical tree taxid -> parent-taxid\n";
open FP1, "$taxPath/nodes.dmp" or die "Error opening $taxPath/nodes.dmp: $!";
while ( <FP1> ) {
	chomp ;

    # Format: 
	#  1 tax_id					-- node id in GenBank taxonomy database
 	#  2 parent tax_id				-- parent node id in GenBank taxonomy database
 	#  3 rank					-- rank of this node (superkingdom, kingdom, ...) 
 	#  4 embl code				-- locus-name prefix; not unique
 	#  5 division id				-- see division.dmp file
 	#  6 inherited div flag  (1 or 0)		-- 1 if node inherits division from parent
 	#  7 genetic code id				-- see gencode.dmp file
 	#  8 inherited GC  flag  (1 or 0)		-- 1 if node inherits genetic code from parent
 	#  9 mitochondrial genetic code id		-- see gencode.dmp file
 	# 10 inherited MGC flag  (1 or 0)		-- 1 if node inherits mitochondrial gencode from parent
 	# 11 GenBank hidden flag (1 or 0)            -- 1 if name is suppressed in GenBank entry lineage
 	# 12 hidden subtree root flag (1 or 0)       -- 1 if this subtree has no sequence data yet
 	# 13 comments				-- free-text comments and citations

    #   tax_id  |   parent tax_id  |    rank      |  embl code | division id | genetic code |   inh GC  |
	my ($tid, undef, $parentTid, undef, $rank ) = split;

	$taxTree{ $tid } = $parentTid ;
	if ( $rank eq "species" ) {
		$species{ $tid } = 1 ;
	}
	elsif ( $rank eq "genus" )
	{
		$genus{ $tid } = 1 ;
	}
}
close FP1 ;

###############################################################################
# Put the sub-species taxonomy id into the corresponding species.
print STDERR "Step 4: Associating sub-species taxonomy ids with the species taxonomy id\n";
for my $tid ( keys %tidToGid )
{
	my $p = $tid ;
	my $found = 0 ;
	while ( 1 )
	{
		last if ( $p <= 1 ) ;
		if ( defined $species{ $p } ) 
		{
			$found = 1 ;
			last ;
		}
        last unless defined $taxTree{ $p };
		$p = $taxTree{ $p } ;
	}

	push(@{ $speciesList{ $p } }, $tid) if $found;
}


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
	if ( $file =~ /\/(\w*?)uid/ )
	{
		$speciesName = ( split /_/, $1 )[0]."_". ( split /_/, $1 )[1] ;
	}
	else
	{
		$speciesName = "" ;
	}
	my $id = ( $speciesId << 32 ) | $genusId ;
	my $header = ">$id $speciesName $genomeSize ".scalar( @subspeciesList ) ;
	print $header, "\n" ;

#return ;
# Build the sequence
	my $seq = "" ;
	if ( $noCompress == 0 )
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
    print STDERR "SYSTEM CALL: ".join(" ",@_)."\n" ;
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


###############################################################################
print STDERR "Step 5: Merging sub-species\n";
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

#Merge the unused files
open FP2, ">>", "tmp_output.fa" ;
my $seq = "" ;
my $genomeSize = 0 ;
my $k = 0 ;
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
		$genomeSize += GetGenomeSize( $i ) ;
		++$k ;
	}
}
$genomeSize /= $k ;
$genomeSize = int( $genomeSize ) ;
print FP2 ">0 unknown_taxno_genomes $genomeSize $k\n$seq\n" ;
close FP2 ;

# Remove the Ns from the file
system_call("perl $bssPath/RemoveN.pl tmp_output.fa > $output") ;

unlink glob("tmp_*") ;
