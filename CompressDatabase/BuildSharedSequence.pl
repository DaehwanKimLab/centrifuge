#!/bin/perl

use strict ;
use Getopt::Long;
use File::Basename;

my $usage = "perl ".basename($0)." file_list [-prefix tmp -kmerSize 53 -kmerPortion 0.01 -nucmerIdy 99 -overlap 250 [-fragment] ]" ;

my @fileNames ; # Assume the genes for each subspecies are concatenated.
my @used ;
my $i ;
my $j ;
my $k ;
my %globalKmer ;
my $kmerSize = 53 ;
my $useKmerPortion = 0.01 ;
my %localKmer ;
my $id = 0 ;
my $kmer ;
my %chroms ;
my %shared ;
my $id = "" ;
my $seq = "" ;
my $index ;
my $prefix = "tmp";
my %sharedKmerCnt ;
my $nucmerIdy = 99 ;
my $overlap = 250 ;
my $fragment = 0 ;
my $jellyfish = "jellyfish";

print `jellyfish --version`;

GetOptions (
	"prefix=s" => \$prefix,
	"kmerSize=s" => \$kmerSize,
	"kmerPortion=s" => \$useKmerPortion,
	"nucmerIdy=s" => \$nucmerIdy,
	"overlap=s" => \$overlap,
	"fragment" => \$fragment,
	"jellyfish=s" => \$jellyfish)
or die("Error in command line arguments. \n\n$usage");

die "$usage\n" if ( scalar( @ARGV ) == 0 ) ;
open FP1, $ARGV[0] ; 

# Create the temporary files, while making sure the header id is unique within each file
my $listSize = 0 ;
while ( <FP1> )
{
	chomp ;
	open FP2, $_ ;
	my $fileName = $prefix."_".$listSize.".fa" ;
	open FPtmp, ">$fileName" ;
	my $chromCnt = 0 ;
	while ( <FP2> )
	{
		if ( /^>/ )
		{
			s/\s/\|$chromCnt / ;
			print FPtmp ;
			++$chromCnt ;
		}
		else
		{
			print FPtmp ;
		}
	}
	++$listSize ;
}

# Read in the file names
#while ( <FP1> )
#{
#	chomp ;
#	push @fileNames, $_ ;
#	push @used, 0 ;
#	++$index ;
#}
#close( FP1 ) ;
for ( my $i = 0 ; $i < $listSize ; ++$i )
{
	my $fileName = $prefix."_".$i.".fa" ;
	push @fileNames, $fileName ;
	push @used, 1 ;
	++$index ;
}



# Count and select kmers
print "Find the kmers for testing\n" ;
system_call("$jellyfish count -o tmp_$prefix.jf -m $kmerSize -s 5000000 -C -t 8 @fileNames") ;
system_call("$jellyfish dump tmp_$prefix.jf > tmp_$prefix.jf_dump") ;

srand( 17 ) ;
open FP1, "tmp_$prefix.jf_dump" ;
open FP2, ">tmp_$prefix.testingKmer" ;
while ( <FP1> )
{
	my $line = $_ ;
	$kmer = <FP1> ;
	#chomp $kmer ;
	
	#++$testingKmer{ $kmer} if ( rand() < $useKmerPortion ) ;
	print FP2 "$line$kmer" if ( rand() < $useKmerPortion ) ;
		
}
close FP1 ;
close FP2 ;

print "Get the kmer profile for each input file\n" ;
for ( $i = 0 ; $i < scalar( @fileNames )  ; ++$i )
{
	my $fileName = $fileNames[ $i ] ;
	system_call("$jellyfish count --if tmp_$prefix.testingKmer -o tmp_$prefix.jf -m $kmerSize -s 5000000 -C -t 8 $fileName") ;
	system_call("$jellyfish dump tmp_$prefix.jf > tmp_$prefix.jf_dump") ;
	open FP1, "tmp_$prefix.jf_dump";

	while ( <FP1> )
	{
		chomp ;
		my $cnt = substr( $_, 1 )  ;
		if ( $cnt eq "0" ) 
		{
			$kmer = <FP1> ;
			next ;
		}
		#print "$cnt\n" ;
		$kmer = <FP1> ;
		chomp $kmer ;
		$localKmer{$i}->{$kmer} = 1 ; 
	}
	close FP1 ;
}

#for ( $i = 0 ; $i < scalar( @fileNames ) ; ++$i )
print "Begin merge files\n" ;
my $maxSharedKmerCnt = -1 ;
while ( 1 ) 
{
	# Find the suitable files
	my @maxPair ;
	my $max = 0 ;
	print "Selecting two genomes to merge.\n" ;
	foreach $i (keys %localKmer )
	{
		foreach $j ( keys %localKmer )
		{
			next if ( $i <= $j ) ;

			my $cnt = 0 ;

			if ( defined $sharedKmerCnt{ "$i $j" } )
			{
				$cnt = $sharedKmerCnt{ "$i $j" }
			}
			else
			{
				foreach $kmer ( keys %{$localKmer{ $i } } )
				{
					#print $kmer, "\n" ;
					if ( defined $localKmer{ $j }->{ $kmer} ) 
					{
						++$cnt ;				
					}	
				}
				$sharedKmerCnt{ "$i $j" } = $cnt ;
			}

			if ( $cnt > $max )
			{
				$max = $cnt ;
				$maxPair[0] = $i ;
				$maxPair[1] = $j ;
			}
		}
	}

	$maxSharedKmerCnt = $max if ( $maxSharedKmerCnt == -1 ) ;
	last if ( $max == 0 || $max < $maxSharedKmerCnt * 0.03 ) ;

	my @commonRegion ;
	my $fileNameA ;
	my $fileNameB ;
	$i = $maxPair[0] ;
	$j = $maxPair[1] ;
	if ( $i < scalar( @fileNames ) )
	{
		$fileNameA = $fileNames[ $i ] ;
	}
	else
	{
		$fileNameA = $prefix."_".$i.".fa" ;
	}
	if ( $j < scalar( @fileNames ) )
	{
		$fileNameB = $fileNames[ $j ] ;
	}
	else
	{
		$fileNameB = $prefix."_".$j.".fa" ;
	}
	my $nucmerC = 3 * $overlap ;
	print "nucmer --maxmatch --coords -l $kmerSize -g 10 -b 10 -c $nucmerC -p nucmer_$prefix $fileNameA $fileNameB\n" ; 
	my $nucRet = system("nucmer --maxmatch --coords -l $kmerSize -g 10 -b 10 -c $nucmerC -p nucmer_$prefix $fileNameA $fileNameB") ; # if the call to nucmer failed, we just not compress at all. 

	open FPCoords, "nucmer_$prefix.coords" ;
	my $line ;
	$line = <FPCoords> ;
	$line = <FPCoords> ;
	$line = <FPCoords> ;
	$line = <FPCoords> ;
	$line = <FPCoords> ;
	
	#1     5195  |        1     5195  |     5195     5195  |    99.98  | gi|385223048|ref|NC_017374.1|	gi|385230889|ref|NC_017381.1|
	print "Merging $fileNameA $fileNameB\n" ;

	my $cnt = 0 ;	
	while ( <FPCoords> )
	{
		chomp ;
		$line = $_ ;
		my @cols = split /\s+/, $line ;

		shift @cols if ( $cols[0] eq "" ) ;

		next if ( $cols[6] <= 3 * $overlap || $cols[9] < $nucmerIdy ) ;
		
		++$cnt ;
		my $ind = scalar( @commonRegion ) ;
		push @{ $commonRegion[$ind] }, ( $cols[11], $cols[0] + $overlap, $cols[1] - $overlap )  ;
		if ( $cols[3] < $cols[4] )
		{
			push @{ $commonRegion[ $ind ] }, ( $cols[12], $cols[3] + $overlap, $cols[4] - $overlap ) ;
		}
		else
		{
			push @{ $commonRegion[ $ind ] }, ( $cols[12], $cols[4] + $overlap, $cols[3] - $overlap ) ;
		}
	}
	last if ( $cnt == 0 ) ;
	close FPCoords ;

	my $outputSeq = "" ;
	my $outputHeader = ">${prefix}_${index}" ;
	# Use fileNameA to represent the shared sequences
	if ( $fragment == 0 )
	{
		# Just a big chunk.
		open FP1, $fileNameA ;
		while ( <FP1> )
		{
			chomp ;
			next if ( /^>/ ) ;
			$outputSeq .= $_ ;
			#print "$_\n$outputSeq\n" ;
		}
		close FP1 ;
	}
	else
	{
		open FP1, $fileNameA ;
		$id = "" ;
		$seq = "" ;
		undef %chroms ;
		undef %shared ;
		while ( <FP1> )
		{
			chomp ;
			if ( /^>/ )
			{
				$chroms{ $id } = $seq if ( $id ne "" ) ;
				$id = ( split /\s+/, substr( $_, 1 ) )[0] ;
				#print $id, "\n" ;
				$seq = "" ;
				@#{ $shared{ $id } } = length( $seq ) + 1 ; 
			}
			else
			{
				$seq .= $_ ;
			}
		}
		$chroms{ $id } = $seq if ( $id ne "" ) ;
		@#{ $shared{ $id } } = length( $seq ) + 1 ; 
		close FP1 ;

		for ( $i = 0 ; $i < scalar( @commonRegion) ; ++$i )
		{
			#print "hi $i $commonRegion[$i]->[1] $commonRegion[$i]->[2]\n" ;	
			my $tmpArray = \@{ $shared{ $commonRegion[$i]->[0] } } ;
			$commonRegion[$i]->[1] -= $overlap if ( $commonRegion[$i]->[1] == $overlap + 1 ) ;
			$commonRegion[$i]->[2] += $overlap if ( $commonRegion[$i]->[2] + $overlap == length( $chroms{ $commonRegion[$i]->[0] } ) ) ;
			for ( $j = $commonRegion[$i]->[1] - 1 ; $j < $commonRegion[$i]->[2] ; ++$j ) # Shift the coordinates
			{
				$tmpArray->[$j] = 1 ;
			}
		}

		# Print the information of genome A, including the shared part
		my $fileName = $prefix."_".$index.".fa" ;
		open fpOut, ">$fileName" ;
		foreach $i (keys %chroms )
		{
			my $tmpArray = \@{ $shared{ $i } } ;
			my $len = length( $chroms{ $i } ) ;
			my $header = ( split /\|Range:/, $i )[0] ;
			my $origStart = ( split /-/, ( ( split /\|Range:/, $i )[1] ) )[0] ;
			for ( $j = 0 ; $j < $len ;  )
			{
				if ( $tmpArray->[$j] == 1 )
				{
					my $start = $j ;
					my $end ;
					for ( $end = $j + 1 ; $end < $len && $tmpArray->[$end] == 1 ; ++$end )
					{
						;
					}
					--$end ;
					my $rangeStart = $origStart + $start ;
					my $rangeEnd = $origStart + $end ;
					print fpOut ">$header|Range:$rangeStart-$rangeEnd shared\n" ;
					print fpOut substr( $chroms{$i}, $start, $end - $start + 1 ), "\n" ;

					$j = $end + 1 ;
				}
				else
				{
					my $start = $j ;
					my $end ;

					for ( $end = $j + 1 ; $end < $len && $tmpArray->[$end] == 0 ; ++$end )
					{
						;
					}
					--$end ;
					my $rangeStart = $origStart + $start ;
					my $rangeEnd = $origStart + $end ;

					print fpOut ">$header|Range:$rangeStart-$rangeEnd non-shared\n" ;
					print fpOut substr( $chroms{$i}, $start, $end - $start + 1 ), "\n" ;

					$j = $end + 1 ;
				}
			}
		}
	} # end if fragment. There might be bugs in the fragment mode

	# Print the sequence from genome B, only including the non-shared part
	open FP1, $fileNameB ;
	$id = "" ;
	$seq = "" ;
	undef %chroms ;
	undef %shared ;
	while ( <FP1> )
	{
		chomp ;
		if ( /^>/ )
		{
			$chroms{ $id } = $seq if ( $id ne "" ) ;
			$id = ( split /\s+/, substr( $_, 1 ) )[0] ;
			$seq = "" ;
			@#{ $shared{ $id } } = length( $seq ) + 1 ; 
		}
		else
		{
			$seq .= $_ ;
		}
	}
	$chroms{ $id } = $seq if ( $id ne "" ) ;
	@#{ $shared{ $id } } = length( $seq ) + 1 ; 
	close FP1 ;

	for ( $i = 0 ; $i < scalar( @commonRegion) ; ++$i )
	{
		my $tmpArray = \@{ $shared{ $commonRegion[$i]->[3] } } ;
		$commonRegion[$i]->[4] -= $overlap if ( $commonRegion[$i]->[4] == $overlap + 1 ) ;
		$commonRegion[$i]->[5] += $overlap if ( $commonRegion[$i]->[5] + $overlap == length( $chroms{ $commonRegion[$i]->[3] } ) ) ;
		for ( $j = $commonRegion[$i]->[4] - 1 ; $j < $commonRegion[$i]->[5] ; ++$j )	
		{
			$tmpArray->[$j] = 1 ;
		}
	}

	foreach $i (keys %chroms )
	{
		my $tmpArray = \@{ $shared{ $i } } ;
		my $len = length( $chroms{ $i } ) ;
		my $header = ( split /\|Range:/, $i )[0] ;
		my $origStart = ( split /-/, ( ( split /\|Range:/, $i )[1] ) )[0] ;
		for ( $j = 0 ; $j < $len ;  )
		{
			if ( $tmpArray->[$j] == 1 ) 
			{
				++$j ;
				next ;
			}

			my $start = $j ;
			my $end ;
			for ( $end = $j + 1 ; $end < $len && $tmpArray->[$end] == 0 ; ++$end )
			{
				;
			}
			--$end ;
			$j = $end + 1 ;
			
			next if ( $end - $start < $overlap ) ;

			my $rangeStart = $origStart + $start ;
			my $rangeEnd = $origStart + $end ;

			if ( $fragment == 0 )
			{
				#print length( $outputSeq ), "\n" ;
				$outputSeq .= substr( $chroms{$i}, $start, $end - $start + 1 ) ; 	
				#print length( $outputSeq ), "\n" ;
			}
			else
			{
				my $fileName = $prefix."_".$index.".fa" ;
				#open fpOut, ">>$fileName" ;
				print fpOut ">$header|Range:$rangeStart-$rangeEnd non-shared\n" ;
				print fpOut substr( $chroms{$i}, $start, $end - $start + 1 ), "\n" ;
				#close fpOut ;
			}
		}
	}

	delete $localKmer{ $maxPair[0] } ;
	delete $localKmer{ $maxPair[1] } ;

	#print defined( $localKmer{ $maxPair[0] }) ;
	unlink glob("$fileNameA") ; #if ( $maxPair[0] >= scalar( @fileNames ) ) ;
	unlink glob("$fileNameB") ; #if ( $maxPair[1] >= scalar( @fileNames ) ) ;

	# Count the kmer for the new genome
	my $fileName = $prefix."_".$index.".fa" ;
	if ( $fragment == 0 )
	{
		open fpOut, ">$fileName" ;
		print fpOut "$outputHeader\n$outputSeq\n" ;
	}
	close fpOut ;

	system_call("$jellyfish count --if tmp_$prefix.testingKmer -o tmp_$prefix.jf -m $kmerSize -s 5000000 -C -t 4 $fileName") ;
	system_call("$jellyfish dump tmp_$prefix.jf > tmp_$prefix.jf_dump") ;
	open FP1, "tmp_$prefix.jf_dump";

	while ( <FP1> )
	{
		chomp ;
		my $cnt = substr( $_, 1 )  ;
		if ( $cnt eq "0" )
		{
			$kmer = <FP1> ;
			next ;
		}
		$kmer = <FP1> ;
		chomp $kmer ;
		$localKmer{$index}->{$kmer} = 1 ;
	}
	close FP1 ;
	
	++$index ;
} # while 1

#foreach $i (keys %localKmer )
#{
#	if ( $i < scalar( @fileNames ) ) 
#	{
#		my $path = $fileNames[ $i ] ;
#		my $fileName = $prefix."_".$i.".fa" ;
#		system_call("cp $path $fileName") ;
#	}
#}

# clean up
unlink glob("tmp_$prefix.jf tmp_$prefix.jf_dump tmp_$prefix.testingKmer") ;
unlink glob("nucmer_$prefix*") ;
print "Finish.\n" ;

sub system_call {
    print STDERR "SYSTEM CALL: ".join(" ",@_)."\n" ;
	system(@_) == 0
	  or die "system @_ failed: $?";
    print STDERR " finished\n";
}


