#!/bin/perl

my $usage = "perl a.pl xxx.oneline.fa" ;

die "$usage\n" if ( scalar( @ARGV ) == 0 ) ;

open FP1, $ARGV[0] ;
while ( <FP1> )
{
	chomp ;
	my $id = $_ ;
	my $seq = <FP1> ;
	chomp $seq ;

	my $i ;
	my $cntN = 0 ;
	my $startN = 0 ;
	my $stretchN = 0 ;
	my $max ;
	my @cols ;
	for ( $i = 0 ; $i < -1 ; ++$i )
	{
		if ( substr( $seq, $i, 1 ) eq "N" ) 
		{
			++$cntN ;
			if ( $startN == 0 )
			{
				$stretchN = 1 ;
				$startN = 1 ;
			}
			else
			{
				++$stretchN ;
			}
		}
		else
		{
			if ( $startN == 1 && $stretchN > $max )
			{
				$max = $stretchN ;
			}
			$startN = 0 ;
		}
	}
	#print "$id $cntN $i ", $cntN * 100.0 / $i, " $stretchN\n" ;

	@cols = split /N/, $seq ;
	my $outputSeq = "" ;
	for ( $i = 0 ; $i < scalar( @cols ) ; ++$i )
	{
		$outputSeq .= $cols[$i] ;
	}
	print "$id\n$outputSeq\n" ;
}
