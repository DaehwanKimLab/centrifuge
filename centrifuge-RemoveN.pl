#!/usr/bin/env perl

use strict ;
use warnings ;

die "usage: a.pl xxx.fa > output.fa" if ( @ARGV == 0 ) ;

my $LINE_WIDTH = 80 ;
my $BUFFER_SIZE = 100000 ;
open FP1, $ARGV[0] ;
my $buffer = "" ;

while ( <FP1> )
{
	if ( /^>/ )
	{
		my $line = $_ ;
		if ( $buffer ne "" )
		{
			my $len = length( $buffer ) ;
			for ( my $i = 0 ; $i < $len ; $i += $LINE_WIDTH )
			{
				print substr( $buffer, $i, $LINE_WIDTH ), "\n" ;
			}
			$buffer = "" ;
		}
		print $line ;
	}
	else
	{
		chomp ;
		tr/nN//d ;
		#my $line = uc( $_ ) ;
		my $line = $_ ;
		$buffer .= $line ;

		if ( length( $buffer ) > $BUFFER_SIZE )
		{
			my $len = length( $buffer ) ;
			my $i = 0 ;
			for ( $i = 0 ; $i + $LINE_WIDTH - 1 < $len ; $i += $LINE_WIDTH )
			{
				print substr( $buffer, $i, $LINE_WIDTH ), "\n" ;
			}
			substr( $buffer, 0, $i, "" ) ;
		}
	}
}
if ( $buffer ne "" )
{
	my $len = length( $buffer ) ;
	for ( my $i = 0 ; $i < $len ; $i += $LINE_WIDTH )
	{
		print substr( $buffer, $i, $LINE_WIDTH ), "\n" ;
	}
	$buffer = "" ;
}
