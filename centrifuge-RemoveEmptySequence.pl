#!/usr/bin/env perl

# remove the headers with empty sequences. possible introduced by dustmask

use strict ;

# die "usage: a.pl < in.fa >out.fa"

my $tag ;
my @lines ;

$lines[0] = <> ;
$tag = 0 ;
while ( <> )
{
	$lines[1-$tag] = $_ ;
	if ( /^>/ && $lines[$tag] =~ /^>/ ) 
	{
		$tag = 1 - $tag ;
		next ;
	}
	print $lines[$tag] ;
	$tag = 1- $tag ;
}
if ( !( $lines[$tag] =~ /^>/ ) )
{
	print $lines[$tag] ;
}
