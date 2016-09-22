#! /usr/bin/env perl
#
# Sort nt file sequences according to their taxonomy ID
#
# Author fbreitwieser <fbreitwieser@sherman>
# Version 0.2
# Copyright (C) 2016 fbreitwieser <fbreitwieser@sherman>
# Modified On 2016-09-22
# Created  2016-02-28 12:56
#
use strict;
use warnings;
use File::Basename;

$#ARGV==1 or die "USAGE: ".basename($0)." <sequence file> <mapping file>\n";
my ($nt_file,$gi_taxid_file) = @ARGV;

my %gi_to_pos;
my %taxid_to_gi;
my %gi_to_taxid;

print STDERR "Reading headers from $nt_file ...\n";
open(my $NT, "<", $nt_file) or die $!;
while (<$NT>) {
    if (/(^>(gi\|[0-9]*).*)/) {
        #print STDERR "\$gi_to_pos{$2} = [$1,".tell($NT)."]\n";
        $gi_to_pos{$2} = [tell($NT),$1];
    }
}

print STDERR "Reading gi to taxid mapping $gi_taxid_file ...\n";
my $FP1;
if ($gi_taxid_file =~ /.zip$/) {
    open($FP1, "-|", "unzip -c '$gi_taxid_file'") or die $!;
} else {
    open($FP1, "<", $gi_taxid_file) or die $!;
}
while ( <$FP1> ) {
	chomp;
	my ($gi,$taxid) = split;
    next unless defined $taxid;
	if ( defined( $gi_to_pos{ $gi } ) )
	{
		push @{ $taxid_to_gi{ $taxid } }, $gi;
		$gi_to_taxid{ $gi } = $taxid;
	}
}
close $FP1;

print STDERR "Outputting sorted FASTA ...\n";
foreach my $taxid (sort {$a <=> $b} keys %taxid_to_gi) {
    my @gis = @{$taxid_to_gi{$taxid}};
    my @sorted_gis = sort { $gi_to_pos{$a}->[0] <=> $gi_to_pos{$b}->[0] } @gis;
    foreach (@sorted_gis) {
        print $gi_to_pos{$_}->[1],"\n";
        seek($NT, $gi_to_pos{$_}->[0], 0);
        while (<$NT>) {
            last if (/^>/);
            print $_;
        }
    }
}
close $NT;
