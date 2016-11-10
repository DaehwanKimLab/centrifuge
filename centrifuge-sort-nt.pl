#! /usr/bin/env perl
#
# Sort nt file sequences according to their taxonomy ID
# Uses the new mappinf file format available at
# ftp://ftp.ncbi.nih.gov/pub/taxonomy/accession2taxid/
#
# Author fbreitwieser <fbreitwieser@sherman>
# Version 0.2
# Copyright (C) 2016 fbreitwieser <fbreitwieser@sherman>
# Modified On 2016-09-26
# Created  2016-02-28 12:56
#
use strict;
use warnings;
use File::Basename;
use Getopt::Long;

my $new_map_file;
my $ac_wo_mapping_file;
my $opt_help;

my $USAGE = 
"USAGE: ".basename($0)." OPTIONS <sequence file> <mapping file> [<mapping file>*]\n

OPTIONS:
  -m str      Output mappings that are present in sequence file to file str
  -a str      Output ACs w/o mapping to file str
  -h          This message
";

GetOptions(
	"m|map=s" => \$new_map_file,
	"a=s" => \$ac_wo_mapping_file,
	"h|help" => \$opt_help) or die "Error in command line arguments";

scalar(@ARGV) >= 2 or die $USAGE;
if (defined $opt_help) {
	print STDERR $USAGE;
	exit(0);
}

my ($nt_file, @ac_taxid_files) = @ARGV;

my %ac_to_pos;
my %taxid_to_ac;
my %ac_to_taxid;

print STDERR "Reading headers from $nt_file ... ";
open(my $NT, "<", $nt_file) or die $!;
while (<$NT>) {
    # get the headers with (!) the version number
    if (/(^>([^ ]*).*)/) {
        # record the position of this AC
        $ac_to_pos{$2} = [tell($NT),$1];
    }
}
print STDERR "found ", scalar(keys %ac_to_pos), " ACs\n";

foreach my $ac_taxid_file (@ac_taxid_files) {
print STDERR "Reading ac to taxid mapping from $ac_taxid_file ...\n";
    my $FP1;
    if ($ac_taxid_file =~ /.gz$/) {
        open($FP1, "-|", "gunzip -c '$ac_taxid_file'") or die $!;
    } else {
        open($FP1, "<", $ac_taxid_file) or die $!;
    }

    # format: accession <tab> accession.version <tab> taxid <tab> gi
    # currently we look for a mapping with the version number
    while ( <$FP1> ) {
    	my (undef, $ac, $taxid) = split;
        next unless defined $taxid;
	    if ( defined( $ac_to_pos{ $ac } ) ) {
            push @{ $taxid_to_ac{ $taxid } }, $ac;
		    $ac_to_taxid{ $ac } = $taxid;
        }
    }
    close $FP1;
}
print STDERR "Got taxonomy mappings for ", scalar(keys %ac_to_taxid), " ACs\n";
if (defined $ac_wo_mapping_file && scalar(keys %ac_to_taxid) < scalar(keys %ac_to_pos)) {
    print STDERR "Writing ACs without taxonomy mapping to $ac_wo_mapping_file\n";
    open(my $FP2, ">", $ac_wo_mapping_file) or die $!;
    foreach my $ac (keys %ac_to_pos) {
        next unless defined $ac_to_taxid{$ac};
        print $FP2 $ac, "\n";
    }
    close($FP2);
}

if (defined $new_map_file) {
print STDERR "Writing taxonomy ID mapping to $new_map_file\n";
open(my $FP3, ">", $new_map_file) or die $!;
foreach my $ac (keys %ac_to_taxid) {
    print $FP3 $ac,"\t",$ac_to_taxid{$ac},"\n";
}
close($FP3);
}


print STDERR "Outputting sorted FASTA ...\n";
foreach my $taxid (sort {$a <=> $b} keys %taxid_to_ac) {
    my @acs = @{$taxid_to_ac{$taxid}};
    my @sorted_acs = sort { $ac_to_pos{$a}->[0] <=> $ac_to_pos{$b}->[0] } @acs;
    foreach (@sorted_acs) {
        print $ac_to_pos{$_}->[1],"\n";
        seek($NT, $ac_to_pos{$_}->[0], 0);
        while (<$NT>) {
            last if (/^>/);
            print $_;
        }
    }
}
close $NT;
