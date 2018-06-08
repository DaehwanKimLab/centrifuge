#!/usr/bin/env perl

# Give a Kraken-style report from a Centrifuge output
#
# Based on kraken-report by Derrick Wood
# Copyright 2013-2016, Derrick Wood <dwood@cs.jhu.edu>
#

use strict;
use warnings;
use Getopt::Long;
use File::Basename;
use Cwd;
use Cwd 'cwd' ;
use Cwd 'abs_path' ;

my ($centrifuge_index, $min_score, $min_length);
my $no_lca = 0;
my $show_zeros = 0;
my $is_cnts_table = 0;
my $PROG = "centrifuge-kreport";
my $CWD = dirname( abs_path( $0 ) ) ;

GetOptions(
  "help" => \&display_help,
  "x=s" => \$centrifuge_index,
  "show-zeros" => \$show_zeros,
  "is-count-table" => \$is_cnts_table,
  "min-score=i" => \$min_score,
  "min-length=i"=> \$min_length,
  "no-lca" => \$no_lca
) or usage();

usage() unless defined $centrifuge_index;
if (!defined $ARGV[0]) {
    print STDERR "Reading centrifuge out file from STDIN ... \n";
}

sub usage {
  my $exit_code = @_ ? shift : 64;
  print STDERR "
Usage: centrifuge-kreport -x <index name> OPTIONS <centrifuge output file(s)>

centrifuge-kreport creates Kraken-style reports from centrifuge out files.

Options:
    -x INDEX            (REQUIRED) Centrifuge index

    --no-lca             Do not report the LCA of multiple assignments, but report count fractions at the taxa.
    --show-zeros         Show clades that have zero reads, too
    --is-count-table     The format of the file is 'taxID<tab>COUNT' instead of the standard
                         Centrifuge output format

    --min-score SCORE    Require a minimum score for reads to be counted
    --min-length LENGTH  Require a minimum alignment length to the read
  
  ";
  exit $exit_code;
}

sub display_help {
  usage(0);
}

my (%child_lists, %name_map, %rank_map, %parent_map);
print STDERR "Loading taxonomy ...\n";
load_taxonomy();

my %taxo_counts;
my $seq_count = 0;
$taxo_counts{0} = 0;
if ($is_cnts_table) {
  while (<>) {
    my ($taxID,$count) = split;
    $taxo_counts{$taxID} = $count;
    $seq_count += $count;
  }
} else {
  chomp(my $header = <>);
  my @cols = split /\t/, $header ;
  my %headerMap ;
  for ( my $i = 0 ; $i < scalar( @cols ) ; ++$i )
  {
  	$headerMap{ $cols[$i] } = $i ;
  }

  my $prevReadID;
  my $prevTaxID;
  while (<>) {
    #my (undef,$seqID,$taxID,$score, undef, $hitLength, $queryLength, $numMatches) = split /\t/;
    my @cols = split /\t/ ;
    my $readID = $cols[ $headerMap{ "readID" } ] ;
    my $seqID = $cols[ $headerMap{ "seqID" } ] ; 
    my $taxID = $cols[ $headerMap{ "taxID" } ] ; 
    my $score = $cols[ $headerMap{ "score" } ] ; 
    my $hitLength = $cols[ $headerMap{ "hitLength" } ] ; 
    my $queryLength = $cols[ $headerMap{ "queryLength" } ] ; 
    my $numMatches = $cols[ $headerMap{ "numMatches" } ] ; 

    $taxID = 1 if ( !isTaxIDInTree( $taxID ) ) ;

    if ($no_lca) {
      next if defined $min_length && $hitLength < $min_length;
      next if defined $min_score && $score < $min_score;
      $taxo_counts{$taxID} += 1/$numMatches;
      $seq_count += 1/$numMatches;
    } else {
      if ( ( defined $prevReadID ) && ( $readID eq $prevReadID ) ) {
        --$taxo_counts{$prevTaxID};
        $prevTaxID = lca($prevTaxID, $taxID);
        ++$taxo_counts{$prevTaxID};
      } else {
        ++$taxo_counts{$taxID};
        ++$seq_count;
        $prevTaxID = $taxID;
      }
    }
    $prevReadID = $readID ;
  }
}
my $classified_count = $seq_count - $taxo_counts{0};

my %clade_counts = %taxo_counts;
dfs_summation(1);

for (keys %name_map) {
  $clade_counts{$_} ||= 0;
}

die "No sequence matches with given settings" unless $seq_count > 0;

printf "%6.2f\t%d\t%d\t%s\t%d\t%s%s\n",
  $clade_counts{0} * 100 / $seq_count,
  $clade_counts{0}, $taxo_counts{0}, "U",
  0, "", "unclassified";
dfs_report(1, 0);

sub dfs_report {
  my $node = shift;
  my $depth = shift;
  if (! $clade_counts{$node} && ! $show_zeros) {
    return;
  }
  printf "%6.2f\t%d\t%d\t%s\t%d\t%s%s\n",
    ($clade_counts{$node} || 0) * 100 / $seq_count,
    ($clade_counts{$node} || 0),
    ($taxo_counts{$node} || 0),
    rank_code($rank_map{$node}),
    $node,
    "  " x $depth,
    $name_map{$node};
  my $children = $child_lists{$node};
  if ($children) {
    my @sorted_children = sort { $clade_counts{$b} <=> $clade_counts{$a} } @$children;
    for my $child (@sorted_children) {
      dfs_report($child, $depth + 1);
    }
  }
}

sub isTaxIDInTree {
	my $a = $_[0] ;

	while ( $a > 1 )
	{
		if ( !defined $parent_map{ $a } )
		{
			print STDERR "Couldn't find parent of taxID $a - directly assigned to root.\n";
			return 0 ;
		}
		last if ( $a eq $parent_map{$a} ) ;
		$a = $parent_map{ $a } ;
	}
	return 1 ;
}

sub lca {
  my ($a, $b) = @_;
  return $b if $a eq 0;
  return $a if $b eq 0;
  return $a if $a eq $b;
  my %a_path;
  while ($a ge 1) {
    $a_path{$a} = 1;
    if (!defined $parent_map{$a}) {
      print STDERR "Couldn't find parent of taxID $a - directly assigned to root.\n";
      last;
    }
    last if $a eq $parent_map{$a};
    $a = $parent_map{$a};
  }
  while ($b > 1) {
    return $b if (defined $a_path{$b});
    if (!defined $parent_map{$b}) {
      print STDERR "Couldn't find parent of taxID $b - directly assigned to root.\n";
      last;
    }
    last if $b eq $parent_map{$b};
    $b = $parent_map{$b};
  }
  return 1;
}

sub rank_code {
  my $rank = shift;
  for ($rank) {
    $_ eq "species" and return "S";
    $_ eq "genus" and return "G";
    $_ eq "family" and return "F";
    $_ eq "order" and return "O";
    $_ eq "class" and return "C";
    $_ eq "phylum" and return "P";
    $_ eq "kingdom" and return "K";
    $_ eq "superkingdom" and return "D";
  }
  return "-";
}

sub dfs_summation {
  my $node = shift;
  my $children = $child_lists{$node};
  if ($children) {
    for my $child (@$children) {
      dfs_summation($child);
      $clade_counts{$node} += ($clade_counts{$child} || 0);
    }
  }
}

sub load_taxonomy {

  print STDERR "Loading names file ...\n";
  open NAMES, "-|", "$CWD/centrifuge-inspect --name-table $centrifuge_index"
    or die "$PROG: can't open names file: $!\n";
  while (<NAMES>) {
    chomp;
    s/\t\|$//;
    my @fields = split /\t/;
    my ($node_id, $name) = @fields[0,1];
    $name_map{$node_id} = $name;
  }
  close NAMES;

  print STDERR "Loading nodes file ...\n";
  open NODES, "-|", "$CWD/centrifuge-inspect --taxonomy-tree $centrifuge_index"
    or die "$PROG: can't open nodes file: $!\n";
  while (<NODES>) {
    chomp;
    my @fields = split /\t\|\t/;
    my ($node_id, $parent_id, $rank) = @fields[0,1,2];
    if ($node_id == 1) {
      $parent_id = 0;
    }
    $child_lists{$parent_id} ||= [];
    push @{ $child_lists{$parent_id} }, $node_id;
    $rank_map{$node_id} = $rank;
    $parent_map{$node_id} = $parent_id;
  }
  close NODES;
}
