#!/bin/bash

set -eu
source `dirname $0`/functions.sh

#########################################################

function download_n_process() {
    set -x
    IFS=$'\t' read -r TAXID FILEPATH <<< "$1"

    NAME=`basename $FILEPATH .gz`
    wget --reject="*index.html*" -nc -qO- $FILEPATH | gunzip -c > $BASEDIR/$DOMAIN/$NAME
    if [[ "$CHANGE_HEADER" == "1" ]]; then
        sed -i "s/^>/>kraken:taxid|$TAXID /" $BASEDIR/$DOMAIN/$NAME
    else 
        sed -n "s/^>\([^ ]\+\)\( .*\)\?/\1 $TAXID/p" $BASEDIR/$DOMAIN/$NAME >> "$BASEDIR/$DOMAIN.map"
    fi

    if [[ "$DO_DUST" == "1" ]]; then
      dustmasker -infmt fasta -in - -level 20 -outfmt fasta | sed '/^>/! s/[^AGCT]/N/g' > $BASEDIR/$DOMAIN/$NAME
    fi

}
export -f download_n_process

function count {
   typeset C=0
   while read L; do
      C=$(( C + 1 ))
      if [ $(( $C % $1 )) -eq 0 ]; then
         echo $C 1>&2
      fi
      echo "$L"
   done
}

## Check if GNU parallel exists
command -v parallel >/dev/null 2>&1 && PARALLEL=1 || PARALLEL=0


ALL_GENOMES="bacteria viral archaea fungi protozoa invertebrate plant vertebrate_mammalian vertebrate_other"
ALL_DATABASES="refseq genbank"
ALL_ASSEMBLY_LEVELS="Complete\ Genome Chromosome Scaffold Contig"

## Option parsing
DATABASE="refseq"
ASSEMBLY_LEVEL="Complete Genome"
BASE_DIR="."

N_PROC=5
CHANGE_HEADER=0
DOWNLOAD_RNA=0
DO_DUST=0

USAGE="
`basename $0` [OPTIONS] DOMAINS

OPTIONS
 -g DATABASE            Either refseq or genbank. Default: '$DATABASE'.
 -a ASSEMBLY_LEVEL      Assembly level. Default: '$ASSEMBLY_LEVEL'.
 -l LIBRARY_DIRECTORY   Directory to which the files are downloaded. Default: '$BASE_DIR'.
 -P NPROC               Number of processes when downloading (uses xargs). Default: '$N_PROC'
 -d                     Mask low-complexity regions using dustmasker. Default: off.
 -t                     Modify header to include taxonomy ID. Default: off.

DOMAINS
  Domains to download. For RefSeq, this is one or more of
  $ALL_GENOMES
"

# arguments: $OPTFIND (current index), $OPTARG (argument for option), $OPTERR (bash-specific)
while getopts "g:a:P:j:rdt" OPT "$@"; do
    case $OPT in
        g) DATABASE=$OPTARG ;;
        a) ASSEMBLY_LEVEL="$OPTARG" ;;
        l) BASE_DIR="$OPTARG" ;;
        P) N_PROC=$OPTARG ;;
        r) DOWNLOAD_RNA=1 ;;
        d) DO_DUST=1 ;;
        t) CHANGE_HEADER=1 ;;
        \?) echo "Invalid option: -$OPTARG" >&2 
            exit 1 
        ;;
        :) echo "Option -$OPTARG requires an argument." >&2
           exit 1
        ;;
    esac
done
shift $((OPTIND-1))

DOMAINS=$*
[[ $DOMAINS != "" ]] || { echo "$USAGE" && exit 1; };

export BASEDIR="$BASE_DIR"
export DO_DUST="$DO_DUST"
export CHANGE_HEADER="$CHANGE_HEADER"

## Fields in the assembly_summary.txt file
TAXID_FIELD=6
SPECIES_TAXID_FIELD=7
VERSION_STATUS_FIELD=11
ASSEMBLY_LEVEL_FIELD=12
FTP_PATH_FIELD=20

echo "Downloading genomes for $DOMAINS at assembly level $ASSEMBLY_LEVEL"
wget -qO- --no-remove-listing ftp://ftp.ncbi.nlm.nih.gov/genomes/$DATABASE/ > /dev/null

for DOMAIN in $DOMAINS; do
    if [[ `grep " $DOMAIN" .listing` ]]; then
        c_echo "Downloading $DOMAIN genomes"
    else
        c_echo "$DOMAIN is not a valid domain - use one of the following:"
        grep '^d' .listing  | sed 's/.* //'
        exit 1
    fi
    export DOMAIN=$DOMAIN
    check_or_mkdir $BASEDIR/$DOMAIN

    ASSEMBLY_SUMMARY_FILE="$BASEDIR/$DOMAIN/assembly_summary.txt"

    echo "Downloading and filtering the assembly_summary.txt file"
    wget -qO- -nc ftp://ftp.ncbi.nlm.nih.gov/genomes/$DATABASE/$DOMAIN/assembly_summary.txt |\
        awk -F "\t" "\$$ASSEMBLY_LEVEL_FIELD==\"$ASSEMBLY_LEVEL\" && \$$VERSION_STATUS_FIELD==\"latest\"" > "$ASSEMBLY_SUMMARY_FILE"

    echo "Download all genomic.fna.gz files"
    cut -f "$TAXID_FIELD,$FTP_PATH_FIELD" "$ASSEMBLY_SUMMARY_FILE" | sed 's#\([^/]*\)$#\1/\1_genomic.fna.gz#' | \
       tr '\n' '\0' | count 100 | xargs -0 -n1 -P $N_PROC bash -c 'download_n_process "$@"' _

    N_EXPECTED=`cat "$ASSEMBLY_SUMMARY_FILE" | wc -l`
    N_GENOMIC=`find $BASEDIR/$DOMAIN -maxdepth 1 -type f -name '*_genomic.fna' | wc -l`
    c_echo "$DOMAIN: expected $N_EXPECTED files, downloaded $N_GENOMIC"

    if [[ "$DOWNLOAD_RNA" == "1" && ! `echo $DOMAIN | egrep 'bacteria|viral|archaea'` ]]; then
    	echo "Downloading all rna.fna.gz files"
        cut -f $TAXID_FIELD,$FTP_PATH_FIELD  "$ASSEMBLY_SUMMARY_FILE"| sed 's#\([^/]*\)$#\1/\1_rna.fna.gz#' |\
            tr '\n' '\0' | count 100 |  xargs -0 -n1 -P $N_PROC bash -c 'download_n_process "$@"' _
        N_RNA=`find $BASEDIR/$DOMAIN -maxdepth 1 -type f -name '*_rna.fna' | wc -l`
        c_echo "$DOMAIN: further downloaded $N_RNA RNAs"
    fi
done

cat $BASEDIR/*.map > seqid_to_taxid.map

