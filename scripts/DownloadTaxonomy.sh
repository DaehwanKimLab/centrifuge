#!/bin/bash

set -eu

source `dirname $0`/functions.sh

FTP="ftp://ftp.ncbi.nih.gov"
BASE_DIR="taxonomy"
DOWNLOAD_GI_MAP=0

USAGE="
`basename $0` [-l DIRECTORY] [-g]

OPTIONS
 -l DIRECTORY   Directory to which the taxonomy files are downloaded. Default: $BASE_DIR
 -g             Download GI to taxonomy ID map (~ 8 GB)
"

# arguments: $OPTFIND (current index), $OPTARG (argument for option), $OPTERR (bash-specific)
while getopts "l:g" OPT "$@"; do
    case $OPT in
        l) BASE_DIR="$OPTARG" ;;
        g) DOWNLOAD_GI_MAP=1 ;;
        \?) echo "Invalid option: -$OPTARG" >&2 
            exit 1 
        ;;
        :) echo "Option -$OPTARG requires an argument." >&2
           exit 1
        ;;
    esac
done
shift $((OPTIND-1))

## Taxonomy
if check_or_mkdir_no_fail "$BASE_DIR"; then
    cd "$BASE_DIR"
    wget $FTP/pub/taxonomy/taxdump.tar.gz
    tar -zxvf taxdump.tar.gz nodes.dmp
    tar -zxvf taxdump.tar.gz names.dmp
    rm taxdump.tar.gz
    if [[ "$DOWNLOAD_GI_MAP" == "1" ]]; then
        wget $FTP/pub/taxonomy/gi_taxid_nucl.dmp.gz
        gunzip gi_taxid_nucl.dmp.gz
    fi
    cd -
fi


