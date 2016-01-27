# Centrifuge
Classifier for metagenomic sequences

[Centrifuge] is a novel microbial classification engine that enables
rapid, accurate and sensitive labeling of reads and quantification of
species on desktop computers.  The system uses a novel indexing scheme
based on the Burrows-Wheeler transform (BWT) and the Ferragina-Manzini
(FM) index, optimized specifically for the metagenomic classification
problem. Centrifuge requires a relatively small index (4.7 GB for all
complete bacterial and viral genomes plus the human genome) and
classifies sequences at very high speed, allowing it to process the
millions of reads from a typical high-throughput DNA sequencing run
within a few minutes.  Together these advances enable timely and
accurate analysis of large metagenomics data sets on conventional
desktop computers

The Centrifuge hompage is  http://www.ccb.jhu.edu/software/centrifuge

The Centrifuge poster is available at http://www.ccb.jhu.edu/people/infphilo/data/Centrifuge-poster.pdf

For more details on installing and running Centrifuge, look at MANUAL

## Quick guide
### Installation from source

    git clone https://github.com/infphilo/centrifuge
    cd centrifuge
    make

### Building indexes

We provide several indexes on the Centrifuge homepage at http://www.ccb.jhu.edu/software/centrifuge.
To create your own, you need to download and index the genomes.

#### Downloading genomes
[Centrifuge] needs the NCBI taxonomy for the naming and hierarchy. To download the taxonomy and 
GI to taxonomy mapping, as well as all complete viral, archaeal and bacterial genomes from RefSeq, 
use the following commands:

    scripts/DownloadTaxonomy.sh 
    scripts/DownloadGenomes.sh -g refseq -a 'Complete Genome' -l library bacteria archaea viral
