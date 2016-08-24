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

The Centrifuge preprint is available at http://biorxiv.org/content/early/2016/05/25/054965.abstract

The Centrifuge poster is available at http://www.ccb.jhu.edu/people/infphilo/data/Centrifuge-poster.pdf

For more details on installing and running Centrifuge, look at MANUAL

## Quick guide
### Installation from source

    git clone https://github.com/infphilo/centrifuge
    cd centrifuge
    make
    sudo make install prefix=/usr/local

### Building indexes

We provide several indexes on the Centrifuge homepage at http://www.ccb.jhu.edu/software/centrifuge.
Centrifuge needs sequence and taxonomy files,  as well as sequence ID to taxonomy ID mapping. 
See the MANUAL files for details. We provide a Makefile that simplifies the building of several
standard and custom indices

    cd indices
    make b+h+v                   # bacterial, human, and viral genomes [~12G]
    make b_compressed            # bacterial genomes compressed at the species level [~4.2G]
    make b_compressed+h+v        # combination of the two above [~8G]
