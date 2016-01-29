#
# Makefile
# fbreitwieser, 2016-01-29 13:00
#

THREADS?=1

all: index-build.complete clean

library-download.complete:
	## Downloading and dusting archaeal, bacterial, and viral genomes
	../../scripts/DownloadGenomes.sh -d -l library archaea bacteria viruses
	## Downloading (and not dusting) human and mouse genomes
	../../scripts/DownloadGenomes.sh -l library -a 'Chromosome' -t 9606,10090 -c 'reference genome' vertebrate_mammalian
	## Downloading contaminant sequences
	TODO
	## TODO: Download additional viral genomes
	touch library-download.complete

taxonomy.complete:
	../../scripts/DownloadTaxonomy.sh -l taxonomy
	touch taxonomy.complete

input-sequences.fna: library-download.complete 
	cat library/*/*.fna > input-sequences.fna

index-build.complete: input-sequences.fna taxonomy.complete
	../../centrifuge-build -p $(THREADS) \
		--conversion-table seqid_to_taxid.map --taxonomy-tree taxonomy/nodes.dmp --name-table taxonomy/names.dmp \
		input-sequences.fna bacteria
	touch index-build.complete

clean:
	# Removing input sequences (not needed anymore)
	rm -rf taxonomy
	rm -rf library
	rm input-sequences.fna

# vim:ft=make
#
