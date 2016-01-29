#
# Makefile
# fbreitwieser, 2016-01-29 13:00
#

THREADS?=1

all: index-build.complete

library-download.complete:
	../../scripts/DownloadGenomes.sh -d -l library archaea bacteria
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
	rm library-download.complete taxonomy.complete

# vim:ft=make
#
