#
# Makefile
# fbreitwieser, 2016-01-29 13:00
#

THREADS?=1

all: index-build.complete

library-download.complete:
	## Download and dust archaeal, bacterial, and viral genomes
	../../scripts/centrifuge-download -m -d "archaea,bacteria,viral" -P $(THREADS) refseq > microbial_seqid2taxid.map
	## Downloading (and not dusting) human and mouse reference genomes and transcript sequences
	../../scripts/centrifuge-download -d "vertebrate_mammalian" -a "Chromosome" -t 9606,10090 -c 'reference genome' -r -P $(THREADS) refseq > host_seqid2taxid.map
	## Downloading contaminant sequences
	../../scripts/centrifuge-download contaminants > contaminant_seqid2taxid.map
	## TODO: Download additional viral genomes, and check no same header is with different taxIDs in the file
	cat microbial_seqid2taxid.map host_seqid2taxid.map contaminant_seqid2taxid.map > seqid2taxid.map
	touch library-download.complete

taxonomy.complete:
	../../scripts/centrifuge-download taxonomy
	touch taxonomy.complete

input-sequences.fna: library-download.complete 
	cat library/*/*.fna > input-sequences.fna

index-build.complete: input-sequences.fna taxonomy.complete
	../../centrifuge-build -p $(THREADS) \
		--conversion-table seqid2taxid.map --taxonomy-tree taxonomy/nodes.dmp --name-table taxonomy/names.dmp \
		input-sequences.fna bacteria
	touch index-build.complete

clean:
	# Removing input sequences (not needed anymore)
	rm -rf taxonomy
	rm -rf library
	rm -f input-sequences.fna
	rm -f library-download.complete taxonomy.complete
	rm -f *.map

# vim:ft=make
#
