#
# Makefile
# fbreitwieser, 2016-01-29 13:00
#

THREADS?=1

all: index-build.complete

seqid-to-taxid.csv: library

library:
	## Download and dust archaeal, bacterial, and viral genomes
	../../scripts/centrifuge-download -m -d "archaea,bacteria,viral" -P $(THREADS) refseq > microbial_seqid-to-taxid.csv
	## Downloading contaminant sequences
	../../scripts/centrifuge-download contaminants > contaminant_seqid-to-taxid.csv
	## Write sequence ID to taxonomy ID map
	cat microbial_seqid-to-taxid.csv contaminant_seqid-to-taxid.csv > seqid-to-taxid.csv

taxonomy:
	../../scripts/centrifuge-download taxonomy

input-sequences.fna: library 
	cat library/*/*.fna > input-sequences.fna

index-build.complete: input-sequences.fna taxonomy seqid-to-taxid.csv
	../../centrifuge-build -p $(THREADS) \
		--conversion-table seqid-to-taxid.csv \
		--taxonomy-tree taxonomy/nodes.dmp --name-table taxonomy/names.dmp \
		input-sequences.fna bacteria
	touch index-build.complete

nkmers: ../../scripts/estimate-nkmers.pl input-sequences.fna
	../../scripts/estimate-nkmers.pl input-sequences.fna

clean:
	# Removing input sequences (not needed anymore)
	rm -rf taxonomy
	rm -rf library
	rm -f input-sequences.fna
	rm -f library-download.complete taxonomy.complete

# vim:ft=make
#
