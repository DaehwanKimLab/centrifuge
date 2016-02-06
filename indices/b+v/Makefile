#
# Makefile
# fbreitwieser, 2016-01-29 13:00
#

THREADS?=1
IDX_NAME := $(shell basename $(shell dirname $(abspath $(lastword $(MAKEFILE_LIST)))))
PATH_OK  := $(shell command -v centrifuge-build 2> /dev/null && command -v centrifuge-download 2> /dev/null )
CF_BASE_DIR := $(shell dirname $(shell dirname $(shell dirname $(abspath $(lastword $(MAKEFILE_LIST))))))

error_msg := centrifuge-build or centrifuge-download is not available - please make sure they are in the path.
define n


endef

ifneq ("$(wildcard $(CF_BASE_DIR)/centrifuge-build)","")
error_msg := $(error_msg)$n$nThe following command may solve this problem:$n  export PATH=$$PATH:"$(CF_BASE_DIR)"$n
endif

all: index

.PHONY: index-name index-size index .path-ok .no-tmp

.path-ok:
ifndef PATH_OK
    $(error $n$(error_msg))
else
	@echo Found centrifuge-download and centrifuge-build.
endif

.no-tmp:
ifneq ($(wildcard tmp),)
	$(error tmp exists in current directory - please remove it.)
endif

index-name:
	echo $(IDX_NAME)

index-size:
	du -csh $(IDX_NAME).[123].cf

index: $(IDX_NAME).1.cf

seqid2taxid.map: library
	cat microbial_seqid2taxid.map contaminant_seqid2taxid.map > seqid2taxid.map

library: .path-ok .no-tmp
	## Download and dust archaeal, bacterial, and viral genomes
	centrifuge-download -f tmp -m -d "archaea,bacteria,viral" -P $(THREADS) refseq > microbial_seqid2taxid.map
	## Downloading contaminant sequences
	centrifuge-download -f tmp contaminants > contaminant_seqid2taxid.map
	mv tmp/library . && rmdir tmp

taxonomy: .path-ok
	centrifuge-download -f tmp taxonomy
	mv tmp/taxonomy . && rmdir tmp

input-sequences.fna: library 
	cat library/*/*.fna > input-sequences.fna

$(IDX_NAME).1.cf: input-sequences.fna taxonomy seqid2taxid.map .path-ok .no-tmp
	mkdir tmp
	/usr/bin/time -v centrifuge-build -p $(THREADS) \
		--conversion-table seqid2taxid.map \
		--taxonomy-tree taxonomy/nodes.dmp --name-table taxonomy/names.dmp \
		input-sequences.fna tmp/$(IDX_NAME) 2>&1 | tee centrifuge-build.log
	mv tmp/$(IDX_NAME).*.cf . && rmdir tmp

nkmers: scripts/estimate-nkmers.pl input-sequences.fna
	scripts/estimate-nkmers.pl input-sequences.fna

clean:
	# Removing input sequences (all required information is in the index)
	rm -rf taxonomy
	rm -rf library
	rm -f input-sequences.fna
	rm -f library-download.complete taxonomy.complete
	rm -f *.map

# vim:ft=make
#
