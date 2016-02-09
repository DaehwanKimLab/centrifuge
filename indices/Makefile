#
# Makefile
# fbreitwieser, 2016-01-29 13:00
#

SHELL := /bin/bash

THREADS?=1

KEEP_FILES?=0

#COMPLETE_GENOMES?=archaea bacteria viral protozoa fungi
#CHROMOSOME_GENOMES?=protozoa fungi
#MAMMALIAN_TAXIDS?=9606 10090

define USAGE

Makefile to create common indices to use with Centrifuge.

  make [OPTIONS] TARGET

OPTIONS:
    THREADS=n          Number of threads for downloading, compression and
                       index building

TARGETS:

    b_compressed        Download all bacteria genomes from RefSeq,
                        and compresses them at the species level

    b_compressed+h+v    b_compressed + human genome and transcripts,
                        contaminant sequences from UniVec and EmVec,
                        and all viral genomes

    b+h+v               As above, but with uncompressed bacterial genomes

endef
export USAGE

.SECONDEXPANSION:

all: 
	@echo "$$USAGE"

get_ref_file_names = $(call get_full_ref_names, $(1), \
	$(addprefix all-,$(COMPLETE_GENOMES)) \
	$(addprefix mammalian-reference-,$(MAMMALIAN_TAXIDS)) \
	$(addprefix all-compressed-,$(COMPLETE_GENOMES_COMPRESSED)) \
	$(if $(INCLUDE_CONTAMINANTS),contaminants))

get_full_ref_names = $(addprefix $(REFERENCE_SEQUENCES_DIR)/,$(addsuffix $(1), $(2) ))

DL_DIR=downloaded-seq
#TMP_DIR:=tmp_$(shell date +%Y-%m-%d_%H-%M-%S)
TMP_DIR?=tmp_$(IDX_NAME)
TAXID_SUFFIX:=.map

REFERENCE_SEQUENCES_DIR:=reference-sequences

REFERENCE_SEQUENCES=$(call get_ref_file_names, .fna)
TAXID_MAPS=$(call get_ref_file_names, $(TAXID_SUFFIX))
IDX_NAME?=$(shell basename $(shell dirname $(abspath $(lastword $(MAKEFILE_LIST)))))

THIS_FILE := $(lastword $(MAKEFILE_LIST))

# Join a list with commas
COMMA:=,
EMPTY:=
SPACE:= $(EMPTY) $(EMPTY)
join_w_comma = $(subst $(SPACE),$(COMMA),$(strip $1))

PATH_OK  := $(shell command -v centrifuge-build 2> /dev/null && command -v centrifuge-download 2> /dev/null )
CF_BASE_DIR := $(shell dirname $(shell dirname $(shell dirname $(abspath $(lastword $(MAKEFILE_LIST))))))

error_msg := centrifuge-build or centrifuge-download is not available - please make sure they are in the path.
define n


endef

ifneq ("$(wildcard $(CF_BASE_DIR)/centrifuge-build)","")
error_msg := $(error_msg)$n$nThe following command may solve this problem:$n  export PATH=$$PATH:"$(CF_BASE_DIR)"$n
endif


.PHONY: index index-name index-size .path-ok .mk-tmp


b+h+v: COMPLETE_GENOMES=archaea bacteria viral
b+h+v: MAMMALIAN_TAXIDS=9606
b+h+v: IDX_NAME=b+h+v
b+h+v: INCLUDE_CONTAMINANTS=
b+h+v: REFERENCE_SEQUENCES=$(call get_ref_file_names, .fna)
b+h+V: TAXID_MAPS=$(call get_ref_file_names, $(TAXID_SUFFIX))

b+h+v: index

b_compressed: COMPLETE_GENOMES=archaea
b_compressed: COMPLETE_GENOMES_COMPRESSED=bacteria
b_compressed: INCLUDE_CONTAMINANTS=0
b_compressed: MAMMALIAN_TAXIDS=
b_compressed: IDX_NAME=b_compressed
b_compressed: index

b_compressed+h+v: COMPLETE_GENOMES=archaea viral
b_compressed+h+v: COMPLETE_GENOMES_COMPRESSED=bacteria
b_compressed+h+v: MAMMALIAN_TAXIDS=9606
b_compressed+h+v: IDX_NAME=b_compressed+h+v
b_compressed+h+v: index



.path-ok:
ifndef PATH_OK
    $(error $n$(error_msg))
else
	@echo Found centrifuge-download and centrifuge-build.
endif


.mk-tmp:
ifneq ($(wildcard $(TMP_DIR)),)
	@echo tmp directory exists - cleaning it.
	rm -rf $(TMP_DIR)
endif
	mkdir $(TMP_DIR)

index-name:
	echo $(IDX_NAME)

index-size:
	du -csh $(IDX_NAME).[123].cf

index: $$(IDX_NAME).1.cf

$(REFERENCE_SEQUENCES_DIR):
	mkdir -p $(REFERENCE_SEQUENCES_DIR)

$(TAXID_MAPS): | $(REFERENCE_SEQUENCES_DIR)
	rm $(patsubst %$(TAXID_SUFFIX),%.fna, $@)
	$(MAKE) -f $(THIS_FILE) $(patsubst %$(TAXID_SUFFIX),%.fna, $@)

$(REFERENCE_SEQUENCES_DIR)/mammalian-reference-%.fna: | $(REFERENCE_SEQUENCES_DIR)
	$(MAKE) -f $(THIS_FILE) TMP_DIR=$(TMP_DIR) .mk-tmp
	centrifuge-download -o $(TMP_DIR) -d "vertebrate_mammalian" -a "Chromosome" -t $* -c 'reference genome' -r -P $(THREADS) refseq > \
		$(TMP_DIR)/$(patsubst %.fna,%$(TAXID_SUFFIX), $(notdir $@))
	cat $(TMP_DIR)/vertebrate_mammalian/*.fna > $@.tmp && mv $@.tmp $@
	mv $(TMP_DIR)/$(patsubst %.fna,%$(TAXID_SUFFIX),$(notdir $@)) $(patsubst %.fna,%$(TAXID_SUFFIX),$@)
ifeq (1,$(KEEP_FILES))
    [[ -d $(TMP_DIR)/verterbrate_mammalian ]] || mkdir -p $(TMP_DIR)/verterbrate_mammalian
	mv $(TMP_DIR)/verterbrate_mammalian/* $(DL_DIR)/vertebrate_mammalian
else
	rm -rf $(TMP_DIR)
endif

$(REFERENCE_SEQUENCES_DIR)/all-compressed-%.fna: | $(REFERENCE_SEQUENCES_DIR)
	#$(MAKE) -f $(THIS_FILE) TMP_DIR=$(TMP_DIR) .mk-tmp
	#centrifuge-download -o $(TMP_DIR) -d "$*" -P $(THREADS) refseq > \
	#	$(TMP_DIR)/all-$*.map
	/usr/bin/time -v centrifuge-BuildDatabaseFasta.pl $(TMP_DIR)/$* taxonomy -map $(TMP_DIR)/all-$*.map \
		-o $@.tmp -t $(THREADS) -maxG 50000000 2>&1 | tee centrifuge-compress-$(IDX_NAME).log
	mv $@.tmp.fa $@
	mv $@.tmp.map $(patsubst %.fna,%$(TAXID_SUFFIX),$@)
ifeq (1,$(KEEP_FILES))
    [[ -d $(TMP_DIR)/$* ]] || mkdir -p $(TMP_DIR)/$*
	mv $(TMP_DIR)/$*/* $(DL_DIR)/$*
else
	rm -rf $(TMP_DIR)
endif

$(REFERENCE_SEQUENCES_DIR)/all-%.fna: | $(REFERENCE_SEQUENCES_DIR)
	$(MAKE) -f $(THIS_FILE) TMP_DIR=$(TMP_DIR) .mk-tmp
	@echo Downloading and dust-masking $*
	centrifuge-download -o $(TMP_DIR) -m -d "$*" -P $(THREADS) refseq > \
		$(TMP_DIR)/$(patsubst %.fna,%$(TAXID_SUFFIX),$(notdir $@))
	cat $(TMP_DIR)/$*/*.fna > $@.tmp && mv $@.tmp $@
	mv $(TMP_DIR)/$(patsubst %.fna,%$(TAXID_SUFFIX),$(notdir $@)) $(patsubst %.fna,%$(TAXID_SUFFIX),$@)
ifeq (1,$(KEEP_FILES))
    [[ -d $(TMP_DIR)/$* ]] || mkdir -p $(TMP_DIR)/$*
	mv $(TMP_DIR)/$*/* $(DL_DIR)/$*
else
	rm -rf $(TMP_DIR)
endif

$(REFERENCE_SEQUENCES_DIR)/contaminants.fna: | $(REFERENCE_SEQUENCES_DIR)
	$(MAKE) -f $(THIS_FILE) TMP_DIR=$(TMP_DIR) .mk-tmp
	centrifuge-download -o $(TMP_DIR) contaminants > $(TMP_DIR)/$(patsubst %.fna,%$(TAXID_SUFFIX),$(notdir $@))
	cat $(TMP_DIR)/contaminants/*.fna > $@.tmp && mv $@.tmp $@
	mv $(TMP_DIR)/$(patsubst %.fna,%$(TAXID_SUFFIX),$(notdir $@)) $(patsubst %.fna,%$(TAXID_SUFFIX),$@)
ifeq (1,$(KEEP_FILES))
    [[ -d $(TMP_DIR)/contaminants ]] || mkdir -p $(TMP_DIR)/contaminants
	mv $(TMP_DIR)/contaminants/* $(DL_DIR)/$*
else
	rm -rf $(TMP_DIR)
endif


taxonomy/names.dmp: | taxonomy
taxonomy/nodes.dmp: | taxonomy

taxonomy: | .path-ok
	$(MAKE) -f $(THIS_FILE) TMP_DIR=$(TMP_DIR) .mk-tmp
	centrifuge-download -o $(TMP_DIR)/taxonomy taxonomy
	mv $(TMP_DIR)/taxonomy . && rmdir $(TMP_DIR)

$(IDX_NAME).1.cf: $$(REFERENCE_SEQUENCES) $$(TAXID_MAPS) taxonomy/nodes.dmp taxonomy/names.dmp | .path-ok
	@echo $<
	$(MAKE) -f $(THIS_FILE) TMP_DIR=$(TMP_DIR) .mk-tmp
	/usr/bin/time -v centrifuge-build -p $(THREADS) \
		--conversion-table <(cat $$(TAXID_MAPS)) \
		--taxonomy-tree taxonomy/nodes.dmp --name-table taxonomy/names.dmp \
		$(call join_w_comma,$(REFERENCE_SEQUENCES)) $(TMP_DIR)/$(IDX_NAME) 2>&1 | tee centrifuge-build-$(IDX_NAME).log
	mv $(TMP_DIR)/$(IDX_NAME).*.cf . && rmdir $(TMP_DIR)

nkmers: scripts/estimate-nkmers.pl input-sequences.fna
	scripts/estimate-nkmers.pl input-sequences.fna

clean:
	# Removing input sequences (all required information is in the index)
	rm -rf taxonomy
	rm -rf $(DL_DIR)
	rm -rf $(TMP_DIR)
	rm -f *.map

# vim:ft=make
#
