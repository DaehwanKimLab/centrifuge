#
# Copyright 2014, Daehwan Kim <infphilo@gmail.com>
#
# This file is part of Centrifuge, which is copied and modified from Makefile in the Bowtie2 package.
#
# Centrifuge is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# Centrifuge is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with Centrifuge.  If not, see <http://www.gnu.org/licenses/>.
#
#
# Makefile for centrifuge-bin, centrifuge-build, centrifuge-inspect
#

INC =
GCC_PREFIX = $(shell dirname `which gcc`)
GCC_SUFFIX =
CC = $(GCC_PREFIX)/gcc$(GCC_SUFFIX)
CPP = $(GCC_PREFIX)/g++$(GCC_SUFFIX)
CXX = $(CPP) #-fdiagnostics-color=always
HEADERS = $(wildcard *.h)
BOWTIE_MM = 1
BOWTIE_SHARED_MEM = 0

# Detect Cygwin or MinGW
WINDOWS = 0
CYGWIN = 0
MINGW = 0
ifneq (,$(findstring CYGWIN,$(shell uname)))
	WINDOWS = 1 
	CYGWIN = 1
	# POSIX memory-mapped files not currently supported on Windows
	BOWTIE_MM = 0
	BOWTIE_SHARED_MEM = 0
else
	ifneq (,$(findstring MINGW,$(shell uname)))
		WINDOWS = 1
		MINGW = 1
		# POSIX memory-mapped files not currently supported on Windows
		BOWTIE_MM = 0
		BOWTIE_SHARED_MEM = 0
	endif
endif

MACOS = 0
ifneq (,$(findstring Darwin,$(shell uname)))
	MACOS = 1
endif

POPCNT_CAPABILITY ?= 1
ifeq (1, $(POPCNT_CAPABILITY))
    EXTRA_FLAGS += -DPOPCNT_CAPABILITY
    INC += -I third_party
endif

MM_DEF = 

ifeq (1,$(BOWTIE_MM))
	MM_DEF = -DBOWTIE_MM
endif

SHMEM_DEF = 

ifeq (1,$(BOWTIE_SHARED_MEM))
	SHMEM_DEF = -DBOWTIE_SHARED_MEM
endif

PTHREAD_PKG =
PTHREAD_LIB = 

ifeq (1,$(MINGW))
	PTHREAD_LIB = 
else
	PTHREAD_LIB = -lpthread
endif

SEARCH_LIBS = 
BUILD_LIBS = 
INSPECT_LIBS =

ifeq (1,$(MINGW))
	BUILD_LIBS = 
	INSPECT_LIBS = 
endif

USE_SRA = 0
SRA_DEF =
SRA_LIB =
SERACH_INC = 
ifeq (1,$(USE_SRA))
	SRA_DEF = -DUSE_SRA
	SRA_LIB = -lncbi-ngs-c++-static -lngs-c++-static -lncbi-vdb-static -ldl
	SEARCH_INC += -I$(NCBI_NGS_DIR)/include -I$(NCBI_VDB_DIR)/include
	SEARCH_LIBS += -L$(NCBI_NGS_DIR)/lib64 -L$(NCBI_VDB_DIR)/lib64
endif

LIBS = $(PTHREAD_LIB)

SHARED_CPPS = ccnt_lut.cpp ref_read.cpp alphabet.cpp shmem.cpp \
	edit.cpp bt2_idx.cpp \
	reference.cpp ds.cpp limit.cpp \
	random_source.cpp tinythread.cpp
SEARCH_CPPS = qual.cpp pat.cpp \
	read_qseq.cpp ref_coord.cpp mask.cpp \
	pe.cpp aligner_seed_policy.cpp \
	scoring.cpp presets.cpp \
	simple_func.cpp random_util.cpp outq.cpp

BUILD_CPPS = diff_sample.cpp

CENTRIFUGE_CPPS_MAIN = $(SEARCH_CPPS) centrifuge_main.cpp
CENTRIFUGE_BUILD_CPPS_MAIN = $(BUILD_CPPS) centrifuge_build_main.cpp
CENTRIFUGE_COMPRESS_CPPS_MAIN = $(BUILD_CPPS) \
	aligner_seed.cpp \
	aligner_sw.cpp \
	aligner_cache.cpp \
	dp_framer.cpp \
	aligner_bt.cpp sse_util.cpp \
	aligner_swsse.cpp \
	aligner_swsse_loc_i16.cpp \
	aligner_swsse_ee_i16.cpp \
	aligner_swsse_loc_u8.cpp \
	aligner_swsse_ee_u8.cpp \
	scoring.cpp \
	mask.cpp \
	qual.cpp

CENTRIFUGE_REPORT_CPPS_MAIN=$(BUILD_CPPS)

SEARCH_FRAGMENTS = $(wildcard search_*_phase*.c)
VERSION = $(shell cat VERSION)
GIT_VERSION = $(VERSION)
#GIT_VERSION = $(shell command -v git 2>&1 > /dev/null && git describe --long --tags --dirty --always --abbrev=10 || cat VERSION)

# Convert BITS=?? to a -m flag
BITS=32
ifeq (x86_64,$(shell uname -m))
BITS=64
endif
# msys will always be 32 bit so look at the cpu arch instead.
ifneq (,$(findstring AMD64,$(PROCESSOR_ARCHITEW6432)))
	ifeq (1,$(MINGW))
		BITS=64
	endif
endif
BITS_FLAG =

ifeq (32,$(BITS))
	BITS_FLAG = -m32
endif

ifeq (64,$(BITS))
	BITS_FLAG = -m64
endif
SSE_FLAG=-msse2

DEBUG_FLAGS    = -O0 -g3 $(BIToS_FLAG) $(SSE_FLAG) -std=c++11
DEBUG_DEFS     = -DCOMPILER_OPTIONS="\"$(DEBUG_FLAGS) $(EXTRA_FLAGS)\""
RELEASE_FLAGS  = -O3 $(BITS_FLAG) $(SSE_FLAG) -funroll-loops -g3 -std=c++11
RELEASE_DEFS   = -DCOMPILER_OPTIONS="\"$(RELEASE_FLAGS) $(EXTRA_FLAGS)\""
NOASSERT_FLAGS = -DNDEBUG
FILE_FLAGS     = -D_LARGEFILE_SOURCE -D_FILE_OFFSET_BITS=64 -D_GNU_SOURCE
CFLAGS         = 
#CFLAGS         = -fdiagnostics-color=always

ifeq (1,$(USE_SRA))
	ifeq (1, $(MACOS))
		DEBUG_FLAGS += -mmacosx-version-min=10.6
		RELEASE_FLAGS += -mmacosx-version-min=10.6
	endif
endif


CENTRIFUGE_BIN_LIST = centrifuge-build-bin \
	centrifuge-class \
	centrifuge-inspect-bin

CENTRIFUGE_BIN_LIST_AUX = centrifuge-build-bin-debug \
	centrifuge-class-debug \
	centrifuge-inspect-bin-debug

CENTRIFUGE_SCRIPT_LIST = 	centrifuge \
	centrifuge-build \
	centrifuge-inspect \
	centrifuge-download \
	centrifuge-kreport \
	$(wildcard centrifuge-*.pl)


GENERAL_LIST = $(wildcard scripts/*.sh) \
	$(wildcard scripts/*.pl) \
	$(wildcard *.py) \
	$(wildcard *.pl) \
	doc/manual.inc.html \
	doc/README \
	doc/style.css \
	$(wildcard example/index/*.cf) \
	$(wildcard example/reads/*.fa) \
	$(wildcard example/reference/*) \
	indices/Makefile \
	$(PTHREAD_PKG) \
	$(CENTRIFUGE_SCRIPT_LIST) \
	AUTHORS \
	LICENSE \
	NEWS \
	MANUAL \
	MANUAL.markdown \
	TUTORIAL \
	VERSION

ifeq (1,$(WINDOWS))
	CENTRIFUGE_BIN_LIST := $(CENTRIFUGE_BIN_LIST) centrifuge.bat centrifuge-build.bat centrifuge-inspect.bat 
endif

# This is helpful on Windows under MinGW/MSYS, where Make might go for
# the Windows FIND tool instead.
FIND=$(shell which find)

SRC_PKG_LIST = $(wildcard *.h) \
	$(wildcard *.hh) \
	$(wildcard *.c) \
	$(wildcard *.cpp) \
	$(wildcard third_party/*.h) \
	$(wildcard third_party/*.cpp) \
	doc/strip_markdown.pl \
	Makefile \
	$(GENERAL_LIST)

BIN_PKG_LIST = $(GENERAL_LIST)

.PHONY: all allall both both-debug

all: $(CENTRIFUGE_BIN_LIST)

allall: $(CENTRIFUGE_BIN_LIST) $(CENTRIFUGE_BIN_LIST_AUX)

both: centrifuge-class centrifuge-build-bin

both-debug: centrifuge-class-debug centrifuge-build-bin-debug

DEFS=-fno-strict-aliasing \
     -DCENTRIFUGE_VERSION="\"$(GIT_VERSION)\"" \
     -DBUILD_HOST="\"`hostname`\"" \
     -DBUILD_TIME="\"`date`\"" \
     -DCOMPILER_VERSION="\"`$(CXX) -v 2>&1 | tail -1`\"" \
     $(FILE_FLAGS) \
	 $(CFLAGS) \
     $(PREF_DEF) \
     $(MM_DEF) \
     $(SHMEM_DEF)

#
# centrifuge targets
#

centrifuge-class: centrifuge.cpp $(SEARCH_CPPS) $(SHARED_CPPS) $(HEADERS) $(SEARCH_FRAGMENTS)
	$(CXX) $(RELEASE_FLAGS) $(RELEASE_DEFS) $(EXTRA_FLAGS) \
	$(DEFS) $(SRA_DEF) -DCENTRIFUGE -DBOWTIE2 -DBOWTIE_64BIT_INDEX $(NOASSERT_FLAGS) -Wall \
	$(INC) $(SEARCH_INC) \
	-o $@ $< \
	$(SHARED_CPPS) $(CENTRIFUGE_CPPS_MAIN) \
	$(LIBS) $(SRA_LIB) $(SEARCH_LIBS)

centrifuge-class-debug: centrifuge.cpp $(SEARCH_CPPS) $(SHARED_CPPS) $(HEADERS) $(SEARCH_FRAGMENTS)
	$(CXX) $(DEBUG_FLAGS) $(DEBUG_DEFS) $(EXTRA_FLAGS) \
	$(DEFS) $(SRA_DEF) -DCENTRIFUGE -DBOWTIE2 -DBOWTIE_64BIT_INDEX -Wall \
	$(INC) $(SRA_LIB) $(SEARCH_INC) \
	-o $@ $< \
	$(SHARED_CPPS) $(CENTRIFUGE_CPPS_MAIN) \
	$(LIBS) $(SRA_LIB) $(SEARCH_LIBS)

centrifuge-build-bin: centrifuge_build.cpp $(SHARED_CPPS) $(HEADERS)
	$(CXX) $(RELEASE_FLAGS) $(RELEASE_DEFS) $(EXTRA_FLAGS) \
	$(DEFS) -DCENTRIFUGE -DBOWTIE2 -DBOWTIE_64BIT_INDEX $(NOASSERT_FLAGS) -Wall \
	$(INC) \
	-o $@ $< \
	$(SHARED_CPPS) $(CENTRIFUGE_BUILD_CPPS_MAIN) \
	$(LIBS) $(BUILD_LIBS)

centrifuge-build-bin-debug: centrifuge_build.cpp $(SHARED_CPPS) $(HEADERS)
	$(CXX) $(DEBUG_FLAGS) $(DEBUG_DEFS) $(EXTRA_FLAGS) \
	$(DEFS) -DCENTRIFUGE -DBOWTIE2 -DBOWTIE_64BIT_INDEX -Wall \
	$(INC) \
	-o $@ $< \
	$(SHARED_CPPS) $(CENTRIFUGE_BUILD_CPPS_MAIN) \
	$(LIBS) $(BUILD_LIBS)

centrifuge-compress-bin: centrifuge_compress.cpp $(SHARED_CPPS) $(CENTRIFUGE_COMPRESS_CPPS_MAIN) $(HEADERS)
	$(CXX) $(RELEASE_FLAGS) $(RELEASE_DEFS) $(EXTRA_FLAGS) \
	$(DEFS) -DCENTRIFUGE -DBOWTIE2 -DBOWTIE_64BIT_INDEX $(NOASSERT_FLAGS) -Wall \
	$(INC) \
	-o $@ $< \
	$(SHARED_CPPS) $(CENTRIFUGE_COMPRESS_CPPS_MAIN) \
	$(LIBS) $(BUILD_LIBS)

centrifuge-compress-bin-debug: centrifuge_compress.cpp $(SHARED_CPPS) $(CENTRIFUGE_COMPRESS_CPPS_MAIN) $(HEADERS)
	$(CXX) $(DEBUG_FLAGS) $(DEBUG_DEFS) $(EXTRA_FLAGS) \
	$(DEFS) -DCENTRIFUGE -DBOWTIE2 -DBOWTIE_64BIT_INDEX -Wall \
	$(INC) \
	-o $@ $< \
	$(SHARED_CPPS) $(CENTRIFUGE_COMPRESS_CPPS_MAIN) \
	$(LIBS) $(BUILD_LIBS)

centrifuge-report-bin: centrifuge_report.cpp $(SHARED_CPPS) $(CENTRIFUGE_REPORT_CPPS_MAIN) $(HEADERS)
	$(CXX) $(RELEASE_FLAGS) $(RELEASE_DEFS) $(EXTRA_FLAGS) \
	$(DEFS) -DCENTRIFUGE -DBOWTIE2 -DBOWTIE_64BIT_INDEX $(NOASSERT_FLAGS) -Wall \
	$(INC) \
	-o $@ $< \
	$(SHARED_CPPS) $(CENTRIFUGE_REPORT_CPPS_MAIN) \
	$(LIBS) $(BUILD_LIBS)

centrifuge-report-bin-debug: centrifuge_report.cpp $(SHARED_CPPS) $(CENTRIFUGE_REPORT_CPPS_MAIN) $(HEADERS)
	$(CXX) $(DEBUG_FLAGS) $(DEBUG_DEFS) $(EXTRA_FLAGS) \
	$(DEFS) -DCENTRIFUGE -DBOWTIE2 -DBOWTIE_64BIT_INDEX -Wall \
	$(INC) \
	-o $@ $< \
	$(SHARED_CPPS) $(CENTRIFUGE_REPORT_CPPS_MAIN) \
	$(LIBS) $(BUILD_LIBS)

#centrifuge-RemoveN: centrifuge-RemoveN.cpp 
#	$(CXX) $(RELEASE_FLAGS) $(RELEASE_DEFS) $(EXTRA_FLAGS) \
#	$(DEFS) -DCENTRIFUGE -DBOWTIE2 -DBOWTIE_64BIT_INDEX $(NOASSERT_FLAGS) -Wall \
#	$(INC) \
#	-o $@ $< 


#
# centrifuge-inspect targets
#

centrifuge-inspect-bin: centrifuge_inspect.cpp $(HEADERS) $(SHARED_CPPS)
	$(CXX) $(RELEASE_FLAGS) \
	$(RELEASE_DEFS) $(EXTRA_FLAGS) \
	$(DEFS) -DCENTRIFUGE -DBOWTIE2 -DBOWTIE_64BIT_INDEX -Wall \
	$(INC) -I . \
	-o $@ $< \
	$(SHARED_CPPS) \
	$(LIBS) $(INSPECT_LIBS)

centrifuge-inspect-bin-debug: centrifuge_inspect.cpp $(HEADERS) $(SHARED_CPPS) 
	$(CXX) $(DEBUG_FLAGS) \
	$(DEBUG_DEFS) $(EXTRA_FLAGS) \
	$(DEFS) -DCENTRIFUGE -DBOWTIE2 -DBOWTIE_64BIT_INDEX -Wall \
	$(INC) -I . \
	-o $@ $< \
	$(SHARED_CPPS) \
	$(LIBS) $(INSPECT_LIBS)


centrifuge: ;

centrifuge.bat:
	echo "@echo off" > centrifuge.bat
	echo "perl %~dp0/centrifuge %*" >> centrifuge.bat

centrifuge-build.bat:
	echo "@echo off" > centrifuge-build.bat
	echo "python %~dp0/centrifuge-build %*" >> centrifuge-build.bat

centrifuge-inspect.bat:
	echo "@echo off" > centrifuge-inspect.bat
	echo "python %~dp0/centrifuge-inspect %*" >> centrifuge-inspect.bat


.PHONY: centrifuge-src
centrifuge-src: $(SRC_PKG_LIST)
	mkdir .src.tmp
	mkdir .src.tmp/centrifuge-$(VERSION)
	zip tmp.zip $(SRC_PKG_LIST)
	mv tmp.zip .src.tmp/centrifuge-$(VERSION)
	cd .src.tmp/centrifuge-$(VERSION) ; unzip tmp.zip ; rm -f tmp.zip
	cd .src.tmp ; zip -r centrifuge-$(VERSION)-source.zip centrifuge-$(VERSION)
	cp .src.tmp/centrifuge-$(VERSION)-source.zip .
	rm -rf .src.tmp

.PHONY: centrifuge-bin
centrifuge-bin: $(BIN_PKG_LIST) $(CENTRIFUGE_BIN_LIST) $(CENTRIFUGE_BIN_LIST_AUX) 
	rm -rf .bin.tmp
	mkdir .bin.tmp
	mkdir .bin.tmp/centrifuge-$(VERSION)
	if [ -f centrifuge.exe ] ; then \
		zip tmp.zip $(BIN_PKG_LIST) $(addsuffix .exe,$(CENTRIFUGE_BIN_LIST) $(CENTRIFUGE_BIN_LIST_AUX)) ; \
	else \
		zip tmp.zip $(BIN_PKG_LIST) $(CENTRIFUGE_BIN_LIST) $(CENTRIFUGE_BIN_LIST_AUX) ; \
	fi
	mv tmp.zip .bin.tmp/centrifuge-$(VERSION)
	cd .bin.tmp/centrifuge-$(VERSION) ; unzip tmp.zip ; rm -f tmp.zip
	cd .bin.tmp ; zip -r centrifuge-$(VERSION)-$(BITS).zip centrifuge-$(VERSION)
	cp .bin.tmp/centrifuge-$(VERSION)-$(BITS).zip .
	rm -rf .bin.tmp

.PHONY: doc
doc: doc/manual.inc.html MANUAL

doc/manual.inc.html: MANUAL.markdown
	pandoc -T "Centrifuge Manual" -o $@ \
	 --from markdown --to HTML --toc $^
	perl -i -ne \
	 '$$w=0 if m|^</body>|;print if $$w;$$w=1 if m|^<body>|;' $@

MANUAL: MANUAL.markdown
	perl doc/strip_markdown.pl < $^ > $@

prefix=/usr/local

.PHONY: install
install: all
	mkdir -p $(prefix)/bin
	mkdir -p $(prefix)/share/centrifuge/indices
	install -m 0644 indices/Makefile $(prefix)/share/centrifuge/indices
	install -d -m 0755 $(prefix)/share/centrifuge/doc
	install -m 0644 doc/* $(prefix)/share/centrifuge/doc
	for file in $(CENTRIFUGE_BIN_LIST) $(CENTRIFUGE_SCRIPT_LIST); do \
		install -m 0755 $$file $(prefix)/bin ; \
	done

.PHONY: uninstall
uninstall: all
	for file in $(CENTRIFUGE_BIN_LIST) $(CENTRIFUGE_SCRIPT_LIST); do \
		rm -v $(prefix)/bin/$$file ; \
		rm -v $(prefix)/share/centrifuge; \
	done


.PHONY: clean
clean:
	rm -f $(CENTRIFUGE_BIN_LIST) $(CENTRIFUGE_BIN_LIST_AUX) \
	$(addsuffix .exe,$(CENTRIFUGE_BIN_LIST) $(CENTRIFUGE_BIN_LIST_AUX)) \
	centrifuge-src.zip centrifuge-bin.zip
	rm -f core.* .tmp.head
	rm -rf *.dSYM
push-doc: doc/manual.inc.html
	scp doc/*.*html igm1:/data1/igm3/www/ccb.jhu.edu/html/software/centrifuge/
