/*
 * Copyright 2015, Daehwan Kim <infphilo@gmail.com>
 *
 * This file is part of Centrifuge.
 *
 * Centrifuge is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * Centrifuge is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with Centrifuge.  If not, see <http://www.gnu.org/licenses/>.
 */

#include <iostream>
#include <fstream>
#include <string>
#include <cassert>
#include <getopt.h>
#include "assert_helpers.h"
#include "endian_swap.h"
#include "bt2_idx.h"
#include "bt2_io.h"
#include "bt2_util.h"
#include "formats.h"
#include "sequence_io.h"
#include "tokenize.h"
#include "timer.h"
#include "ref_read.h"
#include "filebuf.h"
#include "reference.h"
#include "ds.h"
#include "aligner_sw.h"

/**
 * \file Driver for the bowtie-build indexing tool.
 */

// Build parameters
int verbose;
static int sanityCheck;
static int format;
static TIndexOffU bmax;
static TIndexOffU bmaxMultSqrt;
static uint32_t bmaxDivN;
static int dcv;
static int noDc;
static int entireSA;
static int seed;
static int showVersion;
//   Ebwt parameters
static int32_t lineRate;
static int32_t linesPerSide;
static int32_t offRate;
static int32_t ftabChars;
static int32_t localOffRate;
static int32_t localFtabChars;
static int  bigEndian;
static bool nsToAs;
static bool doSaFile;  // make a file with just the suffix array in it
static bool doBwtFile; // make a file with just the BWT string in it
static bool autoMem;
static bool packed;
static bool writeRef;
static bool justRef;
static bool reverseEach;
static string wrapper;
static int across;

static void resetOptions() {
	verbose        = true;  // be talkative (default)
	sanityCheck    = 0;     // do slow sanity checks
	format         = FASTA; // input sequence format
	bmax           = OFF_MASK; // max blockwise SA bucket size
	bmaxMultSqrt   = OFF_MASK; // same, as multplier of sqrt(n)
	bmaxDivN       = 4;          // same, as divisor of n
	dcv            = 1024;  // bwise SA difference-cover sample sz
	noDc           = 0;     // disable difference-cover sample
	entireSA       = 0;     // 1 = disable blockwise SA
	seed           = 0;     // srandom seed
	showVersion    = 0;     // just print version and quit?
	//   Ebwt parameters
	lineRate       = Ebwt<TIndexOffU>::default_lineRate;
	linesPerSide   = 1;  // 1 64-byte line on a side
	offRate        = 4;  // sample 1 out of 16 SA elts
	ftabChars      = 10; // 10 chars in initial lookup table
    localOffRate   = 3;
    localFtabChars = 6;
	bigEndian      = 0;  // little endian
	nsToAs         = false; // convert reference Ns to As prior to indexing
    doSaFile       = false; // make a file with just the suffix array in it
    doBwtFile      = false; // make a file with just the BWT string in it
	autoMem        = true;  // automatically adjust memory usage parameters
	packed         = false; //
	writeRef       = true;  // write compact reference to .3.bt2/.4.bt2
	justRef        = false; // *just* write compact reference, don't index
	reverseEach    = false;
    across         = 60; // number of characters across in FASTA output
    wrapper.clear();
}

// Argument constants for getopts
enum {
	ARG_BMAX = 256,
	ARG_BMAX_MULT,
	ARG_BMAX_DIV,
	ARG_DCV,
	ARG_SEED,
	ARG_CUTOFF,
	ARG_PMAP,
	ARG_NTOA,
	ARG_USAGE,
	ARG_REVERSE_EACH,
    ARG_SA,
	ARG_WRAPPER,
    ARG_LOCAL_OFFRATE,
    ARG_LOCAL_FTABCHARS
};

/**
 * Print a detailed usage message to the provided output stream.
 */
static void printUsage(ostream& out) {
	out << "Centrifuge version " << string(CENTRIFUGE_VERSION).c_str() << " by Daehwan Kim (infphilo@gmail.com, http://www.ccb.jhu.edu/people/infphilo)" << endl;
    
#ifdef BOWTIE_64BIT_INDEX
	string tool_name = "hisat-build-l";
#else
	string tool_name = "hisat-build-s";
#endif
	if(wrapper == "basic-0") {
		tool_name = "hisat-build";
	}
    
	out << "Usage: hisat2-build [options]* <reference_in> <bt2_index_base>" << endl
	    << "    reference_in            comma-separated list of files with ref sequences" << endl
	    << "    hisat_index_base          write " << gEbwt_ext << " data to files with this dir/basename" << endl
        << "Options:" << endl
        << "    -c                      reference sequences given on cmd line (as" << endl
        << "                            <reference_in>)" << endl;
    if(wrapper == "basic-0") {
        out << "    --large-index           force generated index to be 'large', even if ref" << endl
		<< "                            has fewer than 4 billion nucleotides" << endl;
	}
    out << "    -a/--noauto             disable automatic -p/--bmax/--dcv memory-fitting" << endl
	    << "    -p/--packed             use packed strings internally; slower, uses less mem" << endl
	    << "    --bmax <int>            max bucket sz for blockwise suffix-array builder" << endl
	    << "    --bmaxdivn <int>        max bucket sz as divisor of ref len (default: 4)" << endl
	    << "    --dcv <int>             diff-cover period for blockwise (default: 1024)" << endl
	    << "    --nodc                  disable diff-cover (algorithm becomes quadratic)" << endl
	    << "    -r/--noref              don't build .3/.4.bt2 (packed reference) portion" << endl
	    << "    -3/--justref            just build .3/.4.bt2 (packed reference) portion" << endl
	    << "    -o/--offrate <int>      SA is sampled every 2^offRate BWT chars (default: 5)" << endl
	    << "    -t/--ftabchars <int>    # of chars consumed in initial lookup (default: 10)" << endl
        << "    --localoffrate <int>    SA (local) is sampled every 2^offRate BWT chars (default: 3)" << endl
        << "    --localftabchars <int>  # of chars consumed in initial lookup in a local index (default: 6)" << endl
	    << "    --seed <int>            seed for random number generator" << endl
	    << "    -q/--quiet              verbose output (for debugging)" << endl
	    << "    -h/--help               print detailed description of tool and its options" << endl
	    << "    --usage                 print this usage message" << endl
	    << "    --version               print version information and quit" << endl
	    ;
    
    if(wrapper.empty()) {
		cerr << endl
        << "*** Warning ***" << endl
        << "'" << tool_name << "' was run directly.  It is recommended "
        << "that you run the wrapper script 'bowtie2-build' instead."
        << endl << endl;
	}
}

static const char *short_options = "qraph?nscfl:i:o:t:h:3C";

static struct option long_options[] = {
	{(char*)"quiet",          no_argument,       0,            'q'},
	{(char*)"sanity",         no_argument,       0,            's'},
	{(char*)"packed",         no_argument,       0,            'p'},
	{(char*)"little",         no_argument,       &bigEndian,   0},
	{(char*)"big",            no_argument,       &bigEndian,   1},
	{(char*)"bmax",           required_argument, 0,            ARG_BMAX},
	{(char*)"bmaxmultsqrt",   required_argument, 0,            ARG_BMAX_MULT},
	{(char*)"bmaxdivn",       required_argument, 0,            ARG_BMAX_DIV},
	{(char*)"dcv",            required_argument, 0,            ARG_DCV},
	{(char*)"nodc",           no_argument,       &noDc,        1},
	{(char*)"seed",           required_argument, 0,            ARG_SEED},
	{(char*)"entiresa",       no_argument,       &entireSA,    1},
	{(char*)"version",        no_argument,       &showVersion, 1},
	{(char*)"noauto",         no_argument,       0,            'a'},
	{(char*)"noblocks",       required_argument, 0,            'n'},
	{(char*)"linerate",       required_argument, 0,            'l'},
	{(char*)"linesperside",   required_argument, 0,            'i'},
	{(char*)"offrate",        required_argument, 0,            'o'},
	{(char*)"ftabchars",      required_argument, 0,            't'},
    {(char*)"localoffrate",   required_argument, 0,            ARG_LOCAL_OFFRATE},
	{(char*)"localftabchars", required_argument, 0,            ARG_LOCAL_FTABCHARS},
	{(char*)"help",           no_argument,       0,            'h'},
	{(char*)"ntoa",           no_argument,       0,            ARG_NTOA},
	{(char*)"justref",        no_argument,       0,            '3'},
	{(char*)"noref",          no_argument,       0,            'r'},
	{(char*)"color",          no_argument,       0,            'C'},
    {(char*)"sa",             no_argument,       0,            ARG_SA},
	{(char*)"reverse-each",   no_argument,       0,            ARG_REVERSE_EACH},
	{(char*)"usage",          no_argument,       0,            ARG_USAGE},
    {(char*)"wrapper",        required_argument, 0,            ARG_WRAPPER},
	{(char*)0, 0, 0, 0} // terminator
};

/**
 * Parse an int out of optarg and enforce that it be at least 'lower';
 * if it is less than 'lower', then output the given error message and
 * exit with an error and a usage message.
 */
template<typename T>
static int parseNumber(T lower, const char *errmsg) {
	char *endPtr= NULL;
	T t = (T)strtoll(optarg, &endPtr, 10);
	if (endPtr != NULL) {
		if (t < lower) {
			cerr << errmsg << endl;
			printUsage(cerr);
			throw 1;
		}
		return t;
	}
	cerr << errmsg << endl;
	printUsage(cerr);
	throw 1;
	return -1;
}

/**
 * Read command-line arguments
 */
static void parseOptions(int argc, const char **argv) {
	int option_index = 0;
	int next_option;
	do {
		next_option = getopt_long(
			argc, const_cast<char**>(argv),
			short_options, long_options, &option_index);
		switch (next_option) {
            case ARG_WRAPPER:
				wrapper = optarg;
				break;
			case 'f': format = FASTA; break;
			case 'c': format = CMDLINE; break;
			case 'p': packed = true; break;
			case 'C':
				cerr << "Error: -C specified but Bowtie 2 does not support colorspace input." << endl;
				throw 1;
				break;
			case 'l':
				lineRate = parseNumber<int>(3, "-l/--lineRate arg must be at least 3");
				break;
			case 'i':
				linesPerSide = parseNumber<int>(1, "-i/--linesPerSide arg must be at least 1");
				break;
			case 'o':
				offRate = parseNumber<int>(0, "-o/--offRate arg must be at least 0");
				break;
            case ARG_LOCAL_OFFRATE:
                localOffRate = parseNumber<int>(0, "-o/--localoffrate arg must be at least 0");
                break;
			case '3':
				justRef = true;
				break;
			case 't':
				ftabChars = parseNumber<int>(1, "-t/--ftabChars arg must be at least 1");
				break;
            case ARG_LOCAL_FTABCHARS:
				localFtabChars = parseNumber<int>(1, "-t/--localftabchars arg must be at least 1");
				break;
			case 'n':
				// all f-s is used to mean "not set", so put 'e' on end
				bmax = 0xfffffffe;
				break;
			case 'h':
			case ARG_USAGE:
				printUsage(cout);
				throw 0;
				break;
			case ARG_BMAX:
				bmax = parseNumber<TIndexOffU>(1, "--bmax arg must be at least 1");
				bmaxMultSqrt = OFF_MASK; // don't use multSqrt
				bmaxDivN = 0xffffffff;     // don't use multSqrt
				break;
			case ARG_BMAX_MULT:
				bmaxMultSqrt = parseNumber<TIndexOffU>(1, "--bmaxmultsqrt arg must be at least 1");
				bmax = OFF_MASK;     // don't use bmax
				bmaxDivN = 0xffffffff; // don't use multSqrt
				break;
			case ARG_BMAX_DIV:
				bmaxDivN = parseNumber<uint32_t>(1, "--bmaxdivn arg must be at least 1");
				bmax = OFF_MASK;         // don't use bmax
				bmaxMultSqrt = OFF_MASK; // don't use multSqrt
				break;
			case ARG_DCV:
				dcv = parseNumber<int>(3, "--dcv arg must be at least 3");
				break;
			case ARG_SEED:
				seed = parseNumber<int>(0, "--seed arg must be at least 0");
				break;
			case ARG_REVERSE_EACH:
				reverseEach = true;
				break;
            case ARG_SA:
                doSaFile = true;
                break;
			case ARG_NTOA: nsToAs = true; break;
			case 'a': autoMem = false; break;
			case 'q': verbose = false; break;
			case 's': sanityCheck = true; break;
			case 'r': writeRef = false; break;

			case -1: /* Done with options. */
				break;
			case 0:
				if (long_options[option_index].flag != 0)
					break;
			default:
				printUsage(cerr);
				throw 1;
		}
	} while(next_option != -1);
	if(bmax < 40) {
		cerr << "Warning: specified bmax is very small (" << bmax << ").  This can lead to" << endl
		     << "extremely slow performance and memory exhaustion.  Perhaps you meant to specify" << endl
		     << "a small --bmaxdivn?" << endl;
	}
}

EList<string> filesWritten;

static void print_fasta_record(
                               ostream& fout,
                               const string& defline,
                               const SString<char>& seq,
                               size_t len)
{
    fout << ">";
    fout << defline.c_str() << endl;
    
    if(across > 0) {
        size_t i = 0;
        while (i + across < len)
        {
            for(size_t j = 0; j < across; j++) {
                int base = seq.get(i + j);
                assert_lt(base, 4);
                fout << "ACGTN"[base];
            }
            fout << endl;
            i += across;
        }
        if (i < len) {
            for(size_t j = i; j < len; j++) {
                int base = seq.get(j);
                assert_lt(base, 4);
                fout << "ACGTN"[base];
            }
            fout << endl;
        }
    } else {
        for(size_t j = 0; j < len; j++) {
            int base = seq.get(j);
            assert_lt(base, 4);
            fout << "ACGTN"[base];
        }
        fout << endl;
    }
}

struct RegionSimilar {
    bool fw;
    size_t pos;
    uint32_t fw_length; // this includes the base at 'pos'
    uint32_t bw_length;
    uint32_t mismatches;
    uint32_t gaps;
    
    void reset() {
        fw = false;
        pos = 0;
        fw_length = bw_length = 0;
        mismatches = gaps = 0;
    }
    
    bool operator< (const RegionSimilar& o) const {
        return pos < o.pos;
    }
};

struct Region {
    size_t pos;
    size_t fw_length;
    bool   low_complexity;
    EList<RegionSimilar, 4> hits;
    
    void reset() {
        pos = fw_length = 0, low_complexity = false;
        hits.clear();
    }
    
    bool operator< (const Region& o) const {
        return pos < o.pos;
    }
};

struct RegionToMerge {
    bool processed;
    EList<pair<uint32_t, uint32_t> > list;
    
    void reset() {
        processed = false;
        list.clear();
    }
};

/**
 * Drive the index construction process and optionally sanity-check the
 * result.
 */
static void driver(
	const string& fafile,
	const string& safile,
	bool packed,
	int reverse)
{
    EList<FileBuf*> is(MISC_CAT);
	bool bisulfite = false;
	RefReadInParams refparams(false, reverse, nsToAs, bisulfite);
    FILE *f = fopen(fafile.c_str(), "r");
    if (f == NULL) {
        cerr << "Error: could not open "<<fafile.c_str() << endl;
        throw 1;
    }
    FileBuf *fb = new FileBuf(f);
    assert(fb != NULL);
    if(fb->peek() == -1 || fb->eof()) {
        cerr << "Warning: Empty fasta file: '" << fafile.c_str() << "'" << endl;
        throw 1;
    }
    assert(!fb->eof());
    assert(fb->get() == '>');
    ASSERT_ONLY(fb->reset());
    assert(!fb->eof());
    is.push_back(fb);
    if(is.empty()) {
		cerr << "Warning: All fasta inputs were empty" << endl;
		throw 1;
	}
	// Vector for the ordered list of "records" comprising the input
	// sequences.  A record represents a stretch of unambiguous
	// characters in one of the input sequences.
	EList<RefRecord> szs(MISC_CAT);
	BitPairReference::szsFromFasta(is, string(), bigEndian, refparams, szs, sanityCheck);
	assert_gt(szs.size(), 0);
    
    EList<string> refnames;
    
    SString<char> s; EList<uint32_t> sa;
    assert_eq(szs.size(), 1);
    size_t jlen = szs[0].len;
    try {
        Timer _t(cerr, "  (1/5) Time reading reference sequence and suffix array: ", verbose);
        
        s.resize(jlen);
        RefReadInParams rpcp = refparams;
        // For each filebuf
        assert_eq(is.size(), 1);
        FileBuf *fb = is[0];
        assert(!fb->eof());
        // For each *fragment* (not necessary an entire sequence) we
        // can pull out of istream l[i]...
        if(!fb->eof()) {
            // Push a new name onto our vector
            refnames.push_back("");
            TIndexOffU distoff = 0;
            fastaRefReadAppend(*fb, true, s, distoff, rpcp, &refnames.back());
        }
        fb->reset();
        assert(!fb->eof());
        
        ifstream in(safile.c_str(), ios::binary);
        uint64_t sa_size = readIndex<uint64_t>(in, false);
        assert_eq(s.length() + 1, sa_size);
        sa.reserveExact(sa_size); sa.clear();
        for(size_t i = 0; i < sa_size; i++) {
            uint64_t sa_pos = readIndex<uint64_t>(in, false);
            sa.push_back((uint32_t)sa_pos);
        }
        assert_eq(s.length() + 1, sa.size());
        
        // Joined reference sequence now in 's'
    } catch(bad_alloc& e) {
        // If we throw an allocation exception in the try block,
        // that means that the joined version of the reference
        // string itself is too larger to fit in memory.  The only
        // alternatives are to tell the user to give us more memory
        // or to try again with a packed representation of the
        // reference (if we haven't tried that already).
        cerr << "Could not allocate space for a joined string of " << jlen << " elements." << endl;
        // There's no point passing this exception on.  The fact
        // that we couldn't allocate the joined string means that
        // --bmax is irrelevant - the user should re-run with
        // ebwt-build-packed
        if(packed) {
            cerr << "Please try running bowtie-build on a computer with more memory." << endl;
        } else {
            cerr << "Please try running bowtie-build in packed mode (-p/--packed) or in automatic" << endl
            << "mode (-a/--auto), or try again on a computer with more memory." << endl;
        }
        if(sizeof(void*) == 4) {
            cerr << "If this computer has more than 4 GB of memory, try using a 64-bit executable;" << endl
            << "this executable is 32-bit." << endl;
        }
        throw 1;
    }
    // Succesfully obtained joined reference string
    assert_eq(s.length(), jlen);
    assert_eq(s.length() % 2, 0);
    size_t sense_seq_len = s.length() / 2;
    assert_geq(sense_seq_len, 2);
    
    SwAligner sw;

    SimpleFunc scoreMin; scoreMin.init(SIMPLE_FUNC_LINEAR, DEFAULT_MIN_CONST, DEFAULT_MIN_LINEAR);
    SimpleFunc nCeil; nCeil.init(SIMPLE_FUNC_LINEAR, 0.0f, std::numeric_limits<double>::max(), 2.0f, 0.1f);
    const int gGapBarrier = 4;
    Scoring sc(
               DEFAULT_MATCH_BONUS,          // constant reward for match
               DEFAULT_MATCH_BONUS_TYPE,     // how to penalize mismatches
               DEFAULT_MM_PENALTY_MAX,       // max mm pelanty
               DEFAULT_MM_PENALTY_MIN,       // min mm pelanty
               scoreMin,                     // min score as function of read len
               nCeil,                        // max # Ns as function of read len
               DEFAULT_N_PENALTY_TYPE,       // how to penalize Ns in the read
               DEFAULT_N_PENALTY,            // constant if N pelanty is a constant
               DEFAULT_N_CAT_PAIR,           // whether to concat mates before N filtering
               DEFAULT_READ_GAP_CONST,       // constant coeff for read gap cost
               DEFAULT_REF_GAP_CONST,        // constant coeff for ref gap cost
               DEFAULT_READ_GAP_LINEAR,      // linear coeff for read gap cost
               DEFAULT_REF_GAP_LINEAR,       // linear coeff for ref gap cost
               gGapBarrier);                 // # rows at top/bot only entered diagonally
    
    size_t tmp_sense_seq_len = sense_seq_len;
    size_t min_kmer = 0;
    while(tmp_sense_seq_len > 0) {
        tmp_sense_seq_len >>= 2;
        min_kmer++;
    }
    //
    min_kmer += 4;
    
    //
    const size_t min_seed_length = min_kmer * 2;
    
    EList<Region> matches;
    {
        EList<uint16_t> prefix_lengths;
        prefix_lengths.resizeExact(sense_seq_len);
        prefix_lengths.fillZero();
        {
            Timer _t(cerr, "  (2/5) Time finding seeds: ", verbose);
            
            // Compress sequences by removing redundant sub-sequences
            size_t last_i1 = 0;
            for(size_t i1 = 0; i1 < sa.size() - 1; i1++) {
                // daehwan - for debugging purposes
                if((i1 + 1) % 50000000 == 0) {
                    cerr << "\t\t" << (i1 + 1) / 1000000 << " million" << endl;
                }
                   
                size_t pos1 = sa[i1];
                if(pos1 == s.length()) continue;
                if(pos1 + min_seed_length >= sense_seq_len) continue;
                // Check if the sequence (defined by cpos1) had been removed
                if(s[pos1] > 3) continue;
                
                // Compare with the following sequences
                bool expanded = false;
                for(size_t i2 = last_i1 + 1; i2 < sa.size(); i2++) {
                    if(i1 == i2) continue;
                    
                    size_t pos2 = sa[i2];
                    if(pos2 == s.length()) continue;
                    // opos2 is relative pos of pos2 on the other strand
                    size_t opos2 = s.length() - pos2 - 1;
                    // cpos2 is canonical pos on the sense strand
                    size_t cpos2 = min(pos2, opos2);
                    // Check if the sequence (defined by cpos2) had been removed
                    if(s[cpos2] > 3) continue;
                    
                    bool fw = pos2 == cpos2;
                    
                    size_t j1 = 0; // includes the base at 'pos1'
                    while(pos1 + j1 < sense_seq_len && pos2 + j1 < (fw ? sense_seq_len : s.length())) {
                        int base1 = s[pos1 + j1];
                        int base2 = s[pos2 + j1];
                        if(base1 > 3 || base2 > 3) break;
                        if(base1 != base2) break;
                        j1++;
                    }
                    size_t j2 = 0; // doesn't include the base at 'pos1'
                    while(j2 <= pos1 && (fw ? 0 : sense_seq_len) + j2 <= pos2) {
                        int base1 = s[pos1 - j2];
                        int base2 = s[pos2 - j2];
                        if(base1 > 3 || base2 > 3) break;
                        if(base1 != base2) break;
                        j2++;
                    }
                    if(j2 > 0) j2--;
                    
                    size_t j = j1 + j2;
                    
                    if(j1 < min_kmer + 10) {
                        if(i2 > i1) break;
                        else        continue;
                    }
                    
                    // Do not proceed if two sequences are not similar
                    if(j < min_seed_length) continue;
                    
                    assert_leq(pos1 + j1, prefix_lengths.size());
                    if(!expanded && j1 <= prefix_lengths[pos1]) continue;
                    
                    if(!expanded) {
                        matches.expand();
                        matches.back().pos = pos1;
                        matches.back().fw_length = j1;
                        for(size_t k = 0; k < j1; k++) {
                            if(prefix_lengths[pos1 + k] < j1 - k) {
                                prefix_lengths[pos1 + k] = j1 - k;
                            }
                        }
                        expanded = true;
                    }
                    
                    if(matches.back().fw_length < j1) {
                        matches.back().fw_length = j1;
                    }
                    
                    matches.back().hits.expand();
                    matches.back().hits.back().reset();
                    matches.back().hits.back().fw = fw;
                    matches.back().hits.back().pos = cpos2;
                    if(fw) {
                        matches.back().hits.back().fw_length = j1;
                        matches.back().hits.back().bw_length = j2;
                    } else {
                        matches.back().hits.back().fw_length = j1 > 0 ? j2 + 1 : 0;
                        matches.back().hits.back().bw_length = j1 > 0 ? j1 - 1 : 0;
                    }
                    
                    if(matches.back().hits.size() >= 20) break;
                }
                
                last_i1 = i1;
                
                if(expanded) {
                    assert_gt(matches.size(), 0);
                    EList<RegionSimilar, 4>& hits = matches.back().hits;
                    if(hits.size() > 1) {
                        hits.sort();
                        size_t cur_pos = 1;
                        for(size_t i = 1; i < hits.size(); i++) {
                            assert_gt(cur_pos, 0);
                            const RegionSimilar& last_region = hits[cur_pos-1];
                            const RegionSimilar& new_region = hits[i];
                            if(last_region.fw == new_region.fw) {
                                if(last_region.fw) {
                                    if(last_region.pos + last_region.fw_length >= new_region.pos) {
                                        continue;
                                    }
                                } else {
                                    if(last_region.pos + last_region.fw_length >= new_region.pos) {
                                        hits[cur_pos-1] = new_region;
                                        continue;
                                    }
                                }
                            }
                            if(cur_pos != i) {
                                assert_lt(cur_pos, hits.size());
                                hits[cur_pos] = new_region;
                            }
                            cur_pos++;
                        }
                        if(cur_pos < hits.size()) {
                            matches.back().low_complexity = true;
                            hits.resizeExact(cur_pos);
                        }
                    }
                }
            }
        }
        
        sa.resizeExact(0);
        
        {
            Timer _t(cerr, "  (3/5) Time sorting seeds and then removing redundant seeds: ", verbose);
            matches.sort();
            
            if(matches.size() > 1) {
                size_t cur_pos = 1;
                for(size_t i = 1; i < matches.size(); i++) {
                    assert_gt(cur_pos, 0);
                    const Region& last_region = matches[cur_pos-1];
                    const Region& new_region = matches[i];
                    if(last_region.low_complexity && last_region.pos + last_region.fw_length > new_region.pos) continue;
                    if(last_region.pos + last_region.fw_length >= new_region.pos + new_region.fw_length) continue;
                    if(cur_pos != i) {
                        assert_lt(cur_pos, matches.size());
                        matches[cur_pos] = new_region;
                    }
                    cur_pos++;
                }
                matches.resizeExact(cur_pos);
            }
        }
    }

    // Print matches
#if 0
    cout << "no. of matches: " << matches.size() << endl << endl;
    for(size_t i = 0; i < matches.size(); i++) {
        const Region& region = matches[i];
        cout << "At " << region.pos << "\t" << region.fw_length << " bps" << endl;
        for(size_t j = 0; j < region.hits.size(); j++) {
            const RegionSimilar& region2 = region.hits[j];
            cout << "\t" << (region2.fw ? "+" : "-") << "\tat " << region2.pos
            << "\t-" << region2.bw_length
            << "\t+" << region2.fw_length << endl;
        }
        cout << endl << endl;
    }
#endif
    
    // daehwan - for debugging purposes
    // const size_t min_sim_length = min_seed_length * 2;
    const size_t min_sim_length = 100;
    
    {
        Timer _t(cerr, "  (4/5) Time merging seeds and masking sequence: ", verbose);
        
        EList<uint8_t> mask;
        mask.resizeExact(sense_seq_len);
        mask.fillZero();
        
        EList<RegionToMerge> merge_list;
        for(size_t i = 0; i < matches.size(); i++) {
            const Region& region = matches[i];
            if(i == 0) {
                for(size_t j = 0; j < region.hits.size(); j++) {
                    merge_list.expand();
                    merge_list.back().reset();
                    merge_list.back().list.expand();
                    merge_list.back().list.back().first = i;
                    merge_list.back().list.back().second = j;
                }
            } else {
                assert_gt(i, 0);
                for(size_t j = 0; j < region.hits.size(); j++) {
                    const RegionSimilar& cmp_region = region.hits[j];
                    bool added = false;
                    for(size_t k = 0; k < merge_list.size(); k++) {
                        RegionToMerge& merge = merge_list[k];
                        uint32_t region_id1 = merge.list.back().first;
                        if(region_id1 >= i) break;
                        uint32_t region_id2 = merge.list.back().second;
                        assert_lt(region_id1, matches.size());
                        const Region& prev_region = matches[region_id1];
                        assert_lt(region_id2, prev_region.hits.size());
                        
                        assert_lt(prev_region.pos, region.pos);
                        size_t gap = region.pos - prev_region.pos;
                        
                        const RegionSimilar& prev_cmp_region = matches[region_id1].hits[region_id2];
                        if(prev_cmp_region.fw != cmp_region.fw) continue;
                        if(prev_cmp_region.pos + cmp_region.bw_length == cmp_region.pos + prev_cmp_region.bw_length &&
                           prev_cmp_region.pos + prev_cmp_region.fw_length == cmp_region.pos + cmp_region.fw_length)
                            continue;
                        
                        size_t cmp_gap = 0;
                        if(cmp_region.fw) {
                            if(prev_cmp_region.pos >= cmp_region.pos) continue;
                            cmp_gap = cmp_region.pos - prev_cmp_region.pos;
                        } else {
                            if(prev_cmp_region.pos <= cmp_region.pos) continue;
                            cmp_gap = prev_cmp_region.pos - cmp_region.pos;
                        }
                        if(cmp_gap + 10 < gap || gap + 10 < cmp_gap) continue;
                        
                        if(prev_region.fw_length + 200 < gap) continue;
                        if(cmp_region.fw) {
                            if(prev_cmp_region.fw_length + 200 < cmp_gap) continue;
                        } else {
                            if(cmp_region.fw_length + 200 < cmp_gap) continue;
                        }
                        
                        added = true;
                        merge.list.expand();
                        merge.list.back().first = i;
                        merge.list.back().second = j;
                    }
                    
                    if(!added) {
                        added = true;
                        merge_list.expand();
                        merge_list.back().reset();
                        merge_list.back().list.expand();
                        merge_list.back().list.back().first = i;
                        merge_list.back().list.back().second = j;
                    }
                }
            }
            
            for(size_t j = 0; j < merge_list.size(); j++) {
                RegionToMerge& merge = merge_list[j];
                uint32_t region_id1 = merge.list.back().first;
                if(i + 1 < matches.size()) {
                    if(region_id1 == i) continue;
                    assert_lt(region_id1, matches.size());
                    const Region& prev_region = matches[region_id1];
                    if(prev_region.pos + 200 > region.pos) continue;
                }
                merge_list[j].processed = true;
                
#if 0
                bool skip_merge = true;
                for(size_t k = 0; k < merge.list.size(); k++) {
                    uint32_t region_id1 = merge.list[k].first;
                    uint32_t region_id2 = merge.list[k].second;
                    assert_lt(region_id1, matches.size());
                    const Region& region = matches[region_id1];
                    assert_lt(region_id2, region.hits.size());
                    const RegionSimilar& sim_region = region.hits[region_id2];
                    assert_lt(region.pos, mask.size()); assert_lt(sim_region.pos, mask.size());
                    if(mask[region.pos] == 0 || mask[sim_region.pos] == 0) {
                        skip_merge = false;
                        break;
                    }
                }
                if(skip_merge) continue;
#endif
                
#if 1
                bool output_merge = merge.list.size() > 1;
                if(!output_merge) {
                    assert_gt(merge.list.size(), 0);
                    uint32_t region_id1 = merge.list[0].first;
                    uint32_t region_id2 = merge.list[0].second;
                    assert_lt(region_id1, matches.size());
                    const Region& region = matches[region_id1];
                    assert_lt(region_id2, region.hits.size());
                    const RegionSimilar& sim_region = region.hits[region_id2];
                    if(sim_region.bw_length + sim_region.fw_length >= min_sim_length) {
                        output_merge = true;
                    }
                }
                if(output_merge) {
                    cout << endl << ":" << endl;
                    for(size_t k = 0; k < merge.list.size(); k++) {
                        uint32_t region_id1 = merge.list[k].first;
                        uint32_t region_id2 = merge.list[k].second;
                        assert_lt(region_id1, matches.size());
                        const Region& region = matches[region_id1];
                        assert_lt(region_id2, region.hits.size());
                        const RegionSimilar& sim_region = region.hits[region_id2];
                        
                        cout << "\t";
                        cout << k << ") at " << region.pos << "\t" << region.fw_length << " bps\t"
                        << (sim_region.fw ? "+" : "-") << "\tat " << sim_region.pos << "\t-" << sim_region.bw_length << "\t+" << sim_region.fw_length << endl;
                    }
                    cout << endl << endl;
                }
#endif
                
                assert_gt(merge.list.size(), 0);
                for(size_t k = merge.list.size() - 1; k > 0; k--) {
                    uint32_t region_id1 = merge.list[k-1].first;
                    uint32_t region_id2 = merge.list[k].first;
                    
                    assert_lt(region_id1, region_id2);
                    assert_lt(region_id2, matches.size());
                    
                    Region& region1 = matches[region_id1];
                    Region& region2 = matches[region_id2];
                    
                    uint32_t cmp_region_id1 = merge.list[k-1].second;
                    assert_lt(cmp_region_id1, region1.hits.size());
                    RegionSimilar& cmp_region1 = region1.hits[cmp_region_id1];
                    
                    uint32_t cmp_region_id2 = merge.list[k].second;
                    assert_lt(cmp_region_id2, region2.hits.size());
                    const RegionSimilar& cmp_region2 = region2.hits[cmp_region_id2];
                    
                    assert_eq(cmp_region1.fw, cmp_region2.fw);
                    const bool fw = cmp_region1.fw;
                    
                    {
                        // Set up query
                        BTString seq;
                        
                        BTDnaString cmp_seq;
                        BTString cmp_qual;
                        
                        size_t query_len, left = region1.pos, right = region2.pos, cmp_left, cmp_right;
                        if(fw) {
                            assert_lt(cmp_region1.pos, cmp_region2.pos);
                            query_len = cmp_region2.pos - cmp_region1.pos + cmp_region2.fw_length + cmp_region1.bw_length;
                            cmp_left = cmp_region1.pos, cmp_right = cmp_region2.pos;
                            
                            assert_gt(cmp_region1.fw_length, 0);
                            left = left + cmp_region1.fw_length - 1;
                            cmp_left = cmp_left + cmp_region1.fw_length - 1;
                            
                            assert_geq(right, cmp_region2.bw_length);
                            right = right - cmp_region2.bw_length;
                            assert_geq(cmp_right, cmp_region2.bw_length);
                            cmp_right = cmp_right - cmp_region2.bw_length;
                            
                        } else {
                            assert_lt(cmp_region2.pos, cmp_region1.pos);
                            query_len = cmp_region1.pos - cmp_region2.pos + cmp_region1.fw_length + cmp_region2.bw_length;
                            cmp_left = cmp_region2.pos, cmp_right = cmp_region1.pos;
                            
                            left = left + cmp_region1.bw_length;
                            assert_gt(cmp_region2.fw_length, 0);
                            cmp_left = cmp_left + cmp_region2.fw_length - 1;
                            
                            assert_geq(right + 1, cmp_region2.fw_length);
                            right = right + 1 - cmp_region2.fw_length;
                            assert_geq(cmp_right, cmp_region1.bw_length);
                            cmp_right = cmp_right - cmp_region1.bw_length;
                        }
                        
#if 0
                        cout << "query length: " << query_len << endl;
                        cout << "left-right: " << left << "\t" << right << endl;
                        cout << "cmp left-right: " << cmp_left << "\t" << cmp_right << endl;
#endif
                        
                        size_t max_diffs = (query_len + 9) / 10;
                        if(max_diffs > cmp_region1.mismatches + cmp_region1.gaps) {
                            max_diffs -= (cmp_region1.mismatches + cmp_region1.gaps);
                        } else {
                            max_diffs = 0;
                        }
                        
                        bool aligned = false, do_swalign = max_diffs > 0;
                        if(left >= right && cmp_left >= cmp_right) {
                            aligned = true;
                        } else if(left >= right) {
                            assert_lt(cmp_left, cmp_right);
                            size_t gap = cmp_right - cmp_left + 1 + left - right;
                            if(gap <= max_diffs) {
                                aligned = true;
                                cmp_region1.gaps += gap;
                            } else {
                                do_swalign = false;
                            }
                            
                        } else if(cmp_left >= cmp_right) {
                            assert_lt(left, right);
                            size_t gap = right - left + 1 + cmp_left - cmp_right;
                            if(gap <= max_diffs) {
                                aligned = true;
                                cmp_region1.gaps += gap;
                            } else {
                                do_swalign = false;
                            }
                        }
                        
                        if(!aligned && do_swalign) {
                            // daehwan - for debugging purposes
                            left -= 5; right += 5;
                            cmp_left -= 5; cmp_right += 5;
                            
                            assert_lt(region1.pos, region2.pos);
                            for(size_t pos = left; pos <= right; pos++) {
                                assert_lt(pos, s.length());
                                seq.append(1 << s[pos]);
                            }
                            
                            for(size_t pos = cmp_left; pos <= cmp_right; pos++) {
                                assert_lt(pos, s.length());
                                cmp_seq.append(s[pos]);
                            }
                            cmp_qual.resize(cmp_seq.length());
                            cmp_qual.fill('I');
                            if(!fw) {
                                cmp_seq.reverseComp();
                                cmp_qual.reverse();
                            }
                            
                            sw.initRead(cmp_seq, cmp_seq, cmp_qual, cmp_qual, 0, cmp_seq.length(), sc);
                            
                            DPRect rect;
                            rect.refl = rect.refl_pretrim = rect.corel = 0;
                            rect.refr = rect.refr_pretrim = rect.corer = seq.length();
                            rect.triml = rect.trimr = 0;
                            rect.maxgap = 10;
                            
                            TAlScore minsc = -max_diffs * 6;
                            if(minsc < 0) {
                                sw.initRef(
                                           true, // fw
                                           0, // refidx
                                           rect,
                                           const_cast<char *>(seq.toZBuf()),
                                           0,
                                           seq.length(),
                                           seq.length(),
                                           sc,
                                           minsc,
                                           true, // enable8
                                           2000, // cminlen
                                           4, // cpow2
                                           false, // doTri
                                           true); // extend);
                                
                                // Perform dynamic programing
                                RandomSource rnd(seed);
                                TAlScore bestCell = std::numeric_limits<TAlScore>::min();
                                if(seq.length() <= 200) {
                                    aligned = sw.align(rnd, bestCell);
                                }
                                
#if 1
                                if(aligned) {
                                    BTDnaString seqstr;
                                    for(size_t bi = 0; bi < seq.length(); bi++) {
                                        seqstr.append(firsts5[(int)seq[bi]]);
                                    }
                                    cout << seqstr << endl;
                                    cout << cmp_seq << endl;
                                    
                                    SwResult res;
                                    res.reset();
                                    sw.nextAlignment(res, minsc, rnd);
                                    res.alres.ned().reverse();
                                    cout << "Succeeded (" << bestCell << "): "; Edit::print(cout, res.alres.ned()); cout << endl;
                                }
#endif
                            }
                        }
                        
                        if(aligned) {
                            if(cmp_region1.fw) {
                                assert_lt(cmp_region1.pos, cmp_region2.pos);
                                cmp_region1.fw_length = cmp_region2.pos - cmp_region1.pos + cmp_region2.fw_length;
                            } else {
                                assert_lt(cmp_region2.pos, cmp_region1.pos);
                                cmp_region1.bw_length = cmp_region1.pos - cmp_region2.pos + cmp_region2.bw_length;
                            }
                        } else {
                            if(cmp_region2.bw_length + cmp_region2.fw_length - 1 >= min_sim_length) {
                                assert_leq(cmp_region2.bw_length, cmp_region2.pos + 1);
                                assert_leq(cmp_region2.pos + cmp_region2.fw_length, sense_seq_len);
                                for(size_t pos = cmp_region2.pos + 1 - cmp_region2.bw_length; pos < cmp_region2.pos + cmp_region2.fw_length; pos ++) {
                                    assert_lt(pos, mask.size());
                                    mask[pos] = 1;
                                }
                            }
                        }
                    }
                    
                    merge.list.pop_back();
                }
                
                assert_eq(merge.list.size(), 1);
                uint32_t region_id = merge.list[0].first;
                assert_lt(region_id, matches.size());
                const Region& region = matches[region_id];
                uint32_t cmp_region_id = merge.list[0].second;
                assert_lt(cmp_region_id, region.hits.size());
                const RegionSimilar& cmp_region = region.hits[cmp_region_id];
                
                if(cmp_region.bw_length + cmp_region.fw_length - 1 >= min_sim_length) {
                    assert_leq(cmp_region.bw_length, cmp_region.pos + 1);
                    assert_leq(cmp_region.pos + cmp_region.fw_length, sense_seq_len);
                    for(size_t pos = cmp_region.pos + 1 - cmp_region.bw_length; pos < cmp_region.pos + cmp_region.fw_length; pos ++) {
                        assert_lt(pos, mask.size());
                        mask[pos] = 1;
                    }
#if 1
                    cout << (cmp_region.fw_length + cmp_region.bw_length - 1) << " bps is masked" << endl;
#endif
                }
            }
            
            size_t cur_pos = 0;
            for(size_t j = 0; j < merge_list.size(); j++) {
                if(merge_list[j].processed) continue;
                if(j != cur_pos) {
                    merge_list[cur_pos] = merge_list[j];
                }
                cur_pos++;
            }
             merge_list.resize(cur_pos);
        }
        
        assert_eq(merge_list.size(), 0);
        assert_eq(mask.size(), sense_seq_len);
        for(size_t i = 0; i < mask.size(); i++){
            if(mask[i] != 0) {
                s.set(4, i);
            }
        }
    }
    
    // Output compressed sequence
    const size_t min_seq_len = 31;
    size_t cur_pos = 0;
    {
        Timer _t(cerr, "  (5/5) Time outputing compressed sequence: ", verbose);

        size_t cur_seq_len = 0;
        for(size_t i = 0; i < sense_seq_len; i++) {
            int base = s[i];
            assert_leq(base, 4);
            if(base < 4) {
                s.set(base, cur_pos);
                cur_pos++;
                cur_seq_len++;
            } else {
                if(cur_seq_len < min_seq_len) {
                    assert_leq(cur_seq_len, i);
                    assert_leq(cur_seq_len, cur_pos);
                    for(size_t j = i - cur_seq_len; j < i; j++) {
                        assert_lt(s[j], 4);
                        s.set(4, j);
                    }
                    cur_pos -= cur_seq_len;
                }
                cur_seq_len = 0;
            }
        }
        print_fasta_record(cout, refnames[0], s, cur_pos);
    }
    
    cerr << endl;
    cerr << "Compressed: " << sense_seq_len << " to " << cur_pos
         << " bps (" << (sense_seq_len - cur_pos) * 100.0 / sense_seq_len << "%)" << endl;
}

static const char *argv0 = NULL;

/**
 * main function.  Parses command-line arguments.
 */
int centrifuge_compress(int argc, const char **argv) {
	string outfile;
	try {
		// Reset all global state, including getopt state
		opterr = optind = 1;
		resetOptions();

		string fafile;
        string safile;
		
		parseOptions(argc, argv);
		argv0 = argv[0];
		if(showVersion) {
			cout << argv0 << " version " << string(CENTRIFUGE_VERSION).c_str() << endl;
			if(sizeof(void*) == 4) {
				cout << "32-bit" << endl;
			} else if(sizeof(void*) == 8) {
				cout << "64-bit" << endl;
			} else {
				cout << "Neither 32- nor 64-bit: sizeof(void*) = " << sizeof(void*) << endl;
			}
			cout << "Built on " << BUILD_HOST << endl;
			cout << BUILD_TIME << endl;
			cout << "Compiler: " << COMPILER_VERSION << endl;
			cout << "Options: " << COMPILER_OPTIONS << endl;
			cout << "Sizeof {int, long, long long, void*, size_t, off_t}: {"
				 << sizeof(int)
				 << ", " << sizeof(long) << ", " << sizeof(long long)
				 << ", " << sizeof(void *) << ", " << sizeof(size_t)
				 << ", " << sizeof(off_t) << "}" << endl;
			return 0;
		}

		// Get input filename
		if(optind >= argc) {
			cerr << "No input sequence or sequence file specified!" << endl;
			printUsage(cerr);
			return 1;
		}
		fafile = argv[optind++];

		// Get output filename
		if(optind >= argc) {
			cerr << "No output file specified!" << endl;
			printUsage(cerr);
			return 1;
		}
		safile = argv[optind++];

		// Optionally summarize
		if(verbose) {
#if 0
			cout << "Settings:" << endl
				 << "  Output files: \"" << outfile.c_str() << ".*." << gEbwt_ext << "\"" << endl
				 << "  Line rate: " << lineRate << " (line is " << (1<<lineRate) << " bytes)" << endl
				 << "  Lines per side: " << linesPerSide << " (side is " << ((1<<lineRate)*linesPerSide) << " bytes)" << endl
				 << "  Offset rate: " << offRate << " (one in " << (1<<offRate) << ")" << endl
				 << "  FTable chars: " << ftabChars << endl
				 << "  Strings: " << (packed? "packed" : "unpacked") << endl
                 << "  Local offset rate: " << localOffRate << " (one in " << (1<<localOffRate) << ")" << endl
                 << "  Local fTable chars: " << localFtabChars << endl
				 ;
			if(bmax == OFF_MASK) {
				cout << "  Max bucket size: default" << endl;
			} else {
				cout << "  Max bucket size: " << bmax << endl;
			}
			if(bmaxMultSqrt == OFF_MASK) {
				cout << "  Max bucket size, sqrt multiplier: default" << endl;
			} else {
				cout << "  Max bucket size, sqrt multiplier: " << bmaxMultSqrt << endl;
			}
			if(bmaxDivN == 0xffffffff) {
				cout << "  Max bucket size, len divisor: default" << endl;
			} else {
				cout << "  Max bucket size, len divisor: " << bmaxDivN << endl;
			}
			cout << "  Difference-cover sample period: " << dcv << endl;
			cout << "  Endianness: " << (bigEndian? "big":"little") << endl
				 << "  Actual local endianness: " << (currentlyBigEndian()? "big":"little") << endl
				 << "  Sanity checking: " << (sanityCheck? "enabled":"disabled") << endl;
	#ifdef NDEBUG
			cout << "  Assertions: disabled" << endl;
	#else
			cout << "  Assertions: enabled" << endl;
	#endif
			cout << "  Random seed: " << seed << endl;
			cout << "  Sizeofs: void*:" << sizeof(void*) << ", int:" << sizeof(int) << ", long:" << sizeof(long) << ", size_t:" << sizeof(size_t) << endl;
			cout << "Input files DNA, " << file_format_names[format].c_str() << ":" << endl;
			for(size_t i = 0; i < infiles.size(); i++) {
				cout << "  " << infiles[i].c_str() << endl;
			}
#endif
        }
		// Seed random number generator
		srand(seed);
		{
			try {
                driver(fafile, safile, false, REF_READ_FORWARD);
            } catch(bad_alloc& e) {
                if(autoMem) {
                    cerr << "Switching to a packed string representation." << endl;
                    packed = true;
                } else {
                    throw e;
                }
            }
		}
		return 0;
	} catch(std::exception& e) {
		cerr << "Error: Encountered exception: '" << e.what() << "'" << endl;
		cerr << "Command: ";
		for(int i = 0; i < argc; i++) cerr << argv[i] << " ";
		cerr << endl;
		return 1;
	} catch(int e) {
		if(e != 0) {
			cerr << "Error: Encountered internal Bowtie 2 exception (#" << e << ")" << endl;
			cerr << "Command: ";
			for(int i = 0; i < argc; i++) cerr << argv[i] << " ";
			cerr << endl;
		}
		return e;
	}
}

/**
 * bowtie-build main function.  It is placed in a separate source file
 * to make it slightly easier to compile as a library.
 *
 * If the user specifies -A <file> as the first two arguments, main
 * will interpret that file as having one set of command-line arguments
 * per line, and will dispatch each batch of arguments one at a time to
 * bowtie-build.
 */
int main(int argc, const char **argv) {
    if(argc > 2 && strcmp(argv[1], "-A") == 0) {
        const char *file = argv[2];
        ifstream in;
        in.open(file);
        char buf[4096];
        int lastret = -1;
        while(in.getline(buf, 4095)) {
            EList<string> args(MISC_CAT);
            args.push_back(string(argv[0]));
            tokenize(buf, " \t", args);
            const char **myargs = (const char**)malloc(sizeof(char*)*args.size());
            for(size_t i = 0; i < args.size(); i++) {
                myargs[i] = args[i].c_str();
            }
            if(args.size() == 1) continue;
            lastret = centrifuge_compress((int)args.size(), myargs);
            free(myargs);
        }
        if(lastret == -1) {
            cerr << "Warning: No arg strings parsed from " << file << endl;
            return 0;
        }
        return lastret;
    } else {
        return centrifuge_compress(argc, argv);
    }
}
