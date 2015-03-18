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

struct RegionIdentical {
    bool fw;
    size_t pos;
    size_t fw_length; // this includes the base at 'pos'
    size_t bw_length;
};

struct Region {
    size_t pos;
    size_t fw_length;
    EList<RegionIdentical, 4> hits;
    
    bool operator< (const Region& o) const {
        return pos < o.pos;
    }
};

struct RegionToMerge {
    bool processed;
    EList<pair<uint32_t, uint32_t> > list;
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
    
    {
        BTDnaString seq, seqrc;
        BTString qual, qualrc;
        
        // Set up read
        seq.install("TAGACGTCGTCCT", true);
        qual.resize(seq.length());
        for(size_t i = 0; i < seq.length(); i++) {
            qual.set('I', i);
        }
        seqrc = seq;
        seqrc.reverseComp();
        qualrc = qual;
        qualrc.reverse();
        sw.initRead(seq, seqrc, qual, qualrc, 0, seq.length(), sc);
        
        // Set up ref
        BTString ref;
        ref.install("TACACGTACGTCCT");
        for(size_t i = 0; i < ref.length(); i++) {
            int m = asc2dnamask[(int)ref[i]];
            if(m == 15) {
                m = 16; // N
            }
            ref.set(m, i);
        }
        DPRect rect;
        rect.refl = rect.refl_pretrim = rect.corel = 0;
        rect.refr = rect.refr_pretrim = rect.corer = ref.length();
        rect.triml = rect.trimr = 0;
        rect.maxgap = 10;
        
        sw.initRef(
                   true, // fw
                   0, // refidx
                   rect,
                   const_cast<char *>(ref.toZBuf()),
                   0,
                   ref.length(),
                   ref.length(),
                   sc,
                   -18,
                   true, // enable8
                   2000, // cminlen
                   4, // cpow2
                   false, // doTri
                   true); // extend);
        
        // Perform dynamic programing
        RandomSource rnd(seed);
        TAlScore bestCell = std::numeric_limits<TAlScore>::min();
        bool aligned = sw.align(rnd, bestCell);
        if(aligned) {
            int kk = 100;
            kk += 20;
        }
    }
    
    // daehwan - for debugging purposes
    const size_t min_sim_length = 20;
    
    EList<Region> temp_matches;
    EList<uint16_t> prefix_lengths;
    prefix_lengths.resizeExact(sense_seq_len);
    prefix_lengths.fillZero();
    
    // Compress sequences by removing redundant sub-sequences
    size_t last_i1 = 0;
    for(size_t i1 = 0; i1 < sa.size() - 1; i1++) {
        size_t pos1 = sa[i1];
        if(pos1 == s.length()) continue;
        if(pos1 + min_sim_length >= sense_seq_len) continue;
        // Check if the sequence (defined by cpos1) had been removed
        if(s[pos1] > 3) continue;
        
        // daehwan - for debugging purposes
        if(pos1 >= 76 && pos1 <= 98) {
            int kk = 0;
            kk += 20;
        }
        
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
            
            size_t j1 = 0;
            while(pos1 + j1 < sense_seq_len && pos2 + j1 < (fw ? sense_seq_len : s.length())) {
                int base1 = s[pos1 + j1];
                int base2 = s[pos2 + j1];
                if(base1 > 3 || base2 > 3) break;
                if(base1 != base2) break;
                j1++;
            }
            size_t j2 = 0;
            while(j2 <= pos1 && (fw ? 0 : sense_seq_len) + j2 <= pos2) {
                int base1 = s[pos1 - j2];
                int base2 = s[pos2 - j2];
                if(base1 > 3 || base2 > 3) break;
                if(base1 != base2) break;
                j2++;
            }
            size_t j = j1 + j2;
            
            if(j < min<size_t>(min_sim_length, 31)) break;
            
            // Do not proceed if two sequences are not similar
            if(j < min_sim_length) continue;
            
            assert_leq(pos1 + j1, prefix_lengths.size());
            if(!expanded && j1 <= prefix_lengths[pos1]) continue;
            
            if(!expanded) {
                temp_matches.expand();
                temp_matches.back().pos = pos1;
                temp_matches.back().fw_length = j1;
                for(size_t k = 0; k < j1; k++) {
                    if(prefix_lengths[pos1 + k] < j1 - k) {
                        prefix_lengths[pos1 + k] = j1 - k;
                    }
                }
                expanded = true;
            }
            
            temp_matches.back().hits.expand();
            temp_matches.back().hits.back().fw = fw;
            temp_matches.back().hits.back().pos = cpos2;
            if(fw) {
                temp_matches.back().hits.back().fw_length = j1;
                temp_matches.back().hits.back().bw_length = j2 > 0 ? j2 - 1 : 0;
            } else {
                temp_matches.back().hits.back().fw_length = j2;
                temp_matches.back().hits.back().bw_length = j1 > 0 ? j1 - 1 : 0;
            }
            
            if(temp_matches.back().hits.size() >= 20) break;
        }
        
        last_i1 = i1;
    }
    
    temp_matches.sort();
    
    EList<Region> matches;
    if(temp_matches.size() > 0) {
        matches.push_back(temp_matches[0]);
        for(size_t i = 1; i < temp_matches.size(); i++) {
            const Region& last_region = matches.back();
            const Region& new_region = temp_matches[i];
            if(last_region.pos + last_region.fw_length >= new_region.pos + new_region.fw_length) {
                continue;
            }
            matches.push_back(new_region);
        }
    }
    temp_matches.resizeExact(0);
    prefix_lengths.resizeExact(0);

    // Print matches
#if 0
    cout << "no. of matches: " << matches.size() << endl << endl;
    for(size_t i = 0; i < matches.size(); i++) {
        const Region& region = matches[i];
        cout << "At " << region.pos << "\t" << region.fw_length << " bps" << endl;
        for(size_t j = 0; j < region.hits.size(); j++) {
            const RegionIdentical& region2 = region.hits[j];
            cout << "\t" << (region2.fw ? "+" : "-") << "\tat " << region2.pos
            << "\t-" << region2.bw_length
            << "\t+" << region2.fw_length << endl;
        }
        cout << endl << endl;
    }
#endif
    
    EList<RegionToMerge> merge_list;
    for(size_t i = 0; i < matches.size(); i++) {
        const Region& region = matches[i];
        if(i == 0) {
            for(size_t j = 0; j < region.hits.size(); j++) {
                merge_list.expand();
                merge_list.back().processed = false;
                merge_list.back().list.clear();
                merge_list.back().list.expand();
                merge_list.back().list.back().first = i;
                merge_list.back().list.back().second = j;
            }
        } else {
            assert_gt(i, 0);
            for(size_t j = 0; j < region.hits.size(); j++) {
                const RegionIdentical& cmp_region = region.hits[j];
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
                    
                    const RegionIdentical& prev_cmp_region = matches[region_id1].hits[region_id2];
                    if(prev_cmp_region.fw != cmp_region.fw) continue;
                    if(cmp_region.fw) {
                        if(prev_cmp_region.pos >= cmp_region.pos) continue;
                        size_t cmp_gap = cmp_region.pos - prev_cmp_region.pos;
                        if(cmp_gap + 10 < gap || gap + 10 < cmp_gap) continue;
                    } else {
                        if(prev_cmp_region.pos <= cmp_region.pos) continue;
                        size_t cmp_gap = prev_cmp_region.pos - cmp_region.pos;
                        if(cmp_gap + 10 < gap || gap + 10 < cmp_gap) continue;
                    }
                    
                    added = true;
                    merge.list.expand();
                    merge.list.back().first = i;
                    merge.list.back().second = j;
                }
                
                if(!added) {
                    added = true;
                    merge_list.expand();
                    merge_list.back().processed = false;
                    merge_list.back().list.clear();
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
                if(prev_region.pos + 200 < region.pos) continue;
            }
            merge_list[j].processed = true;
            
#if 1
            cout << ":" << endl;
            for(size_t k = 0; k < merge.list.size(); k++) {
                uint32_t region_id1 = merge.list[k].first;
                uint32_t region_id2 = merge.list[k].second;
                assert_lt(region_id1, matches.size());
                const Region& region = matches[region_id1];
                assert_lt(region_id2, region.hits.size());
                const RegionIdentical& sim_region = region.hits[region_id2];
                
                cout << "\t";
                cout << "at " << region.pos << "\t" << region.fw_length << " bps\t"
                << (sim_region.fw ? "+" : "-") << "\tat " << sim_region.pos << "\t-" << sim_region.bw_length << "\t+" << sim_region.fw_length << endl;
            }
            cout << endl << endl;
#endif
            
            assert_gt(merge.list.size(), 0);
            for(size_t k = merge.list.size() - 1; k > 0; k--) {
                uint32_t region_id1 = merge.list[k-1].first;
                uint32_t region_id2 = merge.list[k].first;
                
                assert_lt(region_id1, region_id2);
                assert_lt(region_id2, matches.size());
                
                const Region& region1 = matches[region_id1];
                const Region& region2 = matches[region_id2];
                
                uint32_t cmp_region_id1 = merge.list[k-1].second;
                assert_lt(cmp_region_id1, region1.hits.size());
                const RegionIdentical& cmp_region1 = region1.hits[cmp_region_id1];
                
                uint32_t cmp_region_id2 = merge.list[k].second;
                assert_lt(cmp_region_id2, region2.hits.size());
                const RegionIdentical& cmp_region2 = region2.hits[cmp_region_id2];
                
                {
                    BTDnaString seq, seqrc;
                    BTString qual, qualrc;
                    
                    uint32_t leftext = 10 ,rightext = 10;
                    if(cmp_region1.fw) {
                        if(leftext > cmp_region1.bw_length)  leftext  = cmp_region1.bw_length;
                        if(rightext > cmp_region2.fw_length) rightext = cmp_region2.fw_length;
                    } else {
                        if(leftext > cmp_region2.bw_length)  leftext  = cmp_region2.bw_length;
                        if(rightext > cmp_region1.fw_length) rightext = cmp_region1.fw_length;
                    }
                    
                    // Set up read
                    assert_eq(cmp_region1.fw, cmp_region2.fw);
                    if(cmp_region1.fw) {
                        assert_lt(cmp_region1.pos, cmp_region2.pos);
                        for(size_t pos = cmp_region1.pos - leftext; pos < cmp_region2.pos + rightext; pos++) {
                            assert_lt(pos, s.length());
                            seq.append(s[pos]);
                        }
                    } else {
                        assert_lt(cmp_region2.pos, cmp_region1.pos);
                        for(size_t pos = cmp_region2.pos - leftext; pos < cmp_region1.pos + rightext; pos++) {
                            assert_lt(pos, s.length());
                            seq.append(s[pos]);
                        }
                    }
                    qual.resize(seq.length());
                    qual.fill('I');
                    seqrc = seq;
                    seqrc.reverseComp();
                    qualrc = qual;
                    qualrc.reverse();
                    if(cmp_region1.fw) {
                        sw.initRead(seq, seqrc, qual, qualrc, 0, seq.length(), sc);
                    } else {
                        sw.initRead(seqrc, seq, qualrc, qual, 0, seq.length(), sc);
                    }
                    
                    // Set up ref
                    BTString ref;
                    assert_lt(region1.pos, region2.pos);
                    assert_leq(leftext, region1.pos);
                    for(size_t pos = region1.pos - leftext; pos < region2.pos + rightext; pos++) {
                        assert_lt(pos, s.length());
                        ref.append(1 << s[pos]);
                    }
                    
                    DPRect rect;
                    rect.refl = rect.refl_pretrim = rect.corel = 0;
                    rect.refr = rect.refr_pretrim = rect.corer = ref.length();
                    rect.triml = rect.trimr = 0;
                    rect.maxgap = 10;
                    
                    sw.initRef(
                               true, // fw
                               0, // refidx
                               rect,
                               const_cast<char *>(ref.toZBuf()),
                               0,
                               ref.length(),
                               ref.length(),
                               sc,
                               -18,
                               true, // enable8
                               2000, // cminlen
                               4, // cpow2
                               false, // doTri
                               true); // extend);
                    
                    // Perform dynamic programing
                    RandomSource rnd(seed);
                    TAlScore bestCell = std::numeric_limits<TAlScore>::min();
                    bool aligned = sw.align(rnd, bestCell);
                    if(aligned) {
                        SwResult res;
                        res.reset();
                        sw.nextAlignment(res, -18, rnd);
#if 1
                        BTDnaString refstr;
                        for(size_t bi = 0; bi < ref.length(); bi++) {
                            refstr.append(firsts5[(int)ref[bi]]);
                        }
                        cout << refstr << endl;
                        cout << seq << endl;
                        Edit::print(cout, res.alres.ned()); cout << endl;
#endif
                    }
                }
                
                merge.list.pop_back();
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
    
    exit(1);
    
    // Perform sanity checking
    if(sanityCheck) {
	}
    
    
    // Output compressed sequence
    // daehwan - for debugging purposes
    const size_t min_seq_len = 31;
    size_t cur_pos = 0, cur_seq_len = 0;
    for(size_t i = 0; i < sense_seq_len; i++) {
        int base = s[i];
        assert_lt(base, 4);
        if(base < 4) {
            s.set(base, cur_pos);
            cur_pos++;
            cur_seq_len++;
        } else {
            if(cur_seq_len < min_seq_len) {
                assert_lt(cur_seq_len, i);
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
