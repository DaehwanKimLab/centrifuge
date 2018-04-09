/*
 * Copyright 2016
 *
 * This file is part of Centrifuge and based on code from Bowtie 2.
 *
 * Bowtie 2 is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * Bowtie 2 is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with Bowtie 2.  If not, see <http://www.gnu.org/licenses/>.
 */

#include <string>
#include <iostream>
#include <getopt.h>
#include <stdexcept>
#include <set>

#include "assert_helpers.h"
#include "endian_swap.h"
#include "hier_idx.h"
#include "reference.h"
#include "ds.h"
#include "hyperloglogplus.h"

using namespace std;

static bool showVersion = false; // just print version and quit?
int verbose             = 0;  // be talkative
static int names_only   = 0;  // just print the sequence names in the index
static int summarize_only = 0; // just print summary of index and quit
static int across       = 60; // number of characters across in FASTA output
static bool refFromEbwt = false; // true -> when printing reference, decode it from Ebwt instead of reading it from BitPairReference
static string wrapper;
static const char *short_options = "vhnsea:";
static int conversion_table = 0;
static int taxonomy_tree = 0;
static int name_table = 0;
static int size_table = 0;
static int count_kmers = 0;

//#define TEST_KMER_COUNTING

enum {
	ARG_VERSION = 256,
    ARG_WRAPPER,
	ARG_USAGE,
    ARG_CONVERSION_TABLE,
    ARG_TAXONOMY_TREE,
    ARG_NAME_TABLE,
    ARG_SIZE_TABLE,
	ARG_COUNT_KMERS
};

static struct option long_options[] = {
	{(char*)"verbose",  no_argument,        0, 'v'},
	{(char*)"version",  no_argument,        0, ARG_VERSION},
	{(char*)"usage",    no_argument,        0, ARG_USAGE},
	{(char*)"names",    no_argument,        0, 'n'},
	{(char*)"summary",  no_argument,        0, 's'},
	{(char*)"help",     no_argument,        0, 'h'},
	{(char*)"across",   required_argument,  0, 'a'},
	{(char*)"ebwt-ref", no_argument,        0, 'e'},
    {(char*)"wrapper",  required_argument,  0, ARG_WRAPPER},
    {(char*)"conversion-table", no_argument,  0, ARG_CONVERSION_TABLE},
    {(char*)"taxonomy-tree",    no_argument,  0, ARG_TAXONOMY_TREE},
    {(char*)"name-table",       no_argument,  0, ARG_NAME_TABLE},
    {(char*)"size-table",       no_argument,  0, ARG_SIZE_TABLE},
	{(char*)"estimate-n-kmers",      no_argument,  0, ARG_COUNT_KMERS},
	{(char*)0, 0, 0, 0} // terminator
};

/**
 * Print a summary usage message to the provided output stream.
 */
static void printUsage(ostream& out) {
	out << "Centrifuge version " << string(CENTRIFUGE_VERSION).c_str() << " by Daehwan Kim (infphilo@gmail.com, http://www.ccb.jhu.edu/people/infphilo)" << endl;
	out
	<< "Usage: centrifuge-inspect [options]* <cf_base>" << endl
	<< "  <cf_base>         cf filename minus trailing .1." << gEbwt_ext << "/.2." << gEbwt_ext << "/.3." << gEbwt_ext << endl
	<< endl
	<< "  By default, prints FASTA records of the indexed nucleotide sequences to" << endl
	<< "  standard out.  With -n, just prints names.  With -s, just prints a summary of" << endl
	<< "  the index parameters and sequences.  With -e, preserves colors if applicable." << endl
	<< endl
	<< "Options:" << endl;
    if(wrapper == "basic-0") {
		//out << "  --large-index      force inspection of the 'large' index, even if a" << endl
	        //<< "                     'small' one is present." << endl;
	}
	out << "  -a/--across <int>  Number of characters across in FASTA output (default: 60)" << endl
	<< "  -n/--names         Print reference sequence names only" << endl
	<< "  -s/--summary       Print summary incl. ref names, lengths, index properties" << endl
	<< "  -e/--bt2-ref       Reconstruct reference from ." << gEbwt_ext << " (slow, preserves colors)" << endl
    << "  --conversion-table Print conversion table" << endl
    << "  --taxonomy-tree    Print taxonomy tree" << endl
    << "  --name-table       Print names corresponding to taxonomic IDs" << endl
    << "  --size-table       Print the lengths of the sequences belonging to the same taxonomic ID" << endl
	<< "  -v/--verbose       Verbose output (for debugging)" << endl
	<< "  -h/--help          print detailed description of tool and its options" << endl
	<< "  --help             print this usage message" << endl
	;
    if(wrapper.empty()) {
		cerr << endl
        << "*** Warning ***" << endl
        << "'centrifuge-inspect-bin' was run directly.  It is recommended "
        << "to use the wrapper script instead."
        << endl << endl;
	}
}

/**
 * Parse an int out of optarg and enforce that it be at least 'lower';
 * if it is less than 'lower', than output the given error message and
 * exit with an error and a usage message.
 */
static int parseInt(int lower, const char *errmsg) {
	long l;
	char *endPtr= NULL;
	l = strtol(optarg, &endPtr, 10);
	if (endPtr != NULL) {
		if (l < lower) {
			cerr << errmsg << endl;
			printUsage(cerr);
			throw 1;
		}
		return (int32_t)l;
	}
	cerr << errmsg << endl;
	printUsage(cerr);
	throw 1;
	return -1;
}

/**
 * Read command-line arguments
 */
static void parseOptions(int argc, char **argv) {
	int option_index = 0;
	int next_option;
	do {
		next_option = getopt_long(argc, argv, short_options, long_options, &option_index);
		switch (next_option) {
            case ARG_WRAPPER:
				wrapper = optarg;
				break;
			case ARG_USAGE:
			case 'h':
				printUsage(cout);
				throw 0;
				break;
			case 'v': verbose = true; break;
			case ARG_VERSION: showVersion = true; break;
            case ARG_CONVERSION_TABLE:
                conversion_table = true;
                break;
            case ARG_TAXONOMY_TREE:
                taxonomy_tree = true;
                break;
            case ARG_NAME_TABLE:
                name_table = true;
                break;
            case ARG_SIZE_TABLE:
                size_table = true;
                break;
			case ARG_COUNT_KMERS:
				count_kmers = true;
				break;
			case 'e': refFromEbwt = true; break;
			case 'n': names_only = true; break;
			case 's': summarize_only = true; break;
			case 'a': across = parseInt(-1, "-a/--across arg must be at least 1"); break;
			case -1: break; /* Done with options. */
			case 0:
				if (long_options[option_index].flag != 0)
					break;
			default:
				printUsage(cerr);
				throw 1;
		}
	} while(next_option != -1);
}

static void print_fasta_record(
	ostream& fout,
	const string& defline,
	const string& seq)
{
	fout << ">";
	fout << defline.c_str() << endl;

	if(across > 0) {
		size_t i = 0;
		while (i + across < seq.length())
		{
			fout << seq.substr(i, across).c_str() << endl;
			i += across;
		}
		if (i < seq.length())
			fout << seq.substr(i).c_str() << endl;
	} else {
		fout << seq.c_str() << endl;
	}
}

/**
 * Counts the number of unique k-mers in the reference sequence
 * that's reconstructed from the index
 */
template<typename index_t, typename TStr>
static uint64_t count_idx_kmers ( Ebwt<index_t>& ebwt)
{
	TStr cat_ref;
	ebwt.restore(cat_ref);
	cerr << "Index loaded" << endl;
#ifdef TEST_KMER_COUNTING
	std::set<uint64_t> my_set;
#endif

	HyperLogLogPlusMinus<uint64_t> kmer_counter(16);
	uint64_t word = 0;
	uint64_t curr_length = 0;
	uint8_t k = 32;

	TIndexOffU curr_ref = OFF_MASK;
	TIndexOffU last_text_off = 0;
	size_t orig_len = cat_ref.length();
	TIndexOffU tlen = OFF_MASK;
	bool first = true;

	for(size_t i = 0; i < orig_len; i++) {
		TIndexOffU tidx = OFF_MASK;
		TIndexOffU textoff = OFF_MASK;
		tlen = OFF_MASK;
		bool straddled = false;
		ebwt.joinedToTextOff(1 /* qlen */, (TIndexOffU)i, tidx, textoff, tlen, true, straddled);

		if (tidx != OFF_MASK && textoff < tlen) {
			if (curr_ref != tidx) {
				// End of the sequence - reset word and counter
				curr_ref = tidx;
				word = 0; curr_length = 0;
				last_text_off = 0;
				first = true;
			}

			TIndexOffU textoff_adj = textoff;
			if(first && textoff > 0) textoff_adj++;
			if (textoff_adj - last_text_off > 1) {
				// there's an N - reset word and counter
				word = 0; curr_length = 0;
			}
			// add another char.
            int bp = (int)cat_ref[i];

            // shift the first two bits off the word
            word = word << 2;
            // put the base-pair code from pos at that position
            word |= bp;
			++curr_length;
			//cerr << "[" << i << "; " << curr_length << "; " << word << ":" << kmer_counter.cardinality()  << "]" << endl;
			if (curr_length >= k) {
				kmer_counter.add(word);
#ifdef TEST_KMER_COUNTING
				my_set.insert(word);
				cerr << " " << kmer_counter.cardinality()  << " vs " << my_set.size() << endl;
#endif
			}

			last_text_off = textoff;
			first = false;

		}
	}
	if (curr_length >= k) {
		kmer_counter.add(word);
#ifdef TEST_KMER_COUNTING
		my_set.insert(word);
#endif
	}

#ifdef TEST_KMER_COUNTING
	cerr << "Exact count: " << my_set.size() << endl;
#endif

	return kmer_counter.cardinality();
}

/**
 * Given output stream, BitPairReference, reference index, name and
 * length, print the whole nucleotide reference with the appropriate
 * number of columns.
 */
static void print_ref_sequence(
	ostream& fout,
	BitPairReference& ref,
	const string& name,
	size_t refi,
	size_t len)
{
	bool newlines = across > 0;
	int myacross = across > 0 ? across : 60;
	size_t incr = myacross * 1000;
	uint32_t *buf = new uint32_t[(incr + 128)/4];
	fout << ">" << name.c_str() << "\n";
	ASSERT_ONLY(SStringExpandable<uint32_t> destU32);
	for(size_t i = 0; i < len; i += incr) {
		size_t amt = min(incr, len-i);
		assert_leq(amt, incr);
		int off = ref.getStretch(buf, refi, i, amt ASSERT_ONLY(, destU32));
		uint8_t *cb = ((uint8_t*)buf) + off;
		for(size_t j = 0; j < amt; j++) {
			if(newlines && j > 0 && (j % myacross) == 0) fout << "\n";
			assert_range(0, 4, (int)cb[j]);
			fout << "ACGTN"[(int)cb[j]];
		}
		fout << "\n";
	}
	delete [] buf;
}

/**
 * Create a BitPairReference encapsulating the reference portion of the
 * index at the given basename.  Iterate through the reference
 * sequences, sending each one to print_ref_sequence to print.
 */
static void print_ref_sequences(
	ostream& fout,
	bool color,
	const EList<string>& refnames,
	const TIndexOffU* plen,
	const string& adjustedEbwtFileBase)
{
	BitPairReference ref(
		adjustedEbwtFileBase, // input basename
		color,                // true -> expect colorspace reference
		false,                // sanity-check reference
		NULL,                 // infiles
		NULL,                 // originals
		false,                // infiles are sequences
		false,                // memory-map
		false,                // use shared memory
		false,                // sweep mm-mapped ref
		verbose,              // be talkative
		verbose);             // be talkative at startup
	assert_eq(ref.numRefs(), refnames.size());
	for(size_t i = 0; i < ref.numRefs(); i++) {
		print_ref_sequence(
			fout,
			ref,
			refnames[i],
			i,
			plen[i] + (color ? 1 : 0));
	}
}

/**
 * Given an index, reconstruct the reference by LF mapping through the
 * entire thing.
 */
template<typename index_t, typename TStr>
static void print_index_sequences(ostream& fout, Ebwt<index_t>& ebwt)
{
	EList<string>* refnames = &(ebwt.refnames());

	TStr cat_ref;
	ebwt.restore(cat_ref);

	HyperLogLogPlusMinus<uint64_t> kmer_counter;
	TIndexOffU curr_ref = OFF_MASK;
	string curr_ref_seq = "";
	TIndexOffU curr_ref_len = OFF_MASK;
	TIndexOffU last_text_off = 0;
	size_t orig_len = cat_ref.length();
	TIndexOffU tlen = OFF_MASK;
	bool first = true;
	for(size_t i = 0; i < orig_len; i++) {
		TIndexOffU tidx = OFF_MASK;
		TIndexOffU textoff = OFF_MASK;
		tlen = OFF_MASK;
		bool straddled = false;
		ebwt.joinedToTextOff(1 /* qlen */, (TIndexOffU)i, tidx, textoff, tlen, true, straddled);

		if (tidx != OFF_MASK && textoff < tlen)
		{
			if (curr_ref != tidx)
			{
				if (curr_ref != OFF_MASK)
				{
					// Add trailing gaps, if any exist
					if(curr_ref_seq.length() < curr_ref_len) {
						curr_ref_seq += string(curr_ref_len - curr_ref_seq.length(), 'N');
					}
					print_fasta_record(fout, (*refnames)[curr_ref], curr_ref_seq);
				}
				curr_ref = tidx;
				curr_ref_seq = "";
				curr_ref_len = tlen;
				last_text_off = 0;
				first = true;
			}

			TIndexOffU textoff_adj = textoff;
			if(first && textoff > 0) textoff_adj++;
			if (textoff_adj - last_text_off > 1)
				curr_ref_seq += string(textoff_adj - last_text_off - 1, 'N');

            curr_ref_seq.push_back("ACGT"[int(cat_ref[i])]);			
			last_text_off = textoff;
			first = false;
		}
	}
	if (curr_ref < refnames->size())
	{
		// Add trailing gaps, if any exist
		if(curr_ref_seq.length() < curr_ref_len) {
			curr_ref_seq += string(curr_ref_len - curr_ref_seq.length(), 'N');
		}
		print_fasta_record(fout, (*refnames)[curr_ref], curr_ref_seq);
	}

}

static char *argv0 = NULL;

template <typename index_t>
static void print_index_sequence_names(const string& fname, ostream& fout)
{
	EList<string> p_refnames;
	readEbwtRefnames<index_t>(fname, p_refnames);
	for(size_t i = 0; i < p_refnames.size(); i++) {
		cout << p_refnames[i].c_str() << endl;
	}
}

/**
 * Print a short summary of what's in the index and its flags.
 */
template <typename index_t>
static void print_index_summary(
	const string& fname,
	ostream& fout)
{
	int32_t flags = Ebwt<index_t>::readFlags(fname);
	bool color = readEbwtColor(fname);
	Ebwt<index_t> ebwt(
					   fname,
					   color,                // index is colorspace
					   -1,                   // don't require entire reverse
					   true,                 // index is for the forward direction
					   -1,                   // offrate (-1 = index default)
					   0,                    // offrate-plus (0 = index default)
					   false,                // use memory-mapped IO
					   false,                // use shared memory
					   false,                // sweep memory-mapped memory
					   true,                 // load names?
					   false,                // load SA sample?
					   false,                // load ftab?
					   false,                // load rstarts?
					   verbose,              // be talkative?
					   verbose,              // be talkative at startup?
					   false,                // pass up memory exceptions?
					   false);               // sanity check?
	EList<string> p_refnames;
	readEbwtRefnames<index_t>(fname, p_refnames);
	cout << "Flags" << '\t' << (-flags) << endl;
	cout << "SA-Sample" << "\t1 in " << (1 << ebwt.eh().offRate()) << endl;
	cout << "FTab-Chars" << '\t' << ebwt.eh().ftabChars() << endl;
	assert_eq(ebwt.nPat(), p_refnames.size());
	for(size_t i = 0; i < p_refnames.size(); i++) {
		cout << "Sequence-" << (i+1)
		     << '\t' << p_refnames[i].c_str()
		     << '\t' << (ebwt.plen()[i] + (color ? 1 : 0))
		     << endl;
	}
}

extern void initializeCntLut();

static void driver(
	const string& ebwtFileBase,
	const string& query)
{
    initializeCntLut();
    
	// Adjust
	string adjustedEbwtFileBase = adjustEbwtBase(argv0, ebwtFileBase, verbose);

	if(names_only) {
		print_index_sequence_names<TIndexOffU>(adjustedEbwtFileBase, cout);
	} else if(summarize_only) {
		print_index_summary<TIndexOffU>(adjustedEbwtFileBase, cout);
    } else {
        // Initialize Ebwt object
		bool color = readEbwtColor(adjustedEbwtFileBase);
		HierEbwt<TIndexOffU, uint16_t> ebwt(
                                            adjustedEbwtFileBase,
                                            color,                // index is colorspace
                                            -1,                   // don't care about entire-reverse
                                            true,                 // index is for the forward direction
                                            -1,                   // offrate (-1 = index default)
                                            0,                    // offrate-plus (0 = index default)
                                            false,                // use memory-mapped IO
                                            false,                // use shared memory
                                            false,                // sweep memory-mapped memory
                                            true,                 // load names?
                                            true,                 // load SA sample?
                                            true,                 // load ftab?
                                            true,                 // load rstarts?
                                            false,                // be talkative?
                                            false,                // be talkative at startup?
                                            false,                // pass up memory exceptions?
                                            false);               // sanity check?        
        
        if(conversion_table) {
            const EList<pair<string, uint64_t> >& uid_to_tid = ebwt.uid_to_tid();
            for(size_t i = 0; i < uid_to_tid.size(); i++) {
                uint64_t tid = uid_to_tid[i].second;
                cout << uid_to_tid[i].first << "\t"
                     << (tid & 0xffffffff);
                tid >>= 32;
                if(tid > 0) {
                    cout << "." << tid;
                }
                cout << endl;
            }
        } else if(taxonomy_tree) {
            const map<uint64_t, TaxonomyNode>& tree = ebwt.tree();
            for(map<uint64_t, TaxonomyNode>::const_iterator itr = tree.begin(); itr != tree.end(); itr++) {
                string rank = get_tax_rank_string(itr->second.rank);
                cout << itr->first << "\t|\t" << itr->second.parent_tid << "\t|\t" << rank << endl;
            }
        } else if(name_table) {
            const std::map<uint64_t, string>& name_map = ebwt.name();
            for(std::map<uint64_t, string>::const_iterator itr = name_map.begin(); itr != name_map.end(); itr++) {
                uint64_t tid = itr->first;
                cout << (tid & 0xffffffff);
                tid >>= 32;
                if(tid > 0) {
                    cout << "." << tid;
                }
                cout << "\t" << itr->second << endl;
            }
        } else if(size_table) {
            const std::map<uint64_t, uint64_t>& size_map = ebwt.size();
            for(std::map<uint64_t, uint64_t>::const_iterator itr = size_map.begin(); itr != size_map.end(); itr++) {
                uint64_t tid = itr->first;
                uint64_t size = itr->second;
                cout << (tid & 0xffffffff);
                tid >>= 32;
                if(tid > 0) {
                    cout << "." << tid;
                }
                cout << "\t" << size << endl;
            }
        } else if (count_kmers) {
        	ebwt.loadIntoMemory(
        	                                -1,     // color
        	                                -1,     // need entire reverse
        	                                true,   // load SA sample
        	                                true,   // load ftab
        	                                true,   // load rstarts
        	                                true,   // load names
        	                                verbose);  // verbose
        	uint64_t n_kmers = count_idx_kmers<TIndexOffU, SString<char> >(ebwt);
        	cout << "Approximate number of kmers in the reference sequence: " << n_kmers << endl;

        } else {
            ebwt.loadIntoMemory(
                                -1,     // color
                                -1,     // need entire reverse
                                true,   // load SA sample
                                true,   // load ftab
                                true,   // load rstarts
                                true,   // load names
                                verbose);  // verbose
            
            // Load whole index into memory
            if(refFromEbwt || true) {
                print_index_sequences<TIndexOffU, SString<char> >(cout, ebwt);
            } else {
                EList<string> refnames;
                readEbwtRefnames<TIndexOffU>(adjustedEbwtFileBase, refnames);
                print_ref_sequences(
                                    cout,
                                    readEbwtColor(ebwtFileBase),
                                    refnames,
                                    ebwt.plen(),
                                    adjustedEbwtFileBase);
            }
        }
		// Evict any loaded indexes from memory
		if(ebwt.isInMemory()) {
			ebwt.evictFromMemory();
		}
	}
}

/**
 * main function.  Parses command-line arguments.
 */
int main(int argc, char **argv) {
	try {
		string ebwtFile;  // read serialized Ebwt from this file
		string query;   // read query string(s) from this file
		EList<string> queries;
		string outfile; // write query results to this file
		argv0 = argv[0];
		parseOptions(argc, argv);
		if(showVersion) {
			cout << argv0 << " version " << CENTRIFUGE_VERSION << endl;
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
			cerr << "No index name given!" << endl;
			printUsage(cerr);
			return 1;
		}
		ebwtFile = argv[optind++];

		// Optionally summarize
		if(verbose) {
			cout << "Input ebwt file: \"" << ebwtFile.c_str() << "\"" << endl;
			cout << "Output file: \"" << outfile.c_str() << "\"" << endl;
			cout << "Local endianness: " << (currentlyBigEndian()? "big":"little") << endl;
#ifdef NDEBUG
			cout << "Assertions: disabled" << endl;
#else
			cout << "Assertions: enabled" << endl;
#endif
		}
		driver(ebwtFile, query);
		return 0;
	} catch(std::exception& e) {
		cerr << "Error: Encountered exception: '" << e.what() << "'" << endl;
		cerr << "Command: ";
		for(int i = 0; i < argc; i++) cerr << argv[i] << " ";
		cerr << endl;
		return 1;
	} catch(int e) {
		if(e != 0) {
			cerr << "Error: Encountered internal Centrifuge exception (#" << e << ")" << endl;
			cerr << "Command: ";
			for(int i = 0; i < argc; i++) cerr << argv[i] << " ";
			cerr << endl;
		}
		return e;
	}
}
