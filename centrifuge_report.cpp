/*
 * centrifuge-build.cpp
 *
 *  Created on: Apr 8, 2015
 *      Author: fbreitwieser
 */

#include<iostream>
#include<fstream>
#include<sstream>
#include<map>
#include<vector>
#include "assert_helpers.h"
#include "sstring.h"
#include "ds.h"      // EList
#include "bt2_idx.h" // Ebwt
#include "bt2_io.h"
#include "util.h"

using namespace std;
typedef TIndexOffU index_t;

static bool startVerbose = true; // be talkative at startup
int gVerbose = 1; // be talkative always
static const char *argv0 = NULL;
static string adjIdxBase;
static bool useShmem				= false; // use shared memory to hold the index
static bool useMm					= false; // use memory-mapped files to hold the index
static bool mmSweep					= false; // sweep through memory-mapped files immediately after mapping
static int offRate					= -1;    // keep default offRate
static bool noRefNames				= false; // true -> print reference indexes; not names
static int sanityCheck				= 0;  // enable expensive sanity checks
/**
 * Print a summary usage message to the provided output stream.
 */
static void printUsage(ostream& out) {
	out << "Centrifuge version " << string(CENTRIFUGE_VERSION).c_str() << " by Daehwan Kim (infphilo@gmail.com, www.ccb.jhu.edu/people/infphilo)" << endl;
	string tool_name = "centrifuge-class";

	out << "Usage: " << endl
	    << "  " << tool_name.c_str() << " <bt2-idx> <centrifuge-out>" << endl
	    << endl
		<<     "  <bt2-idx>  Index filename prefix (minus trailing .X." << gEbwt_ext << ")." << endl
		<<     "  <centrifuge-out>  Centrifuge result file." << endl;
}

template <typename T>
class Pair2ndComparator{
public:
     bool operator()(const pair<T,T> &left, const pair<T,T> &right){
    	 return left.second < right.second;
     }
};


template<typename TStr>
static void driver(
	const char * type,
	const string& bt2indexBase,
	const string& cf_out)
{
	if(gVerbose || startVerbose)  {
		cerr << "Entered driver(): "; logTime(cerr, true);
	}

    //initializeCntLut();  // FB: test commenting

	// Vector of the reference sequences; used for sanity-checking
	EList<SString<char> > names, os;
	EList<size_t> nameLens, seqLens;

	// Initialize Ebwt object and read in header
	if(gVerbose || startVerbose) {
		cerr << "About to initialize fw Ebwt: "; logTime(cerr, true);
	}
	adjIdxBase = adjustEbwtBase(argv0, bt2indexBase, gVerbose);
	Ebwt<index_t> ebwt(
		adjIdxBase,
	    0,        // index is colorspace
		-1,       // fw index
	    true,     // index is for the forward direction
	    /* overriding: */ offRate,
		0, // amount to add to index offrate or <= 0 to do nothing
	    useMm,    // whether to use memory-mapped files
	    useShmem, // whether to use shared memory
	    mmSweep,  // sweep memory-mapped files
	    !noRefNames, // load names?
		true,        // load SA sample?
		true,        // load ftab?
		true,        // load rstarts?
	    gVerbose, // whether to be talkative
	    startVerbose, // talkative during initialization
	    false /*passMemExc*/,
	    sanityCheck);
	//Ebwt<index_t>* ebwtBw = NULL;


	EList<size_t> reflens;
	EList<string> refnames;
	readEbwtRefnames<index_t>(adjIdxBase, refnames);
	map<uint32_t,pair<string,uint64_t> > speciesID_to_name_len;
	for(size_t i = 0; i < ebwt.nPat(); i++) {
		// cerr << "Push back to reflens: "<<  refnames[i] << " is so long: " << ebwt.plen()[i] << endl;
		reflens.push_back(ebwt.plen()[i]);

		// extract numeric id from refName
		const string& refName = refnames[i];
		uint64_t id = extractIDFromRefName(refName);
		uint32_t speciesID = (uint32_t)(id >> 32);

		// extract name from refName
		const string& name_part = refName.substr(refName.find_first_of(' '));

		//uint32_t genusID = (uint32_t)(id & 0xffffffff);
		speciesID_to_name_len[speciesID] = pair<string,uint64_t>(name_part,ebwt.plen()[i]);

	}
//	EList<string> refnames;
//	readEbwtRefnames<index_t>(adjIdxBase, refnames);

	// Read Centrifuge output file
	ifstream infile(cf_out.c_str());

	string line;
	map<uint32_t,uint32_t> species_to_score;

	while (getline(infile,line)) {
		string rd_name;
		uint32_t genusID;
		uint32_t speciesID;
		uint32_t score;
		uint32_t secbest_score;

		istringstream iss(line);
		iss >> rd_name >> genusID >> speciesID >> score >> secbest_score;
		// cerr << rd_name << " -> " << genusID << " -> " << speciesID << " -> " << score << " -> " << secbest_score << "\n";
		species_to_score[speciesID] += score;
	}

	// Sort the species by their score
	vector<pair<uint32_t,uint32_t> > species_to_score_v(species_to_score.begin(), species_to_score.end());

	sort(species_to_score_v.begin(),species_to_score_v.end(),Pair2ndComparator<uint32_t>());

	cout << "Name\tTaxonID\tLength\tSummed Score\tNormalized Score\n";
	// Output the summed species scores
	for (vector<pair<uint32_t,uint32_t> >::iterator species_score = species_to_score_v.begin();
			species_score != species_to_score_v.end();
			++species_score) {
		uint32_t speciesID = species_score->first;
		pair<string,uint64_t> name_len = speciesID_to_name_len[speciesID];
		uint64_t slength = name_len.second;
		uint64_t sumscore = species_score->second;

		cout << name_len.first << "\t" <<
				speciesID << "\t" <<
				slength << "\t" <<
				sumscore << "\t" <<
				(float)sumscore/slength << "\n";
	}



}

//int centrifuge_report(int argc, const char **argv) {
int main(int argc, const char **argv) {

	if (argc < 3) {
		cerr << "Number of arguments is " << argc << endl;
		printUsage(cerr);
		exit(1);
	}

	argv0 = argv[0];
	const string bt2index = argv[1];
	const string cf_out = argv[2];
	//static string outfile;        // write SAM output to this file

	cout << "Input bt2 file: \"" << bt2index.c_str() << "\"" << endl;
	cout << "Centrifuge results file: \"" << cf_out.c_str() << "\"" << endl;

	driver<SString<char> >("DNA", bt2index, cf_out);
	return 0;
}

