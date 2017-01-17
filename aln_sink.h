/*
 * Copyright 2011, Ben Langmead <langmea@cs.jhu.edu>
 *
 * This file is part of Bowtie 2.
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

#ifndef ALN_SINK_H_
#define ALN_SINK_H_

#include <limits>
#include <utility>
#include <map>
#include "read.h"
#include "ds.h"
#include "simple_func.h"
#include "outq.h"
#include "aligner_result.h"
#include "hyperloglogplus.h"
#include "timer.h"
#include "taxonomy.h"


// Forward decl
template <typename index_t>
class SeedResults;

enum {
	OUTPUT_SAM = 1
};


struct ReadCounts {
	uint32_t n_unique_reads;
	uint32_t n_reads;
	uint32_t sum_score;
	double summed_hit_len;
	double weighted_reads;

	ReadCounts& operator+=(const ReadCounts& b) {
			n_unique_reads += b.n_unique_reads;
			n_reads += b.n_reads;
			sum_score += b.sum_score;
			summed_hit_len += b.summed_hit_len;
			weighted_reads += b.weighted_reads;
			return *this;
		}
};

ReadCounts operator+(ReadCounts a, const ReadCounts& b) {
	return a += b;
}

/**
 * Metrics summarizing the species level information we have
 */
struct SpeciesMetrics {
    
    //
    struct IDs {
        EList<uint64_t, 5> ids;
        bool operator<(const IDs& o) const {
            if(ids.size() != o.ids.size()) return ids.size() < o.ids.size();
            for(size_t i = 0; i < ids.size(); i++) {
                assert_lt(i, o.ids.size());
                if(ids[i] != o.ids[i]) return ids[i] < o.ids[i];
            }
            return false;
        }
        
        IDs& operator=(const IDs& other) {
            if(this == &other)
                return *this;
            
            ids = other.ids;
            return *this;
        }
    };


	SpeciesMetrics():mutex_m() {
	    reset();
	}

	void reset() {
		species_counts.clear();
		//for(map<uint32_t, HyperLogLogPlusMinus<uint64_t> >::iterator it = this->species_kmers.begin(); it != this->species_kmers.end(); ++it) {
		//	it->second.reset();
		//} //TODO: is this required?
		species_kmers.clear();
        num_non_leaves = 0;
	}

	void init(
              const map<uint64_t, ReadCounts>& species_counts_,
              const map<uint64_t, HyperLogLogPlusMinus<uint64_t> >& species_kmers_,
              const map<IDs, uint64_t>& observed_)
	{
		species_counts = species_counts_;
		species_kmers = species_kmers_;
        observed = observed_;
        num_non_leaves = 0;
    }

	/**
	 * Merge (add) the counters in the given ReportingMetrics object
	 * into this object.  This is the only safe way to update a
	 * ReportingMetrics shared by multiple threads.
	 */
	void merge(const SpeciesMetrics& met, bool getLock = false) {
        ThreadSafe ts(&mutex_m, getLock);

        // update species read count
        for(map<uint64_t, ReadCounts>::const_iterator it = met.species_counts.begin(); it != met.species_counts.end(); ++it) {
        	if (species_counts.find(it->first) == species_counts.end()) {
        		species_counts[it->first] = it->second;
        	} else {
        		species_counts[it->first].n_reads += it->second.n_reads;
        		species_counts[it->first].sum_score += it->second.sum_score;
        		species_counts[it->first].summed_hit_len += it->second.summed_hit_len;
        		species_counts[it->first].weighted_reads += it->second.weighted_reads;
        		species_counts[it->first].n_unique_reads += it->second.n_unique_reads;
        	}
        }

        // update species k-mers
        for(map<uint64_t, HyperLogLogPlusMinus<uint64_t> >::const_iterator it = met.species_kmers.begin(); it != met.species_kmers.end(); ++it) {
        	species_kmers[it->first].merge(&(it->second));
        }

        for(map<IDs, uint64_t>::const_iterator itr = met.observed.begin(); itr != met.observed.end(); itr++) {
            const IDs& ids = itr->first;
            uint64_t count = itr->second;
            
            if(observed.find(ids) == observed.end()) {
                observed[ids] = count;
            } else {
                observed[ids] += count;
            }
        }
    }

	void addSpeciesCounts(
                          uint64_t taxID,
                          int64_t score,
                          int64_t max_score,
                          double summed_hit_len,
                          double weighted_read,
                          uint32_t nresult)
    {
		species_counts[taxID].n_reads += 1;
		species_counts[taxID].sum_score += 1;
		species_counts[taxID].weighted_reads += weighted_read;
		species_counts[taxID].summed_hit_len += summed_hit_len;
		if(nresult == 1) {
			species_counts[taxID].n_unique_reads += 1;
		}

        // Only consider good hits for abundance analysis
        // DK - for debugging purposes
        if(score >= max_score) {
            cur_ids.ids.push_back(taxID);
            if(cur_ids.ids.size() == nresult) {
                cur_ids.ids.sort();
                if(observed.find(cur_ids) == observed.end()) {
                    observed[cur_ids] = 1;
                } else {
                    observed[cur_ids] += 1;
                }
                cur_ids.ids.clear();
            }
        }
	}

	void addAllKmers(
                     uint64_t taxID,
                     const BTDnaString &btdna,
                     size_t begin,
                     size_t len) {
#ifdef FLORIAN_DEBUG
		cerr << "add all kmers for " << taxID << " from " << begin << " for " << len << ": " << string(btdna.toZBuf()).substr(begin,len) << endl;
#endif
		uint64_t kmer = btdna.int_kmer<uint64_t>(begin,begin+len);
		species_kmers[taxID].add(kmer);
		size_t i = begin;
		while (i+32 < len) {
			kmer = btdna.next_kmer(kmer,i);
			species_kmers[taxID].add(kmer);
			++i;
		}
	}

	size_t nDistinctKmers(uint64_t taxID) {
		return(species_kmers[taxID].cardinality());
	}
    
    static void EM(
                   const map<IDs, uint64_t>& observed,
                   const map<uint64_t, EList<uint64_t> >& ancestors,
                   const map<uint64_t, uint64_t>& tid_to_num,
                   const EList<double>& p,
                   EList<double>& p_next,
                   const EList<size_t>& len)
    {
        assert_eq(p.size(), len.size());
        
        // E step
        p_next.fill(0.0);
        // for each assigned read set
        for(map<IDs, uint64_t>::const_iterator itr = observed.begin(); itr != observed.end(); itr++) {
            const EList<uint64_t, 5>& ids = itr->first.ids; // all ids assigned to the read set
            uint64_t count = itr->second; // number of reads in the read set
            double psum = 0.0;
            for(size_t i = 0; i < ids.size(); i++) {
                uint64_t tid = ids[i];
                // Leaves?
                map<uint64_t, uint64_t>::const_iterator id_itr = tid_to_num.find(tid);
                if(id_itr != tid_to_num.end()) {
                    uint64_t num = id_itr->second;
                    assert_lt(num, p.size());
                    psum += p[num];
                } else { // Ancestors
                    map<uint64_t, EList<uint64_t> >::const_iterator a_itr = ancestors.find(tid);
                    if(a_itr == ancestors.end())
                        continue;
                    const EList<uint64_t>& children = a_itr->second;
                    for(size_t c = 0; c < children.size(); c++) {
                        uint64_t c_tid = children[c];
                        map<uint64_t, uint64_t>::const_iterator id_itr = tid_to_num.find(c_tid);
                        if(id_itr == tid_to_num.end())
                            continue;
                        uint64_t c_num = id_itr->second;
                        psum += p[c_num];
                    }
                }
            }

            if(psum == 0.0) continue;
            
            for(size_t i = 0; i < ids.size(); i++) {
                uint64_t tid = ids[i];
                // Leaves?
                map<uint64_t, uint64_t>::const_iterator id_itr = tid_to_num.find(tid);
                if(id_itr != tid_to_num.end()) {
                    uint64_t num = id_itr->second;
                    assert_leq(p[num], psum);
                    p_next[num] += (count * (p[num] / psum));
                } else {
                    map<uint64_t, EList<uint64_t> >::const_iterator a_itr = ancestors.find(tid);
                    if(a_itr == ancestors.end())
                        continue;
                    const EList<uint64_t>& children = a_itr->second;
                    for(size_t c = 0; c < children.size(); c++) {
                        uint64_t c_tid = children[c];
                        map<uint64_t, uint64_t>::const_iterator id_itr = tid_to_num.find(c_tid);
                        if(id_itr == tid_to_num.end())
                            continue;
                        uint64_t c_num = id_itr->second;
                        p_next[c_num] += (count * (p[c_num] / psum));
                    }
                }
            }
        }
        
        // M step
        double sum = 0.0;
        for(size_t i = 0; i < p_next.size(); i++) {
            sum += (p_next[i] / len[i]);
        }
        for(size_t i = 0; i < p_next.size(); i++) {
            p_next[i] = p_next[i] / len[i] / sum;
        }
    }
    
    void calculateAbundance(const Ebwt<uint64_t>& ebwt, uint8_t rank)
    {
        const map<uint64_t, TaxonomyNode>& tree = ebwt.tree();
        
        // Find leaves
        set<uint64_t> leaves;
        for(map<IDs, uint64_t>::iterator itr = observed.begin(); itr != observed.end(); itr++) {
            const IDs& ids = itr->first;
            for(size_t i = 0; i < ids.ids.size(); i++) {
                uint64_t tid = ids.ids[i];
                map<uint64_t, TaxonomyNode>::const_iterator tree_itr = tree.find(tid);
                if(tree_itr == tree.end())
                    continue;
                const TaxonomyNode& node = tree_itr->second;
                if(!node.leaf) {
                    //if(tax_rank_num[node.rank] > tax_rank_num[rank]) {
                        continue;
                    //}
                }
                leaves.insert(tree_itr->first);
            }
        }
        
 
#ifdef DAEHWAN_DEBUG
        cerr << "\t\tnumber of leaves: " << leaves.size() << endl;
#endif
        
        // Find all descendants coming from the same ancestor
        map<uint64_t, EList<uint64_t> > ancestors;
        for(map<IDs, uint64_t>::iterator itr = observed.begin(); itr != observed.end(); itr++) {
            const IDs& ids = itr->first;
            for(size_t i = 0; i < ids.ids.size(); i++) {
                uint64_t tid = ids.ids[i];
                if(leaves.find(tid) != leaves.end())
                    continue;
                if(ancestors.find(tid) != ancestors.end())
                    continue;
                ancestors[tid].clear();
                for(set<uint64_t> ::const_iterator leaf_itr = leaves.begin(); leaf_itr != leaves.end(); leaf_itr++) {
                    uint64_t tid2 = *leaf_itr;
                    assert(tree.find(tid2) != tree.end());
                    assert(tree.find(tid2)->second.leaf);
                    uint64_t temp_tid2 = tid2;
                    while(true) {
                        map<uint64_t, TaxonomyNode>::const_iterator tree_itr = tree.find(temp_tid2);
                        if(tree_itr == tree.end())
                            break;
                        const TaxonomyNode& node = tree_itr->second;
                        if(tid == node.parent_tid) {
                            ancestors[tid].push_back(tid2);
                        }
                        if(temp_tid2 == node.parent_tid)
                            break;
                        temp_tid2 = node.parent_tid;
                    }
                }
                ancestors[tid].sort();
            }
        }
        
#ifdef DAEHWAN_DEBUG
        cerr << "\t\tnumber of ancestors: " << ancestors.size() << endl;
        for(map<uint64_t, EList<uint64_t> >::const_iterator itr = ancestors.begin(); itr != ancestors.end(); itr++) {
            uint64_t tid = itr->first;
            const EList<uint64_t>& children = itr->second;
            if(children.size() <= 0)
                continue;
            map<uint64_t, TaxonomyNode>::const_iterator tree_itr = tree.find(tid);
            if(tree_itr == tree.end())
                continue;
            const TaxonomyNode& node = tree_itr->second;
            cerr << "\t\t\t" << tid << ": " << children.size() << "\t" << get_tax_rank(node.rank) << endl;
            cerr << "\t\t\t\t";
            for(size_t i = 0; i < children.size(); i++) {
                cerr << children[i];
                if(i + 1 < children.size())
                    cerr << ",";
                if(i > 10) {
                    cerr << " ...";
                    break;
                }
            }
            cerr << endl;
        }
        
        uint64_t test_tid = 0, test_tid2 = 0;
#endif
        // Lengths of genomes (or contigs)
        const map<uint64_t, uint64_t>& size_table = ebwt.size();
        
        // Initialize probabilities
        map<uint64_t, uint64_t> tid_to_num; // taxonomic ID to corresponding element of a list
        EList<double> p;
        EList<size_t> len; // genome lengths
        for(map<IDs, uint64_t>::iterator itr = observed.begin(); itr != observed.end(); itr++) {
            const IDs& ids = itr->first;
            uint64_t count = itr->second;
            for(size_t i = 0; i < ids.ids.size(); i++) {
                uint64_t tid = ids.ids[i];
                if(leaves.find(tid) == leaves.end())
                    continue;
                
#ifdef DAEHWAN_DEBUG
                if((tid == test_tid || tid == test_tid2) &&
                   count >= 100) {
                    cerr << tid << ": " << count << "\t";
                    for(size_t j = 0; j < ids.ids.size(); j++) {
                        cerr << ids.ids[j];
                        if(j + 1 < ids.ids.size())
                            cerr << ",";
                    }
                    cerr << endl;
                }
#endif
                
                if(tid_to_num.find(tid) == tid_to_num.end()) {
                    tid_to_num[tid] = p.size();
                    p.push_back(1.0 / ids.ids.size() * count);
                    map<uint64_t, uint64_t>::const_iterator size_itr = size_table.find(tid);
                    if(size_itr != size_table.end()) {
                        len.push_back(size_itr->second);
                    } else {
                        len.push_back(std::numeric_limits<size_t>::max());
                    }
                } else {
                    uint64_t num = tid_to_num[tid];
                    assert_lt(num, p.size());
                    p[num] += (1.0 / ids.ids.size() * count);
                }
            }
        }
        
        assert_eq(p.size(), len.size());
        
        {
            double sum = 0.0;
            for(size_t i = 0; i < p.size(); i++) {
                sum += (p[i] / len[i]);
            }
            for(size_t i = 0; i < p.size(); i++) {
                p[i] = (p[i] / len[i]) / sum;
            }
        }
        
        EList<double> p_next; p_next.resizeExact(p.size());
        EList<double> p_next2; p_next2.resizeExact(p.size());
        EList<double> p_r; p_r.resizeExact(p.size());
        EList<double> p_v; p_v.resizeExact(p.size());
        size_t num_iteration = 0;
        double diff = 0.0;
        while(true) {
#ifdef DAEHWAN_DEBUG
            if(num_iteration % 50 == 0) {
                if(test_tid != 0 || test_tid2 != 0)
                    cerr << "iter " << num_iteration << endl;
                if(test_tid != 0)
                    cerr << "\t" << test_tid << ": " << p[tid_to_num[test_tid]] << endl;
                if(test_tid2 != 0)
                    cerr << "\t" << test_tid2 << ": " << p[tid_to_num[test_tid2]] << endl;
            }
#endif
            
            // Accelerated version of EM - SQUAREM iteration
            //    Varadhan, R. & Roland, C. Scand. J. Stat. 35, 335–353 (2008).
            //    Also, this algorithm is used in Sailfish - http://www.nature.com/nbt/journal/v32/n5/full/nbt.2862.html
#if 1
            EM(observed, ancestors, tid_to_num, p, p_next, len);
            EM(observed, ancestors, tid_to_num, p_next, p_next2, len);
            double sum_squared_r = 0.0, sum_squared_v = 0.0;
            for(size_t i = 0; i < p.size(); i++) {
                p_r[i] = p_next[i] - p[i];
                sum_squared_r += (p_r[i] * p_r[i]);
                p_v[i] = p_next2[i] - p_next[i] - p_r[i];
                sum_squared_v += (p_v[i] * p_v[i]);
            }
            if(sum_squared_v > 0.0) {
                double gamma = -sqrt(sum_squared_r / sum_squared_v);
                for(size_t i = 0; i < p.size(); i++) {
                    p_next2[i] = max(0.0, p[i] - 2 * gamma * p_r[i] + gamma * gamma * p_v[i]);
                }
                EM(observed, ancestors, tid_to_num, p_next2, p_next, len);
            }
            
#else
            EM(observed, ancestors, tid_to_num, p, p_next, len);
#endif
            
            diff = 0.0;
            for(size_t i = 0; i < p.size(); i++) {
                diff += (p[i] > p_next[i] ? p[i] - p_next[i] : p_next[i] - p[i]);
            }
            if(diff < 0.0000000001) break;
            if(++num_iteration >= 10000) break;
            p = p_next;
        }
        
        cerr << "Number of iterations in EM algorithm: " << num_iteration << endl;
        cerr << "Probability diff. (P - P_prev) in the last iteration: " << diff << endl;
        
        {
            // Calculate abundance normalized by genome size
            abundance_len.clear();
            double sum = 0.0;
            for(map<uint64_t, uint64_t>::iterator itr = tid_to_num.begin(); itr != tid_to_num.end(); itr++) {
                uint64_t tid = itr->first;
                uint64_t num = itr->second;
                assert_lt(num, p.size());
                abundance_len[tid] = p[num];
                sum += (p[num] * len[num]);
            }
            
            // Calculate abundance without genome size taken into account
            abundance.clear();
            for(map<uint64_t, uint64_t>::iterator itr = tid_to_num.begin(); itr != tid_to_num.end(); itr++) {
                uint64_t tid = itr->first;
                uint64_t num = itr->second;
                assert_lt(num, p.size());
                abundance[tid] = (p[num] * len[num]) / sum;
            }
        }
    }

	map<uint64_t, ReadCounts> species_counts;                        // read count per species
	map<uint64_t, HyperLogLogPlusMinus<uint64_t> > species_kmers;    // unique k-mer count per species
    
    map<IDs, uint64_t>     observed;
    IDs                    cur_ids;
    uint32_t               num_non_leaves;
    map<uint64_t, double>  abundance;      // abundance without genome size taken into consideration
    map<uint64_t, double>  abundance_len;  // abundance normalized by genome size

	MUTEX_T mutex_m;
};


/**
 * Metrics summarizing the work done by the reporter and summarizing
 * the number of reads that align, that fail to align, and that align
 * non-uniquely.
 */
struct ReportingMetrics {

	ReportingMetrics():mutex_m() {
	    reset();
	}

	void reset() {
		init(0, 0, 0, 0);
	}

	void init(
		uint64_t nread_,
		uint64_t npaired_,
		uint64_t nunpaired_,
		uint64_t nconcord_uni_)
	{
		nread         = nread_;
		
		npaired       = npaired_;
		nunpaired     = nunpaired_;
		
		nconcord_uni  = nconcord_uni_;
    }
	
	/**
	 * Merge (add) the counters in the given ReportingMetrics object
	 * into this object.  This is the only safe way to update a
	 * ReportingMetrics shared by multiple threads.
	 */
	void merge(const ReportingMetrics& met, bool getLock = false) {
        ThreadSafe ts(&mutex_m, getLock);
		nread         += met.nread;

		npaired       += met.npaired;
		nunpaired     += met.nunpaired;

		nconcord_uni  += met.nconcord_uni;
    }

	uint64_t  nread;         // # reads
	uint64_t  npaired;       // # pairs
	uint64_t  nunpaired;     // # unpaired reads
	
	// Paired
	
	// Concordant
	uint64_t  nconcord_uni;  // # pairs with unique concordant alns
		
	MUTEX_T mutex_m;
};

// Type for expression numbers of hits
typedef int64_t THitInt;

/**
 * Parameters affecting reporting of alignments, specifically -k & -a,
 * -m & -M.
 */
struct ReportingParams {

	explicit ReportingParams(THitInt khits_, bool compressed_)
	{
		init(khits_, compressed_);
	}

	void init(THitInt khits_, bool compressed_)
	{
		khits = khits_;     // -k (or high if -a)
        if(compressed_) {
            ihits = max<THitInt>(khits, 5) * 4;
        } else {
            ihits = max<THitInt>(khits, 5) * 40;
        }
	}
	
#ifndef NDEBUG
	/**
	 * Check that reporting parameters are internally consistent.
	 */
	bool repOk() const {
		assert_geq(khits, 1);
		return true;
	}
#endif
	
	inline THitInt mult() const {
		return khits;
	}

	// Number of assignments to report
	THitInt khits;
    
    // Number of internal assignments
    THitInt ihits;
};

/**
 * A state machine keeping track of the number and type of alignments found so
 * far.  Its purpose is to inform the caller as to what stage the alignment is
 * in and what categories of alignment are still of interest.  This information
 * should allow the caller to short-circuit some alignment work.  Another
 * purpose is to tell the AlnSinkWrap how many and what type of alignment to
 * report.
 *
 * TODO: This class does not keep accurate information about what
 * short-circuiting took place.  If a read is identical to a previous read,
 * there should be a way to query this object to determine what work, if any,
 * has to be re-done for the new read.
 */
class ReportingState {

public:

	enum {
		NO_READ = 1,        // haven't got a read yet
		CONCORDANT_PAIRS,   // looking for concordant pairs
		DONE                // finished looking
	};

	// Flags for different ways we can finish out a category of potential
	// alignments.
	
	enum {
		EXIT_DID_NOT_EXIT = 1,        // haven't finished
		EXIT_DID_NOT_ENTER,           // never tried search	
		EXIT_SHORT_CIRCUIT_k,         // -k exceeded
		EXIT_NO_ALIGNMENTS,           // none found
		EXIT_WITH_ALIGNMENTS          // some found
	};
	
	ReportingState(const ReportingParams& p) : p_(p) { reset(); }
	
	/**
	 * Set all state to uninitialized defaults.
	 */
	void reset() {
		state_ = ReportingState::NO_READ;
		paired_ = false;
		nconcord_ = 0;
		doneConcord_ = false;
		exitConcord_ = ReportingState::EXIT_DID_NOT_ENTER;
		done_ = false;
	}
	
	/**
	 * Return true iff this ReportingState has been initialized with a call to
	 * nextRead() since the last time reset() was called.
	 */
	bool inited() const { return state_ != ReportingState::NO_READ; }

	/**
	 * Initialize state machine with a new read.  The state we start in depends
	 * on whether it's paired-end or unpaired.
	 */
	void nextRead(bool paired);

	/**
	 * Caller uses this member function to indicate that one additional
	 * concordant alignment has been found.
	 */
	bool foundConcordant();

	/**
	 * Caller uses this member function to indicate that one additional
	 * discordant alignment has been found.
	 */
	bool foundUnpaired(bool mate1);
	
	/**
	 * Called to indicate that the aligner has finished searching for
	 * alignments.  This gives us a chance to finalize our state.
	 *
	 * TODO: Keep track of short-circuiting information.
	 */
	void finish();
	
	/**
	 * Populate given counters with the number of various kinds of alignments
	 * to report for this read.  Concordant alignments are preferable to (and
	 * mutually exclusive with) discordant alignments, and paired-end
	 * alignments are preferable to unpaired alignments.
	 *
	 * The caller also needs some additional information for the case where a
	 * pair or unpaired read aligns repetitively.  If the read is paired-end
	 * and the paired-end has repetitive concordant alignments, that should be
	 * reported, and 'pairMax' is set to true to indicate this.  If the read is
	 * paired-end, does not have any conordant alignments, but does have
	 * repetitive alignments for one or both mates, then that should be
	 * reported, and 'unpair1Max' and 'unpair2Max' are set accordingly.
	 *
	 * Note that it's possible in the case of a paired-end read for the read to
	 * have repetitive concordant alignments, but for one mate to have a unique
	 * unpaired alignment.
	 */
	void getReport(uint64_t& nconcordAln) const; // # concordant alignments to report

	/**
	 * Return an integer representing the alignment state we're in.
	 */
	inline int state() const { return state_; }
	
	/**
	 * If false, there's no need to solve any more dynamic programming problems
	 * for finding opposite mates.
	 */
	inline bool doneConcordant() const { return doneConcord_; }
	
	/**
	 * Return true iff all alignment stages have been exited.
	 */
	inline bool done() const { return done_; }

	inline uint64_t numConcordant() const { return nconcord_; }

	inline int exitConcordant() const { return exitConcord_; }

	/**
	 * Return ReportingParams object governing this ReportingState.
	 */
	const ReportingParams& params() const {
		return p_;
	}

protected:
	const ReportingParams& p_;  // reporting parameters
	int state_;          // state we're currently in
	bool paired_;        // true iff read we're currently handling is paired
	uint64_t nconcord_;  // # concordants found so far
	bool doneConcord_;   // true iff we're no longner interested in concordants
	int exitConcord_;    // flag indicating how we exited concordant state
	bool done_;          // done with all alignments
};

/**
 * Global hit sink for hits from the MultiSeed aligner.  Encapsulates
 * all aspects of the MultiSeed aligner hitsink that are global to all
 * threads.  This includes aspects relating to:
 *
 * (a) synchronized access to the output stream
 * (b) the policy to be enforced by the per-thread wrapper
 *
 * TODO: Implement splitting up of alignments into separate files
 * according to genomic coordinate.
 */
template <typename index_t>
class AlnSink {

	typedef EList<std::string> StrList;

public:

	explicit AlnSink(
                     OutputQueue& oq,
                     const StrList& refnames,
                     const EList<uint32_t>& tab_fmt_cols,
                     bool quiet) :
    oq_(oq),
    refnames_(refnames),
    tab_fmt_cols_(tab_fmt_cols),
    quiet_(quiet)
	{
	}

	/**
	 * Destroy HitSinkobject;
	 */
	virtual ~AlnSink() { }

	/**
	 * Called when the AlnSink is wrapped by a new AlnSinkWrap.  This helps us
	 * keep track of whether the main lock or any of the per-stream locks will
	 * be contended by multiple threads.
	 */
	void addWrapper() { numWrappers_++; }

	/**
	 * Append a single hit to the given output stream.  If
	 * synchronization is required, append() assumes the caller has
	 * already grabbed the appropriate lock.
	 */
	virtual void append(
		BTString&             o,
		size_t                threadId,
		const Read           *rd1,
		const Read           *rd2,
		const TReadId         rdid,
		AlnRes               *rs1,
		AlnRes               *rs2,
		const AlnSetSumm&     summ,
		const PerReadMetrics& prm,
		SpeciesMetrics& sm,
		bool report2,
		size_t n_results) = 0;

	/**
	 * Report a given batch of hits for the given read or read pair.
	 * Should be called just once per read pair.  Assumes all the
	 * alignments are paired, split between rs1 and rs2.
	 *
	 * The caller hasn't decided which alignments get reported as primary
	 * or secondary; that's up to the routine.  Because the caller might
	 * want to know this, we use the pri1 and pri2 out arguments to
	 * convey this.
	 */
	virtual void reportHits(
		BTString&             o,              // write to this buffer
		size_t                threadId,       // which thread am I?
		const Read           *rd1,            // mate #1
		const Read           *rd2,            // mate #2
		const TReadId         rdid,           // read ID
		const EList<size_t>&  select1,        // random subset of rd1s
		const EList<size_t>*  select2,        // random subset of rd2s
		EList<AlnRes>        *rs1,            // alignments for mate #1
		EList<AlnRes>        *rs2,            // alignments for mate #2
		bool                  maxed,          // true iff -m/-M exceeded
		const AlnSetSumm&     summ,           // summary
		const PerReadMetrics& prm,            // per-read metrics
		SpeciesMetrics& sm,             // species metrics
		bool                  getLock = true) // true iff lock held by caller
	{
		assert(rd1 != NULL || rd2 != NULL);
		assert(rs1 != NULL || rs2 != NULL);

        for(size_t i = 0; i < select1.size(); i++) {
            AlnRes* r1 = ((rs1 != NULL) ? &rs1->get(select1[i]) : NULL);
            AlnRes* r2 = ((rs2 != NULL) ? &rs2->get(select1[i]) : NULL);
            append(o, threadId, rd1, rd2, rdid, r1, r2, summ, prm, sm, true, select1.size());
        }
	}

	/**
	 * Report an unaligned read.  Typically we do nothing, but we might
	 * want to print a placeholder when output is chained.
	 */
	virtual void reportUnaligned(
		BTString&             o,              // write to this string
		size_t                threadId,       // which thread am I?
		const Read           *rd1,            // mate #1
		const Read           *rd2,            // mate #2
		const TReadId         rdid,           // read ID
		const AlnSetSumm&     summ,           // summary
		const PerReadMetrics& prm,            // per-read metrics
		bool                  report2,        // report alns for both mates?
		bool                  getLock = true) // true iff lock held by caller
	{
		// FIXME: reportUnaligned does nothing
		//append(o, threadId, rd1, rd2, rdid, NULL, NULL, summ, prm, NULL,report2);
	}

	/**
	 * Print summary of how many reads aligned, failed to align and aligned
	 * repetitively.  Write it to stderr.  Optionally write Hadoop counter
	 * updates.
	 */
	void printAlSumm(
		const ReportingMetrics& met,
		size_t repThresh, // threshold for uniqueness, or max if no thresh
		bool discord,     // looked for discordant alignments
		bool mixed,       // looked for unpaired alignments where paired failed?
		bool hadoopOut);  // output Hadoop counters?

	/**
	 * Called when all alignments are complete.  It is assumed that no
	 * synchronization is necessary.
	 */
	void finish(
		size_t repThresh,
		bool discord,
		bool mixed,
		bool hadoopOut)
	{
		// Close output streams
		if(!quiet_) {
			printAlSumm(
				met_,
				repThresh,
				discord,
				mixed,
				hadoopOut);
		}
	}

#ifndef NDEBUG
	/**
	 * Check that hit sink is internally consistent.
	 */
	bool repOk() const { return true; }
#endif
	
	//
	// Related to reporting seed hits
	//

	/**
	 * Given a Read and associated, filled-in SeedResults objects,
	 * print a record summarizing the seed hits.
	 */
	void reportSeedSummary(
		BTString&          o,
		const Read&        rd,
		TReadId            rdid,
		size_t             threadId,
		const SeedResults<index_t>& rs,
		bool               getLock = true);

	/**
	 * Given a Read, print an empty record (all 0s).
	 */
	void reportEmptySeedSummary(
		BTString&          o,
		const Read&        rd,
		TReadId            rdid,
		size_t             threadId,
		bool               getLock = true);

	/**
	 * Append a batch of unresolved seed alignment results (i.e. seed
	 * alignments where all we know is the reference sequence aligned
	 * to and its SA range, not where it falls in the reference
	 * sequence) to the given output stream in Bowtie's seed-alignment
	 * verbose-mode format.
	 */
	virtual void appendSeedSummary(
		BTString&     o,
		const Read&   rd,
		const TReadId rdid,
		size_t        seedsTried,
		size_t        nonzero,
		size_t        ranges,
		size_t        elts,
		size_t        seedsTriedFw,
		size_t        nonzeroFw,
		size_t        rangesFw,
		size_t        eltsFw,
		size_t        seedsTriedRc,
		size_t        nonzeroRc,
		size_t        rangesRc,
		size_t        eltsRc);

	/**
	 * Merge given metrics in with ours by summing all individual metrics.
	 */
	void mergeMetrics(const ReportingMetrics& met, bool getLock = true) {
		met_.merge(met, getLock);
	}

	/**
	 * Return mutable reference to the shared OutputQueue.
	 */
	OutputQueue& outq() {
		return oq_;
	}

protected:
    
	OutputQueue&       oq_;           // output queue
	int                numWrappers_;  // # threads owning a wrapper for this HitSink
	const StrList&     refnames_;     // reference names
	const EList<uint32_t>&     tab_fmt_cols_; // Columns that are printed in tabular format
	bool               quiet_;        // true -> don't print alignment stats at the end
	ReportingMetrics   met_;          // global repository of reporting metrics
};

/**
 * Per-thread hit sink "wrapper" for the MultiSeed aligner.  Encapsulates
 * aspects of the MultiSeed aligner hit sink that are per-thread.  This
 * includes aspects relating to:
 *
 * (a) Enforcement of the reporting policy
 * (b) Tallying of results
 * (c) Storing of results for the previous read in case this allows us to
 *     short-circuit some work for the next read (i.e. if it's identical)
 *
 * PHASED ALIGNMENT ASSUMPTION
 *
 * We make some assumptions about how alignment proceeds when we try to
 * short-circuit work for identical reads.  Specifically, we assume that for
 * each read the aligner proceeds in a series of stages (or perhaps just one
 * stage).  In each stage, the aligner either:
 *
 * (a)  Finds no alignments, or
 * (b)  Finds some alignments and short circuits out of the stage with some
 *      random reporting involved (e.g. in -k and/or -M modes), or
 * (c)  Finds all of the alignments in the stage
 *
 * In the event of (a), the aligner proceeds to the next stage and keeps
 * trying; we can skip the stage entirely for the next read if it's identical.
 * In the event of (b), or (c), the aligner stops and does not proceed to
 * further stages.  In the event of (b1), if the next read is identical we
 * would like to tell the aligner to start again at the beginning of the stage
 * that was short-circuited.
 *
 * In any event, the rs1_/rs2_/rs1u_/rs2u_ fields contain the alignments found
 * in the last alignment stage attempted.
 *
 * HANDLING REPORTING LIMITS
 *
 * The user can specify reporting limits, like -k (specifies number of
 * alignments to report out of those found) and -M (specifies a ceiling s.t. if
 * there are more alignments than the ceiling, read is called repetitive and
 * best found is reported).  Enforcing these limits is straightforward for
 * unpaired alignments: if a new alignment causes us to exceed the -M ceiling,
 * we can stop looking.
 *
 * The case where both paired-end and unpaired alignments are possible is
 * trickier.  Once we have a number of unpaired alignments that exceeds the
 * ceiling, we can stop looking *for unpaired alignments* - but we can't
 * necessarily stop looking for paired-end alignments, since there may yet be
 * more to find.  However, if the input read is not a pair, then we can stop at
 * this point.  If the input read is a pair and we have a number of paired
 * aligments that exceeds the -M ceiling, we can stop looking.
 *
 * CONCORDANT & DISCORDANT, PAIRED & UNPAIRED
 *
 * A note on paired-end alignment: Clearly, if an input read is
 * paired-end and we find either concordant or discordant paired-end
 * alignments for the read, then we would like to tally and report
 * those alignments as such (and not as groups of 2 unpaired
 * alignments).  And if we fail to find any paired-end alignments, but
 * we do find some unpaired alignments for one mate or the other, then
 * we should clearly tally and report those alignments as unpaired
 * alignments (if the user so desires).
 *
 * The situation is murkier when there are no paired-end alignments,
 * but there are unpaired alignments for *both* mates.  In this case,
 * we might want to pick out zero or more pairs of mates and classify
 * those pairs as discordant paired-end alignments.  And we might want
 * to classify the remaining alignments as unpaired.  But how do we
 * pick which pairs if any to call discordant?
 *
 * Because the most obvious use for discordant pairs is for identifying
 * large-scale variation, like rearrangements or large indels, we would
 * usually like to be conservative about what we call a discordant
 * alignment.  If there's a good chance that one or the other of the
 * two mates has a good alignment to another place on the genome, this
 * compromises the evidence for the large-scale variant.  For this
 * reason, Bowtie 2's policy is: if there are no paired-end alignments
 * and there is *exactly one alignment each* for both mates, then the
 * two alignments are paired and treated as a discordant paired-end
 * alignment.  Otherwise, all alignments are treated as unpaired
 * alignments.
 *
 * When both paired and unpaired alignments are discovered by the
 * aligner, only the paired alignments are reported by default.  This
 * is sensible considering relative likelihoods: if a good paired-end
 * alignment is found, it is much more likely that the placement of
 * the two mates implied by that paired alignment is correct than any
 * placement implied by an unpaired alignment.
 *
 * 
 */
template <typename index_t>
class AlnSinkWrap {
public:

	AlnSinkWrap(
		AlnSink<index_t>& g,       // AlnSink being wrapped
		const ReportingParams& rp, // Parameters governing reporting
		size_t threadId,           // Thread ID
        bool secondary = false) :  // Secondary alignments
		g_(g),
		rp_(rp),
        threadid_(threadId),
    	secondary_(secondary),
		init_(false),   
		maxed1_(false),       // read is pair and we maxed out mate 1 unp alns
		maxed2_(false),       // read is pair and we maxed out mate 2 unp alns
		maxedOverall_(false), // alignments found so far exceed -m/-M ceiling
		bestPair_(std::numeric_limits<TAlScore>::min()),
		best2Pair_(std::numeric_limits<TAlScore>::min()),
		bestUnp1_(std::numeric_limits<TAlScore>::min()),
		best2Unp1_(std::numeric_limits<TAlScore>::min()),
		bestUnp2_(std::numeric_limits<TAlScore>::min()),
		best2Unp2_(std::numeric_limits<TAlScore>::min()),
        bestSplicedPair_(0),
        best2SplicedPair_(0),
        bestSplicedUnp1_(0),
        best2SplicedUnp1_(0),
        bestSplicedUnp2_(0),
        best2SplicedUnp2_(0),
		rd1_(NULL),    // mate 1
		rd2_(NULL),    // mate 2
		rdid_(std::numeric_limits<TReadId>::max()), // read id
		rs_(),         // mate 1 alignments for paired-end alignments
		select_(),     // for selecting random subsets for mate 1
		st_(rp)        // reporting state - what's left to do?
	{
		assert(rp_.repOk());
	}

	AlnSink<index_t>& getSink() {
		return(g_);
	}

	/**
	 * Initialize the wrapper with a new read pair and return an
	 * integer >= -1 indicating which stage the aligner should start
	 * at.  If -1 is returned, the aligner can skip the read entirely.
	 * at.  If .  Checks if the new read pair is identical to the
	 * previous pair.  If it is, then we return the id of the first
	 * stage to run.
	 */
	int nextRead(
		// One of the other of rd1, rd2 will = NULL if read is unpaired
		const Read* rd1,      // new mate #1
		const Read* rd2,      // new mate #2
		TReadId rdid,         // read ID for new pair
		bool qualitiesMatter);// aln policy distinguishes b/t quals?

	/**
	 * Inform global, shared AlnSink object that we're finished with
	 * this read.  The global AlnSink is responsible for updating
	 * counters, creating the output record, and delivering the record
	 * to the appropriate output stream.
	 */
	void finishRead(
		const SeedResults<index_t> *sr1, // seed alignment results for mate 1
		const SeedResults<index_t> *sr2, // seed alignment results for mate 2
		bool               exhaust1,     // mate 1 exhausted?
		bool               exhaust2,     // mate 2 exhausted?
		bool               nfilt1,       // mate 1 N-filtered?
		bool               nfilt2,       // mate 2 N-filtered?
		bool               scfilt1,      // mate 1 score-filtered?
		bool               scfilt2,      // mate 2 score-filtered?
		bool               lenfilt1,     // mate 1 length-filtered?
		bool               lenfilt2,     // mate 2 length-filtered?
		bool               qcfilt1,      // mate 1 qc-filtered?
		bool               qcfilt2,      // mate 2 qc-filtered?
		bool               sortByScore,  // prioritize alignments by score
		RandomSource&      rnd,          // pseudo-random generator
		ReportingMetrics&  met,          // reporting metrics
		SpeciesMetrics&    smet,         // species metrics
		const PerReadMetrics& prm,       // per-read metrics
		bool suppressSeedSummary = true,
		bool suppressAlignments = false);
	
	/**
	 * Called by the aligner when a new unpaired or paired alignment is
	 * discovered in the given stage.  This function checks whether the
	 * addition of this alignment causes the reporting policy to be
	 * violated (by meeting or exceeding the limits set by -k, -m, -M),
	 * in which case true is returned immediately and the aligner is
	 * short circuited.  Otherwise, the alignment is tallied and false
	 * is returned.
	 */
	bool report(
		int stage,
        const AlnRes* rs);

#ifndef NDEBUG
	/**
	 * Check that hit sink wrapper is internally consistent.
	 */
	bool repOk() const {
		if(init_) {
			assert(rd1_ != NULL);
			assert_neq(std::numeric_limits<TReadId>::max(), rdid_);
		}
		return true;
	}
#endif
	
	/**
	 * Return true iff no alignments have been reported to this wrapper
	 * since the last call to nextRead().
	 */
	bool empty() const {
		return rs_.empty();
	}
	
	/**
	 * Return true iff we have already encountered a number of alignments that
	 * exceeds the -m/-M ceiling.  TODO: how does this distinguish between
	 * pairs and mates?
	 */
	bool maxed() const {
		return maxedOverall_;
	}
	
	/**
	 * Return true if the current read is paired.
	 */
	bool readIsPair() const {
		return rd1_ != NULL && rd2_ != NULL;
	}
	
	/**
	 * Return true iff nextRead() has been called since the last time
	 * finishRead() was called.
	 */
	bool inited() const { return init_; }

	/**
	 * Return a const ref to the ReportingState object associated with the
	 * AlnSinkWrap.
	 */
	const ReportingState& state() const { return st_; }
    
    const ReportingParams& reportingParams() { return rp_;}
	
    SpeciesMetrics& speciesMetrics() { return g_.speciesMetrics(); }
	
	/**
	 * Return true iff at least two alignments have been reported so far for an
	 * unpaired read or mate 1.
	 */
	bool hasSecondBestUnp1() const {
		return best2Unp1_ != std::numeric_limits<TAlScore>::min();
	}

	/**
	 * Return true iff at least two alignments have been reported so far for
	 * mate 2.
	 */
	bool hasSecondBestUnp2() const {
		return best2Unp2_ != std::numeric_limits<TAlScore>::min();
	}

	/**
	 * Return true iff at least two paired-end alignments have been reported so
	 * far.
	 */
	bool hasSecondBestPair() const {
		return best2Pair_ != std::numeric_limits<TAlScore>::min();
	}
	
	/**
	 * Get best score observed so far for an unpaired read or mate 1.
	 */
	TAlScore bestUnp1() const {
		return bestUnp1_;
	}

	/**
	 * Get second-best score observed so far for an unpaired read or mate 1.
	 */
	TAlScore secondBestUnp1() const {
		return best2Unp1_;
	}

	/**
	 * Get best score observed so far for mate 2.
	 */
	TAlScore bestUnp2() const {
		return bestUnp2_;
	}

	/**
	 * Get second-best score observed so far for mate 2.
	 */
	TAlScore secondBestUnp2() const {
		return best2Unp2_;
	}

	/**
	 * Get best score observed so far for paired-end read.
	 */
	TAlScore bestPair() const {
		return bestPair_;
	}

	/**
	 * Get second-best score observed so far for paired-end read.
	 */
	TAlScore secondBestPair() const {
		return best2Pair_;
	}
    
    
    /**
     *
     */
    void getPair(const EList<AlnRes>*& rs) const { rs = &rs_; }

protected:

	/**
	 * Return true iff the read in rd1/rd2 matches the last read handled, which
	 * should still be in rd1_/rd2_.
	 */
	bool sameRead(
		const Read* rd1,
		const Read* rd2,
		bool qualitiesMatter);

	/**
	 * Given that rs is already populated with alignments, consider the
	 * alignment policy and make random selections where necessary.  E.g. if we
	 * found 10 alignments and the policy is -k 2 -m 20, select 2 alignments at
	 * random.  We "select" an alignment by setting the parallel entry in the
	 * 'select' list to true.
	 */
	size_t selectAlnsToReport(
		const EList<AlnRes>& rs,     // alignments to select from
		uint64_t             num,    // number of alignments to select
		EList<size_t>&       select, // list to put results in
		RandomSource&        rnd)
		const;

	/**
	 * rs1 (possibly together with rs2 if reads are paired) are populated with
	 * alignments.  Here we prioritize them according to alignment score, and
	 * some randomness to break ties.  Priorities are returned in the 'select'
	 * list.
	 */
	size_t selectByScore(
		const EList<AlnRes>* rs,    // alignments to select from (mate 1)
		uint64_t             num,    // number of alignments to select
		EList<size_t>&       select, // prioritized list to put results in
		RandomSource&        rnd)
		const;

	AlnSink<index_t>& g_;     // global alignment sink
	ReportingParams   rp_;    // reporting parameters: khits, mhits etc
	size_t            threadid_; // thread ID
    bool              secondary_; // allow for secondary alignments
	bool              init_;  // whether we're initialized w/ read pair
	bool              maxed1_; // true iff # unpaired mate-1 alns reported so far exceeded -m/-M
	bool              maxed2_; // true iff # unpaired mate-2 alns reported so far exceeded -m/-M
	bool              maxedOverall_; // true iff # paired-end alns reported so far exceeded -m/-M
	TAlScore          bestPair_;     // greatest score so far for paired-end
	TAlScore          best2Pair_;    // second-greatest score so far for paired-end
	TAlScore          bestUnp1_;     // greatest score so far for unpaired/mate1
	TAlScore          best2Unp1_;    // second-greatest score so far for unpaired/mate1
	TAlScore          bestUnp2_;     // greatest score so far for mate 2
	TAlScore          best2Unp2_;    // second-greatest score so far for mate 2
    index_t           bestSplicedPair_;
    index_t           best2SplicedPair_;
    index_t           bestSplicedUnp1_;
    index_t           best2SplicedUnp1_;
    index_t           bestSplicedUnp2_;
    index_t           best2SplicedUnp2_;
	const Read*       rd1_;   // mate #1
	const Read*       rd2_;   // mate #2
	TReadId           rdid_;  // read ID (potentially used for ordering)
	EList<AlnRes>     rs_;   // paired alignments for mate #1
	EList<size_t>     select_; // parallel to rs1_/rs2_ - which to report
	ReportingState    st_;      // reporting state - what's left to do?
	
	EList<std::pair<TAlScore, size_t> > selectBuf_;
	BTString obuf_;
};

/**
 * An AlnSink concrete subclass for printing SAM alignments.  The user might
 * want to customize SAM output in various ways.  We encapsulate all these
 * customizations, and some of the key printing routines, in the SamConfig
 * class in sam.h/sam.cpp.
 */
template <typename index_t>
class AlnSinkSam : public AlnSink<index_t> {

	typedef EList<std::string> StrList;

public:

	AlnSinkSam(
               Ebwt<index_t>*   ebwt,
               OutputQueue&     oq,           // output queue
               const StrList&   refnames,     // reference names
			   const EList<uint32_t>&   tab_fmt_cols, // columns to output in the tabular format
               bool             quiet) :
    AlnSink<index_t>(oq,
                     refnames,
					 tab_fmt_cols,
                     quiet),
    ebwt_(ebwt)
    { }
	
	virtual ~AlnSinkSam() { }

	/**
	 * Append a single alignment result, which might be paired or
	 * unpaired, to the given output stream in Bowtie's verbose-mode
	 * format.  If the alignment is paired-end, print mate1's alignment
	 * then mate2's alignment.
	 */
	virtual void append(
		BTString&     o,           // write output to this string
		size_t        threadId,    // which thread am I?
		const Read*   rd1,         // mate #1
		const Read*   rd2,         // mate #2
		const TReadId rdid,        // read ID
		AlnRes* rs1,               // alignments for mate #1
		AlnRes* rs2,               // alignments for mate #2
		const AlnSetSumm& summ,    // summary
		const PerReadMetrics& prm, // per-read metrics
		SpeciesMetrics& sm,  // species metrics
		bool report2,              // report alns for both mates
		size_t n_results)          // number of results for read
	{
		assert(rd1 != NULL || rd2 != NULL);
        appendMate(*ebwt_, o, *rd1, rd2, rdid, rs1, rs2, summ, prm, sm, n_results);
	}

protected:

	/**
	 * Append a single per-mate alignment result to the given output
	 * stream.  If the alignment is part of a pair, information about
	 * the opposite mate and its alignment are given in rdo/rso.
	 */
	void appendMate(
                    Ebwt<index_t>& ebwt,
                    BTString&     o,
                    const Read&   rd,
                    const Read*   rdo,
                    const TReadId rdid,
                    AlnRes* rs,
                    AlnRes* rso,
                    const AlnSetSumm& summ,
                    const PerReadMetrics& prm, // per-read metrics
                    SpeciesMetrics& sm,   // species metrics
                    size_t n_results);


    Ebwt<index_t>*   ebwt_;
	BTDnaString      dseq_;    // buffer for decoded read sequence
	BTString         dqual_;   // buffer for decoded quality sequence
};

static inline std::ostream& printPct(
							  std::ostream& os,
							  uint64_t num,
							  uint64_t denom)
{
	double pct = 0.0f;
	if(denom != 0) { pct = 100.0 * (double)num / (double)denom; }
	os << fixed << setprecision(2) << pct << '%';
	return os;
}

/**
 * Print a friendly summary of:
 *
 *  1. How many reads were aligned and had one or more alignments
 *     reported
 *  2. How many reads exceeded the -m or -M ceiling and therefore had
 *     their alignments suppressed or sampled
 *  3. How many reads failed to align entirely
 *
 * Optionally print a series of Hadoop streaming-style counter updates
 * with similar information.
 */
template <typename index_t>
void AlnSink<index_t>::printAlSumm(
								   const ReportingMetrics& met,
								   size_t repThresh,   // threshold for uniqueness, or max if no thresh
								   bool discord,       // looked for discordant alignments
								   bool mixed,         // looked for unpaired alignments where paired failed?
								   bool hadoopOut)     // output Hadoop counters?
{
	// NOTE: there's a filtering step at the very beginning, so everything
	// being reported here is post filtering
#if 0
	bool canRep = repThresh != MAX_SIZE_T;
	if(hadoopOut) {
		cerr << "reporter:counter:Centrifuge,Reads processed," << met.nread << endl;
	}
	uint64_t totread = met.nread;
	if(totread > 0) {
		cerr << "" << met.nread << " reads (or pairs); of these:" << endl;
	} else {
		assert_eq(0, met.npaired);
		assert_eq(0, met.nunpaired);
		cerr << "" << totread << " reads (or pairs)" << endl;
	}
	if(totread > 0) {
		// Concordants
		cerr << "    " << met.nconcord << " (";
		printPct(cerr, met.nconcord, met.npaired);
		cerr << ") classified 0 times" << endl;
		
        // Print the number that aligned concordantly exactly once
        cerr << "    " << met.nconcord_uni << " (";
        printPct(cerr, met.nconcord_uni, met.npaired);
        cerr << ") classified exactly 1 time" << endl;
        
#if 0
        // Print the number that aligned concordantly more than once
        cerr << "    " << met.nconcord_uni2 << " (";
        printPct(cerr, met.nconcord_uni2, met.npaired);
        cerr << ") classified >1 times" << endl;
#endif
	}
    
#if 0
	uint64_t totunpair = met.nunpaired;
	uint64_t tot_al_cand = totunpair + totpair*2;
	uint64_t tot_al =
	(met.nconcord_uni + met.nconcord_rep)*2 +
	(met.ndiscord)*2 +
	met.nunp_0_uni +
	met.nunp_0_rep + 
	met.nunp_uni +
	met.nunp_rep;
	assert_leq(tot_al, tot_al_cand);
	printPct(cerr, tot_al, tot_al_cand);
#endif
	cerr << " overall classification rate" << endl;
#endif
}

/**
 * Return true iff the read in rd1/rd2 matches the last read handled, which
 * should still be in rd1_/rd2_.
 */
template <typename index_t>
bool AlnSinkWrap<index_t>::sameRead(
									// One of the other of rd1, rd2 will = NULL if read is unpaired
									const Read* rd1,      // new mate #1
									const Read* rd2,      // new mate #2
									bool qualitiesMatter) // aln policy distinguishes b/t quals?
{
	bool same = false;
	if(rd1_ != NULL || rd2_ != NULL) {
		// This is not the first time the sink was initialized with
		// a read.  Check if new read/pair is identical to previous
		// read/pair
		if((rd1_ == NULL) == (rd1 == NULL) &&
		   (rd2_ == NULL) == (rd2 == NULL))
		{
			bool m1same = (rd1 == NULL && rd1_ == NULL);
			if(!m1same) {
				assert(rd1 != NULL);
				assert(rd1_ != NULL);
				m1same = Read::same(
									rd1->patFw,  // new seq
									rd1->qual,   // new quals
									rd1_->patFw, // old seq
									rd1_->qual,  // old quals
									qualitiesMatter);
			}
			if(m1same) {
				bool m2same = (rd2 == NULL && rd2_ == NULL);
				if(!m2same) {
					m2same = Read::same(
										rd2->patFw,  // new seq
										rd2->qual,   // new quals
										rd2_->patFw, // old seq
										rd2_->qual,  // old quals
										qualitiesMatter);
				}
				same = m2same;
			}
		}
	}
	return same;
}

/**
 * Initialize the wrapper with a new read pair and return an integer >= -1
 * indicating which stage the aligner should start at.  If -1 is returned, the
 * aligner can skip the read entirely.  Checks if the new read pair is
 * identical to the previous pair.  If it is, then we return the id of the
 * first stage to run.
 */
template <typename index_t>
int AlnSinkWrap<index_t>::nextRead(
								   // One of the other of rd1, rd2 will = NULL if read is unpaired
								   const Read* rd1,      // new mate #1
								   const Read* rd2,      // new mate #2
								   TReadId rdid,         // read ID for new pair
								   bool qualitiesMatter) // aln policy distinguishes b/t quals?
{
	assert(!init_);
	assert(rd1 != NULL || rd2 != NULL);
	init_ = true;
	// Keep copy of new read, so that we can compare it with the
	// next one
	if(rd1 != NULL) {
		rd1_ = rd1;
	} else rd1_ = NULL;
	if(rd2 != NULL) {
		rd2_ = rd2;
	} else rd2_ = NULL;
	rdid_ = rdid;
	// Caller must now align the read
	maxed1_ = false;
	maxed2_ = false;
	maxedOverall_ = false;
	bestPair_ = best2Pair_ =
	bestUnp1_ = best2Unp1_ =
	bestUnp2_ = best2Unp2_ = std::numeric_limits<THitInt>::min();
    bestSplicedPair_ = best2SplicedPair_ =
    bestSplicedUnp1_ = best2SplicedUnp1_ =
    bestSplicedUnp2_ = best2SplicedUnp2_ = 0;
	rs_.clear();     // clear out paired-end alignments
	st_.nextRead(readIsPair()); // reset state
	assert(empty());
	assert(!maxed());
	// Start from the first stage
	return 0;
}

/**
 * Inform global, shared AlnSink object that we're finished with this read.
 * The global AlnSink is responsible for updating counters, creating the output
 * record, and delivering the record to the appropriate output stream.
 *
 * What gets reported for a paired-end alignment?
 *
 * 1. If there are reportable concordant alignments, report those and stop
 * 2. If there are reportable discordant alignments, report those and stop
 * 3. If unpaired alignments can be reported:
 *    3a. Report 
 #
 * Update metrics.  Only ambiguity is: what if a pair aligns repetitively and
 * one of its mates aligns uniquely?
 *
 * 	uint64_t al;   // # mates w/ >= 1 reported alignment
 *  uint64_t unal; // # mates w/ 0 alignments
 *  uint64_t max;  // # mates withheld for exceeding -M/-m ceiling
 *  uint64_t al_concord;  // # pairs w/ >= 1 concordant alignment
 *  uint64_t al_discord;  // # pairs w/ >= 1 discordant alignment
 *  uint64_t max_concord; // # pairs maxed out
 *  uint64_t unal_pair;   // # pairs where neither mate aligned
 */
template <typename index_t>
void AlnSinkWrap<index_t>::finishRead(
									  const SeedResults<index_t> *sr1, // seed alignment results for mate 1
									  const SeedResults<index_t> *sr2, // seed alignment results for mate 2
									  bool               exhaust1,     // mate 1 exhausted?
									  bool               exhaust2,     // mate 2 exhausted?
									  bool               nfilt1,       // mate 1 N-filtered?
									  bool               nfilt2,       // mate 2 N-filtered?
									  bool               scfilt1,      // mate 1 score-filtered?
									  bool               scfilt2,      // mate 2 score-filtered?
									  bool               lenfilt1,     // mate 1 length-filtered?
									  bool               lenfilt2,     // mate 2 length-filtered?
									  bool               qcfilt1,      // mate 1 qc-filtered?
									  bool               qcfilt2,      // mate 2 qc-filtered?
									  bool               sortByScore,  // prioritize alignments by score
									  RandomSource&      rnd,          // pseudo-random generator
									  ReportingMetrics&  met,          // reporting metrics
									  SpeciesMetrics&    smet,         // species metrics
									  const PerReadMetrics& prm,       // per-read metrics
									  bool suppressSeedSummary,        // = true
									  bool suppressAlignments)         // = false
{
	obuf_.clear();
	OutputQueueMark qqm(g_.outq(), obuf_, rdid_, threadid_);
	assert(init_);
	if(!suppressSeedSummary) {
		if(sr1 != NULL) {
			assert(rd1_ != NULL);
			// Mate exists and has non-empty SeedResults
			g_.reportSeedSummary(obuf_, *rd1_, rdid_, threadid_, *sr1, true);
		} else if(rd1_ != NULL) {
			// Mate exists but has NULL SeedResults
			g_.reportEmptySeedSummary(obuf_, *rd1_, rdid_, true);
		}
		if(sr2 != NULL) {
			assert(rd2_ != NULL);
			// Mate exists and has non-empty SeedResults
			g_.reportSeedSummary(obuf_, *rd2_, rdid_, threadid_, *sr2, true);
		} else if(rd2_ != NULL) {
			// Mate exists but has NULL SeedResults
			g_.reportEmptySeedSummary(obuf_, *rd2_, rdid_, true);
		}
	}

	// TODO FB: Cconsider counting species here, and allow to disable counting

	if(!suppressAlignments) {
		// Ask the ReportingState what to report
		st_.finish();
		uint64_t nconcord = 0;
		bool pairMax = false;
		st_.getReport(nconcord);
		assert_leq(nconcord, rs_.size());
		assert_gt(rp_.khits, 0);
		met.nread++;

		if(readIsPair()) {
			met.npaired++;
		} else {
			met.nunpaired++;
		}
		// Report concordant paired-end alignments if possible
		if(nconcord > 0) {
			AlnSetSumm concordSumm(rd1_, rd2_, &rs_);
			// Possibly select a random subset
			size_t off;
			if(sortByScore) {
				// Sort by score then pick from low to high
				off = selectByScore(&rs_, nconcord, select_, rnd);
			} else {
				// Select subset randomly
				off = selectAlnsToReport(rs_, nconcord, select_, rnd);
			}
			assert_lt(off, rs_.size());
			_unused(off); // make production build happy
			assert(!select_.empty());
			g_.reportHits(
						  obuf_,
						  threadid_,
						  rd1_,
						  rd2_,
						  rdid_,
						  select_,
						  NULL,
						  &rs_,
						  NULL,
						  pairMax,
						  concordSumm,
                          prm,
						  smet);
			if(pairMax) {
				// met.nconcord_rep++;
			} else {
				met.nconcord_uni++;
				assert(!rs_.empty());
				if(rs_.size() == 1) {
					// met.nconcord_uni1++;
				} else {
					// met.nconcord_uni2++;
				}
			}
			init_ = false;
			// write read to file
			//g_.outq().finishRead(obuf_, rdid_, threadid_);
			return;
		}
		
#if 0
		// Update counters given that one mate didn't align
		if(readIsPair()) {
			met.nconcord_0++;
		}
		if(rd1_ != NULL) {
			if(nunpair1 > 0) {
				// Update counters
				if(readIsPair()) {
					if(unpair1Max) met.nunp_0_rep++;
					else {
						met.nunp_0_uni++;
						assert(!rs1u_.empty());
						if(rs1u_.size() == 1) {
							met.nunp_0_uni1++;
						} else {
							met.nunp_0_uni2++;
						}
					}
				} else {
					if(unpair1Max) met.nunp_rep++;
					else {
						met.nunp_uni++;
						assert(!rs1u_.empty());
						if(rs1u_.size() == 1) {
							met.nunp_uni1++;
						} else {
							met.nunp_uni2++;
						}
					}
				}
			} else if(unpair1Max) {
				// Update counters
				if(readIsPair())   met.nunp_0_rep++;
				else               met.nunp_rep++;
			} else {
				// Update counters
				if(readIsPair())   met.nunp_0_0++;
				else               met.nunp_0++;
			}
		}
		if(rd2_ != NULL) {
			if(nunpair2 > 0) {
				// Update counters
				if(readIsPair()) {
					if(unpair2Max) met.nunp_0_rep++;
					else {
						assert(!rs2u_.empty());
						met.nunp_0_uni++;
						if(rs2u_.size() == 1) {
							met.nunp_0_uni1++;
						} else {
							met.nunp_0_uni2++;
						}
					}
				} else {
					if(unpair2Max) met.nunp_rep++;
					else {
						assert(!rs2u_.empty());
						met.nunp_uni++;
						if(rs2u_.size() == 1) {
							met.nunp_uni1++;
						} else {
							met.nunp_uni2++;
						}
					}
				}
			} else if(unpair2Max) {
				// Update counters
				if(readIsPair())   met.nunp_0_rep++;
				else               met.nunp_rep++;
			} else {
				// Update counters
				if(readIsPair())   met.nunp_0_0++;
				else               met.nunp_0++;
			}
		}
        
#endif
	} // if(suppress alignments)
	init_ = false;
	return;
}

/**
 * Called by the aligner when a new unpaired or paired alignment is
 * discovered in the given stage.  This function checks whether the
 * addition of this alignment causes the reporting policy to be
 * violated (by meeting or exceeding the limits set by -k, -m, -M),
 * in which case true is returned immediately and the aligner is
 * short circuited.  Otherwise, the alignment is tallied and false
 * is returned.
 */
template <typename index_t>
bool AlnSinkWrap<index_t>::report(int stage,
								  const AlnRes* rs)
{
	assert(init_);
	assert(rs != NULL);
    st_.foundConcordant();
    rs_.push_back(*rs);

	// Tally overall alignment score
	TAlScore score = rs->score();
	// Update best score so far
    if(score > bestPair_) {
        best2Pair_ = bestPair_;
        bestPair_ = score;
    } else if(score > best2Pair_) {
        best2Pair_ = score;
    }
	return st_.done();
}

/**
 * rs1 (possibly together with rs2 if reads are paired) are populated with
 * alignments.  Here we prioritize them according to alignment score, and
 * some randomness to break ties.  Priorities are returned in the 'select'
 * list.
 */
template <typename index_t>
size_t AlnSinkWrap<index_t>::selectByScore(
										   const EList<AlnRes>* rs,    // alignments to select from (mate 1)
										   uint64_t             num,    // number of alignments to select
										   EList<size_t>&       select, // prioritized list to put results in
										   RandomSource&        rnd)
const
{
	assert(init_);
	assert(repOk());
	assert_gt(num, 0);
	assert(rs != NULL);
	size_t sz = rs->size(); // sz = # alignments found
	assert_leq(num, sz);
	if(sz < num) {
		num = sz;
	}
	// num = # to select
	if(sz < 1) {
		return 0;
	}
	select.resize((size_t)num);
	// Use 'selectBuf_' as a temporary list for sorting purposes
	EList<std::pair<TAlScore, size_t> >& buf =
	const_cast<EList<std::pair<TAlScore, size_t> >& >(selectBuf_);
	buf.resize(sz);
	// Sort by score.  If reads are pairs, sort by sum of mate scores.
	for(size_t i = 0; i < sz; i++) {
        buf[i].first = (*rs)[i].score();
		buf[i].second = i; // original offset
	}
	buf.sort(); buf.reverse(); // sort in descending order by score
	
	// Randomize streaks of alignments that are equal by score
	size_t streak = 0;
	for(size_t i = 1; i < buf.size(); i++) {
		if(buf[i].first == buf[i-1].first) {
			if(streak == 0) { streak = 1; }
			streak++;
		} else {
			if(streak > 1) {
				assert_geq(i, streak);
				buf.shufflePortion(i-streak, streak, rnd);
			}
			streak = 0;
		}
	}
	if(streak > 1) {
		buf.shufflePortion(buf.size() - streak, streak, rnd);
	}
	
	for(size_t i = 0; i < num; i++) { select[i] = buf[i].second; }
    
    if(!secondary_) {
        assert_geq(buf.size(), select.size());
        for(size_t i = 0; i + 1 < select.size(); i++) {
            if(buf[i].first != buf[i+1].first) {
                select.resize(i+1);
                break;
            }
        }
    }
    
	// Returns index of the representative alignment, but in 'select' also
	// returns the indexes of the next best selected alignments in order by
	// score.
	return selectBuf_[0].second;
}

/**
 * Given that rs is already populated with alignments, consider the
 * alignment policy and make random selections where necessary.  E.g. if we
 * found 10 alignments and the policy is -k 2 -m 20, select 2 alignments at
 * random.  We "select" an alignment by setting the parallel entry in the
 * 'select' list to true.
 *
 * Return the "representative" alignment.  This is simply the first one
 * selected.  That will also be what SAM calls the "primary" alignment.
 */
template <typename index_t>
size_t AlnSinkWrap<index_t>::selectAlnsToReport(
												const EList<AlnRes>& rs,     // alignments to select from
												uint64_t             num,    // number of alignments to select
												EList<size_t>&       select, // list to put results in
												RandomSource&        rnd)
const
{
	assert(init_);
	assert(repOk());
	assert_gt(num, 0);
	size_t sz = rs.size();
	if(sz < num) {
		num = sz;
	}
	if(sz < 1) {
		return 0;
	}
	select.resize((size_t)num);
	if(sz == 1) {
		assert_eq(1, num);
		select[0] = 0;
		return 0;
	}
	// Select a random offset into the list of alignments
	uint32_t off = rnd.nextU32() % (uint32_t)sz;
	uint32_t offOrig = off;
	// Now take elements starting at that offset, wrapping around to 0 if
	// necessary.  Leave the rest.
	for(size_t i = 0; i < num; i++) {
		select[i] = off;
		off++;
		if(off == sz) {
			off = 0;
		}
	}
	return offOrig;
}

#define NOT_SUPPRESSED !suppress_[field++]
#define BEGIN_FIELD { \
if(firstfield) firstfield = false; \
else o.append('\t'); \
}
#define WRITE_TAB { \
if(firstfield) firstfield = false; \
else o.append('\t'); \
}
#define WRITE_NUM(o, x) { \
itoa10(x, buf); \
o.append(buf); \
}

#define WRITE_STRING(o, x) { \
o.append(x.c_str()); \
}

template <typename T> inline 
void appendNumber(BTString& o, const T x, char* buf) {
	itoa10<T>(x, buf);
	o.append(buf);
}


/**
 * Print a seed summary to the first output stream in the outs_ list.
 */
template <typename index_t>
void AlnSink<index_t>::reportSeedSummary(
										 BTString&          o,
										 const Read&        rd,
										 TReadId            rdid,
										 size_t             threadId,
										 const SeedResults<index_t>& rs,
										 bool               getLock)
{
#if 0
	appendSeedSummary(
					  o,                     // string to write to
					  rd,                    // read
					  rdid,                  // read id
					  rs.numOffs()*2,        // # seeds tried
					  rs.nonzeroOffsets(),   // # seeds with non-empty results
					  rs.numRanges(),        // # ranges for all seed hits
					  rs.numElts(),          // # elements for all seed hits
					  rs.numOffs(),          // # seeds tried from fw read
					  rs.nonzeroOffsetsFw(), // # seeds with non-empty results from fw read
					  rs.numRangesFw(),      // # ranges for seed hits from fw read
					  rs.numEltsFw(),        // # elements for seed hits from fw read
					  rs.numOffs(),          // # seeds tried from rc read
					  rs.nonzeroOffsetsRc(), // # seeds with non-empty results from fw read
					  rs.numRangesRc(),      // # ranges for seed hits from fw read
					  rs.numEltsRc());       // # elements for seed hits from fw read
#endif
}

/**
 * Print an empty seed summary to the first output stream in the outs_ list.
 */
template <typename index_t>
void AlnSink<index_t>::reportEmptySeedSummary(
											  BTString&          o,
											  const Read&        rd,
											  TReadId            rdid,
											  size_t             threadId,
											  bool               getLock)
{
	appendSeedSummary(
					  o,                     // string to append to
					  rd,                    // read
					  rdid,                  // read id
					  0,                     // # seeds tried
					  0,                     // # seeds with non-empty results
					  0,                     // # ranges for all seed hits
					  0,                     // # elements for all seed hits
					  0,                     // # seeds tried from fw read
					  0,                     // # seeds with non-empty results from fw read
					  0,                     // # ranges for seed hits from fw read
					  0,                     // # elements for seed hits from fw read
					  0,                     // # seeds tried from rc read
					  0,                     // # seeds with non-empty results from fw read
					  0,                     // # ranges for seed hits from fw read
					  0);                    // # elements for seed hits from fw read
}

/**
 * Print the given string.  If ws = true, print only up to and not
 * including the first space or tab.  Useful for printing reference
 * names.
 */
template<typename T>
static inline void printUptoWs(
							   BTString& s,
							   const T& str,
							   bool chopws)
{
	size_t len = str.length();
	for(size_t i = 0; i < len; i++) {
		if(!chopws || (str[i] != ' ' && str[i] != '\t')) {
			s.append(str[i]);
		} else {
			break;
		}
	}
}

/**
 * Append a batch of unresolved seed alignment summary results (i.e.
 * seed alignments where all we know is the reference sequence aligned
 * to and its SA range, not where it falls in the reference
 * sequence) to the given output stream in Bowtie's seed-sumamry
 * verbose-mode format.
 *
 * The seed summary format is:
 *
 *  - One line per read
 *  - A typical line consists of a set of tab-delimited fields:
 *
 *    1. Read name
 *    2. Total number of seeds extracted from the read
 *    3. Total number of seeds that aligned to the reference at
 *       least once (always <= field 2)
 *    4. Total number of distinct BW ranges found in all seed hits
 *       (always >= field 3)
 *    5. Total number of distinct BW elements found in all seed
 *       hits (always >= field 4)
 *    6-9.:   Like 2-5. but just for seeds extracted from the
 *            forward representation of the read
 *    10-13.: Like 2-5. but just for seeds extracted from the
 *            reverse-complement representation of the read
 *
 *    Note that fields 6 and 10 should add to field 2, 7 and 11
 *    should add to 3, etc.
 *
 *  - Lines for reads that are filtered out for any reason (e.g. too
 *    many Ns) have columns 2 through 13 set to 0.
 */
template <typename index_t>
void AlnSink<index_t>::appendSeedSummary(
										 BTString&     o,
										 const Read&   rd,
										 const TReadId rdid,
										 size_t        seedsTried,
										 size_t        nonzero,
										 size_t        ranges,
										 size_t        elts,
										 size_t        seedsTriedFw,
										 size_t        nonzeroFw,
										 size_t        rangesFw,
										 size_t        eltsFw,
										 size_t        seedsTriedRc,
										 size_t        nonzeroRc,
										 size_t        rangesRc,
										 size_t        eltsRc)
{
	char buf[1024];
	bool firstfield = true;
	//
	// Read name
	//
	BEGIN_FIELD;
	printUptoWs(o, rd.name, true);
	
	//
	// Total number of seeds tried
	//
	BEGIN_FIELD;
	WRITE_NUM(o, seedsTried);
	
	//
	// Total number of seeds tried where at least one range was found.
	//
	BEGIN_FIELD;
	WRITE_NUM(o, nonzero);
	
	//
	// Total number of ranges found
	//
	BEGIN_FIELD;
	WRITE_NUM(o, ranges);
	
	//
	// Total number of elements found
	//
	BEGIN_FIELD;
	WRITE_NUM(o, elts);
	
	//
	// The same four numbers, but only for seeds extracted from the
	// forward read representation.
	//
	BEGIN_FIELD;
	WRITE_NUM(o, seedsTriedFw);
	
	BEGIN_FIELD;
	WRITE_NUM(o, nonzeroFw);
	
	BEGIN_FIELD;
	WRITE_NUM(o, rangesFw);
	
	BEGIN_FIELD;
	WRITE_NUM(o, eltsFw);
	
	//
	// The same four numbers, but only for seeds extracted from the
	// reverse complement read representation.
	//
	BEGIN_FIELD;
	WRITE_NUM(o, seedsTriedRc);
	
	BEGIN_FIELD;
	WRITE_NUM(o, nonzeroRc);
	
	BEGIN_FIELD;
	WRITE_NUM(o, rangesRc);
	
	BEGIN_FIELD;
	WRITE_NUM(o, eltsRc);
	
	o.append('\n');
}


inline
void appendReadID(BTString& o, const BTString & rd_name) {
    size_t namelen = rd_name.length();
    if(namelen >= 2 &&
       rd_name[namelen-2] == '/' &&
       (rd_name[namelen-1] == '1' || rd_name[namelen-1] == '2' || rd_name[namelen-1] == '3'))
    {
        namelen -= 2;
    }
    for(size_t i = 0; i < namelen; i++) {
        if(isspace(rd_name[i])) {
            break;
        }
        o.append(rd_name[i]);
    }
}

inline
void appendSeqID(BTString& o, const AlnRes* rs, const std::map<uint64_t, TaxonomyNode>& tree) { 
    bool leaf = true;
    std::map<uint64_t, TaxonomyNode>::const_iterator itr = tree.find(rs->taxID());
    if(itr != tree.end()) {
        const TaxonomyNode& node = itr->second;
        leaf = node.leaf;
    }

    // unique ID
    if(leaf) {
        o.append(rs->uid().c_str());
    } else {
        o.append(get_tax_rank_string(rs->taxRank()));
    }
}

inline 
void appendTaxID(BTString& o, uint64_t tid) {

	char buf[1024];
    uint64_t tid1 = tid & 0xffffffff;
    uint64_t tid2 = tid >> 32;
    itoa10<int64_t>(tid1, buf);
    o.append(buf);

    if(tid2 > 0) {
        o.append(".");
        itoa10<int64_t>(tid2, buf);
        o.append(buf);
    }
}


enum FIELD_DEF {
	PLACEHOLDER=0,
	PLACEHOLDER_STAR,
	PLACEHOLDER_ZERO,
    READ_ID,
	SEQ_ID,
	TAX_ID,
	TAX_RANK,
	TAX_NAME,
	SCORE,
	SCORE2,
	HIT_LENGTH,
	QUERY_LENGTH,
	NUM_MATCHES,
	SEQ,
	SEQ1,
	SEQ2,
	QUAL,
	QUAL1,
	QUAL2
};

/**
 * Append a single hit to the given output stream in Bowtie's
 * verbose-mode format.
 */
template <typename index_t>
void AlnSinkSam<index_t>::appendMate(
                                     Ebwt<index_t>& ebwt,
									 BTString&      o,           // append to this string
									 const Read&    rd,
									 const Read*    rdo,
									 const TReadId  rdid,
									 AlnRes* rs,
									 AlnRes* rso,
									 const AlnSetSumm& summ,
									 const PerReadMetrics& prm,
									 SpeciesMetrics& sm,
									 size_t n_results)
{
	if(rs == NULL) {
		return;
	}

	char buf[1024];
	bool firstfield = true;
	uint64_t taxid =  rs->taxID();
	const basic_string<char> empty_string = "";

	for (int i=0; i < this->tab_fmt_cols_.size(); ++i) {
		BEGIN_FIELD;
		switch (this->tab_fmt_cols_[i]) {
			case READ_ID:      appendReadID(o, rd.name); break;
			case SEQ_ID:       appendSeqID(o, rs, ebwt.tree()); break;
			case SEQ:          o.append((string(rd.patFw.toZBuf()) + 
										(rdo == NULL? "" : "N" + string(rdo->patFw.toZBuf()))).c_str()); break;
			case QUAL:         o.append((string(rd.qual.toZBuf()) + 
										(rdo == NULL? "" : "I" + string(rdo->qual.toZBuf()))).c_str()); break;

			case SEQ1:         o.append(rd.patFw.toZBuf()); break;
			case QUAL1:        o.append(rd.qual.toZBuf()); break;
			case SEQ2:         o.append(rdo == NULL? empty_string.c_str() : rdo->patFw.toZBuf()); break;
			case QUAL2:        o.append(rdo == NULL? empty_string.c_str() : rdo->qual.toZBuf()); break;

			case TAX_ID:       appendTaxID(o, taxid); break;
			case TAX_RANK:     o.append(get_tax_rank_string(rs->taxRank())); break;
			case TAX_NAME:     o.append(find_or_use_default(ebwt.name(), taxid, empty_string).c_str()); break;

			case SCORE:        appendNumber<uint64_t>(o, rs->score(), buf); break;
			case SCORE2:       appendNumber<uint64_t>(o,
									   (summ.secbest().valid()? summ.secbest().score() : 0),
									   buf); break;
			case HIT_LENGTH:   appendNumber<uint64_t>(o, rs->summedHitLen(), buf); break;
			case QUERY_LENGTH: appendNumber<uint64_t>(o, 
									   (rd.patFw.length() + (rdo != NULL ? rdo->patFw.length() : 0)),
									   buf); break;
			case NUM_MATCHES:  appendNumber<uint64_t>(o, n_results, buf); break;
			case PLACEHOLDER:  o.append("");
			case PLACEHOLDER_STAR:  o.append("*");
			case PLACEHOLDER_ZERO:  o.append("0");
			default: ;
		}

	}
	o.append("\n");


	// species counting
	sm.addSpeciesCounts(
                        rs->taxID(),
                        rs->score(),
                        rs->max_score(),
                        rs->summedHitLen(),
                        1.0 / n_results,
                        (uint32_t)n_results);

	// only count k-mers if the read is unique
    if (n_results == 1) {
		for (size_t i = 0; i< rs->nReadPositions(); ++i) {
			sm.addAllKmers(rs->taxID(),
                           rs->isFw()? rd.patFw : rd.patRc,
                           rs->readPositions(i).first,
                           rs->readPositions(i).second);
		}
	}

//    (sc[rs->speciesID_])++;
   
}

// #include <iomanip>

/**
 * Initialize state machine with a new read.  The state we start in depends
 * on whether it's paired-end or unpaired.
 */
void ReportingState::nextRead(bool paired) {
    paired_ = paired;
    state_ = CONCORDANT_PAIRS;
    doneConcord_ = false;
    exitConcord_ = ReportingState::EXIT_DID_NOT_EXIT;
    done_ = false;
    nconcord_ = 0;
}

/**
 * Caller uses this member function to indicate that one additional
 * concordant alignment has been found.
 */
bool ReportingState::foundConcordant() {
    assert_geq(state_, ReportingState::CONCORDANT_PAIRS);
    assert(!doneConcord_);
    nconcord_++;
    if(doneConcord_) {
        // If we're finished looking for concordant alignments, do we have to
        // continue on to search for unpaired alignments?  Only if our exit
        // from the concordant stage is EXIT_SHORT_CIRCUIT_M.  If it's
        // EXIT_SHORT_CIRCUIT_k or EXIT_WITH_ALIGNMENTS, we can skip unpaired.
        assert_neq(ReportingState::EXIT_NO_ALIGNMENTS, exitConcord_);
    }
    return done();
}

/**
 * Caller uses this member function to indicate that one additional unpaired
 * mate alignment has been found for the specified mate.
 */
bool ReportingState::foundUnpaired(bool mate1) {
    return done();
}

/**
 * Called to indicate that the aligner has finished searching for
 * alignments.  This gives us a chance to finalize our state.
 *
 * TODO: Keep track of short-circuiting information.
 */
void ReportingState::finish() {
    if(!doneConcord_) {
        doneConcord_ = true;
        exitConcord_ =
        ((nconcord_ > 0) ?
         ReportingState::EXIT_WITH_ALIGNMENTS :
         ReportingState::EXIT_NO_ALIGNMENTS);
    }
    assert_gt(exitConcord_, EXIT_DID_NOT_EXIT);
    done_ = true;
    assert(done());
}


/**
 * Populate given counters with the number of various kinds of alignments
 * to report for this read.  Concordant alignments are preferable to (and
 * mutually exclusive with) discordant alignments, and paired-end
 * alignments are preferable to unpaired alignments.
 *
 * The caller also needs some additional information for the case where a
 * pair or unpaired read aligns repetitively.  If the read is paired-end
 * and the paired-end has repetitive concordant alignments, that should be
 * reported, and 'pairMax' is set to true to indicate this.  If the read is
 * paired-end, does not have any conordant alignments, but does have
 * repetitive alignments for one or both mates, then that should be
 * reported, and 'unpair1Max' and 'unpair2Max' are set accordingly.
 *
 * Note that it's possible in the case of a paired-end read for the read to
 * have repetitive concordant alignments, but for one mate to have a unique
 * unpaired alignment.
 */
void ReportingState::getReport(uint64_t& nconcordAln) const // # concordant alignments to report
{
    nconcordAln = 0;
    assert_gt(p_.khits, 0);
    // Do we have 1 or more concordant alignments to report?
    if(exitConcord_ == ReportingState::EXIT_SHORT_CIRCUIT_k) {
        // k at random
        assert_geq(nconcord_, (uint64_t)p_.khits);
        nconcordAln = p_.khits;
        return;
    } else if(exitConcord_ == ReportingState::EXIT_WITH_ALIGNMENTS) {
        assert_gt(nconcord_, 0);
        // <= k at random
        nconcordAln = min<uint64_t>(nconcord_, p_.khits);
        return;
    }
}

#if 0
/**
 * Given the number of alignments in a category, check whether we
 * short-circuited out of the category.  Set the done and exit arguments to
 * indicate whether and how we short-circuited.
 */
inline void ReportingState::areDone(
                                    uint64_t cnt,    // # alignments in category
                                    bool& done,      // out: whether we short-circuited out of category
                                    int& exit) const // out: if done, how we short-circuited (-k? -m? etc)
{
    assert(!done);
    // Have we exceeded the -k limit?
    assert_gt(p_.khits, 0);
    assert_gt(p_.mhits, 0);
    if(cnt >= (uint64_t)p_.khits && !p_.mhitsSet()) {
        done = true;
        exit = ReportingState::EXIT_SHORT_CIRCUIT_k;
    }
    // Have we exceeded the -m or -M limit?
    else if(p_.mhitsSet() && cnt > (uint64_t)p_.mhits) {
        done = true;
        assert(p_.msample);
        exit = ReportingState::EXIT_SHORT_CIRCUIT_M;
    }
}
#endif

#endif /*ndef ALN_SINK_H_*/
