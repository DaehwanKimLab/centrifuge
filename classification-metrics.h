/*
 * classification-metrics.h
 *
 *  Created on: Jan 19, 2017
 *      Author: fbreitwieser
 */

#ifndef CLASSIFICATION_METRICS_H_
#define CLASSIFICATION_METRICS_H_

#include "assert_helpers.h"
#include "bt2_idx.h"
#include "bt2_io.h"
#include "hitmap.h"
#include <unordered_set>

typedef TIndexOffU index_t;

struct ReadCounts {
	uint64_t n_reads = 0;
	uint64_t total_score = 0; //summed score
	uint64_t total_hit_len = 0.0;   //summed hit length

	ReadCounts() {}

	ReadCounts(uint64_t _n_reads): n_reads(_n_reads) {}

	ReadCounts& operator+=(const ReadCounts& b) {
		n_reads += b.n_reads;
		total_score += b.total_score;
		total_hit_len += b.total_hit_len;
		return *this;
	}

	ReadCounts& addHit(const HitScore hs) {
		++n_reads;
		total_hit_len += hs[0];
		total_score += hs[1];
		return *this;
	}
};

uint64_t read_counts(uint64_t num) { return num; }
uint64_t read_counts(ReadCounts num) { return num.n_reads; }

ReadCounts operator+(ReadCounts a, const ReadCounts& b) {
	return a += b;
}

typedef set<uint64_t> UId_set;

/**
 * Metrics summarizing the classification information we have
 */
struct ClassificationMetrics {
	/*
    struct IDs {
        //EList<uint64_t, 5> ids;
        set<uint64_t> ids;
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
	 */


	ClassificationMetrics():mutex_m() {
		reset();
	}

	void reset() {
		counts.clear();
		//for(map<uint32_t, HyperLogLogPlusMinus<uint64_t> >::iterator it = this->kmers.begin(); it != this->kmers.end(); ++it) {
		//	it->second.reset();
		//} //TODO: is this required?
		kmers.clear();
		n_unclassified_reads = 0;
	}

	/*
	void init(
              const map<uint64_t, ReadCounts>& counts_,
              const map<uint64_t, HyperLogLogPlusMinus<uint64_t> >& kmers_,
              const map<IDs, uint64_t>& observed_)
	{
		counts = counts_;
		kmers = kmers_;
        observed = observed_;
        num_non_leaves = 0;
    }
	 */

	/**
	 * Merge (add) the counters in the given ReportingMetrics object
	 * into this object.  This is the only safe way to update a
	 * ReportingMetrics shared by multiple threads.
	 */
	void merge(const ClassificationMetrics& met, bool getLock = false) {
		ThreadSafe ts(&mutex_m, getLock);

		// update read count
		for(const auto& c:  met.counts) {
			counts[c.first] += c.second;
		}

		n_unclassified_reads += met.n_unclassified_reads;

		// update k-mers
		for(map<uint64_t, HyperLogLogPlusMinus<uint64_t> >::const_iterator it = met.kmers.begin(); it != met.kmers.end(); ++it) {
			kmers[it->first] += &(it->second);
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
		kmers[taxID].add(kmer);
		size_t i = begin;
		while (i+32 < len) {
			kmer = btdna.next_kmer(kmer,i);
			kmers[taxID].add(kmer);
			++i;
		}
	}

	size_t nDistinctKmers(uint64_t taxID) {
		return(kmers[taxID].cardinality());
	}

	template<typename RC>
	static void TaxID_EM(
			const map<UId_set, RC>& observed,  // Change to ReadCounts?!
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
		for(auto const itr = observed.begin(); itr != observed.end(); itr++) {
			const auto& ids = itr->first; // all ids assigned to the read set
			uint64_t count = read_counts(itr->second); // number of reads in the read set
			double psum = 0.0;
			for(auto tid : ids) {
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

			for(auto tid : ids) {
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

	template<typename RC>
	static void UID_EM(
			const map<UId_set, RC>& counts,
			const map<uint64_t, uint64_t> uid_to_num,
			const EList<double>& p,
			EList<double>& p_next,
			const EList<size_t>& len)
	{
		assert_eq(p.size(), len.size());

		// E step
		p_next.fill(0.0);

		// for each UID set
		for(map<UId_set, ReadCounts>::const_iterator itr = counts.begin(); itr != counts.end(); itr++) {
			const auto& uids = itr->first; // all UIDs in the set
			uint64_t count = read_counts(itr->second); // number of reads for this set of UIDs
			double psum = 0.0;
			for(auto uid : uids) {
				psum += p[uid_to_num.at(uid)];
			}

			if(psum == 0.0) continue;

			for(auto uid : uids) {
				p_next[uid_to_num.at(uid)] += (count * (p[uid_to_num.at(uid)] / psum));
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

	map<UId, array<double,2> > calculateAbundanceOnUIDs(const Ebwt<index_t>& ebwt) const {
		// only used in the end - filled by calculateAbundance

		//const map<uint64_t, TaxonomyNode>& tree = ebwt.tree();

		// The IDs and the lengths are in the UID space. Since not all UIDs will be detected,
		// we map from UIDs to NUMs to be more efficient in the EM algorithm

		// Lengths of sequences for each UID
		const index_t* uid_seq_lengths = ebwt.plen();  // sequence lengths in the UID space
		EList<size_t> seq_lengths; 				   // sequence lengths in the NUM space

		DEBUG_MSG("Calc abundance with " <<  counts.size() << " UID sets" << endl);

		// Initialize probabilities
		map<uint64_t, uint64_t> uid_to_num; // taxonomic ID to corresponding element of a list
		EList<double> p;

		for(const pair<UId_set,ReadCounts>& uid_set_count : counts) {
			const UId_set& uid_set = uid_set_count.first;
			uint64_t count = uid_set_count.second.n_reads;
			for(uint64_t uid : uid_set) {
				if(uid_to_num.find(uid) == uid_to_num.end()) {
					// DEBUG_MSG("Adding UID " << uid << endl);
					uid_to_num[uid] = p.size();
					p.push_back(1.0 / uid_set.size() * count);
					seq_lengths.push_back(uid_seq_lengths[uid]);
				} else {
					uint64_t num = uid_to_num[uid];
					assert_lt(num, p.size());
					p[num] += (1.0 / uid_set.size() * count);
				}
			}
		}

		assert_eq(p.size(), seq_lengths.size());

		{
			double sum = 0.0;
			for(size_t i = 0; i < p.size(); i++) {
				sum += (p[i] / seq_lengths[i]);
			}
			for(size_t i = 0; i < p.size(); i++) {
				p[i] = (p[i] / seq_lengths[i]) / sum;
			}
		}

		EList<double> p_next; p_next.resizeExact(p.size());
		EList<double> p_next2; p_next2.resizeExact(p.size());
		EList<double> p_r; p_r.resizeExact(p.size());
		EList<double> p_v; p_v.resizeExact(p.size());
		size_t num_iteration = 0;
		double diff = 0.0;
		while(true) {

			// Accelerated version of EM - SQUAREM iteration
			//    Varadhan, R. & Roland, C. Scand. J. Stat. 35, 335–353 (2008).
			//    Also, this algorithm is used in Sailfish - http://www.nature.com/nbt/journal/v32/n5/full/nbt.2862.html

			UID_EM(counts, uid_to_num, p, p_next, seq_lengths);
			UID_EM(counts, uid_to_num, p_next, p_next2, seq_lengths);
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
				UID_EM(counts, uid_to_num, p_next2, p_next, seq_lengths);
			}


			diff = 0.0;
			for(size_t i = 0; i < p.size(); i++) {
				diff += (p[i] > p_next[i] ? p[i] - p_next[i] : p_next[i] - p[i]);
			}
			if(diff < 0.0000000001) break;
			if(++num_iteration >= 1000) break;
			p = p_next;
		}

		cerr << "Number of iterations in EM algorithm: " << num_iteration << endl;
		cerr << "Probability diff. (P - P_prev) in the last iteration: " << diff << endl;

		// Return the abundance
		double sum = 0.0;
		// Calculate abundance normalized by genome size (standard)
		//   and without the genome size taken into account
		map<UId, array<double, 2> > abundance;
		for (size_t num = 0; num < p.size(); ++num) {
			sum += (p[num] * seq_lengths[num]);
		}

		for (const auto& un : uid_to_num) {
			assert_lt(un.second, p.size());
			abundance[un.first] = {{ p[un.second], p[un.second]*seq_lengths[un.second]/sum }};
		}
	
		return abundance;
	}


	template <typename T>
	map<TaxId, array<double,2> > calculateAbundanceOnTaxids(const Ebwt<T>& ebwt, uint8_t rank)
	{
		const map<TaxId, TaxonomyNode>& tree = ebwt.tree();

		// Lengths of genomes (or contigs)
		const map<uint64_t, uint64_t>& size_table = ebwt.size();

		// Find leaves
		set<uint64_t> leaves;
		for (auto itr = counts.begin(); itr != counts.end(); itr++) {
			const UId_set& ids = itr->first;
			for(uint64_t tid : ids) {
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
		for(auto itr = counts.begin(); itr != counts.end(); itr++) {
			const UId_set& ids = itr->first;
			for(uint64_t tid : ids) {
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
		for(auto itr = ancestors.begin(); itr != ancestors.end(); itr++) {
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

		// Initialize probabilities
		map<uint64_t, uint64_t> tid_to_num; // taxonomic ID to corresponding element of a list
		EList<double> p;
		EList<size_t> len; // genome lengths
		for(auto itr = counts.begin(); itr != counts.end(); itr++) {
			const UId_set& ids = itr->first;
			uint64_t count = itr->second;
			for(uint64_t tid : ids) {
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
					p.push_back(1.0 / ids.size() * count);
					map<uint64_t, uint64_t>::const_iterator size_itr = size_table.find(tid);
					if(size_itr != size_table.end()) {
						len.push_back(size_itr->second);
					} else {
						len.push_back(numeric_limits<size_t>::max());
					}
				} else {
					uint64_t num = tid_to_num[tid];
					assert_lt(num, p.size());
					p[num] += (1.0 / ids.size() * count);
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
			TaxID_EM(counts, ancestors, tid_to_num, p, p_next, len);
			TaxID_EM(counts, ancestors, tid_to_num, p_next, p_next2, len);
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
				TaxID_EM(counts, ancestors, tid_to_num, p_next2, p_next, len);
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

		// Return the abundance
		double sum = 0.0;
		// Calculate abundance normalized by genome size (standard)
		//   and without the genome size taken into account
		map<TaxId, array<double, 2> > abundance;
		for (size_t num = 0; num < p.size(); ++num) {
			sum += (p[num] * size_table.at(num));
		}

		for (const auto& un : tid_to_num) {
			assert_lt(un.second, p.size());
			abundance[un.first] = {{ p[un.second], p[un.second]*size_table.at(un.second)/sum }};
		}

		return abundance;
	}

	uint64_t n_unclassified_reads;
	map<UId_set, ReadCounts> counts;                        // read counts for each sequence
	map<TaxId, HyperLogLogPlusMinus<uint64_t> > kmers;    // unique k-mer count for each sequence

	MUTEX_T mutex_m;
};


#endif /* CLASSIFICATION_METRICS_H_ */
