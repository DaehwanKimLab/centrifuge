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

#ifndef ALIGNER_RESULT_H_
#define ALIGNER_RESULT_H_

#include <utility>
#include <limits>
#include <vector>
#include "mem_ids.h"
#include "ref_coord.h"
#include "read.h"
#include "filebuf.h"
#include "ds.h"
#include "edit.h"
#include "limit.h"

typedef int64_t TAlScore;

#define VALID_AL_SCORE(x)   ((x).score_ > MIN_I64)
#define VALID_SCORE(x)      ((x) > MIN_I64)
#define INVALIDATE_SCORE(x) ((x) = MIN_I64)

/**
 * A generic score object for an alignment.  Used for accounting during
 * SW and elsewhere.  Encapsulates the score, the number of N positions
 * and the number gaps in the alignment.
 *
 * The scale for 'score' is such that a perfect alignment score is 0
 * and a score with non-zero penalty is less than 0.  So differences
 * between scores work as expected, but interpreting an individual
 * score (larger is better) as a penalty (smaller is better) requires
 * taking the absolute value.
 */
class AlnScore {

public:

	/**
	 * Gapped scores are invalid until proven valid.
	 */
	inline AlnScore() {
		reset();
		invalidate();
		assert(!valid());
	}

	/**
	 * Gapped scores are invalid until proven valid.
	 */
	inline AlnScore(TAlScore score) {
		score_ = score;
	}
	
	/**
	 * Reset the score.
	 */
	void reset() {
		score_ = 0;
	}

	/**
	 * Return an invalid SwScore.
	 */
	inline static AlnScore INVALID() {
		AlnScore s;
		s.invalidate();
		assert(!s.valid());
		return s;
	}
	
	/**
	 * Return true iff this score has a valid value.
	 */
	inline bool valid() const {
		return score_ != MIN_I64;
	}

	/**
	 * Make this score invalid (and therefore <= all other scores).
	 */
	inline void invalidate() {
		score_ = MIN_I64;
		assert(!valid());
	}
	

	/**
	 * Scores are equal iff they're bitwise equal.
	 */
	inline bool operator==(const AlnScore& o) const {
		// Profiling shows cache misses on following line
		return VALID_AL_SCORE(*this) && VALID_AL_SCORE(o) && score_ == o.score_;
	}

	/**
	 * Return true iff the two scores are unequal.
	 */
	inline bool operator!=(const AlnScore& o) const {
		return !(*this == o);
	}

	/**
	 * Return true iff this score is >= score o.
	 */
	inline bool operator>=(const AlnScore& o) const {
		if(!VALID_AL_SCORE(o)) {
			if(!VALID_AL_SCORE(*this)) {
				// both invalid
				return false;
			} else {
				// I'm valid, other is invalid
				return true;
			}
		} else if(!VALID_AL_SCORE(*this)) {
			// I'm invalid, other is valid
			return false;
		}
		return score_ >= o.score_;
	}

	/**
	 * Return true iff this score is < score o.
	 */
	inline bool operator<(const AlnScore& o) const {
		return !operator>=(o);
	}

	/**
	 * Return true iff this score is <= score o.
	 */
	inline bool operator<=(const AlnScore& o) const {
		return operator<(o) || operator==(o);
	}
    
    /**
     * Return true iff this score is < score o.
     */
    inline bool operator>(const AlnScore& o) const {
        return !operator<=(o);
    }

	TAlScore score()   const { return  score_; }

    // Score accumulated so far (penalties are subtracted starting at 0)
	TAlScore score_;
};

static inline ostream& operator<<(ostream& os, const AlnScore& o) {
	os << o.score();
	return os;
}

// Forward declaration
class BitPairReference;


/**
 * Encapsulates an alignment result.  The result comprises:
 *
 * 1. All the nucleotide edits for both mates ('ned').
 * 2. All "edits" where an ambiguous reference char is resolved to an
 *    unambiguous char ('aed').
 * 3. The score for the alginment, including summary information about the
 *    number of gaps and Ns involved.
 * 4. The reference id, strand, and 0-based offset of the leftmost character
 *    involved in the alignment.
 * 5. Information about trimming prior to alignment and whether it was hard or
 *    soft.
 * 6. Information about trimming during alignment and whether it was hard or
 *    soft.  Local-alignment trimming is usually soft when aligning nucleotide
 *    reads.
 *
 * Note that the AlnRes, together with the Read and an AlnSetSumm (*and* the
 * opposite mate's AlnRes and Read in the case of a paired-end alignment),
 * should contain enough information to print an entire alignment record.
 *
 * TRIMMING
 *
 * Accounting for trimming is tricky.  Trimming affects:
 *
 * 1. The values of the trim* and pretrim* fields.
 * 2. The offsets of the Edits in the EList<Edit>s.
 * 3. The read extent, if the trimming is soft.
 * 4. The read extent and the read sequence and length, if trimming is hard.
 *
 * Handling 1. is not too difficult.  2., 3., and 4. are handled in setShape().
 */
class AlnRes {

public:

	AlnRes()
	{
		reset();
	}
    
    AlnRes(const AlnRes& other)
    {
        score_ = other.score_;
        max_score_ = other.max_score_;
        uid_ = other.uid_;
        tid_ = other.tid_;
        isLeaf_ = other.isLeaf_;
        taxRank_ = other.taxRank_;
        summedHitLen_ = other.summedHitLen_;
		readPositions_ = other.readPositions_;
		isFw_ = other.isFw_;
    }
    
    AlnRes& operator=(const AlnRes& other) {
        if(this == &other) return *this;
        score_ = other.score_;
        max_score_ = other.max_score_;
        uid_ = other.uid_;
        tid_ = other.tid_;
        isLeaf_ = other.isLeaf_;
        taxRank_ = other.taxRank_;
        summedHitLen_ = other.summedHitLen_;
		readPositions_ = other.readPositions_;
		isFw_ = other.isFw_;
        return *this;
    }
    
    ~AlnRes() {}

	/**
	 * Clear all contents.
	 */
	void reset() {
        score_ = 0;
        max_score_ = 0;
        uid_ = "";
        tid_ = 0;
        isLeaf_ = true;
        taxRank_ = RANK_UNKNOWN;
        summedHitLen_ = 0.0;
		readPositions_.clear();
    }
    
	/**
	 * Set alignment score for this alignment.
	 */
	void setScore(TAlScore score) {
		score_ = score;
	}

	TAlScore           score()          const { return score_;     }
    TAlScore           max_score()      const { return max_score_; }
    string             uid()            const { return uid_;   }
    uint64_t           taxID()          const { return tid_;   }
    bool               leaf()           const { return isLeaf_; }
    uint8_t            taxRank()        const { return taxRank_; }
    double             summedHitLen()   const { return summedHitLen_; }

	const EList<pair<uint32_t,uint32_t> >& readPositionsPtr() const { return readPositions_; }

	const pair<uint32_t,uint32_t> readPositions(size_t i) const { return readPositions_[i]; }
	size_t nReadPositions() const { return readPositions_.size(); }

	bool               isFw()           const { return isFw_;      }
    
   /**
	 * Print the sequence for the read that aligned using A, C, G and
	 * T.  This will simply print the read sequence (or its reverse
	 * complement).
	 */
 	void printSeq(
		const Read& rd,
		const BTDnaString* dns,
		BTString& o) const
    {
        assert(!rd.patFw.empty());
        ASSERT_ONLY(size_t written = 0);
        // Print decoded nucleotides
        assert(dns != NULL);
        size_t len = dns->length();
        size_t st = 0;
        size_t en = len;
        for(size_t i = st; i < en; i++) {
            int c = dns->get(i);
            assert_range(0, 3, c);
            o.append("ACGT"[c]);
            ASSERT_ONLY(written++);
        }
    }

	/**
	 * Print the quality string for the read that aligned.  This will
	 * simply print the read qualities (or their reverse).
	 */
 	void printQuals(
		const Read& rd,
		const BTString* dqs,
		BTString& o) const
    {
        assert(dqs != NULL);
        size_t len = dqs->length();
        // Print decoded qualities from upstream to downstream Watson
        for(size_t i = 1; i < len-1; i++) {
            o.append(dqs->get(i));
        }
    }
	

	/**
	 * Initialize new AlnRes.
	 */
	void init(
              TAlScore score,           // alignment score
              TAlScore max_score,
              const string& uniqueID,
              uint64_t taxID,
              bool leaf,
              uint8_t taxRank,
			  double summedHitLen,
			  const EList<pair<uint32_t, uint32_t> >& readPositions,
			  bool isFw)
    {
        score_  = score;
        max_score_ = max_score;
        uid_ = uniqueID;
        tid_ = taxID;
        isLeaf_ = leaf;
        taxRank_ = taxRank;
        summedHitLen_ = summedHitLen;
		readPositions_ = readPositions;
		isFw_ = isFw;
    }

protected:
	TAlScore     score_;        //
    TAlScore     max_score_;
    string       uid_;
    uint64_t     tid_;
    bool         isLeaf_;
    uint8_t      taxRank_;
    double       summedHitLen_; // sum of the length of all partial hits, divided by the number of genome matches
	bool         isFw_;
  
	EList<pair<uint32_t, uint32_t> > readPositions_;
};

typedef uint64_t TNumAlns;

/**
 * Encapsulates a concise summary of a set of alignment results for a
 * given pair or mate.  Referring to the fields of this object should
 * provide enough information to print output records for the read.
 */
class AlnSetSumm {

public:

	AlnSetSumm() { reset(); }

	/**
	 * Given an unpaired read (in either rd1 or rd2) or a read pair
	 * (mate 1 in rd1, mate 2 in rd2).
	 */
	explicit AlnSetSumm(
		const Read* rd1,
		const Read* rd2,
		const EList<AlnRes>* rs)
	{
		init(rd1, rd2, rs);
	}

	explicit AlnSetSumm(
		AlnScore best,
		AlnScore secbest)
	{
		init(best, secbest);
	}
	
	/**
	 * Set to uninitialized state.
	 */
	void reset() {
		best_.invalidate();
		secbest_.invalidate();
	}
	
    /**
     * Given all the paired and unpaired results involving mates #1 and #2,
     * calculate best and second-best scores for both mates.  These are
     * used for future MAPQ calculations.
     */
	void init(
		const Read* rd1,
		const Read* rd2,
		const EList<AlnRes>* rs)
    {
        assert(rd1 != NULL || rd2 != NULL);
        assert(rs != NULL);
        AlnScore best, secbest;
        size_t szs = 0;
        best.invalidate();    secbest.invalidate();
        szs = rs->size();
        //assert_gt(szs[j], 0);
        for(size_t i = 0; i < rs->size(); i++) {
            AlnScore sc = (*rs)[i].score();
            if(sc > best) {
                secbest = best;
                best = sc;
                assert(VALID_AL_SCORE(best));
            } else if(sc > secbest) {
                secbest = sc;
                assert(VALID_AL_SCORE(best));
                assert(VALID_AL_SCORE(secbest));
            }
        }
        if(szs > 0) {
            init(best, secbest);
        } else {
            reset();
        }
    }

	
	/**
	 * Initialize given fields.  See constructor for how fields are set.
	 */
	void init(
		AlnScore best,
		AlnScore secbest)
	{
		best_         = best;
		secbest_      = secbest;
		assert(repOk());
	}
	
	/**
	 * Return true iff there is at least a best alignment
	 */
	bool empty() const {
		assert(repOk());
		return !VALID_AL_SCORE(best_);
	}
	
#ifndef NDEBUG
	/**
	 * Check that the summary is internally consistent.
	 */
	bool repOk() const {
		return true;
	}
#endif
	
	AlnScore best()         const { return best_;         }
	AlnScore secbest()      const { return secbest_;      }


protected:
	
	AlnScore best_;         // best full-alignment score found for this read
	AlnScore secbest_;      // second-best
};

#endif
