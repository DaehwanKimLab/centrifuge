/*
 * Copyright 2014, Daehwan Kim <infphilo@gmail.com>
 *
 * This file is part of HISAT.
 *
 * HISAT is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * HISAT is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with HISAT.  If not, see <http://www.gnu.org/licenses/>.
 */

#ifndef HI_ALIGNER_H_
#define HI_ALIGNER_H_

#include <iostream>
#include <utility>
#include <limits>
#include "qual.h"
#include "ds.h"
#include "sstring.h"
#include "alphabet.h"
#include "edit.h"
#include "read.h"
// Threading is necessary to synchronize the classes that dump
// intermediate alignment results to files.  Otherwise, all data herein
// is constant and shared, or per-thread.
#include "threading.h"
#include "aligner_result.h"
#include "scoring.h"
#include "mem_ids.h"
#include "simple_func.h"
#include "group_walk.h"

// Minimum intron length
static const uint32_t minIntronLen = 20;
// Maximum intron length
static const uint32_t maxIntronLen = 500000;
// Maximum insertion length
static const uint32_t maxInsLen = 3;
// Maximum deletion length
static const uint32_t maxDelLen = 3;

// Minimum anchor length required for canonical splice sites
static const uint32_t minAnchorLen = 7;
// Minimum anchor length required for non-canonical splice sites
static const uint32_t minAnchorLen_noncan = 14;

// Allow longer introns for long anchored reads involving canonical splice sites
inline uint32_t MaxIntronLen(uint32_t anchor) {
    uint32_t intronLen = 0;
    if(anchor >= minAnchorLen) {
        assert_geq(anchor, 2);
        uint32_t shift = (anchor << 1) - 4;
        shift = min<uint32_t>(max<uint32_t>(shift, 13), 30);
        intronLen = 1 << shift;
    }
    return intronLen;
}

inline float intronLen_prob(uint32_t anchor, uint32_t intronLen) {
    uint32_t expected_intron_len = maxIntronLen;
    if(anchor < 14) expected_intron_len = 1 << ((anchor << 1) + 4);
    if(expected_intron_len > maxIntronLen) expected_intron_len = maxIntronLen;
    assert_gt(expected_intron_len, 0);
    float result = ((float)intronLen) / ((float)expected_intron_len);
    if(result > 1.0f) result = 1.0f;
    return result;
}

// Allow longer introns for long anchored reads involving non-canonical splice sites
inline uint32_t MaxIntronLen_noncan(uint32_t anchor) {
    uint32_t intronLen = 0;
    if(anchor >= minAnchorLen_noncan) {
        assert_geq(anchor, 5);
        uint32_t shift = (anchor << 1) - 10;
        shift = min<uint32_t>(shift, 30);
        intronLen = 1 << shift;
    }
    return intronLen;
}

inline float intronLen_prob_noncan(uint32_t anchor, uint32_t intronLen) {
    uint32_t expected_intron_len = maxIntronLen;
    if(anchor < 16) expected_intron_len = 1 << (anchor << 1);
    if(expected_intron_len > maxIntronLen) expected_intron_len = maxIntronLen;
    assert_gt(expected_intron_len, 0);
    float result = ((float)intronLen) / ((float)expected_intron_len);
    if(result > 1.0f) result = 1.0f;
    return result;
}

/**
 * Hit types for BWTHit class below
 * Three hit types to anchor a read on the genome
 *
 */
enum {
    CANDIDATE_HIT = 1,
    PSEUDOGENE_HIT,
    ANCHOR_HIT,
};

/**
 * Simple struct for holding a partial alignment for the read
 * The alignment locations are represented by FM offsets [top, bot),
 * and later genomic offsets are calculated when necessary
 */
template <typename index_t>
struct BWTHit {
	
	BWTHit() { reset(); }
	
	void reset() {
		_top = _bot = 0;
		_fw = true;
		_bwoff = (index_t)OFF_MASK;
		_len = 0;
		_coords.clear();
        _anchor_examined = false;
        _hit_type = CANDIDATE_HIT;
	}
	
	void init(
			  index_t top,
			  index_t bot,
  			  bool fw,
			  uint32_t bwoff,
			  uint32_t len,
              index_t hit_type = CANDIDATE_HIT)
	{
		_top = top;
        _bot = bot;
		_fw = fw;
		_bwoff = bwoff;
		_len = len;
        _coords.clear();
        _anchor_examined = false;
        _hit_type = hit_type;
	}
    
    bool hasGenomeCoords() const { return !_coords.empty(); }
	
	/**
	 * Return true iff there is no hit.
	 */
	bool empty() const {
		return _bot <= _top;
	}
	
	/**
	 * Higher score = higher priority.
	 */
	bool operator<(const BWTHit& o) const {
		return _len > o._len;
	}
	
	/**
	 * Return the size of the alignments SA ranges.
	 */
	index_t size() const {
        assert_leq(_top, _bot);
        return _bot - _top;
    }
    
    index_t len() const {
        assert_gt(_len, 0);
        return _len;
    }
	
#ifndef NDEBUG
	/**
	 * Check that hit is sane w/r/t read.
	 */
	bool repOk(const Read& rd) const {
		assert_gt(_bot, _top);
		assert_neq(_bwoff, (index_t)OFF_MASK);
		assert_gt(_len, 0);
		return true;
	}
#endif
	
	index_t         _top;               // start of the range in the FM index
	index_t         _bot;               // end of the range in the FM index
	bool            _fw;                // whether read is forward or reverse complemented
	index_t         _bwoff;             // current base of a read to search from the right end
	index_t         _len;               // read length
	
    EList<Coord>    _coords;            // genomic offsets corresponding to [_top, _bot)
    
    bool            _anchor_examined;   // whether or not this hit is examined
    index_t         _hit_type;          // hit type (anchor hit, pseudogene hit, or candidate hit)
};


/**
 * Simple struct for holding alignments for the read
 * The alignments are represented by chains of BWTHits
 */
template <typename index_t>
struct ReadBWTHit {
	
	ReadBWTHit() { reset(); }
	
	void reset() {
        _fw = true;
		_len = 0;
        _cur = 0;
        _done = false;
        _numPartialSearch = 0;
        _numUniqueSearch = 0;
        _partialHits.clear();
	}

	void init(
			  bool fw,
              index_t len)
	{
        _fw = fw;
        assert_gt(len, 0);
        _len = len;
        _cur = 0;
        _done = false;
        _numPartialSearch = 0;
        _numUniqueSearch = 0;
        _partialHits.clear();
	}
    
    bool done() {
#ifndef NDEBUG
        assert_gt(_len, 0);
        if(_cur >= _len) {
            assert(_done);
        }
#endif
        return _done;
    }
    
    void done(bool done) {
        assert(!_done);
        assert(done);
        _done = done;
    }
    
    index_t len() const { return _len; }
    index_t cur() const { return _cur; }
    
    size_t  offsetSize()             { return _partialHits.size(); }
    size_t  numPartialSearch()       { return _numPartialSearch; }
    size_t  numActualPartialSearch()
    {
        assert_leq(_numUniqueSearch, _numPartialSearch);
        return _numPartialSearch - _numUniqueSearch;
    }
    
    bool width(index_t offset_) {
        assert_lt(offset_, _partialHits.size());
        return _partialHits[offset_].size();
    }
    
    bool hasGenomeCoords(index_t offset_) {
        assert_lt(offset_, _partialHits.size());
        index_t width_ = width(offset_);
        if(width_ == 0) {
            return true;
        } else {
            return _partialHits[offset_].hasGenomeCoords();
        }
    }
    
    bool hasAllGenomeCoords() {
        if(_cur < _len) return false;
        if(_partialHits.size() <= 0) return false;
        for(size_t oi = 0; oi < _partialHits.size(); oi++) {
            if(!_partialHits[oi].hasGenomeCoords())
                return false;
        }
        return true;
    }
    
    /**
     *
     */
    index_t minWidth(index_t& offset) const {
        index_t minWidth_ = (index_t)OFF_MASK;
        index_t minWidthLen_ = 0;
        for(size_t oi = 0; oi < _partialHits.size(); oi++) {
            const BWTHit<index_t>& hit = _partialHits[oi];
            if(hit.empty()) continue;
            // if(!hit.hasGenomeCoords()) continue;
            assert_gt(hit.size(), 0);
            if((minWidth_ > hit.size()) ||
               (minWidth_ == hit.size() && minWidthLen_ < hit.len())) {
                minWidth_ = hit.size();
                minWidthLen_ = hit.len();
                offset = (index_t)oi;
            }
        }
        return minWidth_;
    }
    
    // add policy for calculating a search score
    int64_t searchScore(index_t minK) {
        int64_t score = 0;
        const int64_t penaltyPerOffset = minK * minK;
        for(size_t i = 0; i < _partialHits.size(); i++) {
            index_t len = _partialHits[i]._len;
            score += (len * len);
        }
        
        assert_geq(_numPartialSearch, _partialHits.size());
        index_t actualPartialSearch = numActualPartialSearch();
        score -= (actualPartialSearch * penaltyPerOffset);
        score -= (1 << (actualPartialSearch << 1));
        return score;
    }
    
    BWTHit<index_t>& getPartialHit(index_t offset_) {
        assert_lt(offset_, _partialHits.size());
        return _partialHits[offset_];
    }
    
    bool adjustOffset(index_t minK) {
        assert_gt(_partialHits.size(), 0);
        const BWTHit<index_t>& hit = _partialHits.back();
        if(hit.len() >= minK + 3) {
            return false;
        }
        assert_geq(_cur, hit.len());
        index_t origCur = _cur - hit.len();
        _cur = origCur + max(hit.len(), minK + 1) - minK;
        _partialHits.pop_back();
        return true;
    }
    
    void setOffset(index_t offset) {
        assert_lt(offset, _len);
        _cur = offset;
    }
    
#ifndef NDEBUG
	/**
	 */
	bool repOk() const {
        for(size_t i = 0; i < _partialHits.size(); i++) {
            if(i == 0) {
                assert_geq(_partialHits[i]._bwoff, 0);
            }
            
            if(i + 1 < _partialHits.size()) {
                assert_leq(_partialHits[i]._bwoff + _partialHits[i]._len, _partialHits[i+1]._bwoff);
            } else {
                assert_eq(i+1, _partialHits.size());
                assert_eq(_partialHits[i]._bwoff + _partialHits[i]._len, _cur);
            }
        }
		return true;
	}
#endif
	
	bool     _fw;
	index_t  _len;
    index_t  _cur;
    bool     _done;
    index_t  _numPartialSearch;
    index_t  _numUniqueSearch;
    index_t  _cur_local;
    
    EList<BWTHit<index_t> >       _partialHits;
};


/**
 * this is per-thread data, which are shared by GenomeHit classes
 * the main purpose of this struct is to avoid extensive use of memory related functions
 * such as new and delete - those are really slow and lock based
 */
template <typename index_t>
struct SharedTempVars {
    SStringExpandable<char> raw_refbuf;
    SStringExpandable<char> raw_refbuf2;
    EList<int64_t> temp_scores;
    EList<int64_t> temp_scores2;
    ASSERT_ONLY(SStringExpandable<uint32_t> destU32);
    
    ASSERT_ONLY(BTDnaString editstr);
    ASSERT_ONLY(BTDnaString partialseq);
    ASSERT_ONLY(BTDnaString refstr);
    ASSERT_ONLY(EList<index_t> reflens);
    ASSERT_ONLY(EList<index_t> refoffs);
    
    LinkedEList<EList<Edit> > raw_edits;
};

/**
 * GenomeHit represents read alignment or alignment of a part of a read
 * Two GenomeHits that represents alignments of different parts of a read
 * can be combined together.  Also, GenomeHit can be extended in both directions.
 */
template <typename index_t>
struct GenomeHit {
	
	GenomeHit() :
    _fw(false),
    _rdoff((index_t)OFF_MASK),
    _len((index_t)OFF_MASK),
    _trim5(0),
    _trim3(0),
    _tidx((index_t)OFF_MASK),
    _toff((index_t)OFF_MASK),
    _edits(NULL),
    _score(MIN_I64),
    _hitcount(1),
    _edits_node(NULL),
    _sharedVars(NULL)
    {
    }
    
    GenomeHit(const GenomeHit& otherHit) :
    _fw(false),
    _rdoff((index_t)OFF_MASK),
    _len((index_t)OFF_MASK),
    _trim5(0),
    _trim3(0),
    _tidx((index_t)OFF_MASK),
    _toff((index_t)OFF_MASK),
    _edits(NULL),
    _score(MIN_I64),
    _hitcount(1),
    _edits_node(NULL),
    _sharedVars(NULL)
    {
        init(otherHit._fw,
             otherHit._rdoff,
             otherHit._len,
             otherHit._trim5,
             otherHit._trim3,
             otherHit._tidx,
             otherHit._toff,
             *(otherHit._sharedVars),
             otherHit._edits,
             otherHit._score,
             otherHit._splicescore);
    }
    
    GenomeHit<index_t>& operator=(const GenomeHit<index_t>& otherHit) {
        if(this == &otherHit) return *this;
        init(otherHit._fw,
             otherHit._rdoff,
             otherHit._len,
             otherHit._trim5,
             otherHit._trim3,
             otherHit._tidx,
             otherHit._toff,
             *(otherHit._sharedVars),
             otherHit._edits,
             otherHit._score,
             otherHit._splicescore);
        
        return *this;
    }
    
    ~GenomeHit() {
        if(_edits_node != NULL) {
            assert(_edits != NULL);
            assert(_sharedVars != NULL);
            _sharedVars->raw_edits.delete_node(_edits_node);
            _edits = NULL;
            _edits_node = NULL;
            _sharedVars = NULL;
        }
    }
	
	void init(
              bool                      fw,
			  index_t                   rdoff,
			  index_t                   len,
              index_t                   trim5,
              index_t                   trim3,
              index_t                   tidx,
              index_t                   toff,
              SharedTempVars<index_t>&  sharedVars,
              EList<Edit>*              edits = NULL,
              int64_t                   score = 0,
              double                    splicescore = 0.0)
	{
		_fw = fw;
		_rdoff = rdoff;
		_len = len;
        _trim5 = trim5;
        _trim3 = trim3;
        _tidx = tidx;
        _toff = toff;
		_score = score;
        _splicescore = splicescore;
        
        assert(_sharedVars == NULL || _sharedVars == &sharedVars);
        _sharedVars = &sharedVars;
        if(_edits == NULL) {
            assert(_edits_node == NULL);
            _edits_node = _sharedVars->raw_edits.new_node();
            assert(_edits_node != NULL);
            _edits = &(_edits_node->payload);
        }
        assert(_edits != NULL);
        _edits->clear();
        
        if(edits != NULL) *_edits = *edits;
        _hitcount = 1;
	}
    
    bool inited() const {
        return _len >= 0 && _len < (index_t)OFF_MASK;
    }
    
    index_t rdoff() const { return _rdoff; }
    index_t len()   const { return _len; }
    index_t trim5() const { return _trim5; }
    index_t trim3() const { return _trim3; }
    
    void trim5(index_t trim5) { _trim5 = trim5; }
    void trim3(index_t trim3) { _trim3 = trim3; }
    
    index_t ref()    const { return _tidx; }
    index_t refoff() const { return _toff; }
    index_t fw()     const { return _fw; }
    
    index_t hitcount() const { return _hitcount; }
    
    /**
     * Leftmost coordinate
     */
    Coord coord() const {
        return Coord(_tidx, _toff, _fw);
    }
    
    const EList<Edit>& edits() const { return *_edits; }
    
    bool operator== (const GenomeHit<index_t>& other) const {
        if(_fw != other._fw ||
           _rdoff != other._rdoff ||
           _len != other._len ||
           _tidx != other._tidx ||
           _toff != other._toff ||
           _trim5 != other._trim5 ||
           _trim3 != other._trim3) {
            return false;
        }
        
        if(_edits->size() != other._edits->size()) return false;
        for(index_t i = 0; i < _edits->size(); i++) {
            if(!((*_edits)[i] == (*other._edits)[i])) return false;
        }
        // daehwan - this may not be true when some splice sites are provided from outside
        // assert_eq(_score, other._score);
        return true;
    }
    
    bool contains(const GenomeHit<index_t>& other) const {
        return (*this) == other;
    }


#ifndef NDEBUG
	/**
	 * Check that hit is sane w/r/t read.
	 */
	bool repOk(const Read& rd, const BitPairReference& ref);
#endif
    
public:
	bool            _fw;
	index_t         _rdoff;
	index_t         _len;
    index_t         _trim5;
    index_t         _trim3;
    
    index_t         _tidx;
    index_t         _toff;
	EList<Edit>*    _edits;
    int64_t         _score;
    double          _splicescore;
    
    index_t         _hitcount;  // for selection purposes
    
    LinkedEListNode<EList<Edit> >*  _edits_node;
    SharedTempVars<index_t>* _sharedVars;
};


#ifndef NDEBUG
/**
 * Check that hit is sane w/r/t read.
 */
template <typename index_t>
bool GenomeHit<index_t>::repOk(const Read& rd, const BitPairReference& ref)
{
    assert(_sharedVars != NULL);
    SStringExpandable<char>& raw_refbuf = _sharedVars->raw_refbuf;
    SStringExpandable<uint32_t>& destU32 = _sharedVars->destU32;
    
    BTDnaString& editstr = _sharedVars->editstr;
    BTDnaString& partialseq = _sharedVars->partialseq;
    BTDnaString& refstr = _sharedVars->refstr;
    EList<index_t>& reflens = _sharedVars->reflens;
    EList<index_t>& refoffs = _sharedVars->refoffs;
    
    editstr.clear(); partialseq.clear(); refstr.clear();
    reflens.clear(); refoffs.clear();
    
    const BTDnaString& seq = _fw ? rd.patFw : rd.patRc;
    partialseq.install(seq.buf() + this->_rdoff, (size_t)this->_len);
    Edit::toRef(partialseq, *_edits, editstr);
    
    index_t refallen = 0;
    int64_t reflen = 0;
    int64_t refoff = this->_toff;
    refoffs.push_back(refoff);
    size_t eidx = 0;
    for(size_t i = 0; i < _len; i++, reflen++, refoff++) {
        while(eidx < _edits->size() && (*_edits)[eidx].pos == i) {
            const Edit& edit = (*_edits)[eidx];
            if(edit.isReadGap()) {
                reflen++;
                refoff++;
            } else if(edit.isRefGap()) {
                reflen--;
                refoff--;
            }
            if(edit.isSpliced()) {
                assert_gt(reflen, 0);
                refallen += reflen;
                reflens.push_back((index_t)reflen);
                reflen = 0;
                refoff += edit.splLen;
                assert_gt(refoff, 0);
                refoffs.push_back((index_t)refoff);
            }
            eidx++;
        }
    }
    assert_gt(reflen, 0);
    refallen += (index_t)reflen;
    reflens.push_back(reflen);
    assert_gt(reflens.size(), 0);
    assert_gt(refoffs.size(), 0);
    assert_eq(reflens.size(), refoffs.size());
    refstr.clear();
    for(index_t i = 0; i < reflens.size(); i++) {
        assert_gt(reflens[i], 0);
        if(i > 0) {
            assert_gt(refoffs[i], refoffs[i-1]);
        }
        raw_refbuf.resize(reflens[i] + 16);
        raw_refbuf.clear();
        int off = ref.getStretch(
                                 reinterpret_cast<uint32_t*>(raw_refbuf.wbuf()),
                                 (size_t)this->_tidx,
                                 (size_t)max<TRefOff>(refoffs[i], 0),
                                 reflens[i],
                                 destU32);
        assert_leq(off, 16);
        for(index_t j = 0; j < reflens[i]; j++) {
            char rfc = *(raw_refbuf.buf()+off+j);
            refstr.append(rfc);
        }
    }
    if(refstr != editstr) {
        cerr << "Decoded nucleotides and edits don't match reference:" << endl;
        //cerr << "           score: " << score.score()
        //<< " (" << gaps << " gaps)" << endl;
        cerr << "           edits: ";
        Edit::print(cerr, *_edits);
        cerr << endl;
        cerr << "    decoded nucs: " << partialseq << endl;
        cerr << "     edited nucs: " << editstr << endl;
        cerr << "  reference nucs: " << refstr << endl;
        assert(0);
    }

    return true;
}
#endif


/**
 * Encapsulates counters that measure how much work has been done by
 * hierarchical indexing
 */
struct HIMetrics {
    
	HIMetrics() : mutex_m() {
	    reset();
	}
    
	void reset() {
		anchoratts = 0;
        localatts = 0;
        localindexatts = 0;
        localextatts = 0;
        localsearchrecur = 0;
        globalgenomecoords = 0;
        localgenomecoords = 0;
	}
	
	void init(
              uint64_t localatts_,
              uint64_t anchoratts_,
              uint64_t localindexatts_,
              uint64_t localextatts_,
              uint64_t localsearchrecur_,
              uint64_t globalgenomecoords_,
              uint64_t localgenomecoords_)
	{
        localatts = localatts_;
        anchoratts = anchoratts_;
        localindexatts = localindexatts_;
        localextatts = localextatts_;
        localsearchrecur = localsearchrecur_;
        globalgenomecoords = globalgenomecoords_;
        localgenomecoords = localgenomecoords_;
    }
	
	/**
	 * Merge (add) the counters in the given HIMetrics object into this
	 * object.  This is the only safe way to update a HIMetrics shared
	 * by multiple threads.
	 */
	void merge(const HIMetrics& r, bool getLock = false) {
        ThreadSafe ts(&mutex_m, getLock);
        localatts += r.localatts;
        anchoratts += r.anchoratts;
        localindexatts += r.localindexatts;
        localextatts += r.localextatts;
        localsearchrecur += r.localsearchrecur;
        globalgenomecoords += r.globalgenomecoords;
        localgenomecoords += r.localgenomecoords;
    }
	   
    uint64_t localatts;      // # attempts of local search
    uint64_t anchoratts;     // # attempts of anchor search
    uint64_t localindexatts; // # attempts of local index search
    uint64_t localextatts;   // # attempts of extension search
    uint64_t localsearchrecur;
    uint64_t globalgenomecoords;
    uint64_t localgenomecoords;
	
	MUTEX_T mutex_m;
};

/**
 * With a hierarchical indexing, SplicedAligner provides several alignment strategies
 * , which enable effective alignment of RNA-seq reads
 */
template <typename index_t, typename local_index_t>
class HI_Aligner {

public:
	
	/**
	 * Initialize with index.
	 */
	HI_Aligner(
               const Ebwt<index_t>& ebwt,
               bool secondary = false,
               bool local = false,
               uint64_t threads_rids_mindist = 0,
               bool no_spliced_alignment = false) :
    _secondary(secondary),
    _local(local),
    _gwstate(GW_CAT),
    _gwstate_local(GW_CAT),
    _thread_rids_mindist(threads_rids_mindist),
    _no_spliced_alignment(no_spliced_alignment)
    {
        index_t genomeLen = ebwt.eh().len();
        _minK = 0;
        while(genomeLen > 0) {
            genomeLen >>= 2;
            _minK++;
        }
        _minK_local = 8;
    }
    
    HI_Aligner() {
    }
    
    /**
     */
    void initRead(Read *rd, bool nofw, bool norc, TAlScore minsc, TAlScore maxpen, bool rightendonly = false) {
        assert(rd != NULL);
        _rds[0] = rd;
        _rds[1] = NULL;
		_paired = false;
        _rightendonly = rightendonly;
        _nofw[0] = nofw;
        _nofw[1] = true;
        _norc[0] = norc;
        _norc[1] = true;
        _minsc[0] = minsc;
        _minsc[1] = OFF_MASK;
        _maxpen[0] = maxpen;
        _maxpen[1] = OFF_MASK;
        for(size_t fwi = 0; fwi < 2; fwi++) {
            bool fw = (fwi == 0);
            _hits[0][fwi].init(fw, _rds[0]->length());
        }
        _genomeHits.clear();
        _concordantPairs.clear();
        _hits_searched[0].clear();
        assert(!_paired);
    }
    
    /**
     */
    void initReads(Read *rds[2], bool nofw[2], bool norc[2], TAlScore minsc[2], TAlScore maxpen[2]) {
        assert(rds[0] != NULL && rds[1] != NULL);
		_paired = true;
        _rightendonly = false;
        for(size_t rdi = 0; rdi < 2; rdi++) {
            _rds[rdi] = rds[rdi];
            _nofw[rdi] = nofw[rdi];
            _norc[rdi] = norc[rdi];
            _minsc[rdi] = minsc[rdi];
            _maxpen[rdi] = maxpen[rdi];
            for(size_t fwi = 0; fwi < 2; fwi++) {
                bool fw = (fwi == 0);
		        _hits[rdi][fwi].init(fw, _rds[rdi]->length());
            }
            _hits_searched[rdi].clear();
        }
        _genomeHits.clear();
        _concordantPairs.clear();
        assert(_paired);
        assert(!_rightendonly);
    }
    
    /**
     * Aligns a read or a pair
     * This funcion is called per read or pair
     */
    virtual
    int go(
           const Scoring&           sc,
           const Ebwt<index_t>&     ebwtFw,
           const Ebwt<index_t>&     ebwtBw,
           const BitPairReference&  ref,
           WalkMetrics&             wlm,
           PerReadMetrics&          prm,
           HIMetrics&               him,
           RandomSource&            rnd,
           AlnSinkWrap<index_t>&    sink)
    {
        index_t rdi;
        bool fw;
        bool found[2] = {true, this->_paired};
        // given read and its reverse complement
        //  (and mate and the reverse complement of mate in case of pair alignment),
        // pick up one with best partial alignment
        while(nextBWT(sc, ebwtFw, ebwtBw, ref, rdi, fw, wlm, prm, him, rnd, sink)) {
            // given the partial alignment, try to extend it to full alignments
        	found[rdi] = align(sc, ebwtFw, ebwtBw, ref, rdi, fw, wlm, prm, him, rnd, sink);
            if(!found[0] && !found[1]) {
                break;
            }
            
            // try to combine this alignment with some of mate alignments
            // to produce pair alignment
            if(this->_paired) {
                pairReads(sc, ebwtFw, ebwtBw, ref, wlm, prm, him, rnd, sink);
                // if(sink.bestPair() >= _minsc[0] + _minsc[1]) break;
            }
        }
        
        // if no concordant pair is found, try to use alignment of one-end
        // as an anchor to align the other-end
        if(this->_paired) {
            if(_concordantPairs.size() == 0 &&
               (sink.bestUnp1() >= _minsc[0] || sink.bestUnp2() >= _minsc[1])) {
                bool mate_found = false;
                const EList<AlnRes> *rs[2] = {NULL, NULL};
                sink.getUnp1(rs[0]); assert(rs[0] != NULL);
                sink.getUnp2(rs[1]); assert(rs[1] != NULL);
                index_t rs_size[2] = {rs[0]->size(), rs[1]->size()};
                for(index_t i = 0; i < 2; i++) {
                    for(index_t j = 0; j < rs_size[i]; j++) {
                        const AlnRes& res = (*rs[i])[j];
                        bool fw = (res.orient() == 1);
                        mate_found |= alignMate(
                                                sc,
                                                ebwtFw,
                                                ebwtBw,
                                                ref,
                                                i,
                                                fw,
                                                wlm,
                                                prm,
                                                him,
                                                rnd,
                                                sink,
                                                res.refid(),
                                                res.refoff());
                    }
                }
                
                if(mate_found) {
                    pairReads(sc, ebwtFw, ebwtBw, ref, wlm, prm, him, rnd, sink);
                }
            }
        }
        
        return 0;
    }
    
    /**
     * Given a read or its reverse complement (or mate),
     * align the unmapped portion using the global FM index
     */
    virtual
    bool nextBWT(
                 const Scoring&          sc,
                 const Ebwt<index_t>&    ebwtFw,
                 const Ebwt<index_t>&    ebwtBw,
                 const BitPairReference& ref,
                 index_t&                rdi,
                 bool&                   fw,
                 WalkMetrics&            wlm,
                 PerReadMetrics&         prm,
                 HIMetrics&              him,
                 RandomSource&           rnd,
                 AlnSinkWrap<index_t>&   sink)
    {
        // pick up a candidate from a read or its reverse complement
        // (for pair, also consider mate and its reverse complement)
        while(pickNextReadToSearch(rdi, fw)) {
            size_t mineFw = 0, mineRc = 0;
            index_t fwi = (fw ? 0 : 1);
            ReadBWTHit<index_t>& hit = _hits[rdi][fwi];
            assert(!hit.done());
            bool pseudogeneStop = true, anchorStop = true;
            if(!_secondary) {
                index_t numSearched = hit.numActualPartialSearch();
                int64_t bestScore = 0;
                if(rdi == 0) {
                    bestScore = sink.bestUnp1();
                    if(bestScore >= _minsc[rdi]) {
                        // do not further align this candidate
                        // unless it may be at least as good as the alignment of its reverse complement
                        index_t maxmm = (-bestScore + sc.mmpMax - 1) / sc.mmpMax;
                        if(numSearched > maxmm + sink.bestSplicedUnp1() + 1) {
                            hit.done(true);
                            if(_paired) {
                                if(sink.bestUnp2() >= _minsc[1-rdi] &&
                                   _concordantPairs.size() > 0) return false;
                                else continue;
                            } else {
                                return false;
                            }
                        }
                    }
                } else {
                    assert(_paired);
                    assert_eq(rdi, 1);
                    bestScore = sink.bestUnp2();
                    if(bestScore >= _minsc[rdi]) {
                        // do not further extend this alignment
                        // unless it may be at least as good as the previous alignemnt
                        index_t maxmm = (-bestScore + sc.mmpMax - 1) / sc.mmpMax;
                        if(numSearched > maxmm + sink.bestSplicedUnp2() + 1) {
                            hit.done(true);
                            if(_paired) {
                                if(sink.bestUnp1() >= _minsc[1-rdi] &&
                                   _concordantPairs.size() > 0) return false;
                                else continue;
                            } else {
                                return false;
                            }
                        }
                    }
                }
                
                ReadBWTHit<index_t>& rchit = _hits[rdi][1-fwi];
                if(rchit.done() && bestScore < _minsc[rdi]) {
                    if(numSearched > rchit.numActualPartialSearch() + (anchorStop ? 1 : 0)) {
                        hit.done(true);
                        return false;
                    }
                }
            }

            // align this read beginning from previously stopped base
            // stops when it is uniquelly mapped with at least 28bp or
            // it may involve processed pseudogene
            partialSearch(
                          ebwtFw,
                          *_rds[rdi],
                          sc,
                          fw,
                          0,
                          mineFw,
                          mineRc,
                          hit,
                          rnd,
                          pseudogeneStop,
                          anchorStop);

            assert(hit.repOk());
            if(hit.done()) return true;
            // advance hit._cur by 1
            if(!pseudogeneStop) {
                if(hit._cur + 1 < hit._len) hit._cur++;
            }
            if(anchorStop) {
                hit.done(true);
                return true;
            }
            // hit.adjustOffset(_minK);
        }
        
        return false;
    }
    
    /**
     * Given partial alignments of a read, try to further extend
     * the alignment bidirectionally
     */
    virtual
    bool align(
               const Scoring&                   sc,
               const Ebwt<index_t>&             ebwtFw,
               const Ebwt<index_t>&             ebwtBw,
               const BitPairReference&          ref,
               index_t                          rdi,
               bool                             fw,
               WalkMetrics&                     wlm,
               PerReadMetrics&                  prm,
               HIMetrics&                       him,
               RandomSource&                    rnd,
               AlnSinkWrap<index_t>&            sink);
    
    /**
     * Given the alignment of its mate as an anchor,
     * align the read
     */
    virtual
    bool alignMate(
                   const Scoring&                   sc,
                   const Ebwt<index_t>&             ebwtFw,
                   const Ebwt<index_t>&             ebwtBw,
                   const BitPairReference&          ref,
                   index_t                          rdi,
                   bool                             fw,
                   WalkMetrics&                     wlm,
                   PerReadMetrics&                  prm,
                   HIMetrics&                       him,
                   RandomSource&                    rnd,
                   AlnSinkWrap<index_t>&            sink,
                   index_t                          tidx,
                   index_t                          toff);
    
    /**
     * Given a partial alignment of a read, try to further extend
     * the alignment bidirectionally using a combination of
     * local search, extension, and global search
     */
    virtual
    void hybridSearch(
                      const Scoring&                     sc,
                      const Ebwt<index_t>&               ebwtFw,
                      const Ebwt<index_t>&               ebwtBw,
                      const BitPairReference&            ref,
                      index_t                            rdi,
                      bool                               fw,
                      WalkMetrics&                       wlm,
                      PerReadMetrics&                    prm,
                      HIMetrics&                         him,
                      RandomSource&                      rnd,
                      AlnSinkWrap<index_t>&              sink)
    {}
    
    /**
     * Given a partial alignment of a read, try to further extend
     * the alignment bidirectionally using a combination of
     * local search, extension, and global search
     */
    virtual
    int64_t hybridSearch_recur(
                               const Scoring&                   sc,
                               const Ebwt<index_t>&             ebwtFw,
                               const Ebwt<index_t>&             ebwtBw,
                               const BitPairReference&          ref,
                               index_t                          rdi,
                               const GenomeHit<index_t>&        hit,
                               index_t                          hitoff,
                               index_t                          hitlen,
                               WalkMetrics&                     wlm,
                               PerReadMetrics&                  prm,
                               HIMetrics&                       him,
                               RandomSource&                    rnd,
                               AlnSinkWrap<index_t>&            sink,
                               index_t                          dep = 0)
    { return numeric_limits<int64_t>::min(); }
    
    /**
     * Choose a candidate for alignment from a read or its reverse complement
     * (also from a mate or its reverse complement for pair)
     */
    bool pickNextReadToSearch(index_t& rdi, bool& fw) {
        rdi = 0; fw = true;
        bool picked = false;
        int64_t maxScore = std::numeric_limits<int64_t>::min();
        for(index_t rdi2 = 0; rdi2 < (_paired ? 2 : 1); rdi2++) {
            assert(_rds[rdi2] != NULL);
            for(index_t fwi = 0; fwi < 2; fwi++) {
                if     (fwi == 0 && _nofw[rdi2]) continue;
                else if(fwi == 1 && _norc[rdi2]) continue;

                if(_hits[rdi2][fwi].done()) continue;
                int64_t curScore = _hits[rdi2][fwi].searchScore(_minK);
                if(_hits[rdi2][fwi].cur() == 0) {
                    curScore = std::numeric_limits<int64_t>::max();
                }
                assert_gt(curScore, std::numeric_limits<int64_t>::min());
                if(curScore > maxScore) {
                    maxScore = curScore;
                    rdi = rdi2;
                    fw = (fwi == 0);
                    picked = true;
                }
            }
        }
        
        return picked;
    }

	/**
     * Align a part of a read without any edits
	 */
    size_t partialSearch(
                         const Ebwt<index_t>&    ebwt,    // BWT index
                         const Read&             read,    // read to align
                         const Scoring&          sc,      // scoring scheme
                         bool                    fw,      // don't align forward read
                         size_t                  mineMax, // don't care about edit bounds > this
                         size_t&                 mineFw,  // minimum # edits for forward read
                         size_t&                 mineRc,  // minimum # edits for revcomp read
                         ReadBWTHit<index_t>&    hit,     // holds all the seed hits (and exact hit)
                         RandomSource&           rnd,
                         bool&                   pseudogeneStop,  // stop if mapped to multiple locations due to processed pseudogenes
                         bool&                   anchorStop);
    
    /**
     * Global FM index search
	 */
	size_t globalEbwtSearch(
                            const Ebwt<index_t>& ebwt,  // BWT index
                            const Read&          read,  // read to align
                            const Scoring&       sc,    // scoring scheme
                            bool                 fw,
                            index_t              hitoff,
                            index_t&             hitlen,
                            index_t&             top,
                            index_t&             bot,
                            RandomSource&        rnd,
                            bool&                uniqueStop,
                            index_t              maxHitLen = (index_t)OFF_MASK);
    
    /**
     * Local FM index search
	 */
	size_t localEbwtSearch(
                           const LocalEbwt<local_index_t, index_t>*  ebwtFw,  // BWT index
                           const LocalEbwt<local_index_t, index_t>*  ebwtBw,  // BWT index
                           const Read&                      read,    // read to align
                           const Scoring&                   sc,      // scoring scheme
                           bool                             fw,
                           bool                             searchfw,
                           index_t                          rdoff,
                           index_t&                         hitlen,
                           local_index_t&                   top,
                           local_index_t&                   bot,
                           RandomSource&                    rnd,
                           bool&                            uniqueStop,
                           local_index_t                    minUniqueLen,
                           local_index_t                    maxHitLen = (local_index_t)OFF_MASK);
    
    /**
     * Local FM index search
	 */
	size_t localEbwtSearch_reverse(
                                   const LocalEbwt<local_index_t, index_t>*  ebwtFw,  // BWT index
                                   const LocalEbwt<local_index_t, index_t>*  ebwtBw,  // BWT index
                                   const Read&                      read,    // read to align
                                   const Scoring&                   sc,      // scoring scheme
                                   bool                             fw,
                                   bool                             searchfw,
                                   index_t                          rdoff,
                                   index_t&                         hitlen,
                                   local_index_t&                   top,
                                   local_index_t&                   bot,
                                   RandomSource&                    rnd,
                                   bool&                            uniqueStop,
                                   local_index_t                    minUniqueLen,
                                   local_index_t                    maxHitLen = (local_index_t)OFF_MASK);
    
    /**
     * Convert FM offsets to the corresponding genomic offset (chromosome id, offset)
     **/
    bool getGenomeCoords(
                         const Ebwt<index_t>&       ebwt,
                         const BitPairReference&    ref,
                         RandomSource&              rnd,
                         index_t                    top,
                         index_t                    bot,
                         bool                       fw,
                         index_t                    maxelt,
                         index_t                    rdoff,
                         index_t                    rdlen,
                         EList<Coord>&              coords,
                         WalkMetrics&               met,
                         PerReadMetrics&            prm,
                         HIMetrics&                 him,
                         bool                       rejectStraddle,
                         bool&                      straddled);
    
    /**
     * Convert FM offsets to the corresponding genomic offset (chromosome id, offset)
     **/
    bool getGenomeCoords_local(
                               const Ebwt<local_index_t>&   ebwt,
                               const BitPairReference&      ref,
                               RandomSource&                rnd,
                               local_index_t                top,
                               local_index_t                bot,
                               bool                         fw,
                               index_t                      rdoff,
                               index_t                      rdlen,
                               EList<Coord>&                coords,
                               WalkMetrics&                 met,
                               PerReadMetrics&              prm,
                               HIMetrics&                   him,
                               bool                         rejectStraddle,
                               bool&                        straddled);
    
    /**
     * Given a set of partial alignments for a read,
     * choose some that are longer and mapped to fewer places
     */
    index_t getAnchorHits(
                          const Ebwt<index_t>&              ebwt,
                          const BitPairReference&           ref,
                          RandomSource&                     rnd,
                          index_t                           rdi,
                          bool                              fw,
                          EList<GenomeHit<index_t> >&       genomeHits,
                          index_t                           maxGenomeHitSize,
                          SharedTempVars<index_t>&          sharedVars,
                          WalkMetrics&                      wlm,
                          PerReadMetrics&                   prm,
                          HIMetrics&                        him)
    {
        index_t fwi = (fw ? 0 : 1);
        assert_lt(rdi, 2);
        assert(_rds[rdi] != NULL);
        ReadBWTHit<index_t>& hit = _hits[rdi][fwi];
        assert(hit.done());
        index_t offsetSize = hit.offsetSize();
        assert_gt(offsetSize, 0);
        const index_t max_size = (hit._cur >= hit._len ? maxGenomeHitSize : 1);
        genomeHits.clear();
        for(size_t hi = 0; hi < offsetSize; hi++) {
            index_t hj = 0;
            for(; hj < offsetSize; hj++) {
                BWTHit<index_t>& partialHit_j = hit.getPartialHit(hj);
                if(partialHit_j.empty() ||
                   (partialHit_j._hit_type == CANDIDATE_HIT && partialHit_j.size() > max_size) ||
                   partialHit_j.hasGenomeCoords() ||
                   partialHit_j.len() <= _minK + 2) continue;
                else break;
            }
            if(hj >= offsetSize) break;
            for(index_t hk = hj + 1; hk < offsetSize; hk++) {
                BWTHit<index_t>& partialHit_j = hit.getPartialHit(hj);
                BWTHit<index_t>& partialHit_k = hit.getPartialHit(hk);
                if(partialHit_k.empty() ||
                   (partialHit_k._hit_type == CANDIDATE_HIT && partialHit_k.size() > max_size) ||
                   partialHit_k.hasGenomeCoords() ||
                   partialHit_k.len() <= _minK + 2) continue;
                
                if(partialHit_j._hit_type == partialHit_k._hit_type) {
                    if((partialHit_j.size() > partialHit_k.size()) ||
                       (partialHit_j.size() == partialHit_k.size() && partialHit_j.len() < partialHit_k.len())) {
                        hj = hk;
                    }
                } else {
                    if(partialHit_k._hit_type > partialHit_j._hit_type) {
                        hj = hk;
                    }
                }
            }
            BWTHit<index_t>& partialHit = hit.getPartialHit(hj);
            assert(!partialHit.hasGenomeCoords());
            bool straddled = false;
            getGenomeCoords(
                            ebwt,
                            ref,
                            rnd,
                            partialHit._top,
                            partialHit._bot,
                            fw,
                            partialHit._bot - partialHit._top,
                            hit._len - partialHit._bwoff - partialHit._len,
                            partialHit._len,
                            partialHit._coords,
                            wlm,
                            prm,
                            him,
                            false, // reject straddled
                            straddled);
            if(!partialHit.hasGenomeCoords()) continue;
            EList<Coord>& coords = partialHit._coords;
            assert_gt(coords.size(), 0);
            const index_t genomeHit_size = genomeHits.size();
            if(genomeHit_size + coords.size() > maxGenomeHitSize) {
                coords.shufflePortion(0, coords.size(), rnd);
            }
            for(index_t k = 0; k < coords.size(); k++) {
                const Coord& coord = coords[k];
                index_t len = partialHit._len;
                index_t rdoff = hit._len - partialHit._bwoff - len;
                bool overlapped = false;
                for(index_t l = 0; l < genomeHit_size; l++) {
                    GenomeHit<index_t>& genomeHit = genomeHits[l];
                    if(genomeHit.ref() != (index_t)coord.ref() || genomeHit.fw() != coord.fw()) continue;
                    assert_lt(genomeHit.rdoff(), hit._len);
                    assert_lt(rdoff, hit._len);
                    index_t hitoff = genomeHit.refoff() + hit._len - genomeHit.rdoff();
                    index_t hitoff2 = coord.off() + hit._len - rdoff;
                    if(abs((int64_t)hitoff - (int64_t)hitoff2) <= maxIntronLen) {
                        overlapped = true;
                        genomeHit._hitcount++;
                        break;
                    }
                }
                if(!overlapped) {
                    genomeHits.expand();
                    genomeHits.back().init(
                                           coord.orient(),
                                           rdoff,
                                           straddled ? 1 : len,
                                           0, // trim5
                                           0, // trim3
                                           coord.ref(),
                                           coord.off(),
                                           _sharedVars);
                }
                if(partialHit._hit_type == CANDIDATE_HIT && genomeHits.size() >= maxGenomeHitSize) break;
            }
            if(partialHit._hit_type == CANDIDATE_HIT && genomeHits.size() >= maxGenomeHitSize) break;
        }
        return genomeHits.size();
    }
    
    bool pairReads(
                   const Scoring&          sc,
                   const Ebwt<index_t>&    ebwtFw,
                   const Ebwt<index_t>&    ebwtBw,
                   const BitPairReference& ref,
                   WalkMetrics&            wlm,
                   PerReadMetrics&         prm,
                   HIMetrics&              him,
                   RandomSource&           rnd,
                   AlnSinkWrap<index_t>&   sink);

    /**
     *
     **/
    bool reportHit(
                   const Scoring&                   sc,
                   const Ebwt<index_t>&             ebwt,
                   const BitPairReference&          ref,
                   AlnSinkWrap<index_t>&            sink,
                   index_t                          rdi,
                   const GenomeHit<index_t>&        hit,
                   const GenomeHit<index_t>*        ohit = NULL);
    
    /**
     * check this alignment is already examined
     **/
    bool redundant(
                   AlnSinkWrap<index_t>&    sink,
                   index_t                  rdi,
                   index_t                  tidx,
                   index_t                  toff);
    
    /**
     * check this alignment is already examined
     **/
    bool redundant(
                   AlnSinkWrap<index_t>&            sink,
                   index_t                          rdi,
                   const GenomeHit<index_t>&        hit);
    
    
    /**
     *
     **/
    bool isSearched(
                    const GenomeHit<index_t>&       hit,
                    index_t                         rdi);
    
    /**
     *
     **/
    void addSearched(const GenomeHit<index_t>&       hit,
                     index_t                         rdi);
    
    
protected:
  
    Read *   _rds[2];
    bool     _paired;
    bool     _rightendonly;
    bool     _nofw[2];
    bool     _norc[2];
    TAlScore _minsc[2];
    TAlScore _maxpen[2];
    
    bool     _secondary;  // allow secondary alignments
    bool     _local;      // perform local alignments
    
    ReadBWTHit<index_t> _hits[2][2];
    
    EList<index_t, 16>                                 _offs;
    SARangeWithOffs<EListSlice<index_t, 16> >          _sas;
    GroupWalk2S<index_t, EListSlice<index_t, 16>, 16>  _gws;
    GroupWalkState<index_t>                            _gwstate;
    
    EList<local_index_t, 16>                                       _offs_local;
    SARangeWithOffs<EListSlice<local_index_t, 16> >                _sas_local;
    GroupWalk2S<local_index_t, EListSlice<local_index_t, 16>, 16>  _gws_local;
    GroupWalkState<local_index_t>                                  _gwstate_local;
            
    // temporary and shared variables used for GenomeHit
    // this should be defined before _genomeHits and _hits_searched
    SharedTempVars<index_t> _sharedVars;
    
    // temporary and shared variables for AlnRes
    LinkedEList<EList<Edit> > _rawEdits;
    
    // temporary
    EList<GenomeHit<index_t> >     _genomeHits;
    EList<bool>                    _genomeHits_done;
    ELList<Coord>                  _coords;
    
    EList<pair<index_t, index_t> >  _concordantPairs;
    
    size_t _minK; // log4 of the size of a genome
    size_t _minK_local; // log4 of the size of a local index (8)

    ELList<GenomeHit<index_t> >     _local_genomeHits;
    EList<uint8_t>                  _anchors_added;
    uint64_t max_localindexatts;
    
	uint64_t bwops_;                    // Burrows-Wheeler operations
	uint64_t bwedits_;                  // Burrows-Wheeler edits
    
    //
    EList<GenomeHit<index_t> >     _hits_searched[2];
    
    uint64_t   _thread_rids_mindist;
    bool _no_spliced_alignment;

    // For AlnRes::matchesRef
	ASSERT_ONLY(EList<bool> raw_matches_);
	ASSERT_ONLY(BTDnaString tmp_rf_);
	ASSERT_ONLY(BTDnaString tmp_rdseq_);
	ASSERT_ONLY(BTString tmp_qseq_);
};

#define HIER_INIT_LOCS(top, bot, tloc, bloc, e) { \
	if(bot - top == 1) { \
		tloc.initFromRow(top, (e).eh(), (e).ebwt()); \
		bloc.invalidate(); \
	} else { \
		SideLocus<index_t>::initFromTopBot(top, bot, (e).eh(), (e).ebwt(), tloc, bloc); \
		assert(bloc.valid()); \
	} \
}

#define HIER_SANITY_CHECK_4TUP(t, b, tp, bp) { \
	ASSERT_ONLY(cur_index_t tot = (b[0]-t[0])+(b[1]-t[1])+(b[2]-t[2])+(b[3]-t[3])); \
	ASSERT_ONLY(cur_index_t totp = (bp[0]-tp[0])+(bp[1]-tp[1])+(bp[2]-tp[2])+(bp[3]-tp[3])); \
	assert_eq(tot, totp); \
}

#define LOCAL_INIT_LOCS(top, bot, tloc, bloc, e) { \
    if(bot - top == 1) { \
        tloc.initFromRow(top, (e).eh(), (e).ebwt()); \
        bloc.invalidate(); \
    } else { \
        SideLocus<local_index_t>::initFromTopBot(top, bot, (e).eh(), (e).ebwt(), tloc, bloc); \
        assert(bloc.valid()); \
    } \
}

/**
 * Given partial alignments of a read, try to further extend
 * the alignment bidirectionally
 */
template <typename index_t, typename local_index_t>
bool HI_Aligner<index_t, local_index_t>::align(
                                               const Scoring&                   sc,
                                               const Ebwt<index_t>&             ebwtFw,
                                               const Ebwt<index_t>&             ebwtBw,
                                               const BitPairReference&          ref,
                                               index_t                          rdi,
                                               bool                             fw,
                                               WalkMetrics&                     wlm,
                                               PerReadMetrics&                  prm,
                                               HIMetrics&                       him,
                                               RandomSource&                    rnd,
                                               AlnSinkWrap<index_t>&            sink)
{
    const ReportingParams& rp = sink.reportingParams();
    index_t fwi = (fw ? 0 : 1);
    assert_lt(rdi, 2);
    assert(_rds[rdi] != NULL);
    ReadBWTHit<index_t>& hit = _hits[rdi][fwi];
    assert(hit.done());
    index_t minOff = 0;
    if(hit.minWidth(minOff) > (index_t)(rp.khits * 2)) return false;
    
    // do not try to align if the potential alignment for this read might be
    // worse than the best alignment of its reverse complement
    int64_t bestScore = (rdi == 0 ? sink.bestUnp1() : sink.bestUnp2());
    index_t num_spliced = (rdi == 0 ? sink.bestSplicedUnp1() : sink.bestSplicedUnp2());
    if(bestScore < _minsc[rdi]) bestScore = _minsc[rdi];
    index_t maxmm = (-bestScore + sc.mmpMax - 1) / sc.mmpMax;
    index_t numActualPartialSearch = hit.numActualPartialSearch();
    if(!_secondary && numActualPartialSearch > maxmm + num_spliced + 1) return true;
    
    // choose candidate partial alignments for further alignment
    const index_t maxsize = rp.khits;
    index_t numHits = getAnchorHits(
                                    ebwtFw,
                                    ref,
                                    rnd,
                                    rdi,
                                    fw,
                                    _genomeHits,
                                    maxsize,
                                    _sharedVars,
                                    wlm,
                                    prm,
                                    him);
    if(numHits <= 0) return false;
   
    // limit the number of local index searches used for alignment of the read
    uint64_t add = 0;
    if(_secondary) add = (-_minsc[rdi] / sc.mmpMax) * numHits * 2;
    else           add = (-_minsc[rdi] / sc.mmpMax) * numHits;
    max_localindexatts = him.localindexatts + max<uint64_t>(10, add);
    // extend the partial alignments bidirectionally using
    // local search, extension, and (less often) global search
    hybridSearch(
                 sc,
                 ebwtFw,
                 ebwtBw,
                 ref,
                 rdi,
                 fw,
                 wlm,
                 prm,
                 him,
                 rnd,
                 sink);
    
    return true;
}


/**
 * Given the alignment of its mate as an anchor,
 * align the read
 */
template <typename index_t, typename local_index_t>
bool HI_Aligner<index_t, local_index_t>::alignMate(
                                                   const Scoring&                   sc,
                                                   const Ebwt<index_t>&             ebwtFw,
                                                   const Ebwt<index_t>&             ebwtBw,
                                                   const BitPairReference&          ref,
                                                   index_t                          rdi,
                                                   bool                             fw,
                                                   WalkMetrics&                     wlm,
                                                   PerReadMetrics&                  prm,
                                                   HIMetrics&                       him,
                                                   RandomSource&                    rnd,
                                                   AlnSinkWrap<index_t>&            sink,
                                                   index_t                          tidx,
                                                   index_t                          toff)
{
    assert_lt(rdi, 2);
    index_t ordi = 1 - rdi;
    bool ofw = (fw == gMate2fw ? gMate1fw : gMate2fw);
    assert(_rds[ordi] != NULL);
    const Read& ord = *_rds[ordi];
    index_t rdlen = ord.length();
    assert_gt(rdlen, 0);
    
    _genomeHits.clear();
    if(_coords.size() == 0) {
        _coords.expand();
    }
    EList<Coord>& coords = _coords.front();
    
    // local search to find anchors
    const HierEbwt<index_t, local_index_t>* hierEbwt = (const HierEbwt<index_t, local_index_t>*)(&ebwtFw);
    const LocalEbwt<local_index_t, index_t>* localEbwt = hierEbwt->getLocalEbwt(tidx, toff);
    bool success = false, first = true;
    index_t count = 0;
    index_t max_hitlen = 0;
    while(!success && count++ < 2) {
        if(first) {
            first = false;
        } else {
            localEbwt = hierEbwt->prevLocalEbwt(localEbwt);
            if(localEbwt == NULL || localEbwt->empty()) break;
        }
        index_t hitoff = rdlen - 1;
        while(hitoff >= _minK_local - 1) {
            index_t hitlen = 0;
            local_index_t top = (local_index_t)OFF_MASK, bot = (local_index_t)OFF_MASK;
            bool uniqueStop = false;
            index_t nelt = localEbwtSearch(
                                           localEbwt,   // BWT index
                                           NULL,        // BWT index
                                           ord,          // read to align
                                           sc,          // scoring scheme
                                           ofw,
                                           false,       // searchfw,
                                           hitoff,
                                           hitlen,
                                           top,
                                           bot,
                                           rnd,
                                           uniqueStop,
                                           _minK_local);
            assert_leq(top, bot);
            assert_eq(nelt, (index_t)(bot - top));
            assert_leq(hitlen, hitoff + 1);
            if(nelt > 0 && nelt <= 5 && hitlen > max_hitlen) {
                coords.clear();
                bool straddled = false;
                getGenomeCoords_local(
                                      *localEbwt,
                                      ref,
                                      rnd,
                                      top,
                                      bot,
                                      ofw,
                                      hitoff - hitlen + 1,
                                      hitlen,
                                      coords,
                                      wlm,
                                      prm,
                                      him,
                                      true, // reject straddled?
                                      straddled);
                assert_leq(coords.size(), nelt);
                _genomeHits.clear();
                for(index_t ri = 0; ri < coords.size(); ri++) {
                    const Coord& coord = coords[ri];
                    _genomeHits.expand();
                    _genomeHits.back().init(
                                            coord.orient(),
                                            hitoff - hitlen + 1,
                                            hitlen,
                                            0, // trim5
                                            0, // trim3
                                            coord.ref(),
                                            coord.off(),
                                            _sharedVars);
                }
                max_hitlen = hitlen;
            }
            
            assert_leq(hitlen, hitoff + 1);
            hitoff -= (hitlen - 1);
            if(hitoff > 0) hitoff -= 1;
        } // while(hitoff >= _minK_local - 1)
    } // while(!success && count++ < 2)
    
    if(max_hitlen < _minK_local) return false;
    
    // randomly select
    const index_t maxsize = 5;
    if(_genomeHits.size() > maxsize) {
        _genomeHits.shufflePortion(0, _genomeHits.size(), rnd);
        _genomeHits.resize(maxsize);
    }
    
    // local search using the anchor
    for(index_t hi = 0; hi < _genomeHits.size(); hi++) {
        him.anchoratts++;
        GenomeHit<index_t>& genomeHit = _genomeHits[hi];
        hybridSearch_recur(
                           sc,
                           ebwtFw,
                           ebwtBw,
                           ref,
                           ordi,
                           genomeHit,
                           genomeHit.rdoff(),
                           genomeHit.len(),
                           wlm,
                           prm,
                           him,
                           rnd,
                           sink);
    }
    
    return true;
}


/**
 * convert FM offsets to the corresponding genomic offset (chromosome id, offset)
 **/
template <typename index_t, typename local_index_t>
bool HI_Aligner<index_t, local_index_t>::getGenomeCoords(
                                                         const Ebwt<index_t>&       ebwt,
                                                         const BitPairReference&    ref,
                                                         RandomSource&              rnd,
                                                         index_t                    top,
                                                         index_t                    bot,
                                                         bool                       fw,
                                                         index_t                    maxelt,
                                                         index_t                    rdoff,
                                                         index_t                    rdlen,
                                                         EList<Coord>&              coords,
                                                         WalkMetrics&               met,
                                                         PerReadMetrics&            prm,
                                                         HIMetrics&                 him,
                                                         bool                       rejectStraddle,
                                                         bool&                      straddled)
{
    straddled = false;
    assert_gt(bot, top);
    index_t nelt = bot - top;
    nelt = min<index_t>(nelt, maxelt);
    coords.clear();
    him.globalgenomecoords += (bot - top);
    _offs.resize(nelt);
    _offs.fill(std::numeric_limits<index_t>::max());
    _sas.init(top, rdlen, EListSlice<index_t, 16>(_offs, 0, nelt));
    _gws.init(ebwt, ref, _sas, rnd, met);
    
    for(index_t off = 0; off < nelt; off++) {
        WalkResult<index_t> wr;
        index_t tidx = 0, toff = 0, tlen = 0;
        _gws.advanceElement(
                            off,
                            ebwt,         // forward Bowtie index for walking left
                            ref,          // bitpair-encoded reference
                            _sas,         // SA range with offsets
                            _gwstate,     // GroupWalk state; scratch space
                            wr,           // put the result here
                            met,          // metrics
                            prm);         // per-read metrics
        assert_neq(wr.toff, (index_t)OFF_MASK);
        bool straddled2 = false;
        ebwt.joinedToTextOff(
                             wr.elt.len,
                             wr.toff,
                             tidx,
                             toff,
                             tlen,
                             rejectStraddle,        // reject straddlers?
                             straddled2);  // straddled?
        
        straddled |= straddled2;
        
        if(tidx == (index_t)OFF_MASK) {
            // The seed hit straddled a reference boundary so the seed
            // hit isn't valid
            return false;
        }
        index_t global_toff = toff, global_tidx = tidx;
        if(global_toff < rdoff) continue;
        
        // Coordinate of the seed hit w/r/t the pasted reference string
        coords.expand();
        coords.back().init(global_tidx, (int64_t)global_toff, fw);
    }
    
    return true;
}

/**
 * convert FM offsets to the corresponding genomic offset (chromosome id, offset)
 **/
template <typename index_t, typename local_index_t>
bool HI_Aligner<index_t, local_index_t>::getGenomeCoords_local(
                                                               const Ebwt<local_index_t>&   ebwt,
                                                               const BitPairReference&      ref,
                                                               RandomSource&                rnd,
                                                               local_index_t                top,
                                                               local_index_t                bot,
                                                               bool                         fw,
                                                               index_t                      rdoff,
                                                               index_t                      rdlen,
                                                               EList<Coord>&                coords,
                                                               WalkMetrics&                 met,
                                                               PerReadMetrics&              prm,
                                                               HIMetrics&                   him,
                                                               bool                         rejectStraddle,
                                                               bool&                        straddled)
{
    straddled = false;
    assert_gt(bot, top);
    index_t nelt = bot - top;
    coords.clear();
    him.localgenomecoords += (bot - top);
    _offs_local.resize(nelt);
    _offs_local.fill(std::numeric_limits<local_index_t>::max());
    _sas_local.init(top, rdlen, EListSlice<local_index_t, 16>(_offs_local, 0, nelt));
    _gws_local.init(ebwt, ref, _sas_local, rnd, met);
    
    for(local_index_t off = 0; off < nelt; off++) {
        WalkResult<local_index_t> wr;
        local_index_t tidx = 0, toff = 0, tlen = 0;
        _gws_local.advanceElement(
                                  off,
                                  ebwt,         // forward Bowtie index for walking left
                                  ref,          // bitpair-encoded reference
                                  _sas_local,   // SA range with offsets
                                  _gwstate_local, // GroupWalk state; scratch space
                                  wr,           // put the result here
                                  met,          // metrics
                                  prm);         // per-read metrics
        assert_neq(wr.toff, (local_index_t)OFF_MASK);
        bool straddled2 = false;
        ebwt.joinedToTextOff(
                             wr.elt.len,
                             wr.toff,
                             tidx,
                             toff,
                             tlen,
                             rejectStraddle,        // reject straddlers?
                             straddled2);  // straddled?
        
        straddled |= straddled2;
        
        if(tidx == (local_index_t)OFF_MASK) {
            // The seed hit straddled a reference boundary so the seed
            // hit isn't valid
            return false;
        }
        index_t global_toff = toff, global_tidx = tidx;
        LocalEbwt<local_index_t, index_t>* localEbwt = (LocalEbwt<local_index_t, index_t>*)&ebwt;
        global_tidx = localEbwt->_tidx, global_toff = toff + localEbwt->_localOffset;
        if(global_toff < rdoff) continue;
        
        // Coordinate of the seed hit w/r/t the pasted reference string
        coords.expand();
        coords.back().init(global_tidx, (int64_t)global_toff, fw);
    }
    
    return true;
}


/**
 * examine alignments of left and right reads to produce concordant pair alignment
 **/
template <typename index_t, typename local_index_t>
bool HI_Aligner<index_t, local_index_t>::pairReads(
                                                   const Scoring&          sc,
                                                   const Ebwt<index_t>&    ebwtFw,
                                                   const Ebwt<index_t>&    ebwtBw,
                                                   const BitPairReference& ref,
                                                   WalkMetrics&            wlm,
                                                   PerReadMetrics&         prm,
                                                   HIMetrics&              him,
                                                   RandomSource&           rnd,
                                                   AlnSinkWrap<index_t>&   sink)
{
    assert(_paired);
    const EList<AlnRes> *rs1 = NULL, *rs2 = NULL;
    sink.getUnp1(rs1); assert(rs1 != NULL);
    sink.getUnp2(rs2); assert(rs2 != NULL);
    for(index_t i = 0; i < rs1->size(); i++) {
        for(index_t j = 0; j < rs2->size(); j++) {
            bool exists = false;
            for(index_t k = 0; k < _concordantPairs.size(); k++) {
                const pair<index_t, index_t>& p = _concordantPairs[k];
                if(i == p.first && j == p.second) {
                    exists = true;
                    break;
                }
            }
            if(exists) continue;
	    if(sink.state().doneConcordant()) return true;
            const AlnRes& r1 = (*rs1)[i];
            Coord left = r1.refcoord(), right = r1.refcoord_right();
            assert_eq(left.ref(), right.ref());
            const AlnRes& r2 = (*rs2)[j];
            Coord left2 = r2.refcoord(), right2 = r2.refcoord_right();
            assert_eq(left2.ref(), right2.ref());
            if(left.ref() != left2.ref()) continue;
            assert_eq(left.orient(), right.orient());
            assert_eq(left2.orient(), right2.orient());
            if(left.orient() == gMate1fw) {
                if(left2.orient() != gMate2fw) continue;
            } else {
                if(left2.orient() == gMate2fw) continue;
                Coord temp = left; left = left2; left2 = temp;
                temp = right; right = right2; right2 = temp;
            }
            if(left.off() > left2.off()) continue;
            if(right.off() > right2.off()) continue;
            if(right.off() + 500000 < left2.off()) continue;
            assert_geq(r1.score().score(), _minsc[0]);
            assert_geq(r2.score().score(), _minsc[1]);
            if(r1.score().score() + r2.score().score() >= sink.bestPair() || _secondary) {
                sink.report(0, &r1, &r2);
                _concordantPairs.expand();
                _concordantPairs.back().first = i;
                _concordantPairs.back().second = j;
            }
        }
    }
    return true;
}


/**
 * report read (or pair) alignment
 **/
template <typename index_t, typename local_index_t>
bool HI_Aligner<index_t, local_index_t>::reportHit(
                                                   const Scoring&                   sc,
                                                   const Ebwt<index_t>&             ebwt,
                                                   const BitPairReference&          ref,
                                                   AlnSinkWrap<index_t>&            sink,
                                                   index_t                          rdi,
                                                   const GenomeHit<index_t>&        hit,
                                                   const GenomeHit<index_t>*        ohit)
{
    assert_lt(rdi, 2);
    assert(_rds[rdi] != NULL);
    const Read& rd = *_rds[rdi];
    index_t rdlen = rd.length();
    if(hit.rdoff() - hit.trim5() > 0 || hit.len() + hit.trim5() + hit.trim3() < rdlen) return false;
    if(hit.score() < _minsc[rdi]) return false;
    // else if(hit.spliced() && !hit.spliced_consistently())
    
    // Edits are represented from 5' end of read to 3' end, not an alignment of read
    EList<Edit>& edits = const_cast<EList<Edit>&>(hit.edits());
    if(hit.trim5() > 0) {
        for(size_t i = 0; i < edits.size(); i++) {
            edits[i].pos += hit.trim5();
        }
    }
    if(!hit.fw()) {
        Edit::invertPoss(edits, rdlen, false);
    }
    AlnScore asc(
                 hit.score(),  // numeric score
                 hit.ns(),     // # Ns
                 hit.ngaps(),  // # gaps
                 hit.splicescore()); // splice score
    bool softTrim = hit.trim5() > 0 || hit.trim3() > 0;
    AlnRes rs;
    rs.init(
            rdlen,                      // # chars after hard trimming
            asc,                        // alignment score
            &hit.edits(),               // nucleotide edits array
            0,                          // nucleotide edits first pos
            hit.edits().size(),         // nucleotide edits last pos
            NULL,                       // ambig base array
            0,                          // ambig base first pos
            0,                          // ambig base last pos
            hit.coord(),                // coord of leftmost aligned char in ref
            ebwt.plen()[hit.ref()],     // length of reference aligned to
            &_rawEdits,
            -1,                         // # seed mms allowed
            -1,                         // seed length
            -1,                         // seed interval
            0,                          // minimum score for valid alignment (daehwan)
            -1,                         // nuc5p (for colorspace)
            -1,                         // nuc3p (for colorspace)
            false,                      // soft pre-trimming?
            0,                          // 5p pre-trimming
            0,                          // 3p pre-trimming
            softTrim,                   // soft trimming?
            hit.fw() ? hit.trim5() : hit.trim3(),  // 5p trimming
            hit.fw() ? hit.trim3() : hit.trim5()); // 3p trimming
    if(!hit.fw()) {
        Edit::invertPoss(edits, rdlen, false);
    }
    if(hit.trim5() > 0) {
        for(size_t i = 0; i < edits.size(); i++) {
            edits[i].pos -= hit.trim5();
        }
    }
    //rs.setRefNs(nrefn);
    assert(rs.matchesRef(
                         rd,
                         ref,
                         tmp_rf_,
                         tmp_rdseq_,
                         tmp_qseq_,
                         _sharedVars.raw_refbuf,
                         _sharedVars.destU32,
                         raw_matches_,
                         _sharedVars.raw_refbuf2,
                         _sharedVars.reflens,
                         _sharedVars.refoffs));
    if(ohit == NULL) {
        bool done;
        if(rdi == 0 && !_rightendonly) {
            done = sink.report(0, &rs, NULL);
        } else {
            done = sink.report(0, NULL, &rs);
        }
        return done;
    }
    
    assert(ohit != NULL);
    const Read& ord = *_rds[1-rdi];
    index_t ordlen = ord.length();
    if(ohit->rdoff() - ohit->trim5() > 0 || ohit->len() + ohit->trim5() + ohit->trim3() < ordlen) return false;
    if(ohit->score() < _minsc[1-rdi]) return false;
    EList<Edit>& oedits = const_cast<EList<Edit>&>(ohit->edits());
    if(ohit->trim5() > 0) {
        for(size_t i = 0; i < oedits.size(); i++) {
            oedits[i].pos += ohit->trim5();
        }
    }
    if(!ohit->fw()) {
        Edit::invertPoss(oedits, ordlen, false);
    }
    AlnScore oasc(
                  ohit->score(),  // numeric score
                  ohit->ns(),     // # Ns
                  ohit->ngaps()); // # gaps
    bool osoftTrim = ohit->trim5() > 0 || ohit->trim3() > 0;
    AlnRes ors;
    ors.init(
             ordlen,                     // # chars after hard trimming
             oasc,                       // alignment score
             &ohit->edits(),             // nucleotide edits array
             0,                          // nucleotide edits first pos
             ohit->edits().size(),       // nucleotide edits last pos
             NULL,                       // ambig base array
             0,                          // ambig base first pos
             0,                          // ambig base last pos
             ohit->coord(),              // coord of leftmost aligned char in ref
             ebwt.plen()[ohit->ref()],   // length of reference aligned to
             &_rawEdits,
             -1,                         // # seed mms allowed
             -1,                         // seed length
             -1,                         // seed interval
             0,                          // minimum score for valid alignment (daehwan)
             -1,                         // nuc5p (for colorspace)
             -1,                         // nuc3p (for colorspace)
             false,                      // soft pre-trimming?
             0,                          // 5p pre-trimming
             0,                          // 3p pre-trimming
             osoftTrim,                  // soft trimming?
             ohit->fw() ? ohit->trim5() : ohit->trim3(),  // 5p trimming
             ohit->fw() ? ohit->trim3() : ohit->trim5()); // 3p trimming
    if(!ohit->fw()) {
        Edit::invertPoss(oedits, ordlen, false);
    }
    if(ohit->trim5() > 0) {
        for(size_t i = 0; i < oedits.size(); i++) {
            oedits[i].pos -= ohit->trim5();
        }
    }
    //rs.setRefNs(nrefn);
    assert(ors.matchesRef(
                          ord,
                          ref,
                          tmp_rf_,
                          tmp_rdseq_,
                          tmp_qseq_,
                          _sharedVars.raw_refbuf,
                          _sharedVars.destU32,
                          raw_matches_,
                          _sharedVars.raw_refbuf2,
                          _sharedVars.reflens,
                          _sharedVars.refoffs));
    
    bool done;
    if(rdi == 0) {
        done = sink.report(0, &rs, &ors);
    } else {
        done = sink.report(0, &ors, &rs);
    }
    return done;
}

/**
 * check this alignment is already examined
 **/
template <typename index_t, typename local_index_t>
bool HI_Aligner<index_t, local_index_t>::redundant(
                                                   AlnSinkWrap<index_t>&    sink,
                                                   index_t                  rdi,
                                                   index_t                  tidx,
                                                   index_t                  toff)
{
    assert_lt(rdi, 2);
    const EList<AlnRes>* rs = NULL;
    if(rdi == 0) sink.getUnp1(rs);
    else         sink.getUnp2(rs);
    assert(rs != NULL);
    for(index_t i = 0; i < rs->size(); i++) {
        Coord coord_left = (*rs)[i].refcoord(), coord_right = (*rs)[i].refcoord_right();
        assert_eq(coord_left.ref(), coord_right.ref());
        assert_lt(coord_left.off(), coord_right.off());
        assert_eq(coord_left.orient(), coord_right.orient());
        
        if(tidx != coord_left.ref()) continue;
        if(toff >= coord_left.off() && toff <= coord_right.off()) return true;
    }
    
    return false;
}


/**
 * check this alignment is already examined
 **/
template <typename index_t, typename local_index_t>
bool HI_Aligner<index_t, local_index_t>::redundant(
                                                   AlnSinkWrap<index_t>&            sink,
                                                   index_t                          rdi,
                                                   const GenomeHit<index_t>&        hit)
{
    assert_lt(rdi, 2);
    assert(_rds[rdi] != NULL);
    index_t rdlen = _rds[rdi]->length();
    const EList<AlnRes>* rs = NULL;
    if(rdi == 0) sink.getUnp1(rs);
    else         sink.getUnp2(rs);
    assert(rs != NULL);
    for(index_t i = 0; i < rs->size(); i++) {
        const AlnRes& rsi = (*rs)[i];
        if(rsi.refcoord() == hit.coord()) {
            const EList<Edit>& editsi = rsi.ned();
            const EList<Edit>& edits = hit.edits();
            if(editsi.size() == edits.size()) {
                size_t eidx = 0;
                if(!hit.fw()) {
                    Edit::invertPoss(const_cast<EList<Edit>&>(edits), rdlen, false);
                }
                for(; eidx < editsi.size(); eidx++) {
                    if(!(editsi[eidx] == edits[eidx])) {
                        break;
                    }
                }
                if(!hit.fw()) {
                    Edit::invertPoss(const_cast<EList<Edit>&>(edits), rdlen, false);
                }
                if(eidx >= editsi.size()) {
                    assert_eq(eidx, editsi.size());
                    return true;
                }
            }
        }
    }
    
    return false;
}


/**
 * Sweep right-to-left and left-to-right using exact matching.  Remember all
 * the SA ranges encountered along the way.  Report exact matches if there are
 * any.  Calculate a lower bound on the number of edits in an end-to-end
 * alignment.
 */
template <typename index_t, typename local_index_t>
size_t HI_Aligner<index_t, local_index_t>::partialSearch(
                                                         const Ebwt<index_t>&      ebwt,    // BWT index
                                                         const Read&               read,    // read to align
                                                         const Scoring&            sc,      // scoring scheme
                                                         bool                      fw,
                                                         size_t                    mineMax, // don't care about edit bounds > this
                                                         size_t&                   mineFw,  // minimum # edits for forward read
                                                         size_t&                   mineRc,  // minimum # edits for revcomp read
                                                         ReadBWTHit<index_t>&      hit,     // holds all the seed hits (and exact hit)
                                                         RandomSource&             rnd,     // pseudo-random source
                                                         bool&                     pseudogeneStop,
                                                         bool&                     anchorStop)
{
    bool pseudogeneStop_ = pseudogeneStop, anchorStop_ = anchorStop;
    pseudogeneStop = anchorStop = false;
	const index_t ftabLen = ebwt.eh().ftabChars();
	SideLocus<index_t> tloc, bloc;
	const index_t len = (index_t)read.length();
    const BTDnaString& seq = fw ? read.patFw : read.patRc;
    assert(!seq.empty());
    
    size_t nelt = 0;
    EList<BWTHit<index_t> >& partialHits = hit._partialHits;
    index_t& cur = hit._cur;
    assert_lt(cur, hit._len);
    
    hit._numPartialSearch++;
    
    index_t offset = cur;
    index_t dep = offset;
    index_t top = 0, bot = 0;
    index_t topTemp = 0, botTemp = 0;
    index_t left = len - dep;
    assert_gt(left, 0);
    if(left < ftabLen) {
        cur = hit._len;
        partialHits.expand();
        partialHits.back().init((index_t)OFF_MASK,
                                (index_t)OFF_MASK,
                                fw,
                                (index_t)offset,
                                (index_t)(cur - offset));
        hit.done(true);
		return 0;
    }
    // Does N interfere with use of Ftab?
    for(index_t i = 0; i < ftabLen; i++) {
        int c = seq[len-dep-1-i];
        if(c > 3) {
            cur += (i+1);
            partialHits.expand();
            partialHits.back().init((index_t)OFF_MASK,
                                    (index_t)OFF_MASK,
                                    fw,
                                    (index_t)offset,
                                    (index_t)(cur - offset));
            if(cur >= hit._len) {
                hit.done(true);
            }
			return 0;
        }
    }
    
    // Use ftab
    ebwt.ftabLoHi(seq, len - dep - ftabLen, false, top, bot);
    dep += ftabLen;
    if(bot <= top) {
        cur = dep;
        partialHits.expand();
        partialHits.back().init((index_t)OFF_MASK,
                                (index_t)OFF_MASK,
                                fw,
                                (index_t)offset,
                                (index_t)(cur - offset));
        if(cur >= hit._len) {
            hit.done(true);
        }
        return 0;
    }
    index_t same_range = 0, similar_range = 0;
    HIER_INIT_LOCS(top, bot, tloc, bloc, ebwt);
    // Keep going
    while(dep < len) {
        int c = seq[len-dep-1];
        if(c > 3) {
            topTemp = botTemp = 0;
        } else {
            if(bloc.valid()) {
                bwops_ += 2;
                topTemp = ebwt.mapLF(tloc, c);
                botTemp = ebwt.mapLF(bloc, c);
            } else {
                bwops_++;
                topTemp = ebwt.mapLF1(top, tloc, c);
                if(topTemp == (index_t)OFF_MASK) {
                    topTemp = botTemp = 0;
                } else {
                    botTemp = topTemp + 1;
                }
            }
        }
        if(botTemp <= topTemp) {
            break;
        }

        if(pseudogeneStop_) {
            if(botTemp - topTemp < bot - top && bot - top <= 5) {
                static const index_t minLenForPseudogene = _minK + 6;
                if(dep - offset >= minLenForPseudogene && similar_range >= 5) {
                    hit._numUniqueSearch++;
                    pseudogeneStop = true;
                    break;
                }
            }
            if(botTemp - topTemp != 1) {
                if(botTemp - topTemp + 2 >= bot - top) similar_range++;
                else if(botTemp - topTemp + 4 < bot - top) similar_range = 0;
            } else {
                pseudogeneStop_ = false;
            }
        }
        
        if(anchorStop_) {
            if(botTemp - topTemp != 1 && bot - top == botTemp - topTemp) {
                same_range++;
                if(same_range >= 5) {
                    anchorStop_ = false;
                }
            } else {
                same_range = 0;
            }
        
            if(dep - offset >= _minK + 8 && botTemp - topTemp >= 4) {
                anchorStop_ = false;
            }
        }
        
        top = topTemp;
        bot = botTemp;
        dep++;

        if(anchorStop_) {
            if(dep - offset >= _minK + 12 && bot - top == 1) {
                hit._numUniqueSearch++;
                anchorStop = true;
                break;
            }
        }
        
        HIER_INIT_LOCS(top, bot, tloc, bloc, ebwt);
    }
    
    // Done
    if(bot > top) {
        // This is an exact hit
        assert_gt(dep, offset);
        assert_leq(dep, len);
        partialHits.expand();
        index_t hit_type = CANDIDATE_HIT;
        if(anchorStop) hit_type = ANCHOR_HIT;
        else if(pseudogeneStop) hit_type = PSEUDOGENE_HIT;
        partialHits.back().init(top,
                                bot,
                                fw,
                                (index_t)offset,
                                (index_t)(dep - offset),
                                hit_type);
        
        nelt += (bot - top);
        cur = dep;
        if(cur >= hit._len) {
            if(hit_type == CANDIDATE_HIT) hit._numUniqueSearch++;
            hit.done(true);
        }
    }
    return nelt;
}


/**
 */
template <typename index_t, typename local_index_t>
size_t HI_Aligner<index_t, local_index_t>::globalEbwtSearch(
                                                            const Ebwt<index_t>& ebwt,  // BWT index
                                                            const Read&          read,  // read to align
                                                            const Scoring&       sc,    // scoring scheme
                                                            bool                 fw,
                                                            index_t              hitoff,
                                                            index_t&             hitlen,
                                                            index_t&             top,
                                                            index_t&             bot,
                                                            RandomSource&        rnd,
                                                            bool&                uniqueStop,
                                                            index_t              maxHitLen)
{
    bool uniqueStop_ = uniqueStop;
    uniqueStop = false;
    const index_t ftabLen = ebwt.eh().ftabChars();
	SideLocus<index_t> tloc, bloc;
	const index_t len = (index_t)read.length();
    
	size_t nelt = 0;
    const BTDnaString& seq = fw ? read.patFw : read.patRc;
    assert(!seq.empty());
    
    index_t offset = len - hitoff - 1;
    index_t dep = offset;
    top = 0, bot = 0;
    index_t topTemp = 0, botTemp = 0;
    index_t left = len - dep;
    assert_gt(left, 0);
    if(left < ftabLen) {
#if 1
        hitlen = left;
        return 0;
#else
        // Use fchr
        int c = seq[len-dep-1];
        if(c < 4) {
            top = ebwt.fchr()[c];
            bot = ebwt.fchr()[c+1];
        } else {
            hitlen = left;
            return 0;
        }
        dep++;
#endif
    } else {
        // Does N interfere with use of Ftab?
        for(index_t i = 0; i < ftabLen; i++) {
            int c = seq[len-dep-1-i];
            if(c > 3) {
                hitlen = (i+1);
                return 0;
            }
        }
        
        // Use ftab
        ebwt.ftabLoHi(seq, len - dep - ftabLen, false, top, bot);
        dep += ftabLen;
        if(bot <= top) {
            hitlen = ftabLen;
            return 0;
        }
    }
    
    HIER_INIT_LOCS(top, bot, tloc, bloc, ebwt);
    // Keep going
    while(dep < len) {
        int c = seq[len-dep-1];
        if(c > 3) {
            topTemp = botTemp = 0;
        } else {
            if(bloc.valid()) {
                bwops_ += 2;
                topTemp = ebwt.mapLF(tloc, c);
                botTemp = ebwt.mapLF(bloc, c);
            } else {
                bwops_++;
                topTemp = ebwt.mapLF1(top, tloc, c);
                if(topTemp == (index_t)OFF_MASK) {
                    topTemp = botTemp = 0;
                } else {
                    botTemp = topTemp + 1;
                }
            }
        }
        if(botTemp <= topTemp) {
            break;
        }
        
        top = topTemp;
        bot = botTemp;
        dep++;
        
        if(uniqueStop_) {
            if(bot - top == 1 && dep - offset >= _minK) {
                uniqueStop = true;
                break;
            }
        }
        
        HIER_INIT_LOCS(top, bot, tloc, bloc, ebwt);
    }
    
    // Done
    if(bot > top) {
        assert_gt(dep, offset);
        assert_leq(dep, len);
        nelt += (bot - top);
        hitlen = dep - offset;
    }
    return nelt;
}


/**
 *
 **/
template <typename index_t, typename local_index_t>
size_t HI_Aligner<index_t, local_index_t>::localEbwtSearch(
                                                           const LocalEbwt<local_index_t, index_t>*  ebwtFw,  // BWT index
                                                           const LocalEbwt<local_index_t, index_t>*  ebwtBw,  // BWT index
                                                           const Read&                      read,    // read to align
                                                           const Scoring&                   sc,      // scoring scheme
                                                           bool                             fw,
                                                           bool                             searchfw,
                                                           index_t                          rdoff,
                                                           index_t&                         hitlen,
                                                           local_index_t&                   top,
                                                           local_index_t&                   bot,
                                                           RandomSource&                    rnd,
                                                           bool&                         	uniqueStop,
                                                           local_index_t                    minUniqueLen,
                                                           local_index_t                    maxHitLen)
{
#ifndef NDEBUG
    if(searchfw) {
        assert(ebwtBw != NULL);
    } else {
        assert(ebwtFw != NULL);
    }
#endif
    bool uniqueStop_ = uniqueStop;
    uniqueStop = false;
    const LocalEbwt<local_index_t, index_t>& ebwt = *(searchfw ? ebwtBw : ebwtFw);
	const local_index_t ftabLen = (local_index_t)ebwt.eh().ftabChars();
	SideLocus<local_index_t> tloc, bloc;
	const local_index_t len = (local_index_t)read.length();
	size_t nelt = 0;
    
    const BTDnaString& seq = fw ? read.patFw : read.patRc;
    assert(!seq.empty());
    
    local_index_t offset = searchfw ? rdoff : len - rdoff - 1;
    local_index_t dep = offset;
    top = 0, bot = 0;
    local_index_t topTemp = 0, botTemp = 0;
    local_index_t left = len - dep;
    assert_gt(left, 0);
    if(left < ftabLen) {
        hitlen = left;
		return 0;
    }
    // Does N interfere with use of Ftab?
    for(local_index_t i = 0; i < ftabLen; i++) {
        int c = searchfw ? seq[dep+i] : seq[len-dep-1-i];
        if(c > 3) {
            hitlen = i + 1;
			return 0;
        }
    }
    
    // Use ftab
    if(searchfw) {
        ebwt.ftabLoHi(seq, dep, false, top, bot);
    } else {
        ebwt.ftabLoHi(seq, len - dep - ftabLen, false, top, bot);
    }
    dep += ftabLen;
    if(bot <= top) {
        hitlen = ftabLen;
        return 0;
    }
    LOCAL_INIT_LOCS(top, bot, tloc, bloc, ebwt);
    // Keep going
    while(dep < len) {
        int c = searchfw ? seq[dep] : seq[len-dep-1];
        if(c > 3) {
            topTemp = botTemp = 0;
        } else {
            if(bloc.valid()) {
                bwops_ += 2;
                topTemp = ebwt.mapLF(tloc, c);
                botTemp = ebwt.mapLF(bloc, c);
            } else {
                bwops_++;
                topTemp = ebwt.mapLF1(top, tloc, c);
                if(topTemp == (local_index_t)OFF_MASK) {
                    topTemp = botTemp = 0;
                } else {
                    botTemp = topTemp + 1;
                }
            }
        }
        if(botTemp <= topTemp) {
            break;
        }
        top = topTemp;
        bot = botTemp;
        LOCAL_INIT_LOCS(top, bot, tloc, bloc, ebwt);
        dep++;

        if(uniqueStop_) {
            if(bot - top == 1 && dep - offset >= minUniqueLen) {
                uniqueStop = true;
                break;
            }
        }
        
        if(dep - offset >= maxHitLen) break;
    }
    
    // Done
    if(bot > top) {
        assert_gt(dep, offset);
        assert_leq(dep, len);
        nelt += (bot - top);
        hitlen = dep - offset;
    }

    return nelt;
}

/**
 *
 **/
template <typename index_t, typename local_index_t>
size_t HI_Aligner<index_t, local_index_t>::localEbwtSearch_reverse(
                                                                   const LocalEbwt<local_index_t, index_t>*  ebwtFw,  // BWT index
                                                                   const LocalEbwt<local_index_t, index_t>*  ebwtBw,  // BWT index
                                                                   const Read&                      read,    // read to align
                                                                   const Scoring&                   sc,      // scoring scheme
                                                                   bool                             fw,
                                                                   bool                             searchfw,
                                                                   index_t                          rdoff,
                                                                   index_t&                         hitlen,
                                                                   local_index_t&                   top,
                                                                   local_index_t&                   bot,
                                                                   RandomSource&                    rnd,
                                                                   bool&                         	uniqueStop,
                                                                   local_index_t                    minUniqueLen,
                                                                   local_index_t                    maxHitLen)
{
#ifndef NDEBUG
    if(searchfw) {
        assert(ebwtBw != NULL);
    } else {
        assert(ebwtFw != NULL);
    }
#endif
    bool uniqueStop_ = uniqueStop;
    uniqueStop = false;
    const LocalEbwt<local_index_t, index_t>& ebwt = *(searchfw ? ebwtBw : ebwtFw);
	const local_index_t ftabLen = (local_index_t)ebwt.eh().ftabChars();
	SideLocus<local_index_t> tloc, bloc;
	const local_index_t len = (local_index_t)read.length();
	size_t nelt = 0;
    
    const BTDnaString& seq = fw ? read.patFw : read.patRc;
    assert(!seq.empty());
    
    local_index_t offset = searchfw ? len - rdoff - 1 : rdoff;
    local_index_t dep = offset;
    top = 0, bot = 0;
    local_index_t topTemp = 0, botTemp = 0;
    local_index_t left = len - dep;
    assert_gt(left, 0);
    if(left < ftabLen) {
        hitlen = left;
		return 0;
    }
    // Does N interfere with use of Ftab?
    for(local_index_t i = 0; i < ftabLen; i++) {
        int c = searchfw ? seq[len-dep-1-i] : seq[dep+i];
        if(c > 3) {
            hitlen = i + 1;
			return 0;
        }
    }
    
    // Use ftab
    if(searchfw) {
        ebwt.ftabLoHi(seq, len - dep - ftabLen, false, top, bot);
    } else {
        ebwt.ftabLoHi(seq, dep, false, top, bot);
    }
    dep += ftabLen;
    if(bot <= top) {
        hitlen = ftabLen;
        return 0;
    }
    LOCAL_INIT_LOCS(top, bot, tloc, bloc, ebwt);
    // Keep going
    while(dep < len) {
        int c = searchfw ? seq[len-dep-1] : seq[dep];
        if(c > 3) {
            topTemp = botTemp = 0;
        } else {
            if(bloc.valid()) {
                bwops_ += 2;
                topTemp = ebwt.mapLF(tloc, c);
                botTemp = ebwt.mapLF(bloc, c);
            } else {
                bwops_++;
                topTemp = ebwt.mapLF1(top, tloc, c);
                if(topTemp == (local_index_t)OFF_MASK) {
                    topTemp = botTemp = 0;
                } else {
                    botTemp = topTemp + 1;
                }
            }
        }
        if(botTemp <= topTemp) {
            break;
        }
        top = topTemp;
        bot = botTemp;
        LOCAL_INIT_LOCS(top, bot, tloc, bloc, ebwt);
        dep++;
        
        if(uniqueStop_) {
            if(bot - top == 1 && dep - offset >= minUniqueLen) {
                uniqueStop = true;
                break;
            }
        }
        
        if(dep - offset >= maxHitLen) break;
    }
    
    // Done
    if(bot > top) {
        assert_gt(dep, offset);
        assert_leq(dep, len);
        nelt += (bot - top);
        hitlen = dep - offset;
    }
    
    return nelt;
}

/**
 *
 **/
template <typename index_t, typename local_index_t>
bool HI_Aligner<index_t, local_index_t>::isSearched(
                                                    const GenomeHit<index_t>&   hit,
                                                    index_t                     rdi)
{
    assert_lt(rdi, 2);
    EList<GenomeHit<index_t> >& searchedHits = _hits_searched[rdi];
    for(index_t i = 0; i < searchedHits.size(); i++) {
        if(searchedHits[i].contains(hit)) return true;
    }
    return false;
}

/**
 *
 **/
template <typename index_t, typename local_index_t>
void HI_Aligner<index_t, local_index_t>::addSearched(
                                                     const GenomeHit<index_t>&   hit,
                                                     index_t                     rdi)
{
    assert_lt(rdi, 2);
    assert(!isSearched(hit, rdi));
    EList<GenomeHit<index_t> >& searchedHits = _hits_searched[rdi];
    searchedHits.push_back(hit);
}

#endif /*HI_ALIGNER_H_*/
