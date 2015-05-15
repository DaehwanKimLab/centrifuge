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
        // assert(!_done);
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
        //assert_lt(offset, _len); //FIXME: assertion fails as offset == _len
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
    
    EList<BWTHit<index_t> >  _partialHits;
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
		   SpeciesMetrics&          spm,
           RandomSource&            rnd,
           AlnSinkWrap<index_t>&    sink) = 0;
    
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

#endif /*HI_ALIGNER_H_*/
