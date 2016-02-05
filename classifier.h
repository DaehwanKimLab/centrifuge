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

#ifndef CLASSIFIER_H_
#define CLASSIFIER_H_

//#define LI_DEBUG

#include <algorithm>
#include <vector>
#include "hi_aligner.h"
#include "util.h"

template<typename index_t>
struct HitCount {
    uint64_t uniqueID;
    uint64_t taxID;
    uint32_t count;
    uint32_t score;
    uint32_t scores[2][2];      // scores[rdi][fwi]
    double summedHitLen;
    double summedHitLens[2][2]; // summedHitLens[rdi][fwi]
    uint32_t timeStamp;
    EList<pair<uint32_t,uint32_t> > readPositions;
    
    uint32_t rank;
    EList<uint64_t> path;
    
    void reset() {
        uniqueID = taxID = count = score = timeStamp = 0;
        scores[0][0] = scores[0][1] = scores[1][0] = scores[1][1] = 0;
        summedHitLen = 0.0;
        summedHitLens[0][0] = summedHitLens[0][1] = summedHitLens[1][0] = summedHitLens[1][1] = 0.0;
        readPositions.clear();
    }
    
    HitCount& operator=(const HitCount& o) {
        if(this == &o)
            return *this;
        
        uniqueID = o.uniqueID;
        taxID = o.taxID;
        count = o.count;
        score = o.score;
        scores[0][0] = o.scores[0][0];
        scores[0][1] = o.scores[0][1];
        scores[1][0] = o.scores[1][0];
        scores[1][1] = o.scores[1][1];
        summedHitLen = o.summedHitLen;
        summedHitLens[0][0] = o.summedHitLens[0][0];
        summedHitLens[0][1] = o.summedHitLens[0][1];
        summedHitLens[1][0] = o.summedHitLens[1][0];
        summedHitLens[1][1] = o.summedHitLens[1][1];
        timeStamp = o.timeStamp;
        readPositions = o.readPositions;
        rank = o.rank;
        path = o.path;
        
        return *this;
    }

    void finalize(
                  const Ebwt<index_t>& ebwt,
                  bool paired,
                  bool mate1fw,
                  bool mate2fw) {
        rank = 0; path.clear();
        const TaxonomyPathTable& pathTable = ebwt.paths();
        pathTable.getPath(taxID, path);
        for(; rank < path.size(); rank++) {
            if(path[rank] == taxID) {
                break;
            }
        }
        if(paired) {
            uint32_t score1 = 0, score2 = 0;
            double summedHitLen1 = 0.0, summedHitLen2 = 0.0;
            if(mate1fw == mate2fw) {
                score1 = scores[0][0] + scores[1][0];
                score2 = scores[0][1] + scores[1][1];
                summedHitLen1 = summedHitLens[0][0] + summedHitLens[1][0];
                summedHitLen2 = summedHitLens[0][1] + summedHitLens[1][1];
            } else {
                score1 = scores[0][0] + scores[1][1];
                score2 = scores[0][1] + scores[1][0];
                summedHitLen1 = summedHitLens[0][0] + summedHitLens[1][1];
                summedHitLen2 = summedHitLens[0][1] + summedHitLens[1][0];
            }
            if(score1 >= score2) {
                score = score1;
                summedHitLen = summedHitLen1;
            } else {
                score = score2;
                summedHitLen = summedHitLen2;
            }
        } else {
            score = max<uint32_t>(scores[0][0], scores[0][1]);
            summedHitLen = max<double>(summedHitLens[0][0], summedHitLens[0][1]);
        }
    }
};

/**
 * With a hierarchical indexing, SplicedAligner provides several alignment strategies
 * , which enable effective alignment of RNA-seq reads
 */
template <typename index_t, typename local_index_t>
class Classifier : public HI_Aligner<index_t, local_index_t> {
    
public:
    
    /**
     * Initialize with index.
     */
    Classifier(const Ebwt<index_t>& ebwt,
               const EList<string>& refnames,
               const EList<uint32_t>& hostGenomes,
               bool mate1fw,
               bool mate2fw,
               index_t minHitLen = 22) :
    HI_Aligner<index_t, local_index_t>(
                                       ebwt,
                                       0,    // don't make use of splice sites found by earlier reads
                                       true), // no spliced alignment
    _refnames(refnames),
    _minHitLen(minHitLen),
    _hostGenomes(hostGenomes),
    _mate1fw(mate1fw),
    _mate2fw(mate2fw)
    {
    }
    
    ~Classifier() {
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
           AlnSinkWrap<index_t>&    sink)
    {
        _hitMap.clear();

        const index_t increment = (2 * _minHitLen <= 33) ? 10 : (2 * _minHitLen - 33);
        const ReportingParams& rp = sink.reportingParams();
        index_t maxGenomeHitSize = rp.khits;
		bool isFw = false;
        
        //
        uint32_t ts = 0; // time stamp
        uint32_t human_speciesID = 0, human_genusID = 0;
        // for each mate. only called once for unpaired data
        for(int rdi = 0; rdi < (this->_paired ? 2 : 1); rdi++) {
            assert(this->_rds[rdi] != NULL);
            
            // search for partial hits on the forward and reverse strand (saved in this->_hits[rdi])
            searchForwardAndReverse(rdi, ebwtFw, sc, rnd, rp, increment);
            
            // get forward or reverse hits for this read from this->_hits[rdi]
            //  the strand is chosen based on higher average hit length in either direction
            pair<int, int> fwp = getForwardOrReverseHit(rdi);
            for(int fwi = fwp.first; fwi < fwp.second; fwi++) {
                ReadBWTHit<index_t>& hit = this->_hits[rdi][fwi];
                assert(hit.done());
                isFw = hit._fw;  // TODO: Sync between mates!
                
                // choose candidate partial alignments for further alignment
                index_t offsetSize = hit.offsetSize();
                this->_genomeHits.clear();
                
                // sort partial hits by size (number of genome positions), ascending, and then length, descending
                for(size_t hi = 0; hi < offsetSize; hi++) {
                    const BWTHit<index_t> partialHit = hit.getPartialHit(hi);
#ifdef LI_DEBUG
                    cout << partialHit.len() << " " << partialHit.size() << endl;
#endif
                    if(partialHit.len() >= _minHitLen && partialHit.size() > maxGenomeHitSize) {
                        maxGenomeHitSize = partialHit.size();
                    }
                }
                
                if(maxGenomeHitSize > (index_t)rp.khits) {
                    maxGenomeHitSize += rp.khits;
                }
                
                hit._partialHits.sort(compareBWTHits());
                size_t usedPortion = 0;
                size_t genomeHitCnt = 0;
                for(size_t hi = 0; hi < offsetSize; hi++, ts++) {
                    const BWTHit<index_t>& partialHit = hit.getPartialHit(hi);
                    size_t partialHitLen = partialHit.len();
                    if(partialHitLen <= _minHitLen) continue;
                    if(partialHit.size() == 0) continue;
                    
                    // only keep this partial hit if it is equal to or bigger than minHitLen (default: 22 bp)
                    // TODO: consider not requiring minHitLen when we have already hits to the same genome
                    bool considerOnlyIfPreviouslyObserved = partialHitLen < _minHitLen;
                    
                    // get all coordinates of the hit
                    EList<Coord>& coords = getCoords(
                                                     hit,
                                                     hi,
                                                     ebwtFw,
                                                     ref,
                                                     rnd,
                                                     maxGenomeHitSize,
                                                     wlm,
                                                     prm,
                                                     him);
                    if(coords.empty())
                        continue;
                    
                    usedPortion += partialHitLen;
                    assert_gt(coords.size(), 0);
                    
                    // the maximum number of hits per read is maxGenomeHitSize (change with parameter -k)
                    size_t nHitsToConsider = coords.size();
                    
                    // daehwan - for debugging purposes
#if 1
                    if(coords.size() > rp.ihits) {
                        continue;
                    }
#else
                    if(genomeHitCnt + coords.size() > maxGenomeHitSize) {
                        coords.shufflePortion(0, coords.size(), rnd);
                        nHitsToConsider = maxGenomeHitSize - genomeHitCnt;
                    }
#endif
                    
                    // find the genome id for all coordinates, and count the number of genomes
                    EList<pair<uint64_t, uint64_t> > coord_ids;
                    for(index_t k = 0; k < nHitsToConsider; k++, genomeHitCnt++) {
                        const Coord& coord = coords[k];
                        assert_lt(coord.ref(), _refnames.size()); // gives a warning - coord.ref() is signed integer. why?
                        
                        // extract numeric id from refName
                        const EList<pair<string, uint64_t> >& uid_to_tid = ebwtFw.uid_to_tid();
                        assert_lt(coord.ref(), uid_to_tid.size());
                        uint64_t taxID = uid_to_tid[coord.ref()].second;
                        bool found = false;
                        for(index_t k2 = 0; k2 < coord_ids.size(); k2++) {
                            // count the genome if it is not in coord_ids, yet
                            if(coord_ids[k2].first == (uint64_t)coord.ref()) {
                                found = true;
                                break;
                            }
                        }
                        if(found) continue;
                        // add to coord_ids
                        coord_ids.expand();
                        coord_ids.back().first = coord.ref();
                        coord_ids.back().second = taxID;
                    }
                    
                    ASSERT_ONLY(size_t n_genomes = coord_ids.size());
                    // scoring function: calculate the weight of this partial hit
                    assert_gt(partialHitLen, 15);
                    assert_gt(n_genomes, 0);
                    uint32_t partialHitScore = (uint32_t)((partialHitLen - 15) * (partialHitLen - 15)) ; // / n_genomes;
                    double weightedHitLen = double(partialHitLen) ; // / double(n_genomes) ;
                    
                    // go through all coordinates reported for partial hit
                    for(index_t k = 0; k < coord_ids.size(); ++k) {
                        uint64_t uniqueID = coord_ids[k].first;
                        uint64_t taxID = coord_ids[k].second;
                        uint32_t speciesID = (uint32_t)taxID, genusID = (uint32_t)taxID;
                        const map<uint64_t, TaxonomyNode>& tree = ebwtFw.tree();
                        uint64_t tmp_taxID = taxID;
                        while(tree.find(tmp_taxID) != tree.end()) {
                            const TaxonomyNode& node = tree.find(tmp_taxID)->second;
                            if(node.rank == RANK_SPECIES) {
                                speciesID = (uint32_t)tmp_taxID;
                            } else if(node.rank == RANK_GENUS) {
                                genusID = (uint32_t)tmp_taxID;
                                break;
                            }
                            if(tmp_taxID == node.parent_tid) break;
                            tmp_taxID = node.parent_tid;
                        }
                        
                        if(_hostGenomes.size() > 0 &&
                           _hostGenomes.back() == speciesID &&
                           partialHit.len() >= 31) {
                            human_speciesID = speciesID;
                            human_genusID = genusID;
                        }
                        
                        // add hit to genus map and get new index in the map
                        size_t idx = addHitToHitMap(
                                                    _hitMap,
                                                    rdi,
                                                    fwi,
                                                    uniqueID,
                                                    taxID,
                                                    ts,
                                                    partialHitScore,
                                                    weightedHitLen,
                                                    considerOnlyIfPreviouslyObserved,
                                                    partialHit._bwoff,
                                                    partialHit.len());
                        
                        //if considerOnlyIfPreviouslyObserved and it was not found, genus Idx size is equal to the genus Map size
                        if(idx >= _hitMap.size()) {
                            continue;
                        }
                        
#ifndef NDEBUG //FB
                        std::cerr << speciesID << ';';
#endif
                    }
                    
                    if(genomeHitCnt >= maxGenomeHitSize)
                        break;
                    
#ifndef NDEBUG //FB
                    std::cerr << "  partialHits-done";
#endif
                } // partialHits
            } // fwi
            
#ifndef NDEBUG //FB
            std::cerr << "  rdi-done" << endl;
#endif
        } // rdi
        
        for(size_t i = 0; i < _hitMap.size(); i++) {
            _hitMap[i].finalize(ebwtFw,
                                this->_paired,
                                this->_mate1fw,
                                this->_mate2fw);
        }
        
        // If the number of hits is more than -k,
        //   traverse up the taxonomy tree to reduce the number
        if(_hitMap.size() > 0) {
            uint32_t score = _hitMap[0].score;
            size_t count = 1;
            for(size_t i = 1; i < _hitMap.size(); i++) {
                if(score < _hitMap[i].score) {
                    score = _hitMap[i].score;
                    count = 1;
                } else if(score == _hitMap[i].score) {
                    count += 1;
                }
            }
            
            // daehwan - in the middle of implementation
#if 0
            if(count > (index_t)rp.khits) {
                return 0;
            }
#else
            uint8_t rank = 0;
            while(count > (index_t)rp.khits) {
                _hitTaxCount.clear();
                for(size_t i = 0; i < _hitMap.size(); i++) {
                    assert_geq(_hitMap[i].rank, rank);
                    if(_hitMap[i].rank != rank) continue;
                    if(rank + 1 >= _hitMap[i].path.size()) continue;
                    uint64_t parent_taxID = _hitMap[i].path[rank + 1];
                    if(parent_taxID == 0) continue;
                    
                    size_t j = 0;
                    for(; j < _hitTaxCount.size(); j++) {
                        if(_hitTaxCount[j].second == parent_taxID) {
                            _hitTaxCount[j].first += 1;
                            break;
                        }
                    }
                    if(j == _hitTaxCount.size()) {
                        _hitTaxCount.expand();
                        _hitTaxCount.back().first = 1;
                        _hitTaxCount.back().second = parent_taxID;
                    }
                }
                _hitTaxCount.sort();
                if(_hitTaxCount.size() >= 1) {
                    size_t j = _hitTaxCount.size();
                    while(j-- > 0) {
                        uint64_t parent_taxID = _hitTaxCount[j].second;
                        int64_t max_score = 0;
                        for(size_t i = 0; i < _hitMap.size(); i++) {
                            assert_geq(_hitMap[i].rank, rank);
                            if(_hitMap[i].rank != rank) continue;
                            if(rank + 1 >= _hitMap[i].path.size()) continue;
                            if(parent_taxID == _hitMap[i].path[rank + 1]) {
                                _hitMap[i].uniqueID = 0;
                                _hitMap[i].taxID = parent_taxID;
                                _hitMap[i].rank = rank + 1;
                            }
                            if(parent_taxID == _hitMap[i].taxID) {
                                if(_hitMap[i].score > max_score) {
                                    max_score = _hitMap[i].score;
                                }
                            }
                        }
                        
                        bool first = true;
                        for(size_t i = 0; i < _hitMap.size(); i++) {
                            if(parent_taxID == _hitMap[i].taxID) {
                                bool del = _hitMap[i].score < max_score;
                                del |= (!first && _hitMap[i].score == max_score);
                                if(del) {
                                    if(i + 1 < _hitMap.size()) {
                                        _hitMap.back() = _hitMap[i];
                                    }
                                    _hitMap.pop_back();
                                    i--;
                                } else {
                                    first = false;
                                }
                            }
                        }
                        
                        int kk = 20;
                    }
                }
            }
#endif
        }
       
#if 0
       	// boost up the score if the assignment is unique
        if(_hitMap.size() == 1) {
            HitCount& hitCount = _hitMap[0];
            hitCount.score = (hitCount.summedHitLen - 15) * (hitCount.summedHitLen - 15);
        }
#endif
        
        index_t rdlen = this->_rds[0]->length();
        int64_t max_score = (rdlen > 15 ? (rdlen - 15) * (rdlen - 15) : 0);
        if(this->_paired) {
            rdlen = this->_rds[1]->length();
            max_score += (rdlen > 15 ? (rdlen - 15) * (rdlen - 15) : 0);
        }
        
        for(size_t gi = 0; gi < _hitMap.size(); gi++) {
            assert_gt(_hitMap[gi].score, 0);
            HitCount<index_t>& hitCount = _hitMap[gi];
            if(human_genusID != 0 && human_genusID != hitCount.taxID) continue;
            const EList<pair<string, uint64_t> >& uid_to_tid = ebwtFw.uid_to_tid();
            assert_lt(hitCount.uniqueID, uid_to_tid.size());
            // report
            AlnRes rs;
            rs.init(
                    hitCount.score,
                    max_score,
                    uid_to_tid[hitCount.uniqueID].first,
                    hitCount.taxID,
                    hitCount.summedHitLen,
                    hitCount.readPositions,
                    isFw);
            sink.report(0, &rs);
        }
        return 0;
    }
    
    bool getGenomeIdx(
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
        this->_offs.resize(nelt);
        this->_offs.fill(std::numeric_limits<index_t>::max());
        this->_sas.init(top, rdlen, EListSlice<index_t, 16>(this->_offs, 0, nelt));
        this->_gws.init(ebwt, ref, this->_sas, rnd, met);
        for(index_t off = 0; off < nelt; off++) {
            WalkResult<index_t> wr;
            this->_gws.advanceElement(
                                off,
                                ebwt,         // forward Bowtie index for walking left
                                ref,          // bitpair-encoded reference
                                this->_sas,   // SA range with offsets
                                this->_gwstate,     // GroupWalk state; scratch space
                                wr,           // put the result here
                                met,          // metrics
                                prm);         // per-read metrics
            // Coordinate of the seed hit w/r/t the pasted reference string
            coords.expand();
            coords.back().init(wr.toff, 0, fw);
        }
        
        return true;
    }
private:
    EList<string>                _refnames;
    EList<HitCount<index_t> >    _hitMap;
    index_t                      _minHitLen;
    EList<uint16_t>              _tempTies;
    EList<uint32_t>              _hostGenomes;
    bool                         _mate1fw;
    bool                         _mate2fw;
    
    // Temporary variables
    ReadBWTHit<index_t>          _tempHit;
    EList<pair<uint32_t, uint64_t> > _hitTaxCount;  // pair of count and taxID

    void searchForwardAndReverse(
                                 index_t rdi,
                                 const Ebwt<index_t>& ebwtFw,
                                 const Scoring& sc,
                                 RandomSource& rnd,
                                 const ReportingParams& rp,
                                 const index_t increment)
    {
        const Read& rd = *(this->_rds[rdi]);

        bool done[2] = {false, false};
        size_t cur[2] = {0, 0} ;
        
        index_t rdlen = rd.length();
        //const size_t maxDiff = (rdlen / 2 > 2 * _minHitLen) ? rdlen / 2 : (2 * _minHitLen);
        size_t sum[2] = {0, 0} ;
        
        // search for partial hits on the forward and reverse strand
        while(!done[0] || !done[1]) {
            for(index_t fwi = 0; fwi < 2; fwi++) {
                if(done[fwi])
                    continue;
                
                size_t mineFw = 0, mineRc = 0;
                bool fw = (fwi == 0);
                ReadBWTHit<index_t>& hit = this->_hits[rdi][fwi];
                this->partialSearch(
                                    ebwtFw,
                                    rd,
                                    sc,
                                    fw,
                                    0,
                                    mineFw,
                                    mineRc,
                                    hit,
                                    rnd);
                
                BWTHit<index_t>& lastHit = hit.getPartialHit(hit.offsetSize() - 1);
                if(hit.done()) {
                    done[fwi] = true;
                    cur[fwi] = rdlen;
                    if(lastHit.len() >= _minHitLen) {
                        sum[fwi] += lastHit.len();
                        if(0) //lastHit.len() < 31 && rdlen > 31 && lastHit.size() == 1 )
                        {
                            ReadBWTHit<index_t> testHit ;
                            testHit.init( fw, rdlen ) ;
                            testHit.setOffset(hit.cur() - 1 - 31 + 1);
                            this->partialSearch(ebwtFw,
                                                rd,
                                                sc,
                                                fw,
                                                0,
                                                mineFw,
                                                mineRc,
                                                testHit,
                                                rnd);
                            index_t tmpLen = testHit.getPartialHit( testHit.offsetSize() - 1 ).len();
#ifdef LI_DEBUG
                            cout << "(adjust: " << tmpLen << ")";
#endif
                            if(tmpLen >= 31) {
                                lastHit._len = tmpLen;
                            }
                        }
                    }
                    
                    continue;
                }
                
                cur[fwi] = hit.cur();
#ifdef LI_DEBUG
                cout << fwi << ":" << lastHit.len() << " " << cur[fwi] << " ";
#endif
                if(lastHit.len() >= _minHitLen)
                    sum[fwi] += lastHit.len();
                
                if(lastHit.len() > increment) {
                    if(lastHit.len() < _minHitLen) {
                        // daehwan - for debugging purposes
#if 1
                        hit.setOffset(hit.cur() + 1);
#else
                        hit.setOffset(hit.cur() - increment);
#endif
                    } else {
                        hit.setOffset(hit.cur() + 1);
                        if(0) //lastHit.len() < 31 && hit.cur() >= 31 && lastHit.size() == 1 )
                        {
                            ReadBWTHit<index_t> testHit;
                            testHit.init(fw, rdlen);
                            testHit.setOffset(hit.cur() - 1 - 31); // why not hit.cur() - 1 - 31 + 1? because we "+1" before the if!
                            
                            this->partialSearch(ebwtFw,
                                                rd,
                                                sc,
                                                fw,
                                                0,
                                                mineFw,
                                                mineRc,
                                                testHit,
                                                rnd);
                            index_t tmpLen = testHit.getPartialHit(testHit.offsetSize() - 1 ).len();
#ifdef LI_DEBUG
                            cout << "(adjust: " << tmpLen << ")";
#endif
                            if(tmpLen >= 31) {
                                lastHit._len = tmpLen;
                            }
                        }
                    }
                }
                if(hit.cur() + _minHitLen >= rdlen) {
                    hit.done(true);
                    done[fwi] = true;
                    continue;
                }

                if(lastHit.len() <= 3) {
                    // This happens most likely due to the Ns in the read
                    --fwi ; // Repeat this strand again.
                }
            }
#ifdef LI_DEBUG
            cout << endl;
#endif
            if(sum[0] > sum[1] + (rdlen - cur[1] + 1)) {
                // daehwan - for debugging purposes
#if 0
                this->_hits[rdi][1].done(true);
                done[1] = true;
#endif
            } else if(sum[1] > sum[0] + (rdlen - cur[0] + 1)) {
                // daehwan - for debugging purposes
#if 0
                this->_hits[rdi][0].done(true);
                done[0] = true;
#endif
            }
        }
        
        // Extend partial hits
        if(sum[0] >= _minHitLen && sum[1] >= _minHitLen) {
            ReadBWTHit<index_t>& hits = this->_hits[rdi][0];
            ReadBWTHit<index_t>& rchits = this->_hits[rdi][1];
            for(size_t i = 0; i < hits.offsetSize(); i++) {
                BWTHit<index_t>& hit = hits.getPartialHit(i);
                index_t len = hit.len();
                //if(len < _minHitLen) continue;
                index_t l = hit._bwoff;
                index_t r = hit._bwoff + len;
                for(size_t j = 0; j < rchits.offsetSize(); j++) {
                    BWTHit<index_t>& rchit = rchits.getPartialHit(j);
                    index_t rclen = rchit.len();
                    if(len < _minHitLen && rclen < _minHitLen) continue;
                    index_t rc_l = rdlen - rchit._bwoff - rchit._len;
                    index_t rc_r = rc_l + rclen;
                    if(r <= rc_l) continue;
                    if(rc_r <= l) continue;
                    if(l == rc_l && r == rc_r) continue;
                    if(l < rc_l && r > rc_r) continue;
                    if(l > rc_l && r < rc_r) continue;
                    if(l > rc_l) {
                        _tempHit.init(true /* fw */, rdlen);
                        _tempHit.setOffset(rc_l);
                        size_t mineFw = 0, mineRc = 0;
                        this->partialSearch(ebwtFw,
                                            rd,
                                            sc,
                                            true, // fw
                                            0,
                                            mineFw,
                                            mineRc,
                                            _tempHit,
                                            rnd);
                        BWTHit<index_t>& tmphit = _tempHit.getPartialHit(0);
                        if(tmphit.len() == len + l - rc_l) {
                            hit = tmphit;
                        }
                    }
                    if(r > rc_r) {
                        _tempHit.init(false /* fw */, rdlen);
                        _tempHit.setOffset(rdlen - r);
                        size_t mineFw = 0, mineRc = 0;
                        this->partialSearch(ebwtFw,
                                            rd,
                                            sc,
                                            false, // fw
                                            0,
                                            mineFw,
                                            mineRc,
                                            _tempHit,
                                            rnd);
                        BWTHit<index_t>& tmphit = _tempHit.getPartialHit(0);
                        if(tmphit.len() == rclen + r - rc_r) {
                            rchit = tmphit;
                        }
                    }
                }
            }
            
            // Remove partial hits that are mapped more than user-specified number
            for(size_t i = 0; i < hits.offsetSize(); i++) {
                BWTHit<index_t>& hit = hits.getPartialHit(i);
                index_t len = hit.len();
                index_t l = hit._bwoff;
                index_t r = hit._bwoff + len;
                for(size_t j = 0; j < rchits.offsetSize(); j++) {
                    BWTHit<index_t>& rchit = rchits.getPartialHit(j);
                    index_t rclen = rchit.len();
                    index_t rc_l = rdlen - rchit._bwoff - rchit._len;
                    index_t rc_r = rc_l + rclen;
                    if(rc_l < l) break;
                    if(len != rclen) continue;
                    if(l == rc_l &&
                       r == rc_r &&
                       hit.size() + rchit.size() > rp.ihits) {
                        hit.reset();
                        rchit.reset();
                        break;
                    }
                }
            }
        }
        
        // Trim partial hits
        for(int fwi = 0; fwi < 2; fwi++) {
            ReadBWTHit<index_t>& hits = this->_hits[rdi][fwi];
            if(hits.offsetSize() < 2) continue;
            for(size_t i = 0; i < hits.offsetSize() - 1; i++) {
                BWTHit<index_t>& hit = hits.getPartialHit(i);
                for(size_t j = i + 1; j < hits.offsetSize(); j++) {
                    BWTHit<index_t>& hit2 = hits.getPartialHit(j);
                    if(hit._bwoff >= hit2._bwoff) {
                        hit._len = 0;
                        break;
                    }
                    if(hit._bwoff + hit._len <= hit2._bwoff) break;
                    if(hit._len >= hit2._len) {
                        index_t hit2_end = hit2._bwoff + hit2._len;
                        hit2._bwoff = hit._bwoff + hit._len;
                        hit2._len = hit2_end - hit2._bwoff;
                    } else {
                        hit._len = hit2._bwoff - hit._bwoff;
                    }
                }
            }
        }
    }
    
    pair<int, int> getForwardOrReverseHit(index_t rdi) {
        index_t avgHitLength[2] = {0, 0};
        index_t hitSize[2] = {0, 0} ;
        index_t maxHitLength[2] = {0, 0} ;
        for(index_t fwi = 0; fwi < 2; fwi++) {
            ReadBWTHit<index_t>& hit = this->_hits[rdi][fwi];
            index_t numHits = 0;
            index_t totalHitLength = 0;
#ifdef LI_DEBUG
            cout << fwi << ": ";
#endif
            for(size_t i = 0; i < hit.offsetSize(); i++) {
                index_t len = hit.getPartialHit(i).len();
#ifdef LI_DEBUG
                cout << len << " ";
#endif
                
                if(len < _minHitLen) continue;
                totalHitLength += (len - 15) * (len - 15);
                hitSize[fwi] += hit.getPartialHit(i).size();
                if(len > maxHitLength[fwi])
                    maxHitLength[fwi] = len;
                numHits++;
            }
#ifdef LI_DEBUG
            cout << endl;
#endif
            if(numHits > 0) {
                avgHitLength[fwi] = totalHitLength ; /// numHits;
            }
        }
        
        // choose read direction with a higher average hit length
        //cout<<"strand choosing: "<<avgHitLength[0]<<" "<<avgHitLength[1]<<endl ;
        index_t fwi;//= (avgHitLength[0] > avgHitLength[1])? 0 : 1;
        if(avgHitLength[0] != avgHitLength[1])
            fwi = (avgHitLength[0] > avgHitLength[1]) ? 0 : 1;
        else if(maxHitLength[0] != maxHitLength[1])
            fwi = (maxHitLength[0] > maxHitLength[1])? 0 : 1;
        // daehwan - for debugging purposes
#if 1
        else
            return pair<int, int>(0, 2);
#else
        else if(hitSize[0] != hitSize[1])
            fwi = (hitSize[0] > hitSize[1]) ? 1 : 0;
        else
            return 0;//just randomly pick one
#endif
        
        //return this->_hits[rdi][fwi];
        return pair<int, int>((int)fwi, (int)fwi + 1);
    }
    
    EList<Coord>& getCoords(
                            ReadBWTHit<index_t>& hit,
                            size_t hi,
                            const Ebwt<index_t>& ebwtFw,
                            const BitPairReference& ref,
                            RandomSource& rnd,
                            const index_t maxGenomeHitSize,
                            WalkMetrics& wlm,
                            PerReadMetrics& prm,
                            HIMetrics& him)
    {
        BWTHit<index_t>& partialHit = hit.getPartialHit(hi);
        assert(!partialHit.hasGenomeCoords());
        bool straddled = false;
        this->getGenomeIdx(
                           ebwtFw,     // FB: Why is it called ...FW here?
                           ref,
                           rnd,
                           partialHit._top,
                           partialHit._bot,
                           hit._fw == 0, // FIXME: fwi and hit._fw are defined differently
                           maxGenomeHitSize - this->_genomeHits.size(),
                           hit._len - partialHit._bwoff - partialHit._len,
                           partialHit._len,
                           partialHit._coords,
                           wlm,       // why is it called wlm here?
                           prm,
                           him,
                           false, // reject straddled
                           straddled);
#ifndef NDEBUG //FB
        std::cerr <<  partialHit.len() << ':';
#endif
        // get all coordinates of the hit
        return partialHit._coords;
    }


    // append a hit to genus map or update entry
    size_t addHitToHitMap(
                          EList<HitCount<index_t> >& hitMap,
                          int rdi,
                          int fwi,
                          uint64_t uniqueID,
                          uint64_t taxID,
                          size_t hi,
                          uint32_t partialHitScore,
                          double weightedHitLen,
                          bool considerOnlyIfPreviouslyObserved,
                          size_t offset,
                          size_t length)
    {
        size_t idx = 0;
#ifdef LI_DEBUG
        cout << "Add " << taxID << " " << partialHitScore << " " << weightedHitLen << endl;
#endif
        for(; idx < hitMap.size(); ++idx) {
            if(hitMap[idx].uniqueID == uniqueID) {
                if(hitMap[idx].timeStamp != hi) {
                    hitMap[idx].taxID = taxID;
                    hitMap[idx].count += 1;
                    hitMap[idx].scores[rdi][fwi] += partialHitScore;
                    hitMap[idx].summedHitLens[rdi][fwi] += weightedHitLen;
                    hitMap[idx].timeStamp = (uint32_t)hi;
                    hitMap[idx].readPositions.push_back(make_pair(offset, length));
                }
                break;
            }
        }
        
        if(idx >= hitMap.size() && !considerOnlyIfPreviouslyObserved) {
            hitMap.expand();
            hitMap.back().reset();
            hitMap.back().uniqueID = uniqueID;
            hitMap.back().taxID = taxID;
            hitMap.back().count = 1;
            hitMap[idx].scores[rdi][fwi] = partialHitScore;
            hitMap[idx].summedHitLens[rdi][fwi] = weightedHitLen;
            hitMap.back().timeStamp = (uint32_t)hi;
            hitMap.back().readPositions.clear();
            hitMap.back().readPositions.push_back(make_pair(offset, length));
        }

        //if considerOnlyIfPreviouslyObserved and it was not found, genus Idx size is equal to the genus Map size
        //assert_lt(genusIdx, genusMap.size());
        return idx;
    }

    // compare BWTHits by size, ascending, first, then by length, descending
    //   TODO: move this operator into BWTHits if that is the standard way we would like to sort
    //   TODO: this ordering does not necessarily give the best results
    struct compareBWTHits {
        bool operator()(const BWTHit<index_t>& a, const BWTHit<index_t>& b) const {
            if(a.len() >= 22 || b.len() >= 22) {
                if(a.len() >= 22 && b.len() >= 22) {
                    // sort ascending by size
                    if (a.size() < b.size()) return true;
                    if (a.size() > b.size()) return false;
                }
                
                // sort descending by length
                if (b.len() < a.len()) return true;
                if (b.len() > a.len()) return false;
            }
            
            // sort by the weighted len
            if(b.len() * a.size() < a.len() * b.size()) return true;
            if(b.len() * a.size() > a.len() * b.size()) return false;
            
            // sort ascending by size
            if(a.size() < b.size()) return true;
            if(a.size() > b.size()) return false;
            
            // sort descending by length
            if(b.len() < a.len()) return true;
            if(b.len() > a.len()) return false;
            
            return false;
        }
    };
};


#endif /*CLASSIFIER_H_*/
