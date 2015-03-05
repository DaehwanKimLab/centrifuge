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

#include "hi_aligner.h"

struct SpeciesCount {
	uint32_t id;
	uint32_t count;
	uint32_t weightedCount;
    uint32_t timeStamp;
    
    void reset() {
        id = count = weightedCount = timeStamp = 0;
    }
};

struct GenusCount {
    uint32_t id;
	uint32_t count;
	uint32_t weightedCount;
    uint32_t timeStamp;
    EList<SpeciesCount> speciesMap;
    
    void reset() {
        id = count = weightedCount = timeStamp = 0;
        speciesMap.clear();
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
               index_t minHitLen = 22) :
    HI_Aligner<index_t, local_index_t>(
                                       ebwt,
                                       0,    // don't make use of splice sites found by earlier reads
                                       true), // no spliced alignment
    _refnames(refnames),
    _minHitLen(minHitLen)
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
           RandomSource&            rnd,
           AlnSinkWrap<index_t>&    sink)
    {
        _genusMap.clear();

        const index_t increment = 10;
        size_t bestScore = 0, secondBestScore = 0;
        const ReportingParams& rp = sink.reportingParams();
        const index_t maxGenomeHitSize = rp.khits;


        // for each mate. only called once for unpaired data
        for(index_t rdi = 0; rdi < (this->_paired ? 2 : 1); rdi++) {
            assert(this->_rds[rdi] != NULL);

            // search for partial hits on the forward and reverse strand (saved in this->_hits[rdi])
            searchForwardAndReverse(rdi, ebwtFw, sc, rnd, increment);
            
            // get forward or reverse hits for this read from this->_hits[rdi]
            //  the strand is chosen based on higher average hit length in either direction
    		index_t totalHitLength[2] = {0, 0};
			ReadBWTHit<index_t>& hit = getForwardOrReverseHit(rdi, totalHitLength);
            assert(hit.done());

            // choose candidate partial alignments for further alignment
            index_t offsetSize = hit.offsetSize();
            this->_genomeHits.clear();
            
            int hitLen[100];
            int hitSize[100];
            size_t hiMap[100];
            for(size_t hi = 0; hi < offsetSize; hi++) {
                hitLen[hi] = hit.getPartialHit(hi).len();
                hitSize[hi] = hit.getPartialHit(hi).size();
                hiMap[hi] = hi ;
            }
            
            // Change to quicksort in future
            for(size_t hi = 0; hi < offsetSize; hi++) {
                for(size_t hj = hi + 1; hj < offsetSize; hj++) {
                    if(hitSize[ hiMap[hi] ] > hitSize[ hiMap[hj] ]) { // When use size()
                        size_t tmp = hiMap[hi];
                        hiMap[hi] = hiMap[hj];
                        hiMap[hj] = tmp;
                    } else if(hitSize[ hiMap[hi] ] == hitSize[ hiMap[hj] ] && hitLen[ hiMap[hi] ] < hitLen[ hiMap[hj] ]) {
                        size_t tmp = hiMap[hi];
                        hiMap[hi] = hiMap[hj];
                        hiMap[hj] = tmp;
                    }
                }
            }
            
            size_t usedPortion = 0 ;
            size_t genomeHitCnt = 0 ;
            for(size_t hi = 0; hi < offsetSize; hi++) {
                BWTHit<index_t>& partialHit = hit.getPartialHit(hiMap[ hi ]);
                if(partialHit.len() < _minHitLen)
                    continue;
                assert(!partialHit.hasGenomeCoords());
                usedPortion += partialHit.len();
                bool straddled = false;
                this->getGenomeIdx(
                                   ebwtFw,
                                   ref,
                                   rnd,
                                   partialHit._top,
                                   partialHit._bot,
                                   hit._fw == 0,
                                   maxGenomeHitSize - this->_genomeHits.size(),
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
                if(genomeHitCnt + coords.size() > maxGenomeHitSize) {
                    coords.shufflePortion(0, coords.size(), rnd);
                }
                
                for(index_t k = 0; k < coords.size() && genomeHitCnt < maxGenomeHitSize; k++, genomeHitCnt++) {
                    const Coord& coord = coords[k] ;
                    assert_lt(coord.ref(), _refnames.size());
                    const string& refName = _refnames[coord.ref()];
                    uint64_t id = 0;
                    for(size_t ni = 0; ni < refName.length(); ni++) {
                        if(refName[ni] < '0' || refName[ni] > '9') break;
                        id *= 10;
                        id += (refName[ni] - '0');
                    }
                
                    uint32_t speciesID = (uint32_t)(id >> 32);
                    uint32_t genusID = (uint32_t)(id & 0xffffffff);

		    // florian - scoring function
                    uint32_t addWeight = (uint32_t)((partialHit.len() - 15) * (partialHit.len() - 15));
                    
                    uint32_t newScore = 0;
                    size_t genusIdx = 0;
                    for(; genusIdx < _genusMap.size(); genusIdx++) {
                        if(_genusMap[genusIdx].id == genusID) {
                            if(_genusMap[genusIdx].timeStamp != hi) {
                                _genusMap[genusIdx].count += 1;
                                _genusMap[genusIdx].weightedCount += addWeight;
                                _genusMap[genusIdx].timeStamp = hi;
                            }
                            break;
                        }
                    }
                    
                    if(genusIdx >= _genusMap.size()) {
                        _genusMap.expand();
                        _genusMap.back().reset();
                        _genusMap.back().id = genusID;
                        _genusMap.back().count = 1;
                        _genusMap.back().weightedCount = addWeight;
                        _genusMap.back().timeStamp = hi;
                    }
                    
                    assert_lt(genusIdx, _genusMap.size());
                    GenusCount& genusCount = _genusMap[genusIdx];
                    bool found = false;
                    for(size_t mi = 0; mi < genusCount.speciesMap.size(); mi++) {
                        if(genusCount.speciesMap[mi].id == speciesID) {
                            found = true;
                            if(genusCount.speciesMap[mi].timeStamp != hi) {
                                genusCount.speciesMap[mi].count += 1;
                                genusCount.speciesMap[mi].weightedCount += addWeight;
                                genusCount.speciesMap[mi].timeStamp = hi;
                                newScore = _genusMap[genusIdx].weightedCount;
                            }
                            break;
                        }
                    }
                    
                    if(!found) {
                        genusCount.speciesMap.expand();
                        genusCount.speciesMap.back().reset();
                        genusCount.speciesMap.back().id = speciesID;
                        genusCount.speciesMap.back().count = 1;
                        newScore = genusCount.speciesMap.back().weightedCount = addWeight;
                        genusCount.speciesMap.back().timeStamp = hi;
                    }
                
                    if(newScore > bestScore) {
                        secondBestScore = bestScore;
                        bestScore = newScore;
                    } else if(newScore > secondBestScore) {
                        secondBestScore = newScore;
                    }
                }
                
                if(genomeHitCnt >= maxGenomeHitSize)
                    break;
                
                // FIXME FB: Benchmark the effect of this statement
                bool last_round = rdi == (this->_paired ? 1 : 0);
                if (last_round && bestScore > secondBestScore + (totalHitLength[hit._fw == 0] - usedPortion - 15) * (totalHitLength[hit._fw == 0] - usedPortion - 15)) {
                    break ;
                }
            } // offsetSize
        } // rdi
        
#if 1
        for(size_t gi = 0; gi < _genusMap.size(); gi++) {
            assert_gt(_genusMap[gi].weightedCount, 0);
            GenusCount& genusCount = _genusMap[gi];
            uint32_t speciesWeightedCount = 0;
            uint32_t speciesID = 0xffffffff;
            for(size_t mi = 0; mi < genusCount.speciesMap.size(); mi++) {
                speciesID = genusCount.speciesMap[mi].id; assert_neq(speciesID, 0xffffffff);
                speciesWeightedCount = genusCount.speciesMap[mi].weightedCount;
                
                // report
                AlnScore asc(genusCount.weightedCount + speciesWeightedCount);
                AlnRes rs;
                rs.init(asc,
                        speciesID,
                        genusCount.id);
                sink.report(0, &rs);
            }
        }
#else
        _tempTies.clear();
        uint32_t genusWeightedCount = 0;
        for(size_t mi = 0; mi < _genusMap.size(); mi++) {
            assert_gt(_genusMap[mi].weightedCount, 0);
            if(_genusMap[mi].weightedCount > genusWeightedCount) {
                genusWeightedCount = _genusMap[mi].weightedCount;
                _tempTies.clear();
                _tempTies.push_back((uint16_t)mi);
            } else if(_genusMap[mi].weightedCount == genusWeightedCount) {
                _tempTies.push_back((uint16_t)mi);
            }
        }
        
        if(_tempTies.size() > 0) {
            for(uint16_t gi = 0; gi < _tempTies.size(); gi++) {
                assert_lt(_tempTies[gi], _genusMap.size());
                GenusCount& genusCount = _genusMap[_tempTies[gi]];
                uint32_t speciesWeightedCount = 0;
                uint32_t speciesID = 0xffffffff;
                for(size_t mi = 0; mi < genusCount.speciesMap.size(); mi++) {
                    if(genusCount.speciesMap[mi].weightedCount > speciesWeightedCount) {
                        speciesID = genusCount.speciesMap[mi].id;
                        speciesWeightedCount = genusCount.speciesMap[mi].weightedCount;
                    }
                }
                
                assert_neq(speciesID, 0xffffffff);
                
                // report
                AlnScore asc(genusCount.weightedCount + speciesWeightedCount);
                AlnRes rs;
                rs.init(asc,
                        speciesID,
                        genusCount.id);
                sink.report(0, &rs);
            }
        }
#endif
        
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
    EList<string>      _refnames;
    EList<GenusCount>  _genusMap;
    index_t            _minHitLen;
    
    EList<uint16_t>    _tempTies;

	void searchForwardAndReverse(index_t rdi,
			const Ebwt<index_t>& ebwtFw, const Scoring& sc,
			RandomSource& rnd, const index_t increment) {

		const Read& rd = *(this->_rds[rdi]);

		bool done[2] = {false, false};
        size_t cur[2] = {0, 0} ;

        index_t rdlen = rd.length();
        const size_t maxDiff = (rdlen / 2 > 2 * _minHitLen) ? rdlen / 2 : (2 * _minHitLen);

		// search for partial hits on the forward and reverse strand
		while (!done[0] || !done[1]) {
			for (index_t fwi = 0; fwi < 2; fwi++) {
				if (done[fwi])
					continue;

				size_t mineFw = 0, mineRc = 0;
				bool fw = (fwi == 0);
				ReadBWTHit<index_t>& hit = this->_hits[rdi][fwi];
				bool pseudogeneStop = false, anchorStop = false;
				this->partialSearch(ebwtFw, rd, sc, fw, 0, mineFw, mineRc, hit,
						rnd, pseudogeneStop, anchorStop);
				if (hit.done()) {
					done[fwi] = true;
					cur[fwi] = rdlen;
					continue;
				}
				cur[fwi] = hit.cur();
				BWTHit<index_t>& lastHit = hit.getPartialHit(
						hit.offsetSize() - 1);
				if (lastHit.len() > increment) {
					if (lastHit.len() < _minHitLen)
						hit.setOffset(hit.cur() - increment);
					else
						hit.setOffset(hit.cur() + 1);
				}
				if (hit.cur() + _minHitLen >= rdlen) {
					hit.done(true);
					done[fwi] = true;
					continue;
				}
			}
			if (cur[0] > cur[1] + maxDiff) {
				this->_hits[rdi][1].done(true);
				done[1] = true;
			} else if (cur[1] > cur[0] + maxDiff) {
				this->_hits[rdi][0].done(true);
				done[0] = true;
			}
		}
	}

	ReadBWTHit<index_t>& getForwardOrReverseHit(index_t rdi, index_t (&totalHitLength)[2]) {
		index_t avgHitLength[2] = {0, 0};
		for(index_t fwi = 0; fwi < 2; fwi++) {
			ReadBWTHit<index_t>& hit = this->_hits[rdi][fwi];
			index_t numHits = 0;
			for(size_t i = 0; i < hit.offsetSize(); i++) {
				if(hit.getPartialHit(i).len() < _minHitLen) continue;
				totalHitLength[fwi] += hit.getPartialHit(i).len();
				numHits++;
			}

			if(numHits > 0) {
				avgHitLength[fwi] = totalHitLength[fwi] / numHits;
			}
		}

		// choose read direction with a higher average hit length
		index_t fwi = (avgHitLength[0] > avgHitLength[1])? 0 : 1;
		return this->_hits[rdi][fwi];
	}
};


#endif /*CLASSIFIER_H_*/
