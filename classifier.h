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

struct SpeciesCount {
    uint32_t id;
    uint32_t count;
    uint32_t score;
    double summedHitLen;  // sum of the length of all partial hits, divided by the number of genome matches
    uint32_t timeStamp;

	vector<pair<uint32_t,uint32_t> > readPositions;
    
    void reset() {
        id = count = score = timeStamp = 0;
        summedHitLen = 0.0;
	readPositions.clear();
    }
};

struct GenusCount {
    uint32_t id;
    uint32_t count;
    uint32_t score;
    double summedHitLen;
    uint32_t timeStamp;
    EList<SpeciesCount> speciesMap;
    
    void reset() {
        id = count = score = timeStamp = 0;
        speciesMap.clear();
        summedHitLen = 0.0;
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
		   SpeciesMetrics&          spm,
           RandomSource&            rnd,
           AlnSinkWrap<index_t>&    sink)
    {
        _genusMap.clear();

        const index_t increment = ( 2 * _minHitLen <= 33 ) ? 10 : ( 2 * _minHitLen - 33 ) ;
        //const index_t increment = 10 ;
        size_t bestScore = 0, secondBestScore = 0;
        const ReportingParams& rp = sink.reportingParams();
        index_t maxGenomeHitSize = rp.khits;
		bool isFw = false;


        // for each mate. only called once for unpaired data
        for(index_t rdi = 0; rdi < (this->_paired ? 2 : 1); rdi++) {
            assert(this->_rds[rdi] != NULL);

            // search for partial hits on the forward and reverse strand (saved in this->_hits[rdi])
            searchForwardAndReverse(rdi, ebwtFw, sc, rnd, increment);

            // get forward or reverse hits for this read from this->_hits[rdi]
            //  the strand is chosen based on higher average hit length in either direction
	    int fwi = getForwardOrReverseHit( rdi ) ;
	    if ( fwi == -1 )
	    	return 0 ;
            ReadBWTHit<index_t>& hit = this->_hits[rdi][fwi] ;
            assert(hit.done());
	    isFw = hit._fw;  // TODO: Sync between mates!

            // choose candidate partial alignments for further alignment
            index_t offsetSize = hit.offsetSize();
            this->_genomeHits.clear();
            
            // sort partial hits by size (number of genome positions), ascending, and then length, descending
            for(size_t hi = 0; hi < offsetSize; hi++) 
	    {
		const BWTHit<index_t> partialHit = hit.getPartialHit(hi);
#ifdef LI_DEBUG
		cout<<partialHit.len()<<" "<<partialHit.size()<<endl ;
#endif
		if ( partialHit.len() >= _minHitLen && partialHit.size() > maxGenomeHitSize )
		{
			maxGenomeHitSize = partialHit.size() ;
		}
	    } 
	    if ( maxGenomeHitSize > (index_t)rp.khits )
	    	maxGenomeHitSize += rp.khits ;

            hit._partialHits.sort(compareBWTHits());

            size_t usedPortion = 0 ;
            size_t genomeHitCnt = 0 ;
            for(size_t hi = 0; hi < offsetSize; hi++) {
				const BWTHit<index_t> partialHit = hit.getPartialHit(hi);
                size_t partialHitLen = partialHit.len();

                // only keep this partial hit if it is equal to or bigger than minHitLen (default: 22 bp)
                // TODO: consider not requiring minHitLen when we have already hits to the same genome
                bool considerOnlyIfPreviouslyObserved = partialHitLen < _minHitLen;

                // get all coordinates of the hit
                EList<Coord>& coords = getCoords(hit, hi, ebwtFw, ref, rnd,
                        maxGenomeHitSize, wlm, prm, him);
                if (coords.empty()) continue;

                usedPortion += partialHitLen;

                assert_gt(coords.size(), 0);

				// the maximum number of hits per read is maxGenomeHitSize (change with parameter -k)
				size_t nHitsToConsider = coords.size();

				if (genomeHitCnt + coords.size() > maxGenomeHitSize) {
                    coords.shufflePortion(0, coords.size(), rnd);
					nHitsToConsider = maxGenomeHitSize - genomeHitCnt;
                }
                
                // find the genome id for all coordinates, and count the number of genomes
                size_t n_genomes = coords.size();
                vector<uint64_t> coord_ids(nHitsToConsider);
                for(index_t k = 0; k < nHitsToConsider; k++, genomeHitCnt++) {
                    const Coord& coord = coords[k] ;

                    assert_lt(coord.ref(), _refnames.size()); // gives a warning - coord.ref() is signed integer. why?

                    // extract numeric id from refName
                    const string& refName = _refnames[coord.ref()];
                    uint64_t id = extractIDFromRefName(refName);

                    // count the genome if it is not in coord_ids, yet
                    if (k > 0 && std::find(coord_ids.begin(), coord_ids.begin()+k+1, id)!=coord_ids.begin()+k+1) {
                    	--n_genomes;
                    }

                    // add to coord_ids
                    coord_ids[k] = id;
                }

                // scoring function: calculate the weight of this partial hit
                assert_gt(partialHitLen, 15);
                assert_gt(n_genomes, 0);
                uint32_t partialHitScore = (uint32_t)((partialHitLen - 15) * (partialHitLen - 15)) / n_genomes ;
                double weightedHitLen = double(partialHitLen) / double(n_genomes) ;

                // go through all coordinates reported for partial hit
                for(index_t k = 0; k < nHitsToConsider; ++k) {
                    uint64_t id = coord_ids[k];

                    uint32_t speciesID = (uint32_t)(id >> 32);
                    uint32_t genusID = (uint32_t)(id & 0xffffffff);

                    // add hit to genus map and get new index in the map
                    size_t genusIdx = addHitToGenusMap(_genusMap, genusID, hi, partialHitScore, weightedHitLen, considerOnlyIfPreviouslyObserved);

                    //if considerOnlyIfPreviouslyObserved and it was not found, genus Idx size is equal to the genus Map size
                    if (genusIdx >= _genusMap.size()) {
                    	continue;
                    }

                    // add hit to species map and get new score for the species
                    uint32_t newScore = addHitToSpeciesMap(_genusMap[genusIdx],speciesID,
                    		hi, partialHitScore,weightedHitLen, partialHit._bwoff, partialHit.len());

#ifndef NDEBUG //FB
                    std::cerr << speciesID << ';';
#endif

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
#if 0
                bool last_round = rdi == (this->_paired ? 1 : 0);
                if (last_round && bestScore > secondBestScore + (totalHitLength[hit._fw == 0] - usedPortion - 15) * (totalHitLength[hit._fw == 0] - usedPortion - 15)) {
                    break ;
                }
#endif

#ifndef NDEBUG //FB
            std::cerr << "  partialHits-done";
#endif
            } // partialHits
#ifndef NDEBUG //FB
            std::cerr << "  rdi-done" << endl;
#endif
        } // rdi
        
#if 1

        for(size_t gi = 0; gi < _genusMap.size(); gi++) {
            assert_gt(_genusMap[gi].score, 0);
            GenusCount& genusCount = _genusMap[gi];
            uint32_t speciesID = 0xffffffff;
            for(size_t mi = 0; mi < genusCount.speciesMap.size(); mi++) {
                speciesID = genusCount.speciesMap[mi].id; assert_neq(speciesID, 0xffffffff);
                uint32_t speciesScore = genusCount.speciesMap[mi].score;
                
                // report
                AlnScore asc(genusCount.score + speciesScore);
                AlnRes rs;
                rs.init(asc,
                        speciesID,
                        genusCount.id,
			genusCount.speciesMap[mi].summedHitLen,
			genusCount.speciesMap[mi].readPositions,
			isFw);
                sink.report(0, &rs);
            }
        }
#else
        _tempTies.clear();
        uint32_t genusScore = 0;
        for(size_t mi = 0; mi < _genusMap.size(); mi++) {
            assert_gt(_genusMap[mi].score, 0);
            if(_genusMap[mi].score > genusScore) {
                genusScore = _genusMap[mi].score;
                _tempTies.clear();
                _tempTies.push_back((uint16_t)mi);
            } else if(_genusMap[mi].score == genusScore) {
                _tempTies.push_back((uint16_t)mi);
            }
        }
        
        if(_tempTies.size() > 0) {
            for(uint16_t gi = 0; gi < _tempTies.size(); gi++) {
                assert_lt(_tempTies[gi], _genusMap.size());
                GenusCount& genusCount = _genusMap[_tempTies[gi]];
                uint32_t speciesScore = 0;
                uint32_t speciesID = 0xffffffff;
                for(size_t mi = 0; mi < genusCount.speciesMap.size(); mi++) {
                    if(genusCount.speciesMap[mi].score > speciesScore) {
                        speciesID = genusCount.speciesMap[mi].id;
                        speciesScore = genusCount.speciesMap[mi].score;
                    }
                }
                assert_neq(speciesID, 0xffffffff);
                
                // report
                AlnScore asc(genusCount.score + speciesScore);
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
        //const size_t maxDiff = (rdlen / 2 > 2 * _minHitLen) ? rdlen / 2 : (2 * _minHitLen);
	size_t sum[2] = {0, 0} ;

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
                
		BWTHit<index_t>& lastHit = hit.getPartialHit(
                        hit.offsetSize() - 1);
                if (hit.done()) {
                    done[fwi] = true;
                    cur[fwi] = rdlen;
		    if ( lastHit.len() >= _minHitLen )
		    {
			    sum[ fwi ] += lastHit.len() ;
			    if ( 0 ) //lastHit.len() < 31 && rdlen > 31 && lastHit.size() == 1 )
			    {
				    ReadBWTHit<index_t> testHit ;
				    testHit.init( fw, rdlen ) ;
				    testHit.setOffset( hit.cur() - 1 - 31 + 1 ) ; 

				    this->partialSearch(ebwtFw, rd, sc, fw, 0, mineFw, mineRc, testHit,
						    rnd, pseudogeneStop, anchorStop) ;
				    index_t tmpLen = testHit.getPartialHit( testHit.offsetSize() - 1 ).len() ;
#ifdef LI_DEBUG
				    cout<<"(adjust: "<<tmpLen<<")" ;
#endif
				    if ( tmpLen >= 31 )
				    {
					    lastHit._len = tmpLen ;
				    }
			    }
		    }

                    continue;
                }

		cur[fwi] = hit.cur();
#ifdef LI_DEBUG
		cout<<fwi<<":"<<lastHit.len()<<" "<<cur[fwi]<<" " ;
#endif

		if ( lastHit.len() >= _minHitLen )
			sum[ fwi ] += lastHit.len() ;

		if (lastHit.len() > increment) {
                    if ( lastHit.len() < _minHitLen)
                        hit.setOffset(hit.cur() - increment);
                    else 
		    {
                        hit.setOffset(hit.cur() + 1);
			if ( 0 ) //lastHit.len() < 31 && hit.cur() >= 31 && lastHit.size() == 1 )
			{
				ReadBWTHit<index_t> testHit ;
				testHit.init( fw, rdlen ) ;
				testHit.setOffset( hit.cur() - 1 - 31 ) ; // why not hit.cur() - 1 - 31 + 1? because we "+1" before the if!

				this->partialSearch(ebwtFw, rd, sc, fw, 0, mineFw, mineRc, testHit,
						rnd, pseudogeneStop, anchorStop) ;
				index_t tmpLen = testHit.getPartialHit( testHit.offsetSize() - 1 ).len() ;
#ifdef LI_DEBUG
				cout<<"(adjust: "<<tmpLen<<")" ;
#endif
				if ( tmpLen >= 31 )
				{
					lastHit._len = tmpLen ;
				}
			}
		    }
                }
                if (hit.cur() + _minHitLen >= rdlen) {
                    hit.done(true);
                    done[fwi] = true;
                    continue;
                }

		if ( lastHit.len() <= 3 )
		{
			// This happens most likely due to the Ns in the read
			--fwi ; // Repeat this strand again.
		}
            }
#ifdef LI_DEBUG
	    cout<<endl ;
#endif	    
            if ( sum[0] > sum[1] + ( rdlen - cur[1] + 1 ) ) {
                this->_hits[rdi][1].done(true);
                done[1] = true;
            } else if ( sum[1] > sum[0] + ( rdlen - cur[0] + 1) ) {
                this->_hits[rdi][0].done(true);
                done[0] = true;
            }
        }
    }

    int getForwardOrReverseHit(index_t rdi ) {
        index_t avgHitLength[2] = {0, 0};
	index_t hitSize[2] = {0, 0} ;
	index_t maxHitLength[2] = {0, 0} ;
        for(index_t fwi = 0; fwi < 2; fwi++) {
            ReadBWTHit<index_t>& hit = this->_hits[rdi][fwi];
            index_t numHits = 0;
	    index_t totalHitLength = 0 ;
#ifdef LI_DEBUG
	    cout<<fwi<<": " ;
#endif	    
            for(size_t i = 0; i < hit.offsetSize(); i++) {
	    	index_t len = hit.getPartialHit(i).len() ;
#ifdef LI_DEBUG
	    	cout<<hit.getPartialHit(i).len()<<" " ;
#endif

                if(len < _minHitLen) continue;
                totalHitLength += ( len - 15 ) * ( len - 15 ) ;
		hitSize[fwi] += hit.getPartialHit(i).size() ;
		if ( len > maxHitLength[ fwi ] )
			maxHitLength[fwi] = len ;
                numHits++;
            }
#ifdef LI_DEBUG
	    cout<<endl ;
#endif
            if(numHits > 0) {
                avgHitLength[fwi] = totalHitLength ;/// numHits;
            }
        }

        // choose read direction with a higher average hit length
	//cout<<"strand choosing: "<<avgHitLength[0]<<" "<<avgHitLength[1]<<endl ;
        index_t fwi ;//= (avgHitLength[0] > avgHitLength[1])? 0 : 1;
	if ( avgHitLength[0] != avgHitLength[1] )
        	fwi = (avgHitLength[0] > avgHitLength[1])? 0 : 1;
	else if ( maxHitLength[0] != maxHitLength[1] )
        	fwi = (maxHitLength[0] > maxHitLength[1])? 0 : 1;
	else if ( hitSize[0] != hitSize[1] )
		fwi = (hitSize[0] > hitSize[1])? 1 : 0;
	else
		return 0 ;//just randomly pick one
		
        //return this->_hits[rdi][fwi];
	return fwi ;
    }

    EList<Coord>& getCoords(ReadBWTHit<index_t>& hit, size_t hi,
            const Ebwt<index_t>& ebwtFw, const BitPairReference& ref, RandomSource& rnd, const index_t maxGenomeHitSize,
            WalkMetrics& wlm, PerReadMetrics& prm, HIMetrics& him) {

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
    size_t addHitToGenusMap(EList<GenusCount> &genusMap, uint32_t genusID, size_t hi, uint32_t partialHitScore, double weightedHitLen, bool considerOnlyIfPreviouslyObserved) {
        size_t genusIdx = 0;
#ifdef LI_DEBUG
	cout<<"Add "<<genusID<<" "<<weightedHitLen<<endl ; 
#endif
        for(; genusIdx < genusMap.size(); ++genusIdx) {
            if(genusMap[genusIdx].id == genusID) {
                if(genusMap[genusIdx].timeStamp != hi) {
                    genusMap[genusIdx].count += 1;
                    genusMap[genusIdx].score += partialHitScore;
                    //genusMap[genusIdx].summedHitLen += weightedHitLen;
                    genusMap[genusIdx].timeStamp = hi;
                }
                break;
            }
        }

        if(genusIdx >= genusMap.size() && !considerOnlyIfPreviouslyObserved) {
            genusMap.expand();
            genusMap.back().reset();
            genusMap.back().id = genusID;
            genusMap.back().count = 1;
            genusMap.back().score = partialHitScore;
            genusMap[genusIdx].summedHitLen = 0 ; //weightedHitLen;
            genusMap.back().timeStamp = hi;
        }

        //if considerOnlyIfPreviouslyObserved and it was not found, genus Idx size is equal to the genus Map size
        //assert_lt(genusIdx, genusMap.size());
        return genusIdx;
    }

    // append a hit to species map or update entry
    uint32_t addHitToSpeciesMap(GenusCount& genusCount, uint32_t speciesID, size_t hi,
    		uint32_t partialHitScore, double weightedHitLen, size_t offset, size_t length) {

        uint32_t newScore = 0;
        bool found = false;

        for(size_t mi = 0; mi < genusCount.speciesMap.size(); mi++) {
            if(genusCount.speciesMap[mi].id == speciesID) {
                found = true;
                if(genusCount.speciesMap[mi].timeStamp != hi) {
                    genusCount.speciesMap[mi].count += 1;
                    genusCount.speciesMap[mi].score += partialHitScore;
                    genusCount.speciesMap[mi].summedHitLen += weightedHitLen;
                    genusCount.speciesMap[mi].timeStamp = hi;
                    genusCount.speciesMap[mi].readPositions.push_back(make_pair(offset,length));
                    newScore = genusCount.speciesMap[mi].score;

		    genusCount.summedHitLen += weightedHitLen ;
                } else {
    #ifndef NDEBUG //FB
                        std::cerr <<  "same hit!!";
    #endif
                }
                break;
            }
        }

        if(!found) {
            genusCount.speciesMap.expand();
            genusCount.speciesMap.back().reset();
            genusCount.speciesMap.back().id = speciesID;
            genusCount.speciesMap.back().count = 1;
            newScore = genusCount.speciesMap.back().score = partialHitScore;
            genusCount.speciesMap.back().summedHitLen = weightedHitLen;
            genusCount.speciesMap.back().readPositions = vector<pair<uint32_t,uint32_t> >(1);
            genusCount.speciesMap.back().readPositions[0] = make_pair(offset,length);
            genusCount.speciesMap.back().timeStamp = hi;

	    genusCount.summedHitLen += weightedHitLen ;
       }
        return newScore;
    }



    // compare BWTHits by size, ascending, first, then by length, descending
    //   TODO: move this operator into BWTHits if that is the standard way we would like to sort
    //   TODO: this ordering does not necessarily give the best results
        struct compareBWTHits {
        bool operator()(const BWTHit<index_t>& a, const BWTHit<index_t>& b) const {
	    if ( a.len() >= 22 || b.len() >= 22 )
	    {
	    	    if ( a.len() >= 22 && b.len() >= 22 )
		    {
			    // sort ascending by size
			    if (a.size() < b.size()) return true;
			    if (a.size() > b.size()) return false;
		    }
		    
		    // sort descending by length
		    if (b.len() < a.len()) return true;
		    if (b.len() > a.len()) return false;
	    }

	    // sort by the weighted len
	    if (b.len() * a.size() < a.len() * b.size() ) return true;
            if (b.len() * a.size() > a.len() * b.size() ) return false;

            // sort ascending by size
            if (a.size() < b.size()) return true;
            if (a.size() > b.size()) return false;

	    // sort descending by length
	    if (b.len() < a.len()) return true;
            if (b.len() > a.len()) return false;

                        return false;
        }
    };



};


#endif /*CLASSIFIER_H_*/
