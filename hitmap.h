/*
 * hitmap.h
 *
 *  Created on: Jan 19, 2017
 *      Author: fbreitwieser
 */

#ifndef HITMAP_H_
#define HITMAP_H_

#include "taxonomy.h"

typedef std::array<uint32_t,2> HitScore;  // score for one hit: has hit length and score

template<typename index_t>
struct HitCount {
    uint64_t uniqueID;
    TaxId taxID;

    uint32_t count;
    HitScore best_score;
    HitScore scores[2][2]; // scores[rdi][fwi]
    uint32_t timeStamp;
    EList<pair<uint32_t,uint32_t> > readPositions;
    bool     leaf;
    uint32_t num_leaves;

    uint8_t rank;
    EList<uint64_t> path;

    void reset() {
    	best_score = {0,0};
        uniqueID = taxID = count = timeStamp = 0;
        scores[0][0] = scores[0][1] = scores[1][0] = scores[1][1] = {0,0};
        readPositions.clear();
        rank = 0;
        path.clear();
        leaf = true;
        num_leaves = 1;
    }

    HitCount& operator=(const HitCount& o) {
        if(this == &o)
            return *this;

        uniqueID = o.uniqueID;
        taxID = o.taxID;
        count = o.count;
        best_score = o.best_score;
        scores[0][0] = o.scores[0][0];
        scores[0][1] = o.scores[0][1];
        scores[1][0] = o.scores[1][0];
        scores[1][1] = o.scores[1][1];
        timeStamp = o.timeStamp;
        readPositions = o.readPositions;
        leaf = o.leaf;
        num_leaves = o.num_leaves;
        rank = o.rank;
        path = o.path;

        return *this;
    }

    HitScore finalize(
                  bool paired,
                  bool mate1fw,
                  bool mate2fw) {
        if(paired) {
#if 1
            best_score = max(scores[0][0], scores[0][1]) + max(scores[1][0], scores[1][1]);
#else
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
#endif
        } else {
            best_score = max(scores[0][0], scores[0][1]);
        }
        return best_score;
    }
};



#endif /* HITMAP_H_ */
