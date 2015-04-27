/*
 * hyperloglogplus.h
 *
 * Implementation of HyperLogLog++ algorithm described by Stefan Heule et al.
 *
 *  Created on: Apr 25, 2015
 *      Author: fbreitwieser
 */

#ifndef HYPERLOGLOGPLUS_H_
#define HYPERLOGLOGPLUS_H_

#include<set>
#include<list>
#include<vector>
#include<stdexcept>
#include<iostream>
#include<fstream>
#include<math.h>   //log
#include<bitset>

#include "hyperloglogbias.h"
#include "third_party/MurmurHash3.h"
using namespace std;

#define NDEBUG
#define NDEBUG2
#define arr_len(a) (a + sizeof a / sizeof a[0])

// experimentally determined threshold values for  p - 4
const uint32_t threshold[] = {10, 20, 40, 80, 220, 400, 900, 1800, 3100,
							  6500, 11500, 20000, 50000, 120000, 350000};


///////////////////////

//
/**
 * gives the estimated cardinality for m bins, v of which are non-zero
 * @param m
 * @param v
 * @return
 */
double linearCounting(uint32_t m, uint32_t v) {
#ifdef DEBUG1
	cerr << "linear counting: m=" << m << ", v=" << v << endl;
#endif
	if (v > m) {
		cerr << "v greater then m!!" << endl;
	}
	double fm = double(m);
	return fm * log(fm/double(v));
}

/**
 * Bias correction factors for specific m's
 * @param m
 * @return
 */
double getAlpha(uint32_t m)  {
	switch (m) {
	case 16: return 0.673;
	case 32: return 0.697;
	case 64: return 0.709;
	}

	// m >= 128
	return 0.7213 / (1 + 1.079/double(m));
}

/**
 * calculate the raw estimate as harmonic mean of the ranks in the register
 * @param s
 * @return
 */
double calculateEstimate(vector<uint8_t> M, uint32_t m) {
	double sum = 0.0;
	for (size_t i = 0; i <= M.size(); ++i) {
		sum += 1.0 / (double(uint32_t(1)<<M[i]));
	}

	double dm = double(M.size());
	return getAlpha(m) * dm * dm / sum;
}

uint32_t countZeros(vector<uint8_t> s) {
	return count(s.begin(), s.end(), 0);
}

/**
 * Extract bits (from uint32_t or uint64_t) using LSB 0 numbering from hi to lo, including lo
 * @param bits
 * @param hi
 * @param lo
 * @return
 */
template<typename T>
T extractBits(T bits, uint8_t hi, uint8_t lo) {

    // create a bitmask:
    //            (T(1) << (hi - lo)                 a 1 at the position (hi - lo)
    //           ((T(1) << (hi - lo) - 1)              1's from position 0 to position (hi-lo-1)
    //          (((T(1) << (hi - lo)) - 1) << lo)      1's from position lo to position hi

	// The T(1) is required to not cause overflow on 32bit machines
	// TODO: consider creating a bitmask only once in the beginning
	T bitmask = (((T(1) << (hi - lo)) - 1) << lo);
	
	return ((bits & bitmask) >> lo);
}

template<typename T>
T extractBits(T bits, uint8_t hi) {
    // create a bitmask for first hi bits (LSB 0 numbering)
	T bitmask = (((T(1) << hi) - 1));

	return (bits & bitmask);
}

// functions for counting the number of leading 0-bits (clz)
//           and counting the number of trailing 0-bits (ctz)
//#ifdef __GNUC__

// TODO: switch between builtin clz and 64_clz based on architecture
//#define clz(x) __builtin_clz(x)
static int clz_manual(uint64_t x)
{
  // This uses a binary search (counting down) algorithm from Hacker's Delight.
   uint64_t y;
   int n = 64;
   y = x >>32;  if (y != 0) {n -= 32;  x = y;}
   y = x >>16;  if (y != 0) {n -= 16;  x = y;}
   y = x >> 8;  if (y != 0) {n -=  8;  x = y;}
   y = x >> 4;  if (y != 0) {n -=  4;  x = y;}
   y = x >> 2;  if (y != 0) {n -=  2;  x = y;}
   y = x >> 1;  if (y != 0) return n - 2;
   return n - x;
}


uint32_t clz(const uint32_t x) {
	return __builtin_clz(x);
}

uint32_t clz(const uint64_t x) {
    uint32_t u32 = (x >> 32);
    uint32_t result = u32 ? __builtin_clz(u32) : 32;
    if (result == 32) {
        u32 = x & 0xFFFFFFFFUL;
        result += (u32 ? __builtin_clz(u32) : 32);
    }
    return result;
}
//#else

uint32_t clz_log2(const uint64_t w) {
	return 63 - floor(log2(w));
}
//#endif


// TODO: the sparse list may be encoded with variable length encoding
//   see Heule et al., section 5.3.2
// Also, using sets might give a larger overhead as each insertion costs more
//  consider using vector and sort/unique when merging.
typedef set<uint32_t> SparseListType;
typedef set<uint32_t> TmpSetType;
typedef uint64_t HashSize;

/**
 * HyperLogLogPlus class
 * typename T corresponds to the hash size - usually either uint32_t or uint64_t (implemented for uint64_t)
 */
template <typename T>
class HyperLogLogPlus {

	vector<uint8_t> M;  // registers (M) of size m
	uint8_t p;            // precision
	uint32_t m;           // number of registers
	bool sparse;          // sparse representation of the data?
	TmpSetType tmpSet;
	SparseListType sparseList; // TODO: use a compressed list instead

	vector<vector<double> > rawEstimateData;
	vector<vector<double> > biasData;

	// sparse versions of p and m
	static const uint8_t  pPrime = 25; // precision when using a sparse representation
	static const uint32_t mPrime = 1 << (pPrime -1); // 2^pPrime

public:

	/**
	 * Create new HyperLogLogPlus counter
	 * @param precision
	 * @param sparse
	 */
	HyperLogLogPlus(uint8_t precision, bool sparse=true) {
		if (precision > 18 || precision < 4) {
	        throw std::invalid_argument("precision (number of register = 2^precision) must be between 4 and 18");
		}

		this->p = precision;
		this->m = 1 << precision;
		this->M = vector<uint8_t>(m);
		this->sparse = sparse;
		this->tmpSet = TmpSetType();
		this->sparseList = SparseListType(); // TODO: if SparseListType is changed, initialize with appropriate size
	}

	/**
	 * Aggregation of items. Adds a new item to the counter.
	 * @param item
	 */
	void add(T item) {
		add(item, sizeof(item));
	}

	/**
	 * Aggregation of items. Adds a new item to the counter.
	 * @param item
	 * @param size  size of item
	 */
	void add(T item, size_t size) {

		// compute hash for item
	    HashSize out[2];
	    uint32_t seed = 0;
	    MurmurHash3_x64_128((void *)item, size, seed, &out);
	    HashSize x = out[0];
	
		uint32_t y = encodeHash(x);
		idx_n_rank ir = decodeHash(y);
#ifdef DEBUG2
		cerr << bitset<64>(x) << " ~> " << bitset<32>(y) << " --> ["<< uint32_t(ir.idx) << "]["<< uint32_t(ir.rank)<<":"<<bitset<8>(ir.rank)<<"]" << endl;
#endif

		if (sparse) {
			// sparse mode: put the encoded hash into tmpSet
			tmpSet.insert(encodeHash(x));

			// if the temporary set is too large, merge it with the sparseList
			if (tmpSet.size()*100 > m) {
				mergeSparse();

				// if the sparseList is too large, switch to normal (register) representation
				if (sparseList.size() > m) {
					toNormal();
				}
			}
		} else {
			// normal mode
			uint32_t idx = (uint32_t)extractBits(x, 64, 64-p); // get index: {x63,...,x64-p}
			uint8_t rank = (uint8_t)clz(extractBits(x, 64-p, 1)) - p;         // get rank of w:  {x63-p,...,x0}

#ifdef DEBUG2
			 cerr << bitset<64>(x) << " --> ["<< uint32_t(idx) << "]["<< uint32_t(rank)<<"]";
#endif

			if (rank > M[idx]) {
#ifdef DEBUG2
				 cerr << "i" << endl;
#endif
				M[idx] = rank;
			} else {
#ifdef DEBUG2
				 cerr << "N" << endl;
#endif
			}
		}
		//cerr << "Added item, new cardinality: " << cardinality() << endl;
	}


	/**
	 *	Merge tmpSet and sparseList in the sparse representation.
	 *	  sparseList is compressed using a variable length and difference encoding
	 *	  tmpSet is a list
	 *
	 *	Updates sparseList to contain all the values of sparseList and tmpSet,
	 *	 for entries that were both in sparseList and tmpSet, the
	 *	 entry with the higher nlz is taken
	 *
	 */
	void mergeSparse() {
		sparseList.insert(tmpSet.begin(),tmpSet.end());
		tmpSet = TmpSetType();

	}

	/**
	 * Reset to its initial state.
	 */
	void reset() {
		this->sparse = true;
		this->tmpSet.clear();      //  TODO: if these types are changes, initialize with appropriate size
		this->sparseList.clear();  // 
		this->M.clear();
	}

	/**
	 * Convert from sparse representation (using tmpSet and sparseList) to normal (using register)
	 */
	void toNormal() {
		if (tmpSet.size() > 0) {
			mergeSparse();
		}

		this->M = vector<uint8_t>(this->m);
		for (SparseListType::const_iterator it = sparseList.begin(); it != sparseList.end(); ++it) {
			idx_n_rank ir = decodeHash(*it);
			if (this->M[ir.idx] < ir.rank) {
				this->M[ir.idx] = ir.rank;
			}
		}

		this->sparse = false;
		this->tmpSet.clear();
		this->sparseList.clear();
	}


	/**
	 * Merge another HyperLogLogPlus into this. Converts to normal representation
	 * @param other
	 */
	void merge(HyperLogLogPlus* other) {
		if (this->p != other->p) {
			throw std::invalid_argument("precisions must be equal");
		}

		if (sparse) {
			toNormal();
		}

		if (other->sparse) {
			for (TmpSetType::const_iterator k = other->tmpSet.begin(); k != other->tmpSet.end(); ++k) {
				idx_n_rank ir = other->decodeHash(*k);
				if (M[ir.idx] < ir.rank) {
					M[ir.idx] = ir.rank;
				}
			}

			for (SparseListType::const_iterator iter = other->sparseList.begin();
					iter != other->sparseList.end(); ++iter) {
				idx_n_rank ir = other->decodeHash(*iter);
				if (M[ir.idx] < ir.rank) {
					M[ir.idx] = ir.rank;
				}
			}
		} else {
			for (size_t i = 0; i <= other->M.size(); ++i) {
				if (other->M[i] > M[i]) {
					M[i] = other->M[i];
				}
			}
		}
	}

	/**
	 *
	 * @return cardinality estimate
	 */
	uint64_t cardinality() {
		if (sparse) {
			mergeSparse();

			// if we are still 'sparse', then use linear counting, which is more
			//  accurate for low cardinalities
			return uint64_t(linearCounting(mPrime, mPrime-uint32_t(sparseList.size())));
		}

		if (rawEstimateData.empty()) {
			initRawEstimateData();
		}
		if (biasData.empty()) {
			initBiasData();
		}

		// normal
		double est = calculateEstimate(M, m);
#ifdef DEBUG1
		cerr << "est is " << est << endl;
#endif
		if (est <= double(m)*5.0) {
#ifdef DEBUG1
			cerr << "subtracting bias " << getEstimateBias(est) << endl;
#endif
			est -= getEstimateBias(est);
		}

		uint32_t v = countZeros(M);
#ifdef DEBUG1
		cerr << "m=" << m << "; v=" << v << "; M.size=" << M.size() << endl;
#endif
		if (v != 0) {
			// calculate linear counting estimate
			double lc = linearCounting(m, v);
#ifdef DEBUG1
			cerr << "lc=" << lc << endl;
#endif
			if (lc <= double(threshold[p-4])) {
				if (lc < 0) {
					cerr << "lc smaller then 0 - do something" << endl;
				}
#ifdef DEBUG1
				cerr << "returning lc " << lc << endl;
#endif
				// use it is it is smaller than the threshold
				return uint64_t(lc);
			}
		}
#ifdef DEBUG1
		cerr << "returning est " << est << endl;
#endif
		return uint64_t(est);
	}

private:

    uint8_t rank(HashSize x, uint8_t b) {
        uint8_t v = 1;
        while (v <= b && !(x & 0x80000000)) {
            v++;
            x <<= 1;
        }
        return v;
    }


	void initRawEstimateData() {
	    rawEstimateData = vector<vector<double> >();

	    rawEstimateData.push_back(vector<double>(rawEstimateData_precision4,arr_len(rawEstimateData_precision4)));
	    rawEstimateData.push_back(vector<double>(rawEstimateData_precision5,arr_len(rawEstimateData_precision5)));
	    rawEstimateData.push_back(vector<double>(rawEstimateData_precision6,arr_len(rawEstimateData_precision6)));
	    rawEstimateData.push_back(vector<double>(rawEstimateData_precision7,arr_len(rawEstimateData_precision7)));
	    rawEstimateData.push_back(vector<double>(rawEstimateData_precision8,arr_len(rawEstimateData_precision8)));
	    rawEstimateData.push_back(vector<double>(rawEstimateData_precision9,arr_len(rawEstimateData_precision9)));
	    rawEstimateData.push_back(vector<double>(rawEstimateData_precision10,arr_len(rawEstimateData_precision10)));
	    rawEstimateData.push_back(vector<double>(rawEstimateData_precision11,arr_len(rawEstimateData_precision11)));
	    rawEstimateData.push_back(vector<double>(rawEstimateData_precision12,arr_len(rawEstimateData_precision12)));
	    rawEstimateData.push_back(vector<double>(rawEstimateData_precision13,arr_len(rawEstimateData_precision13)));
	    rawEstimateData.push_back(vector<double>(rawEstimateData_precision14,arr_len(rawEstimateData_precision14)));
	    rawEstimateData.push_back(vector<double>(rawEstimateData_precision15,arr_len(rawEstimateData_precision15)));
	    rawEstimateData.push_back(vector<double>(rawEstimateData_precision16,arr_len(rawEstimateData_precision16)));
	    rawEstimateData.push_back(vector<double>(rawEstimateData_precision17,arr_len(rawEstimateData_precision17)));
	    rawEstimateData.push_back(vector<double>(rawEstimateData_precision18,arr_len(rawEstimateData_precision18)));

	}

	void initBiasData() {
		biasData = vector<vector<double> >();

		biasData.push_back(vector<double>(biasData_precision4,arr_len(biasData_precision4)));
		biasData.push_back(vector<double>(biasData_precision5,arr_len(biasData_precision5)));
		biasData.push_back(vector<double>(biasData_precision6,arr_len(biasData_precision6)));
		biasData.push_back(vector<double>(biasData_precision7,arr_len(biasData_precision7)));
		biasData.push_back(vector<double>(biasData_precision8,arr_len(biasData_precision8)));
		biasData.push_back(vector<double>(biasData_precision9,arr_len(biasData_precision9)));
		biasData.push_back(vector<double>(biasData_precision10,arr_len(biasData_precision10)));
		biasData.push_back(vector<double>(biasData_precision11,arr_len(biasData_precision11)));
		biasData.push_back(vector<double>(biasData_precision12,arr_len(biasData_precision12)));
		biasData.push_back(vector<double>(biasData_precision13,arr_len(biasData_precision13)));
		biasData.push_back(vector<double>(biasData_precision14,arr_len(biasData_precision14)));
		biasData.push_back(vector<double>(biasData_precision15,arr_len(biasData_precision15)));
		biasData.push_back(vector<double>(biasData_precision16,arr_len(biasData_precision16)));
		biasData.push_back(vector<double>(biasData_precision17,arr_len(biasData_precision17)));
		biasData.push_back(vector<double>(biasData_precision18,arr_len(biasData_precision18)));
	}

	/**
	 * Estimate the bias using empirically determined values.
	 * Uses weighted average of the two cells between which the estimate falls.
	 * TODO: Check if nearest neighbor average gives better values, as proposed in the paper
	 * @param est
	 * @return correction value for
	 */
	double getEstimateBias(double estimate) {
		vector<double> rawEstimateTable = rawEstimateData[p-4];
		vector<double> biasTable = biasData[p-4];
	
		// check if estimate is lower than first entry, or larger than last
		if (rawEstimateTable.front() >= estimate) { return rawEstimateTable.front() - biasTable.front(); }
		if (rawEstimateTable.back()  <= estimate) { return rawEstimateTable.back() - biasTable.back(); }
	
		// get iterator to first element that is not smaller than estimate
		iterator it = lower_bound(rawEstimateTable.begin(),rawEstimateTable.end(),estimate);
		size_t pos = it - rawEstimateTable.begin();

		double e1 = rawEstimateTable[pos-1];
		double e2 = rawEstimateTable[pos];
	
		double c = (estimate - e1) / (e2 - e1);

		return biasTable[pos-1]*(1-c) + biasTable[pos]*c;
	}
	

	/**
	 * Encode the hash code x as an integer, to be used in the sparse representation.
	 * see section 5.3 in Heule et al.
	 * TODO: I do not understand yet why we can use the 32bit representation here
	 * @param x the hash bits
	 * @return encoded hash value
	 */
	uint32_t encodeHash(uint64_t x) {
		// get the index (the first pPrime bits)
		uint32_t idx = (uint32_t)extractBits(x,64,64-pPrime);

		// maybe replace with: extractBits(idx,this->p-this->pPrime);
		// are the bits {63-p, ..., 63-p'} all 0?
		if (extractBits(x, 64-this->p, 64-pPrime) == 0) {
			// save the rank (and set the last bit to 1)
			uint8_t rank = clz(extractBits(x, 64-pPrime));
			return idx<<7 | uint32_t(rank<<1) | 1;
		} else {
			// else, return the idx, only (left-shifted, last bit = 0)
			return idx << 1;
		}
	}


	/**
	 * struct holding the index and rank/rho of an entry
	 */
	struct idx_n_rank {
		uint32_t idx;
		uint8_t rank;
		idx_n_rank(uint32_t index, uint8_t rank) : idx(idx), rank(rank) {}
	};

	//
	//
	/**
	 * Decode a hash from the sparse representation.
	 * Returns the index and number of leading zeros (nlz) with precision p stored in k
	 * @param k the hash bits
	 * @return index and rank in non-sparse format
	 */
	idx_n_rank decodeHash(uint32_t k)  {

		// check if the last bit is 1
		if (k&1 == 1) {
			// if yes: the hash was stored with higher precision, bits p to pPrime were 0
			uint8_t pp = p + pPrime;
			return(idx_n_rank((uint32_t)extractBits(k, 32, 32 - this->p), // first p bits are the index
							  (uint8_t) extractBits(k, 7, 1) + pp         // bits 7:1 are the rank (minus pp)
			));
		} else {
			// if no: just the pPrime long index was stored, from which both the
			//   p bits long index, as well as the rank can be restored
			return(idx_n_rank((uint32_t)extractBits(k, this->p + 1, 1), //
							  (uint8_t) clz(extractBits(k,this->pPrime - this->p, 1))
			));
		}
	}



};




#endif /* HYPERLOGLOGPLUS_H_ */
