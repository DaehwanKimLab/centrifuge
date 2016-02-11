/*
 * taxonomy.h
 *
 *  Created on: Feb 10, 2016
 *      Author: fbreitwieser
 */

#ifndef TAXONOMY_H_
#define TAXONOMY_H_

#include<map>
#include<utility>
#include<string>

enum {
    RANK_UNKNOWN = 0,
    RANK_STRAIN,
    RANK_SPECIES,
    RANK_GENUS,
    RANK_FAMILY,
    RANK_ORDER,
    RANK_CLASS,
    RANK_PHYLUM,
    RANK_KINGDOM,
    RANK_DOMAIN,
    RANK_FORMA,
    RANK_INFRA_CLASS,
    RANK_INFRA_ORDER,
    RANK_PARV_ORDER,
    RANK_SUB_CLASS,
    RANK_SUB_FAMILY,
    RANK_SUB_GENUS,
    RANK_SUB_KINGDOM,
    RANK_SUB_ORDER,
    RANK_SUB_PHYLUM,
    RANK_SUB_SPECIES,
    RANK_SUB_TRIBE,
    RANK_SUPER_CLASS,
    RANK_SUPER_FAMILY,
    RANK_SUPER_KINGDOM,
    RANK_SUPER_ORDER,
    RANK_SUPER_PHYLUM,
    RANK_TRIBE,
    RANK_VARIETAS,
};

struct TaxonomyNode {
    uint64_t parent_tid;
    uint8_t  rank;
    uint8_t  leaf;

    TaxonomyNode(uint64_t _parent_tid, uint8_t  _rank, uint8_t _leaf):
    	parent_tid(_parent_tid), rank(_rank), leaf(_leaf) {};

    TaxonomyNode(): parent_tid(0), rank(RANK_UNKNOWN), leaf(false) {};
};

typedef std::map<uint64_t, TaxonomyNode> TaxonomyTree;

struct TaxonomyPathTable {
    static const size_t nranks = 7;

    map<uint64_t, uint32_t> tid_to_pid;  // from taxonomic ID to path ID
    ELList<uint64_t> paths;

    void buildPaths(const EList<pair<string, uint64_t> >& uid_to_tid,
                    const TaxonomyTree& tree)
    {
        map<uint32_t, uint32_t> rank_map;
        rank_map[RANK_STRAIN]      = 0;
        rank_map[RANK_SUB_SPECIES] = 0;
        rank_map[RANK_SPECIES]     = 1;
        rank_map[RANK_GENUS]       = 2;
        rank_map[RANK_FAMILY]      = 3;
        rank_map[RANK_ORDER]       = 4;
        rank_map[RANK_CLASS]       = 5;
        rank_map[RANK_PHYLUM]      = 6;

        tid_to_pid.clear();
        paths.clear();
        for(size_t i = 0; i < uid_to_tid.size(); i++) {
            uint64_t tid = uid_to_tid[i].second;
            if(tid_to_pid.find(tid) != tid_to_pid.end())
                continue;
            if(tree.find(tid) == tree.end())
                continue;

            tid_to_pid[tid] = (uint32_t)paths.size();
            paths.expand();
            EList<uint64_t>& path = paths.back();
            path.resizeExact(nranks);
            path.fillZero();
            bool first = true;
            while(true) {
                TaxonomyTree::const_iterator itr = tree.find(tid);
                if(itr == tree.end()) {
                    break;
                }
                const TaxonomyNode& node = itr->second;
                uint32_t rank = std::numeric_limits<uint32_t>::max();
                if(first && node.rank == RANK_UNKNOWN) {
                    rank = rank_map[RANK_STRAIN];
                } else if(rank_map.find(node.rank) != rank_map.end()) {
                    rank = rank_map[node.rank];
                }
                if(rank < path.size()) {
                    path[rank] = tid;
                }

                first = false;
                if(node.parent_tid == tid) {
                    break;
                }
                tid = node.parent_tid;
            }
        }
    }

    void getPath(uint64_t tid, EList<uint64_t>& path) const {
        map<uint64_t, uint32_t>::const_iterator itr = tid_to_pid.find(tid);
        if(itr != tid_to_pid.end()) {
            uint32_t pid = itr->second;
            assert_lt(pid, paths.size());
            path = paths[pid];
        } else {
            path.clear();
        }
    }
};


inline static const char* getRankName(int rank) {
    switch(rank) {
        case RANK_STRAIN:        return "strain";
        case RANK_SPECIES:       return "species";
        case RANK_GENUS:         return "genus";
        case RANK_FAMILY:        return "family";
        case RANK_ORDER:         return "order";
        case RANK_CLASS:         return "class";
        case RANK_PHYLUM:        return "phylum";
        case RANK_KINGDOM:       return "kingdom";
        case RANK_FORMA:         return "forma";
        case RANK_INFRA_CLASS:   return "infraclass";
        case RANK_INFRA_ORDER:   return "infraorder";
        case RANK_PARV_ORDER:    return "parvorder";
        case RANK_SUB_CLASS:     return "subclass";
        case RANK_SUB_FAMILY:    return "subfamily";
        case RANK_SUB_GENUS:     return "subgenus";
        case RANK_SUB_KINGDOM:   return "subkingdom";
        case RANK_SUB_ORDER:     return "suborder";
        case RANK_SUB_PHYLUM:    return "subphylum";
        case RANK_SUB_SPECIES:   return "subspecies";
        case RANK_SUB_TRIBE:     return "subtribe";
        case RANK_SUPER_CLASS:   return "superclass";
        case RANK_SUPER_FAMILY:  return "superfamily";
        case RANK_SUPER_KINGDOM: return "superkingdom";
        case RANK_SUPER_ORDER:   return "superorder";
        case RANK_SUPER_PHYLUM:  return "superphylum";
        case RANK_TRIBE:         return "tribe";
        case RANK_VARIETAS:      return "varietas";
        default:                 return "no rank";
    };
    return "no rank";
}

inline static uint8_t getNumericRank(string rank_name) {
	uint8_t numeric_rank = RANK_UNKNOWN;

    if(rank_name == "strain") {
        numeric_rank = RANK_STRAIN;
    } else if(rank_name == "species") {
        numeric_rank = RANK_SPECIES;
    } else if(rank_name == "genus") {
        numeric_rank = RANK_GENUS;
    } else if(rank_name == "family") {
        numeric_rank = RANK_FAMILY;
    } else if(rank_name == "order") {
        numeric_rank = RANK_ORDER;
    } else if(rank_name == "class") {
        numeric_rank = RANK_CLASS;
    } else if(rank_name == "phylum") {
        numeric_rank = RANK_PHYLUM;
    } else if(rank_name == "kingdom") {
        numeric_rank = RANK_KINGDOM;
    } else if(rank_name == "forma") {
        numeric_rank = RANK_FORMA;
    } else if(rank_name == "infraclass") {
        numeric_rank = RANK_INFRA_CLASS;
    } else if(rank_name == "infraorder") {
        numeric_rank = RANK_INFRA_ORDER;
    } else if(rank_name == "parvorder") {
        numeric_rank = RANK_PARV_ORDER;
    } else if(rank_name == "subclass") {
        numeric_rank = RANK_SUB_CLASS;
    } else if(rank_name == "subfamily") {
        numeric_rank = RANK_SUB_FAMILY;
    } else if(rank_name == "subgenus") {
        numeric_rank = RANK_SUB_GENUS;
    } else if(rank_name == "subkingdom") {
        numeric_rank = RANK_SUB_KINGDOM;
    } else if(rank_name == "suborder") {
        numeric_rank = RANK_SUB_ORDER;
    } else if(rank_name == "subphylum") {
        numeric_rank = RANK_SUB_PHYLUM;
    } else if(rank_name == "subspecies") {
        numeric_rank = RANK_SUB_SPECIES;
    } else if(rank_name == "subtribe") {
        numeric_rank = RANK_SUB_TRIBE;
    } else if(rank_name == "superclass") {
        numeric_rank = RANK_SUPER_CLASS;
    } else if(rank_name == "superfamily") {
        numeric_rank = RANK_SUPER_FAMILY;
    } else if(rank_name == "superkingdom") {
        numeric_rank = RANK_SUPER_KINGDOM;
    } else if(rank_name == "superorder") {
        numeric_rank = RANK_SUPER_ORDER;
    } else if(rank_name == "superphylum") {
        numeric_rank = RANK_SUPER_PHYLUM;
    } else if(rank_name == "tribe") {
        numeric_rank = RANK_TRIBE;
    } else if(rank_name == "varietas") {
        numeric_rank = RANK_VARIETAS;
    }

    return numeric_rank;
}

inline static uint64_t getParentTaxIDAtRank(const TaxonomyTree& tree, uint64_t taxid, uint8_t at_rank) {
	while (true) {
		TaxonomyTree::const_iterator itr = tree.find(taxid);
		if(itr == tree.end()) {
			break;
		}
		const TaxonomyNode& node = itr->second;

		if (node.rank == at_rank) {
			return taxid;
		} else if (node.rank > at_rank || node.parent_tid == taxid) {
			return 0;
		}

		taxid = node.parent_tid;
	}
	return 0;
}

inline static TaxonomyTree readTaxonomyTree(string taxonomy_fname) {
	TaxonomyTree tree;
	ifstream taxonomy_file(taxonomy_fname.c_str(), ios::in);
	if(taxonomy_file.is_open()) {
		char line[1024];
		while(!taxonomy_file.eof()) {
			line[0] = 0;
			taxonomy_file.getline(line, sizeof(line));
			if(line[0] == 0 || line[0] == '#') continue;
			istringstream cline(line);
			uint64_t tid, parent_tid;
			char dummy; string rank;
			cline >> tid >> dummy >> parent_tid >> dummy >> rank;
			if(tree.find(tid) != tree.end()) {
				cerr << "Warning: " << tid << " already has a parent!" << endl;
				continue;
			}

			tree[tid] = TaxonomyNode(parent_tid, getNumericRank(rank), false);
		}
		taxonomy_file.close();
	} else {
		cerr << "Error: " << taxonomy_fname << " doesn't exist!" << endl;
		throw 1;
	}
	return tree;
}


#endif /* TAXONOMY_H_ */
