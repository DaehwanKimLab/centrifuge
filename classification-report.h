/*
 * classification-report.h
 *
 *  Created on: Jan 24, 2017
 *      Author: fbreitwieser
 */

#ifndef CLASSIFICATION_REPORT_H_
#define CLASSIFICATION_REPORT_H_

class ClassificationReport {
private:

	enum class FORMAT : uint8_t {
		KRAKEN_FORMAT,
		EXT_KRAKEN_FORMAT,
		FLAT
	};

	/*
	 * TaxCounts give the number of reads per clade, number of reads that stay and abundance for each taxon
	 */
	struct TaxCounts : ReadCounts {
		uint64_t n_reads_clade = 0;
		array<double,2> abundance = {{0.0, 0.0}};

		TaxCounts& operator+=(const TaxCounts& b) {
			n_reads += b.n_reads;
			n_reads_clade += b.n_reads_clade;
			total_hit_len += b.total_hit_len;
			total_score += b.total_score;
			abundance += b.abundance;
			return *this;
		}

		TaxCounts& operator+=(const ReadCounts& b) {
			n_reads += b.n_reads;
			n_reads_clade += b.n_reads;
			total_hit_len += b.total_hit_len;
			total_score += b.total_score;
			return *this;
		}


		TaxCounts& addChild(const TaxCounts& child) {
			total_score += child.total_score;
			total_hit_len += child.total_hit_len;
			// add reads of child that come from itself to clade reads
			n_reads_clade += child.n_reads_clade;
			abundance += child.abundance;
			return *this;
		}

		TaxCounts& addChild(const ReadCounts& child) {
			total_score += child.total_score;
			total_hit_len += child.total_hit_len;
			// add reads of child that come from itself to clade reads
			n_reads_clade += child.n_reads;
			return *this;
		}



		TaxCounts() {}

		TaxCounts(const ReadCounts& b, const array<double,2>& abundance_) : ReadCounts(b), n_reads_clade(n_reads), abundance(abundance_) { }
		/*
	TaxCounts& operator+=(const ReadCounts& b) {
		n_unique_reads += b.n_unique_reads;
		n_reads += b.n_reads;
		total_score += b.total_score;
		total_hit_len += b.total_hit_len;
		weighted_reads += b.weighted_reads;
		return *this;
	}
		 */
	};

	struct TaxInfo {
		char rank;
		string name;
		TaxId parent;
		vector<TaxId> children;
		vector<pair<UId, ReadCounts> > uids;
		TaxCounts counts;

		TaxInfo() {}
		TaxInfo(TaxId parent_, char rank_, string name_) : rank(rank_), name(name_), parent(parent_) {}
	};
    const map<UId, std::array<double,2> > _uid_abundance;

	//=== these ELists are accessed by UIds
	const EList<pair<string, uint64_t> >& _uid_to_tid;
	const EList<string >& _uid_refnames;

	const map<TaxId, TaxInfo> _taxinfo;
	//map <TaxId, vector<TaxId> > parents;

	const EList<REPORTCOLS> _report_cols;
	const bool _show_zeros;
	ofstream _reportOfb;
	FORMAT _format;
	double _total_n_reads = 0;

	inline bool set_to_parent(const map<TaxId, TaxonomyNode>& tree, TaxId& taxid) {
		if (taxid == 0 || taxid == TaxId(-1))
			return false;

		auto it = tree.find(taxid);
		if (it == tree.end()) {
			DEBUG_MSG("Taxonomy ID " << taxid << " has no taxonomy mapping in database" << endl);
			taxid = -1;
			return false;
		}
		if (it->second.parent_tid != taxid) {
			taxid = it->second.parent_tid;
			return true;
		}
		return false;
	}

	// return lowest common ancestor of a set of taxid
	TaxId uid_lca(const map<TaxId, TaxonomyNode>& tree, const UId_set& uids) {
		if (uids.size() == 0)
			return 0; // unidentified read

		auto it = uids.begin();
		TaxId taxidA = _uid_to_tid[*it].second;
		if (uids.size() == 1)
			return taxidA;

		//TaxId current_lca = taxidA;
		vector<TaxId> parentsA = { taxidA };
		while (set_to_parent(tree, taxidA)) {
			parentsA.push_back(taxidA);
		}

		size_t lca_pos = 0;
		for (advance(it,1); it != uids.end(); advance(it,1)) {
			TaxId current_taxid = _uid_to_tid[*it].second;
			bool found_it = false;
			do {
				if (current_taxid == parentsA[lca_pos])
					break;

				for (size_t pos = lca_pos; pos < parentsA.size(); ++pos) {
					if (parentsA[pos] == current_taxid) {
						lca_pos = pos;
						found_it = true;
						break;
					}
				}

				if (found_it)
					break;
			} while (set_to_parent(tree, current_taxid));
		}
		return parentsA[lca_pos];
	}



	// return lowest common ancestor of a set of taxid
	TaxId lca(const map<TaxId, TaxonomyNode>& tree, const set<TaxId>& taxids) {
		if (taxids.size() == 0)
			throw invalid_argument("no taxIDs given to lca!");
		if (taxids.size() == 1)
			return *taxids.begin();

		TaxId current_lca;
		TaxId taxidA = *taxids.begin();
		set<TaxId> parentsA = { taxidA };
		while (set_to_parent(tree, taxidA)) {
			parentsA.insert(taxidA);
		}

		int lca_pos = -1;
		for (TaxId current_taxid : taxids) {
			do {
				auto it = find(parentsA.begin(), parentsA.end(), current_taxid);
				if (it != parentsA.end()) {
					//ptrdiff_t pos = it - parentsA.begin(); // doesn't work, no operator- ?
					int pos = std::distance(parentsA.begin(), it);
					if (pos > lca_pos) {
						current_lca = *it;
						lca_pos = pos;
					}
					break;
				}
			} while (set_to_parent(tree, current_taxid));
		}
		return current_lca;
	}

	map<TaxId, TaxInfo> get_taxinfo(const Ebwt<index_t>& ebwt, const ClassificationMetrics& spm) {

		const map<TaxId, TaxonomyNode>& tree = ebwt.tree();
		const map<TaxId, string>& name_map = ebwt.name();

		cerr << "uid_quant1" << endl;
		map<TaxId, TaxInfo > taxinfo = {
			{ 0, { 0, 'U', "unclassified"}},
			{TaxId(-1), {TaxId(-1), 'N', "uncategorized"}}
		};

		cerr << "taxinfo1" << endl;
		// fill rank, name, and children
		for (auto& tree_node : tree) {
			TaxId tax_id = tree_node.first;
			TaxId parent_tax_id = tree_node.second.parent_tid;
			taxinfo[parent_tax_id].children.push_back(tax_id);
			taxinfo[tax_id].parent = parent_tax_id;
			taxinfo[tax_id].rank = get_tax_rank_char(tree_node.second.rank);
			taxinfo[tax_id].name = name_map.at(tax_id);
		}

		cerr << "taxinfo2" << endl;
		// fill uid map and initialize counts
		for (const pair<UId_set, ReadCounts> & uid_counts : spm.counts) {
			UId_set uid_set = uid_counts.first;
			const ReadCounts& rc = uid_counts.second;
			uint64_t lca_taxid = uid_lca(tree, uid_set);

			taxinfo[lca_taxid].counts += rc;

			// If there is only one UId in this set, add the counts to it
			if (uid_set.size() == 1) 
				taxinfo[lca_taxid].uids.push_back({ *uid_set.begin(), rc });

			while (set_to_parent(tree,lca_taxid)) {
				taxinfo[lca_taxid].counts.addChild(rc);
			}
		}

		// Now add the abundances
		for (const pair<UId, array<double,2> >& ua: _uid_abundance) {
			TaxId taxid = _uid_to_tid[ua.first].second;

			taxinfo[taxid].counts.abundance += ua.second;
			while (set_to_parent(tree,taxid)) {
				taxinfo[taxid].counts.abundance += ua.second;
			}		
		}
		cerr << "taxinfoend" << endl;

		return taxinfo;
	}

	void print_line(bool is_tax_id, TaxId tax_id,
			string name, uint64_t depth, TaxCounts counts ) {
		for (auto& col : _report_cols) {
			switch (col) {
			case REPORTCOLS::NAME:        _reportOfb << name ; break;
			case REPORTCOLS::SPACED_NAME:       _reportOfb << string(2*depth, ' ') + name; break;
			case REPORTCOLS::TAX_ID:     _reportOfb << (tax_id == (uint64_t)-1? -1 : (int64_t)tax_id); break;
			case REPORTCOLS::DEPTH:     _reportOfb << depth; break;
			case REPORTCOLS::PERCENTAGE:  _reportOfb << 100*counts.n_reads_clade/_total_n_reads; break;
			case REPORTCOLS::ABUNDANCE:  _reportOfb << 100*counts.abundance[0]; break;
			case REPORTCOLS::ABUNDANCE_LEN:  _reportOfb << 100*counts.abundance[1]; break;
			case REPORTCOLS::NUM_READS_CLADE:  _reportOfb << counts.n_reads_clade; break;
			case REPORTCOLS::NUM_READS:  _reportOfb << counts.n_reads; break;
			//case REPORTCOLS::NUM_UNIQUE_KMERS: _reportOfb << counts.n_kmers; break;
			//case REPORTCOLS::GENOME_SIZE: ; break;
			//case REPORTCOLS::NUM_WEIGHTED_READS: ; break;
			//case REPORTCOLS::SUM_SCORE: ; break;
			case REPORTCOLS::TAX_RANK: _reportOfb << (is_tax_id? _taxinfo.at(tax_id).rank : 'X') ; break;
			default: _reportOfb << "NA";
			}
			if (&col == &_report_cols.back()) {
				_reportOfb << '\n';
			} else {
				_reportOfb << '\t';
			}

		}
	}

public:


	void print_report(string format, string rank) {
		_total_n_reads = 
			_taxinfo.at(0).counts.n_reads_clade + 
			_taxinfo.at(1).counts.n_reads_clade + 
			_taxinfo.at(-1).counts.n_reads_clade;
		if (format == "kraken") {
			// A: print number of unidentified reads
			print_report(0,0);
			// B: print normal results
			print_report(1,0);
			// C: Print Unclassified stuff
			print_report(-1,0);
		} else {
			// print stuff at a certain level ..
			//_uid_abundance;
			//_taxinfo
				
		}
	}

	void print_report(TaxId tax_id, uint64_t depth) {

		if (_show_zeros || _taxinfo.at(tax_id).counts.n_reads_clade > 0) {
			print_line(true, tax_id, _taxinfo.at(tax_id).name, depth, _taxinfo.at(tax_id).counts);

			// If there are UIDs with quant data for this TaxIDs, print them out - but order them first
			// TODO: order by abundance

			// The taxinfo uids must always have counts, as it is initialized from there
			vector<pair<UId, ReadCounts> > uids = _taxinfo.at(tax_id).uids;

			//sort(uids.begin(), uids.end(), [&](auto a, auto b) -> bool { 
			//		return( a.second > b.second );
			//		});
			for (auto& ur : uids) {
				string uid_name = " >" + _uid_refnames[ur.first];
				//string uid_name = "->" + _uid_to_tid[uid].first;
				if (_show_zeros || ur.second.n_reads > 0) 
					print_line(false, _uid_to_tid[ur.first].second, uid_name, depth, 
							TaxCounts(ur.second, _uid_abundance.at(ur.first)));
			}

			// Go to the children of that TaxID
			// TODO: Order by abundance
			vector<UId> children = _taxinfo.at(tax_id).children;

			//sort(children.begin(), children.end(), [&](TaxId a, TaxId b) -> bool { 
			//		return( _taxinfo.at(a).counts.abundance >
			//					_taxinfo.at(b).counts.abundance); });

			for (auto child : children) {
				if (child == tax_id) continue;  // Should only happen for root / taxId 1
				print_report(child, depth+1);
			}
		}

	}

	ClassificationReport(
			string file_name,
			const Ebwt<index_t>& ebwt,
			const ClassificationMetrics& spm,
			EList<REPORTCOLS>& report_cols,
			bool show_zeros) :
				_uid_abundance { spm.calculateAbundance2(ebwt) },
				_uid_to_tid { ebwt.uid_to_tid() },
				_uid_refnames { ebwt.refnames() },
				_taxinfo { get_taxinfo(ebwt, spm) },
				_report_cols { report_cols },
				_show_zeros { show_zeros },
				_reportOfb(file_name.c_str())
				{
					//EList<string> p_refnames;
					//readEbwtRefnames<index_t>(fname, p_refnames);

				}
	
	~ClassificationReport() {
		_reportOfb.close();
	}

				//================
};





#endif /* CLASSIFICATION_REPORT_H_ */
