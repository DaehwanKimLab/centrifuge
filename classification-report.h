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

	/*
	 * TaxCounts give the number of reads per clade, number of reads that stay and abundance for each taxon
	 */
	struct TaxCounts : ReadCounts {
		uint32_t n_kmers = 0;
		uint32_t reads_stay = 0;
		uint32_t reads_clade = 0;
		double abundance = 0.0;
		double abundance_len = 0.0;

		TaxCounts& operator+=(const TaxCounts& b) {
			n_unique_reads += b.n_unique_reads;
			n_reads += b.n_reads;
			sum_score += b.sum_score;
			n_kmers += b.n_kmers;
			reads_stay += b.reads_stay;
			reads_clade += b.reads_clade;
			summed_hit_len += b.summed_hit_len;
			weighted_reads += b.weighted_reads;
			abundance += b.abundance;
			abundance_len += b.abundance_len;
			return *this;
		}

		TaxCounts& addChild(const TaxCounts& child) {
			n_unique_reads += child.n_unique_reads;
			n_reads += child.n_reads;
			sum_score += child.sum_score;
			n_kmers += child.n_kmers;

			// reads_stay is not influenced by child
			reads_clade += child.reads_clade;

			summed_hit_len += child.summed_hit_len;
			weighted_reads += child.weighted_reads;
			abundance += child.abundance;
			abundance_len += child.abundance_len;

			return *this;
		}


		TaxCounts() : ReadCounts() { }

		TaxCounts(const ReadCounts& b, const double abundance_, const double abundance_len_) {
			n_unique_reads = b.n_unique_reads;
			n_reads = b.n_reads;
			sum_score = b.sum_score;
			summed_hit_len = b.summed_hit_len;
			weighted_reads = b.weighted_reads;
			abundance = abundance_;
			abundance_len = abundance_len_;
			reads_clade = n_unique_reads;
			reads_stay = n_unique_reads;

		}
		/*
	TaxCounts& operator+=(const ReadCounts& b) {
		n_unique_reads += b.n_unique_reads;
		n_reads += b.n_reads;
		sum_score += b.sum_score;
		summed_hit_len += b.summed_hit_len;
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
		vector<UId> uids;
		TaxCounts quant_data;

		//TaxIdInfo(char _rank, string _name) : rank{_rank}, name{_name} {}
	};

	//=== for quantitation
	const EList<pair<string, uint64_t> >& _uid_to_tid;
	const EList<string >& _uid_refnames;
	const map<UId, TaxCounts> _uid_quant_data;

	const map<TaxId, TaxInfo> _taxinfo;
	//map <TaxId, vector<TaxId> > parents;

	const EList<REPORTCOLS> _report_cols;
	const bool _show_zeros;
	ofstream& _reportOfb;

	inline TaxId set_to_parent(const map<TaxId, TaxonomyNode>& tree, TaxId& taxid) {
		auto it = tree.find(taxid);
		if (it == tree.end()) {
			cerr << "Ignoring taxonomy ID " << taxid << " which has no taxonomy mapping in database" << endl;
		}
		if (it->second.parent_tid != taxid) {
			taxid = it->second.parent_tid;
			return true;
		}
		return false;
	}

	// return lowest common ancestor of a set of taxid
	TaxId uid_lca(const map<TaxId, TaxonomyNode>& tree, const UId_set& uids) {
		cerr << "uid_lca" << endl;
		if (uids.size() == 0)
			return 0; // unidentified read
		if (uids.size() == 1)
			return _uid_to_tid[*uids.begin()].second;

		TaxId taxidA = _uid_to_tid[*uids.begin()].second;
		//TaxId current_lca = taxidA;
		vector<TaxId> parentsA = { taxidA };
		while (set_to_parent(tree, taxidA)) {
			parentsA.push_back(taxidA);
		}
		size_t lca_pos = 0;

		for (UId current_uid : uids) {
			TaxId current_taxid = _uid_to_tid[current_uid].second;
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
		// requires _uid_quant_data to be set!

		const map<TaxId, TaxonomyNode>& tree = ebwt.tree();
		const map<TaxId, string>& name_map = ebwt.name();
		const EList<pair<string, TaxId>>& uid_to_tid = ebwt.uid_to_tid();
		map<TaxId, TaxInfo > taxinfo;

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
		// fill uid map and initialize quant_data
		for (const pair<UId, TaxCounts> & uid_counts : _uid_quant_data) {
			UId uid = uid_counts.first;
			const TaxCounts& tc = uid_counts.second;

			uint64_t taxid = uid_to_tid[uid].second;
			taxinfo[taxid].uids.push_back(uid);

			DEBUG_MSG("Adding to uid " << uid << " taxid " << taxid << " " << tc.n_reads << endl);
			taxinfo[taxid].quant_data.addChild(tc);
			do {
				taxid = taxinfo[taxid].parent;
				taxinfo[taxid].quant_data.addChild(tc);
			} while (taxid != taxinfo[taxid].parent);
		}

		cerr << "taxinfo3" << endl;
		// Add reads that were assigned to multiple UIds to the LCA
		for (const pair<UId_set, uint64_t> & uid_counts : spm.observed_uid_sets) {
			TaxId taxid = uid_lca(tree, uid_counts.first);
			if (uid_counts.first.size() > 1) {
				taxinfo[taxid].quant_data.reads_stay += uid_counts.second;
				do {
					taxid = taxinfo[taxid].parent;
					taxinfo[taxid].quant_data.reads_clade += uid_counts.second;
				} while (taxid != taxinfo[taxid].parent);
			}
		}

		cerr << "taxinfoend" << endl;

		return taxinfo;
	}

	void print_line(bool is_tax_id, TaxId tax_id,
			string name, uint64_t depth, TaxCounts counts ) {
		DEBUG_MSG("print_line");
		for (auto& col : _report_cols) {
			switch (col) {
			case REPORTCOLS::NAME:        _reportOfb << name ; break;
			case REPORTCOLS::SPACED_NAME:       _reportOfb << string(2*depth, ' ') + name; break;
			case REPORTCOLS::TAX_ID:     _reportOfb << tax_id; break;
			case REPORTCOLS::DEPTH:     _reportOfb << depth; break;
			case REPORTCOLS::ABUNDANCE:  _reportOfb << 100*counts.abundance; break;
			case REPORTCOLS::ABUNDANCE_LEN:  _reportOfb << 100*counts.abundance_len; break;
			case REPORTCOLS::NUM_READS:  _reportOfb << counts.reads_clade; break;
			case REPORTCOLS::NUM_READS_STAY:  _reportOfb << counts.reads_stay; break;
			case REPORTCOLS::NUM_WEIGHTED_READS:  _reportOfb << counts.weighted_reads; break;
			case REPORTCOLS::NUM_UNIQUE_READS: _reportOfb << counts.n_unique_reads; break;
			case REPORTCOLS::NUM_UNIQUE_KMERS: _reportOfb << counts.n_kmers; break;
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

	map<UId, TaxCounts> get_uid_quant_data(const Ebwt<index_t>& ebwt, const ClassificationMetrics& spm) {

		cerr << "uid_quant1" << endl;
		const index_t* uid_seq_lengths = ebwt.plen();  // sequence lengths in the UID space
		const map<UId, double> abundance = spm.calculateAbundance2(ebwt);
		double sum = 0.0;
		for (const auto & it : abundance) {
			sum += it.second*uid_seq_lengths[it.first];
		}

		cerr << "uid_quant2" << endl;
		map<UId, TaxCounts> res;
		for (const pair<UId,ReadCounts>& ur : spm.species_counts) {
			res[ur.first] = TaxCounts(ur.second, abundance.at(ur.first), abundance.at(ur.first)*uid_seq_lengths[ur.first]/sum);
		}
		cerr << "uid_quantend" << endl;
		return res;
	}

public:


	void print_report(TaxId tax_id=1, uint64_t depth=0) {

		if (_show_zeros || _taxinfo.at(tax_id).quant_data.reads_clade > 0) {
			print_line(true, tax_id, _taxinfo.at(tax_id).name, depth, _taxinfo.at(tax_id).quant_data);

			// If there are UIDs with quant data for this TaxIDs, print them out
			for (const auto & uid : _taxinfo.at(tax_id).uids) {
				string uid_name = "->" + _uid_refnames[uid];
				//string uid_name = "->" + _uid_to_tid[uid].first;
				if (_uid_quant_data.find(uid) != _uid_quant_data.end() &&
						(_show_zeros || _uid_quant_data.at(uid).reads_clade > 0))
					print_line(false, tax_id, uid_name, depth, _uid_quant_data.at(uid));
			}

			// Go to the children of that TaxID
			for (auto& child : _taxinfo.at(tax_id).children) {
				if (child == tax_id) continue;  // Should only happen for root / taxId 1
				print_report(child, depth+1);
			}
		}
	}

	ClassificationReport(
			ofstream & reportOfb,
			const Ebwt<index_t>& ebwt,
			const ClassificationMetrics& spm,
			EList<REPORTCOLS>& report_cols,
			bool show_zeros) :
				_uid_to_tid { ebwt.uid_to_tid() },
				_uid_refnames { ebwt.refnames() },
				_uid_quant_data { get_uid_quant_data(ebwt, spm) },
				_taxinfo { get_taxinfo(ebwt, spm) },
				_report_cols { report_cols },
				_show_zeros { show_zeros },
				_reportOfb { reportOfb }
				{
					//EList<string> p_refnames;
					//readEbwtRefnames<index_t>(fname, p_refnames);

				}

				//================
};





#endif /* CLASSIFICATION_REPORT_H_ */
