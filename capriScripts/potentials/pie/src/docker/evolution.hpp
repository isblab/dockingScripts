/* PIEDOCK                                                            *
 * Authors: D. V. S. Ravikant                                         *
 * (C) 2011 D. V. S. Ravikant and Ron Elber                           */

class Alignment{
	public:
		//vector<string> seqids;
		string aseq1, aseq2;
		float size, score, eval, identities, positives;
		// identities and positives are percentages
		int qstart, qend, sstart, send;
		long index;
		
		float positives_interface, identities_interface, gaps_interface;

		Alignment(string,string,long);
		void compute_identities();
		void print_details(ostream *);
};

class Sequence{
	public:
		string id, id_lowercase; 
		string tag, database, organism;
		char chain;
		string *aasequence;
		int malignment_index;
		
		Alignment *seq_alignment, *struct_alignment, *loopp_alignment;
		
	Sequence(string);
};

class Species {
	public:
		short *num_sequences;
		Sequence ***sequences;
		bool *non_redundant_pairs;
		short ***homologues_aatypes;
};
