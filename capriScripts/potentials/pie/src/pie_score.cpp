/* PIEDOCK                                                            *
 * Authors: D. V. S. Ravikant                                         *
 * (C) 2011 D. V. S. Ravikant and Ron Elber
 * Modified by Shruthi Viswanath, 2012  for performance 			  */

#include "object.hpp"
#include <sys/stat.h>

ostream *out;
char scratch_dir[512], command[8192], buf[8192];
extern char *tmp_dir;

/*
 * Arguments: model_pdb_file receptor-chains ligand-chains param-file
 */
int main(int argc, char *argv[]) {
	out = &cout;

	read_molecule_config();

	string pdbfile = string(argv[1]);
	string rchains = string(argv[2]); // enter multiple chains as one single string. e.g. AB
	string lchains = string(argv[3]); // enter multiple chains as one single string. e.g. AB
	string paramfile = string(argv[4]); // location of the PIE parameters file

	int nchains = (rchains + lchains).length();

	Complex *p[2], *c;
	Transformation *tr;

	struct stat filestat;
	string ext = ".pdb";
	string pdbBaseName = pdbfile.substr(0, pdbfile.size() - ext.size());

	if (stat((pdbfile).c_str(), &filestat) == 0) {
		string allchains = rchains + lchains;
		p[0] = new Complex(pdbBaseName, rchains, PDB);
		p[1] = new Complex(pdbBaseName, lchains, PDB);

		c = new Complex(pdbBaseName, allchains, PDB);

		tr = new Transformation(new Vector(0, 0, 0), new Vector(1, 0, 0),
				new Vector(0, 1, 0), 1.0, 0, 0);
		tr->vmetrics = new VerificationMetrics();
		tr->vmetrics->rmsd = tr->vmetrics->lrmsd = tr->vmetrics->irmsd = 100;
		ProtProtDetailedScoringMetrics *details =
				new ProtProtDetailedScoringMetrics();
		tr->detailed_scores = details;

		// compute residue contacts
		for (int i = 0; i < NUM_RESIDUE_TYPES; i++)
			for (int j = 0; j < NUM_RESIDUE_TYPES; j++)
				details->residue_contacts_core[i][j] = 0;

		for (int i = 0; i < NUM_ATOM_TYPES; i++)
			for (int j = 0; j < NUM_ATOM_TYPES; j++)
				if (NUM_ATOM_DISTANCE_DIVISIONS <= 1)
					details->atom_contacts[i][j] = 0;
				else
					for (int k = 0; k < NUM_ATOM_DISTANCE_DIVISIONS; k++)
						details->atom_dcontacts[i][j][k] = 0;

		int num_clashes = 0, num_bbclashes = 0;
		tr->eVdw = 0;
		tr->eVdw_repulsion = 0;
		for (int laindex = 0; laindex < c->num_aminoacids; laindex++) {
			Aminoacid *la = c->aminoacid[laindex];
			if (la->centroid != NULL)
				for (int raindex = 0; raindex < c->num_aminoacids; raindex++) {
					Aminoacid *ra = c->aminoacid[raindex];

					if (lchains.find(la->chain) != string::npos
							&& rchains.find(ra->chain) != string::npos) {

						if (ra->type >= 0 && ra->centroid != NULL) {
							float d2 = Vector::distance_squared(
									*(la->centroid), *(ra->centroid));
							float factor = 0;
#ifdef 	STEP_POTENTIAL 			
							if(d2 < AA_CUTOFF_SQUARED)
#endif
#ifdef LINEAR_SPLINE_POTENTIAL										 			
							if(d2 < AA_SMTHP_CUTOFF_SQUARED)
#endif												
							{
#ifdef 	STEP_POTENTIAL 	
								factor = 1.0;
#endif
#ifdef LINEAR_SPLINE_POTENTIAL	 									
								if(d2 < AA_CUTOFF_SQUARED)
								factor = 1.0;
								else {float d=sqrt(d2); factor = AA_SMTHP_FACTOR(d);}
#endif						
							}
							short ratype = ra->type, latype = la->type;
							if (ratype >= 0 && latype >= 0)
								if (ratype <= latype)
									details->residue_contacts_core[ratype][latype]
											+= factor;
								else
									details->residue_contacts_core[latype][ratype]
											+= factor;
						}

						// backbone centroid and backbone backbone contacts
						{
							Vector *vaa[3];
							vaa[0] = la->centroid;
							vaa[1] = (la->amide_nitrogen == NULL) ? NULL
									: la->amide_nitrogen->position;
							vaa[2] = (la->carbonyl_oxygen == NULL) ? NULL
									: la->carbonyl_oxygen->position;
							for (int aapi = 0; aapi < 3; aapi++)
								if (vaa[aapi] != NULL) {
									Vector vl = Vector(*(vaa[aapi]));
									if (ra->type >= 0)
										for (int aapj = 0; aapj < 3; aapj++) {
											float d, d2, factor = 0;
											if ((aapi == 0 && aapj != 0
													&& la->type >= 0) || (aapi
													!= 0 && aapj == 0
													&& ra->centroid != NULL)) {
												bool computefactor = false;
												if (aapi == 0 && aapj == 1
														&& ra->amide_nitrogen
																!= NULL) {
													computefactor = true;
													d2
															= Vector::distance_squared(
																	vl,
																	*(ra->amide_nitrogen->position));
												}
												if (aapi == 0 && aapj == 2
														&& ra->carbonyl_oxygen
																!= NULL) {
													computefactor = true;
													d2
															= Vector::distance_squared(
																	vl,
																	*(ra->carbonyl_oxygen->position));
												}
												if (aapi != 0) {
													computefactor = true;
													d2
															= Vector::distance_squared(
																	vl,
																	*(ra->centroid));
												}

												if (computefactor)
#ifdef 	STEP_POTENTIAL 			
												if(d2 < BS_CUTOFF_SQUARED)
#endif
#ifdef LINEAR_SPLINE_POTENTIAL										 			
												if(d2 < BS_SMTHP_CUTOFF_SQUARED)
#endif												
												{
#ifdef 	STEP_POTENTIAL 	
													factor = 1.0;
#endif
#ifdef LINEAR_SPLINE_POTENTIAL	 									
													if(d2 < BS_CUTOFF_SQUARED) factor = 1.0;
													else {float d=sqrt(d2); factor = BS_SMTHP_FACTOR(d);}
#endif						
												}
												if (aapi == 0 && aapj != 0
														&& la->type >= 0) {
													details->residue_contacts_core[la->type][19
															+ aapj] += factor;
													details->residue_contacts_core[19
															+ aapj][la->type]
															+= factor;
												} else {
													details->residue_contacts_core[ra->type][19
															+ aapi] += factor;
													details->residue_contacts_core[19
															+ aapi][ra->type]
															+= factor;
												}
											}

											if (aapi > 0 && aapj > 0) {
												bool computefactor = false;
												if (aapj == 1
														&& ra->amide_nitrogen
																!= NULL) {
													computefactor = true;
													d2
															= Vector::distance_squared(
																	vl,
																	*(ra->amide_nitrogen->position));
												}
												if (aapj == 2
														&& ra->carbonyl_oxygen
																!= NULL) {
													computefactor = true;
													d2
															= Vector::distance_squared(
																	vl,
																	*(ra->carbonyl_oxygen->position));
												}

												if (computefactor)
#ifdef 	STEP_POTENTIAL 			
												if(d2 < BB_CUTOFF_SQUARED)
#endif
#ifdef LINEAR_SPLINE_POTENTIAL										 			
												if(d2 < BB_SMTHP_CUTOFF_SQUARED)
#endif												
												{
#ifdef 	STEP_POTENTIAL 	
													factor = 1.0;
#endif
#ifdef LINEAR_SPLINE_POTENTIAL	 									
													if(d2 < BB_CUTOFF_SQUARED) factor = 1.0;
													else {float d=sqrt(d2); factor = BB_SMTHP_FACTOR(d);}
#endif		
												}
												details->residue_contacts_core[19
														+ aapi][19 + aapj]
														+= factor;
												if (aapi != aapj)
													details->residue_contacts_core[19
															+ aapj][19 + aapi]
															+= factor;
											}
										} // aapj
								}
						}

						for (hash_map<const char*, Atom*, hash<const char*> ,
								eqstr>::iterator laitr = la->atom.begin(); laitr
								!= la->atom.end(); laitr++) {
							Atom *al = (Atom *) laitr->second;
							for (hash_map<const char*, Atom*,
									hash<const char*> , eqstr>::iterator raitr =
									ra->atom.begin(); raitr != ra->atom.end(); raitr++) {
								Atom *ar = (Atom *) raitr->second;
								float d2 = Vector::distance_squared(
										al->position, ar->position);

								if (al->atom_type >= 0 && ar->atom_type >= 0)
									if (NUM_ATOM_DISTANCE_DIVISIONS <= 1) {
										float factor = 0;

#ifdef 	STEP_POTENTIAL 			
										if(d2 < ATOM_STEPP_DCUTOFF_SQUARED)// && (d2 >= ATOM_STEPP_DMIN_SQUARED))
#endif
#ifdef LINEAR_SPLINE_POTENTIAL										 			
										if(d2 < ATOM_SMTHP_DCUTOFF_SQUARED)
#endif												
										{
#ifdef 	STEP_POTENTIAL 	
											factor = 1.0;
#endif
#ifdef LINEAR_SPLINE_POTENTIAL				 									
											if(d2 < ATOM_STEPP_DCUTOFF_SQUARED) factor = 1.0;
											else {float d=sqrt(d2); factor = ATOM_SMTHP_FACTOR(d);}
#endif			

											if (al->atom_type < ar->atom_type)
												details->atom_contacts[al->atom_type][ar->atom_type]
														+= factor;
											else
												details->atom_contacts[ar->atom_type][al->atom_type]
														+= factor;
										}
									} else {
										short div = ATOM_DISTANCE_DIVISION(d2);
										if (div >= 0)
											if (al->atom_type < ar->atom_type)
												details->atom_dcontacts[al->atom_type][ar->atom_type][div]
														+= 1.0;
											else
												details->atom_dcontacts[ar->atom_type][al->atom_type][div]
														+= 1.0;
									}

								if (d2 < ATOM_SMTHP_DCUTOFF_SQUARED) {
									float factor = 0;

									float d_scaled = sqrt(d2
											/ (3.2 * ar->sigma));
									factor = VDW_ATTR(d_scaled);
									if (factor > 0) {
										// 6-12 factor = 4 * a->sqrt_eps * ar->sigma_cubed * factor;
										factor = 4 * ar->sqrt_eps * factor;
										tr->eVdw += factor * al->sqrt_eps;
									} else {
										factor = VDW_REPUL(d_scaled);
										factor = 4 * ar->sqrt_eps * factor;
										tr->eVdw_repulsion += factor
												* al->sqrt_eps;
									}
								}
							}
						}
					}
				}
		}

		tr->num_clashes = num_clashes;
		tr->num_bbclashes = num_bbclashes;
		tr->sEvolution_interface = tr->eVdw_repulsion;
		//tr->print_details(&fdetails,TN_PROTPROT_DETAILED_SCORE_VERIFY);
		tr->print_pie_score(TN_PROTPROT_DETAILED_SCORE_VERIFY,paramfile);

		// clean up
		delete c;
		for (int mi = 0; mi < 2; mi++)
			delete p[mi];
		delete tr;
	}

	else 
		 cout << "PDB file " << pdbfile.c_str() <<  " failed to open."<< endl ;

}
