/* PIEDOCK                                                            *
 * Authors: D. V. S. Ravikant                                         *
 * (C) 2011 D. V. S. Ravikant and Ron Elber                           */

#include "object.hpp"

int sorted_cell[2*UNIT_GRID_PHI_DIVISIONS*UNIT_GRID_PHI_DIVISIONS][2*UNIT_GRID_PHI_DIVISIONS*UNIT_GRID_PHI_DIVISIONS];
float distance_bound[2*UNIT_GRID_PHI_DIVISIONS*UNIT_GRID_PHI_DIVISIONS][2*UNIT_GRID_PHI_DIVISIONS*UNIT_GRID_PHI_DIVISIONS];

float avgexposedarea[] = {42.965149484236 ,106.05113602302 ,75.274478813878 ,76.339517320637 ,33.9697265625 ,88.035641630727 ,96.197290821968 ,37.783932731095 ,74.183006535948 ,49.232292967124,
		51.897620553744, 106.50223725286 ,62.944319344933 ,58.839937920331 ,63.79711087035 ,53.602121549552 ,57.558885822476 ,68.416231018173 ,70.406265344201 ,45.507441482885,1,1 };

bool **aacontact_core, **aacontact_rim, *aarcontact, *aalcontact;
extern char scratch_dir[];
float **sas_b, **sas_bprime, **sas_fraclost;

Transformation **localneighbors;

#define MAX_NUM_CONSERVED_RESIDUES 81

extern float ***atom32_dpotential;
extern float vdw_weight;
extern string piedock_home;

Object::Object(Complex *c, Complex *cH){
	this->c = c;
	this->cH = cH;

	c->compute_motions();
	// c->compute_volume_sasa(true,true);
	if(c->filetype != PROCESSED)
		for(int i = 0; i < c->num_atoms; i++)
			c->atom[i]->is_buried = (c->atom[i]->sasa <= 0.1*PI*c->atom[i]->radius*c->atom[i]->radius) || c->atom[i]->is_buried;

	c->compute_aacontacts();
	// c->compute_electrostatics();

	if(c != cH){
		cH->compute_motions();
		// cH->compute_electrostatics();
	}

	// compute the percentile for conservation of residues on surface
	{
		float sum=0, sum_squares=0;
		int count=0;
		for(int i =0; i < c->num_aminoacids; i++){
			Aminoacid* aa = c->aminoacid[i];
			float entropy = aa->entropy_sequence;
			//*out << aa->sa << " "; aa->print_details(out);
			if(entropy != 0 && aa->rsa > 0.1){
				sum += entropy;
				sum_squares += entropy*entropy;
				count++;
			}
		}
		float sigma,avg = sum/count;
		sigma = sqrt(sum_squares/count);
		*out << "num surface residues " << count << endl;

		for(int i =0; i < c->num_aminoacids; i++){
			Aminoacid* aa = c->aminoacid[i];
			if(aa->entropy_sequence != 0 && aa->rsa > 0.1){
				int numlower=0;
				for(int j =0; j < c->num_aminoacids;j++){
					Aminoacid* aa2 = c->aminoacid[j];
					if(aa2->entropy_sequence != 0 && aa2->rsa > 0.1 && aa2->entropy_sequence <= aa->entropy_sequence)
						numlower++;
				}
				// for large proteins avoid labeling one third of surface as highly conserved

				if(count - numlower < MAX_NUM_CONSERVED_RESIDUES){
					if(count <= MAX_NUM_CONSERVED_RESIDUES)
						aa->entropy_sequence_percentile = (numlower*3.0)/count;
					else
						aa->entropy_sequence_percentile = ((numlower - count + MAX_NUM_CONSERVED_RESIDUES)*3.0)/MAX_NUM_CONSERVED_RESIDUES;
					*out << aa->index << " " << aa->rsa << " " << count-numlower << " " << aa->entropy_sequence_percentile << endl;
				}
			}
		}
	}

	vector<Cluster*> faces = c->clusters;
	num_faces = faces.size();
	//*out << "#faces " << num_faces << endl;

	face = (Cluster**) (malloc((num_faces+1)*sizeof(Cluster*)));
	crstype = new int[num_faces+1];
	curvature = new float[num_faces+1];

	int index = 0;
	for(vector<Cluster*>::iterator itr = faces.begin(); itr != faces.end(); itr++){
		Cluster *f = (Cluster*) *itr;
		//f->print_details();
		Vector *v = f->point;
		if(v != NULL){
			face[index] = f;

			/*for(vector<int>::iterator ritr = (f->face_residue_indices).begin(); ritr != (f->face_residue_indices).end(); ritr++){
				residue_indices[index].push_back(*ritr);
			}*/

			index++;
		}
	}
	num_faces = index;
	*out << "#faces " << num_faces << endl;

	transformation_index = 0;
	match_size=NULL;
	match_domain=NULL;
	match_range=NULL;

	transformed_atom_position = NULL;
	partner_cell_of_atom = NULL;
	tr_details = new ProtProtDetailedScoringMetrics();

	num_atoms_in_contact=0;
	atoms_in_contact = (unsigned int *) malloc(sizeof(unsigned int)*c->num_atoms);
}

int Object::find_close_points(Vector *point, int *close_point){
	Vector v = *point - *grid_origin;
	int vx = (int) (v.x/GRID_SPACING);
	int vy = (int) (v.y/GRID_SPACING);
	int vz = (int) (v.z/GRID_SPACING);

	if(vx >= 0 && vx < grid_num_xdivisions && vy >= 0 && vy < grid_num_ydivisions && vz >= 0 && vz < grid_num_zdivisions){
		unsigned int index = (vx*grid_num_ydivisions + vy)*grid_num_zdivisions + vz;

		GridCell *cell;
		// find the close points if there is a contact
		vector<int> cpoints;
		if((cell = grid[index]) != NULL && ((cell->type == INTERIOR) || (cell->type == INTERMEDIATE) || (cell->type == SURFACE) || (cell->curvature != 0))){
			int spread = (int) (RES_ATOM_CUTOFF/GRID_SPACING + 1);
			int minx,maxx,miny,maxy,minz,maxz;
			minx = vx - spread;
			miny = vy - spread;
			minz = vz - spread;

			maxx = vx + spread;
			maxy = vy + spread;
			maxz = vz + spread;

			if(minx <= 0) minx = 0;
			if(miny <= 0) miny = 0;
			if(minz <= 0) minz = 0;

			if(maxx >= grid_num_xdivisions - 1) maxx = grid_num_xdivisions - 1;
			if(maxy >= grid_num_ydivisions - 1) maxy = grid_num_ydivisions - 1;
			if(maxz >= grid_num_zdivisions - 1) maxz = grid_num_zdivisions - 1;

			for(int x = minx; x <= maxx; x++){
				for(int y = miny; y <= maxy; y++){
					for(int z = minz; z <= maxz; z++){
						index = (x*grid_num_ydivisions + y)*grid_num_zdivisions + z;
						if(grid[index] != NULL && grid[index]->points != NULL){
							int num_points = grid[index]->num_points;
							unsigned int* points = grid[index]->points;
							for(int pi = 0; pi < num_points; pi++) {
								int point_index = points[pi];
								if(point_index > num_faces){
									int acindex = point_index - num_faces;
									Vector v = *(c->atom[acindex]->position);	
									float distance = Vector::distance(v,point);
									if(distance <= RES_ATOM_CUTOFF)
										cpoints.push_back(point_index);
								} else {
									Vector v = *(this->face[point_index]->point);
									float distance = Vector::distance(v,point);
									if(distance <= RES_ATOM_CUTOFF){
										cpoints.push_back(point_index);
									}
								}
							}
						}
					}
				}
			}

			if(cpoints.size() > 0){
				int i = 0;
				for(vector<int>::iterator itr = cpoints.begin(); itr != cpoints.end(); itr++){
					close_point[i++] = (int) *itr;
				}
				return cpoints.size();
			}
		}
	}

	return -1;
}

Reference_Frame::Reference_Frame(){
};

Reference_Frame::Reference_Frame(Vector *t, Vector *ex, Vector *ey, float scale,long frame_number){
	translation = t;
	this->ex = ex;
	this->ey = ey;
	ez = new Vector(ex->cross(ey));
	this->scale = scale;
	this->frame_number = frame_number;
}

Reference_Frame::Reference_Frame(Reference_Frame* rf){
	translation = rf->translation;
	ex = rf->ex;
	ey = rf->ey;
	ez = new Vector(ex->cross(ey));
	scale = rf->scale;
}

Vector Reference_Frame::rotate(Vector p){
	Vector v = p * (1.0/scale);
	float vx = v.dot(ex);
	float vy = v.dot(ey);
	float vz = v.dot(ez);
	return Vector(vx,vy,vz);
}

Vector Reference_Frame::inverse_rotate(Vector p){
	float x = p.dot(Vector(ex->x,ey->x,ez->x));
	float y = p.dot(Vector(ex->y,ey->y,ez->y));
	float z = p.dot(Vector(ex->z,ey->z,ez->z));

	Vector v = Vector(x,y,z);
	v = v*scale;
	return v;
}

Vector Reference_Frame::transform(Vector p){
	//*out << "(" << p.x << "," << p.y << "," << p.z << ")" << endl;
	Vector v = p - *translation;
	v = v * (1.0/scale);
	//*out << "(" << v.x << "," << v.y << "," << v.z << ") " << scale << endl;
	float vx = v.dot(ex);
	float vy = v.dot(ey);
	float vz = v.dot(ez);
	return Vector(vx,vy,vz);
}

Vector Reference_Frame::inverse_transform(Vector p){
	float x = p.dot(Vector(ex->x,ey->x,ez->x));
	float y = p.dot(Vector(ex->y,ey->y,ez->y));
	float z = p.dot(Vector(ex->z,ey->z,ez->z));

	Vector v = Vector(x,y,z);
	v = v*scale;
	v = v + *translation;

	return v;
}

// compute t1 o t2
Reference_Frame* Reference_Frame::compose(Reference_Frame *rf1,Reference_Frame *rf2){
	Reference_Frame *res = new Reference_Frame();
	Vector vx = Vector(rf2->ex->x,rf2->ey->x,rf2->ez->x);
	Vector vy = Vector(rf2->ex->y,rf2->ey->y,rf2->ez->y);
	Vector vz = Vector(rf2->ex->z,rf2->ey->z,rf2->ez->z);

	res->ex = new Vector(rf1->ex->dot(vx),rf1->ex->dot(vy),rf1->ex->dot(vz));
	res->ey = new Vector(rf1->ey->dot(vx),rf1->ey->dot(vy),rf1->ey->dot(vz));
	res->ez = new Vector(rf1->ez->dot(vx),rf1->ez->dot(vy),rf1->ez->dot(vz));

	Vector p = Vector(rf1->translation);
	float x = p.dot(Vector(rf2->ex->x,rf2->ey->x,rf2->ez->x));
	float y = p.dot(Vector(rf2->ex->y,rf2->ey->y,rf2->ez->y));
	float z = p.dot(Vector(rf2->ex->z,rf2->ey->z,rf2->ez->z));
	res->translation = new Vector(*(rf2->translation) + Vector(x,y,z));

	res->scale = rf1->scale * rf2->scale;
	return res;
}

Reference_Frame* Reference_Frame::invert(Reference_Frame *rf){
	Reference_Frame *res = new Reference_Frame();

	//*out  << rf->ex << " " << rf->ey << " " << rf->ez << " " << rf->translation << endl; out->flush();
	res->ex = new Vector(rf->ex->x,rf->ey->x,rf->ez->x);
	//*out  << res->ex->x << " " << res->ex->y << " " << res->ex->z << endl; out->flush();
	res->ey = new Vector(rf->ex->y,rf->ey->y,rf->ez->y);
	//*out  << res->ey->x << " " << res->ey->y << " " << res->ey->z << endl; out->flush();
	res->ez = new Vector(rf->ex->z,rf->ey->z,rf->ez->z);
	//*out  << res->ez->x << " " << res->ez->y << " " << res->ez->z << endl; out->flush();

	Vector v = rf->transform(Vector(0,0,0));
	res->translation = new Vector(v);
	res->scale = 1.0/rf->scale;
	return res;
}

Transformation::Transformation(Vector *t, Vector *ex, Vector *ey, float scale, float votes, long transformation_index)
: Reference_Frame(t, ex, ey, scale, transformation_index) {
	this->votes = votes;
	num_neighbors = 0;
	cmr = NULL;
	vmetrics = NULL;
	details = NULL;
	detailed_scores = NULL;
}

/*
 * Convert from euler angles such that the application of rotation corresponds to tr->transform()
 */
Transformation::Transformation(float alpha,float beta,float gamma,long tid){
	scale = 1.0;
	frame_number = tid;
	float sin_alpha = sin(alpha);
	float cos_alpha = cos(alpha);
	float sin_beta = sin(beta);
	float cos_beta = cos(beta);
	float sin_gamma = sin(gamma);
	float cos_gamma = cos(gamma);
	ex = new Vector( cos_alpha*cos_gamma - sin_alpha*cos_beta*sin_gamma , -cos_alpha*sin_gamma - sin_alpha*cos_beta*cos_gamma , sin_beta*sin_alpha );
	ey = new Vector( sin_alpha*cos_gamma + cos_alpha*cos_beta*sin_gamma , -sin_alpha*sin_gamma + cos_alpha*cos_beta*cos_gamma , -sin_beta*cos_alpha );
	ez = new Vector( sin_beta*sin_gamma , sin_beta*cos_gamma , cos_beta );
	translation = new Vector(0,0,0);
	//*out << "tr from euler check " << ex->norm_squared() << " " << ey->norm_squared() << " " << ez->norm_squared() 
	//	<< ex->dot(ey) << " " << ey->dot(ez) << " " << ez->dot(ex) << " " << (ex->cross(ey) - *ez).norm_squared() << endl;
	cmr = NULL;
	vmetrics = NULL;
	details = NULL;
	detailed_scores = NULL;
}

void Reference_Frame::eulerangles(float *alpha, float *beta, float *gamma){
	*beta = acos(ez->z);
	float sinbeta=1.0;
	if( ez->y != 0 || ez->x != 0)
		*gamma = atan(ez->x/ez->y);
	else 
		sinbeta=0;
	if( ey->z != 0 || ex->z != 0)
		*alpha = atan(-ex->z/ey->z);
	else 
		sinbeta=0;

	if(sinbeta==0){
		*gamma = 0;
		*alpha = acos(ex->x);
		if(ey->x < 0)
			*alpha -= PI;
	} else {
		//sinbeta = sin(*beta);
		if(ey->z > 0)	*alpha += PI;
		if(*alpha > PI) *alpha -= 2*PI;
		if(ez->y < 0)	*gamma += PI;
		if(*gamma > PI)	*gamma -= 2*PI;
	}
}

Reference_Frame::~Reference_Frame() {
	//*out << translation << " " << ex << " " << ey << " " << ez << endl; out->flush();
	if(translation != NULL)	delete translation;
	if(ex != NULL)	delete ex;
	if(ey != NULL)	delete ey;
	if(ez != NULL)	delete ez;
}

Transformation::~Transformation() {
	//*out << translation << " " << ex << " " << ey << " " << ez << " " << cmr << " " << vmetrics << " " 
	//	<< details << " " << detailed_scores << endl; out->flush();
	/*if(translation != NULL)	delete translation;
	if(ex != NULL)	delete ex;
	if(ey != NULL)	delete ey;
	if(ez != NULL)	delete ez;*/
	if(cmr != NULL)	delete cmr;
	if(vmetrics != NULL)	delete vmetrics;
	if(details != NULL)	delete details;
	if(detailed_scores != NULL) delete detailed_scores;
}

Transformation::Transformation(char* buf, unsigned short type){
	char *current = buf;
	memcpy(&frame_number,current,sizeof(long));
	current += sizeof(long);

	//*out << "frame_number " << frame_number << endl;
	float x,y,z;
	memcpy(&x,current,sizeof(float));
	current += sizeof(float);
	memcpy(&y,current,sizeof(float));
	current += sizeof(float);
	memcpy(&z,current,sizeof(float));
	current += sizeof(float);
	translation = new Vector(x,y,z);

	memcpy(&x,current,sizeof(float));
	current += sizeof(float);
	memcpy(&y,current,sizeof(float));
	current += sizeof(float);
	memcpy(&z,current,sizeof(float));
	current += sizeof(float);
	ex = new Vector(x,y,z);

	memcpy(&x,current,sizeof(float));
	current += sizeof(float);
	memcpy(&y,current,sizeof(float));
	current += sizeof(float);
	memcpy(&z,current,sizeof(float));
	current += sizeof(float);
	ey = new Vector(x,y,z);

	ez = new Vector(ex->cross(ey));

	scale = 1.0;

	memcpy(&votes,current,sizeof(float));
	current += sizeof(float);
	memcpy(&curvature_score,current,sizeof(float));
	current += sizeof(float);

	memcpy(&eVdw,current,sizeof(float));
	current += sizeof(float);

	memcpy(&num_contacts,current,sizeof(unsigned short));
	current += sizeof(unsigned short);
	memcpy(&eResiduepair,current,sizeof(float));
	current += sizeof(float);
	memcpy(&eElectrostatic,current,sizeof(float));
	current += sizeof(float);
	memcpy(&sEvolution_interface,current,sizeof(float));
	current += sizeof(float);
	memcpy(&eSolvation,current,sizeof(float));
	current += sizeof(float);
	memcpy(&num_neighbors,current,sizeof(unsigned short));
	current += sizeof(unsigned short);

	//*out << "check " << (current - buf) << " " << basic_byte_size << endl; out->flush();

	if(type == TN_VERIFY || type == TN_PROTPROT_DETAILED_SCORE_VERIFY || type == TN_PROTRNA_DETAILED_SCORE_VERIFY){
		vmetrics = new VerificationMetrics();
		memcpy(&(vmetrics->rmsd),current,sizeof(float));
		current += sizeof(float);
		memcpy(&(vmetrics->lrmsd),current,sizeof(float));
		current += sizeof(float);
		memcpy(&(vmetrics->irmsd),current,sizeof(float));
		current += sizeof(float);
		memcpy(&(vmetrics->frac_native_contacts_predicted),current,sizeof(float));
		current += sizeof(float);
		memcpy(&(vmetrics->frac_nonnative_contacts_predicted),current,sizeof(float));
		current += sizeof(float);
		memcpy(&(vmetrics->delta_r),current,sizeof(float));
		current += sizeof(float);
		memcpy(&(vmetrics->delta_U),current,sizeof(float));
		current += sizeof(float);
	} else
		vmetrics = NULL;

	if(type == TN_PROTPROT_DETAILED_SCORE_VERIFY){
		memcpy(&num_clashes,current,sizeof(unsigned short));
		current += sizeof(unsigned short);
		memcpy(&num_bbclashes,current,sizeof(unsigned short));
		current += sizeof(unsigned short);

		memcpy(&delta_sasa,current,sizeof(float));
		current += sizeof(float);

		memcpy(&freeE_rigidmotion,current,sizeof(float));
		current += sizeof(float);

		ProtProtDetailedScoringMetrics *detailed_scores = new ProtProtDetailedScoringMetrics();
		this->detailed_scores = detailed_scores;
		memcpy(&(detailed_scores->points_vs_distance[0]),current,sizeof(unsigned short)*NUM_DTRANSFORM_DIVISIONS);
		current += sizeof(unsigned short)*NUM_DTRANSFORM_DIVISIONS;
		memcpy(&(detailed_scores->atoms_vs_distance[0]),current,sizeof(unsigned short)*NUM_DTRANSFORM_DIVISIONS);
		current += sizeof(unsigned short)*NUM_DTRANSFORM_DIVISIONS;
		for(int i = 0 ; i < NUM_DTRANSFORM_DIVISIONS ; i++){
			memcpy(&(detailed_scores->grid_contacts[i][0]),current,sizeof(unsigned short)*NUM_DTRANSFORM_DIVISIONS);
			current += sizeof(unsigned short)*NUM_DTRANSFORM_DIVISIONS;
		}

		memcpy(&(detailed_scores->delta_vdwsa[0]),current,sizeof(float)*NUM_RESIDUE_TYPES);
		current += sizeof(float)*NUM_RESIDUE_TYPES;
		memcpy(&(detailed_scores->delta_sasa[0]),current,sizeof(float)*NUM_RESIDUE_TYPES);
		current += sizeof(float)*NUM_RESIDUE_TYPES;
		/*memcpy(&(detailed_scores->delta_vol[0]),current,sizeof(float)*NUM_RESIDUE_TYPES);
		current += sizeof(float)*NUM_RESIDUE_TYPES;
		memcpy(&(detailed_scores->delta_sevol[0]),current,sizeof(float)*NUM_RESIDUE_TYPES);
		current += sizeof(float)*NUM_RESIDUE_TYPES;*/

		int num_residue_types = NUM_RESIDUE_TYPES;
		/*{
			for(int i = 0 ; i < num_residue_types ; i++){
				memcpy(&(detailed_scores->residue_contacts_core[i][0]),current,sizeof(unsigned short)*num_residue_types);
				current += sizeof(unsigned short)*num_residue_types;
			}
			for(int i = 0 ; i < num_residue_types ; i++){
				memcpy(&(detailed_scores->residue_contacts_rim[i][0]),current,sizeof(unsigned short)*num_residue_types);
				current += sizeof(unsigned short)*num_residue_types;
			}
		}*/
		{
			for(int i = 0 ; i < num_residue_types ; i++){
				memcpy(&(detailed_scores->residue_contacts_core[i][0]),current,sizeof(float)*num_residue_types);
				current += sizeof(float)*num_residue_types;
			}
			for(int i = 0 ; i < num_residue_types ; i++){
				memcpy(&(detailed_scores->residue_contacts_rim[i][0]),current,sizeof(float)*num_residue_types);
				current += sizeof(float)*num_residue_types;
			}
		}

		num_residue_types = NUM_ENTROPY_DIVISIONS;
		for(int i = 0 ; i < num_residue_types ; i++){
			memcpy(&(detailed_scores->conserv_contacts[i][0]),current,sizeof(unsigned short)*num_residue_types);
			current += sizeof(unsigned short)*num_residue_types;
		}

		memcpy(&(detailed_scores->elec[0]),current,sizeof(float)*NUM_ELEC_FEATURES);
		current += sizeof(float)*NUM_ELEC_FEATURES;

		num_residue_types = NUM_COARSE_RTYPES;
		for(int i = 0; i < num_residue_types; i++)
			for(int j = 0; j < num_residue_types; j++){
				memcpy(&(detailed_scores->threebody_contacts[i][j][0]),current,sizeof(unsigned short)*num_residue_types);
				current += sizeof(unsigned short)*num_residue_types;
			}

		memcpy(&(detailed_scores->delta_vdwsa_atom[0]),current,sizeof(float)*NUM_ATOM_TYPES);
		current += sizeof(float)*NUM_ATOM_TYPES;
		memcpy(&(detailed_scores->delta_sasa_atom[0]),current,sizeof(float)*NUM_ATOM_TYPES);
		current += sizeof(float)*NUM_ATOM_TYPES;
		/*memcpy(&(detailed_scores->delta_vol_atom[0]),current,sizeof(float)*NUM_ATOM_TYPES);
		current += sizeof(float)*NUM_ATOM_TYPES;
		memcpy(&(detailed_scores->delta_sevol_atom[0]),current,sizeof(float)*NUM_ATOM_TYPES);
		current += sizeof(float)*NUM_ATOM_TYPES;*/

		for(int i = 0; i < NUM_ATOM_TYPES; i++){
			if(NUM_ATOM_DISTANCE_DIVISIONS <= 1){
				memcpy(&(detailed_scores->atom_contacts[i][0]),current,sizeof(float)*NUM_ATOM_TYPES);
				current += sizeof(float)*NUM_ATOM_TYPES;
			} else {
				for(int j = 0; j < NUM_ATOM_TYPES; j++){
					memcpy(&(detailed_scores->atom_dcontacts[i][j][0]),current,sizeof(float)*NUM_ATOM_DISTANCE_DIVISIONS);
					current += sizeof(float)*NUM_ATOM_DISTANCE_DIVISIONS;
				}
			}
		}

		memcpy(&(detailed_scores->atom_vdw_repulsion[0]),current,sizeof(unsigned short)*NUM_ATOM_VDWREPULSION_DIVISIONS);
		current += sizeof(unsigned short)*NUM_ATOM_VDWREPULSION_DIVISIONS;
	} else if(type == TN_PROTRNA_DETAILED_SCORE_VERIFY){
		memcpy(&num_clashes,current,sizeof(unsigned short));
		current += sizeof(unsigned short);
		memcpy(&num_bbclashes,current,sizeof(unsigned short));
		current += sizeof(unsigned short);

		memcpy(&delta_sasa,current,sizeof(float));
		current += sizeof(float);

		memcpy(&freeE_rigidmotion,current,sizeof(float));
		current += sizeof(float);

		ProtRnaDetailedScoringMetrics *detailed_scores = new ProtRnaDetailedScoringMetrics();
		this->detailed_scores = detailed_scores;
		for(int i = 0; i < NUM_RESIDUE_TYPES; i++){
			memcpy(&(detailed_scores->coarse_contacts[i][0]),current,sizeof(float)*NUM_RNA_PARTICLE_TYPES);
			current += sizeof(float)*NUM_RNA_PARTICLE_TYPES;
		}

		/*for(int i = 0; i < NUM_RESIDUE_TYPES; i++)
			for(int j = 0; j < NUM_RNA_PARTICLE_TYPES; j++)
				cout << frame_number << " init " << i << " " << j << " " << detailed_scores->coarse_contacts[i][j] << endl;*/
	} else
		detailed_scores = NULL;

	cmr = NULL;
	details = NULL;
}

void Reference_Frame::distance(float *d, float *Ud, Reference_Frame *tr, Vector center){
	*d = Vector::distance(*translation+inverse_rotate(center),*(tr->translation)+tr->inverse_rotate(center));
	distance(Ud,tr);
}

void Reference_Frame::distance(float *d, float *Ud, Reference_Frame *tr){
	*d = Vector::distance(translation,tr->translation);
	distance(Ud,tr);
}

void Reference_Frame::distance(float *Ud, Reference_Frame *tr){
	float ud = 0;
	ud += (ex->dot(tr->ex) -1)*(ex->dot(tr->ex) -1);
	ud += (ex->dot(tr->ey))*(ex->dot(tr->ey));
	ud += (ex->dot(tr->ez))*(ex->dot(tr->ez));

	ud += (ey->dot(tr->ex))*(ey->dot(tr->ex));
	ud += (ey->dot(tr->ey) -1)*(ey->dot(tr->ey) -1);
	ud += (ey->dot(tr->ez))*(ey->dot(tr->ez));

	ud += (ez->dot(tr->ex))*(ez->dot(tr->ex));
	ud += (ez->dot(tr->ey))*(ez->dot(tr->ey));
	ud += (ez->dot(tr->ez) -1)*(ez->dot(tr->ez) -1);

	*Ud = sqrt(ud);
}

MultiTransformation::MultiTransformation(string id, Transformation *tr){
	this->id = id;
	transformations.clear();
	transformations.push_back(tr);
	contacts = NULL;
	delta_r = NULL;
	zscore = 0;

	rmsd = 0;
	lrmsd = 0;
	irmsd = 0;
	frac_native_contacts_predicted =0;
}

MultiTransformation::MultiTransformation(string id, MultiTransformation *cn, Transformation *tr){
	this->id = id;
	vector<Transformation*>::iterator titr = (cn->transformations).begin(), tend = (cn->transformations).end();
	while(titr != tend){
		transformations.push_back(*titr);
		titr++;
	}
	transformations.push_back(tr);
	contacts = NULL;
	delta_r = NULL;
	zscore = 0;

	rmsd = 0;
	lrmsd = 0;
	irmsd = 0;
	frac_native_contacts_predicted =0;
}

void MultiTransformation::compute_details(hash_map<long,float,hash<long>,eqlong> *residue_energy,
		hash_map<int, short, hash<int>,eqint> *lresidue_vdw_energy,
		hash_map<int, short, hash<int>,eqint> *rresidue_vdw_energy){
	for(vector<Transformation*>::iterator titr = transformations.begin(); titr != transformations.end(); titr++){
		Transformation* tr = *titr;

		for(hash_map<long,float,hash<long>,eqlong>::iterator rditr = tr->details->contact_energy_by_residue.begin();
				rditr != tr->details->contact_energy_by_residue.end(); rditr++){
			long index = rditr->first;
			if(residue_energy->count(index) == 0)
				(*residue_energy)[index] = rditr->second;
			else
				(*residue_energy)[index] = minimum((*residue_energy)[index],rditr->second);
			// higher score is better in current metric
		}

		for(hash_map<int, float, hash<int>,eqint>::iterator ritr = tr->details->lresidue_vdw_energy.begin();
				ritr != tr->details->lresidue_vdw_energy.end(); ritr++){
			int residue = ritr->first;
			if(lresidue_vdw_energy->count(residue) == 0)
				(*lresidue_vdw_energy)[residue] = ritr->second;
			else
				(*lresidue_vdw_energy)[residue] = maximum(ritr->second,(*lresidue_vdw_energy)[residue]);
		}

		for(hash_map<int, float, hash<int>,eqint>::iterator ritr = tr->details->rresidue_vdw_energy.begin();
				ritr != tr->details->rresidue_vdw_energy.end(); ritr++){
			int residue = ritr->first;
			if(rresidue_vdw_energy->count(residue) == 0)
				(*rresidue_vdw_energy)[residue] = ritr->second;
			else
				(*rresidue_vdw_energy)[residue] = maximum(ritr->second,(*rresidue_vdw_energy)[residue]);
		}
	}
}

void MultiTransformation::score(){
	hash_map<long,float,hash<long>,eqlong> residue_energy1;
	hash_map<int, short, hash<int>,eqint> lresidue_vdw_energy1;
	hash_map<int, short, hash<int>,eqint> rresidue_vdw_energy1;
	compute_details(&residue_energy1, &lresidue_vdw_energy1 , &rresidue_vdw_energy1);

	num_contacts = residue_energy1.size();
	eResiduepair = 0;
	for(hash_map<long,float,hash<long>,eqlong>::iterator rditr = residue_energy1.begin(); rditr != residue_energy1.end(); rditr++){
		eResiduepair += rditr->second;
	}

	vdw_energy = 0.0;
	residues_in_surface = 0;
	residues_in_intermediate = 0;
	residues_in_interior = 0;
	for(hash_map<int, short, hash<int>,eqint>::iterator ritr = lresidue_vdw_energy1.begin();
			ritr != lresidue_vdw_energy1.end(); ritr++){
		short e = lresidue_vdw_energy1[ritr->first];;
		vdw_energy += e;
		if(e <= SURFACE)
			residues_in_surface++;
		else if (e < INTERIOR)
			residues_in_intermediate++;
		else if (e >= INTERIOR)
			residues_in_interior++;
	}

	for(hash_map<int, short, hash<int>,eqint>::iterator ritr = rresidue_vdw_energy1.begin();
			ritr != rresidue_vdw_energy1.end(); ritr++){
		short e = rresidue_vdw_energy1[ritr->first];
		vdw_energy += e;
		if(e <= SURFACE)
			residues_in_surface++;
		else if (e < INTERIOR)
			residues_in_intermediate++;
		else if (e >= INTERIOR)
			residues_in_interior++;
	}

	hash_set<long, hash<long>, eqlong> neighbors;
	eVdw = 0;
	curvature_score = 0;
	votes = 0;
	for(vector<Transformation*>::iterator titr = transformations.begin(); titr != transformations.end(); titr++){
		Transformation *tr = *titr;
		for(vector<Transformation*>::iterator nnitr = (tr->details->nearest_neighbors).begin(); nnitr != (tr->details->nearest_neighbors).end(); nnitr++){
			Transformation *nntr = *nnitr;
			long nntindex = nntr->frame_number;
			neighbors.insert(nntindex);
		}
		eVdw += tr->eVdw;
		curvature_score += tr->curvature_score;
		votes += tr->votes;
	}
	num_neighbors = neighbors.size();
	votes /= transformations.size();
	eVdw /= transformations.size();
	curvature_score /= transformations.size();
}

/*
 * Compute the minimum distance between transformations in the combinations
 */
void MultiTransformation::distance(float *d, float *ud, MultiTransformation *mtr){
	int count = 0;
	for(vector<Transformation*>::iterator titr = transformations.begin(); titr != transformations.end(); titr++){
		Transformation* tr = *titr;
		for(vector<Transformation*>::iterator mtr_titr = (mtr->transformations).begin(); mtr_titr != (mtr->transformations).end(); mtr_titr++){
			Transformation* mtr_tr = *mtr_titr;
			if(mtr_tr->frame_number != tr->frame_number){
				float pd, pud;
				tr->distance(&pd,&pud,mtr_tr);
				if(count++ == 0){
					*d = pd;
					*ud = pud;
				} else {
					(*d) = minimum((*d),pd);
					(*ud) = minimum((*ud),pud);
				}
			}
		}
	}
}

/*
 * Assume that the details of this transformation have already been computed
 * compute the hamming distance between the contact maps
 */
void MultiTransformation::hamming_distance(int *hamming_distance, int *num_mtr_contacts, MultiTransformation *mtr, int num_lresidues, int num_rresidues){
	*num_mtr_contacts = 0;
	*hamming_distance = 0;

	/*int num_possible_contacts = num_lresidues * num_rresidues;
	int num_bytes = (num_possible_contacts/sizeof(unsigned char))+1;

	if(contacts == NULL){
		contacts = new unsigned char[num_bytes];
		for(int i = 0 ; i < num_bytes; i++)
			contacts[i] = 0;
		for(vector<Transformation*>::iterator titr = transformations.begin(); titr != transformations.end(); titr++){
			Transformation* tr = *titr;	
			for(int i = 0 ; i < num_bytes ; i++)
				contacts[i] |= (tr->contacts)[i];
		}
	}

	unsigned char mtr_contacts[num_bytes];
	for(int i = 0 ; i < num_bytes; i++)
		mtr_contacts[i] = 0;

	for(vector<Transformation*>::iterator titr = (mtr->transformations).begin(); titr != (mtr->transformations).end(); titr++){
		Transformation* tr = *titr;

		for(int i = 0 ; i < num_bytes ; i++)
			mtr_contacts[i] |= (tr->contacts)[i];
	}

	for(int i = 0 ; i < num_bytes ; i++){
		unsigned char symm_diff = mtr_contacts[i] ^ contacts[i];
		if(mtr_contacts[i] || symm_diff){
			for(unsigned char mask = 1; mask != 0 ; mask = (mask << 1)){
				//*out << "mask " << (int) mask << " " << (int) mtr_contacts[i] << " " << (int) contacts[i] << " " << (int) symm_diff << endl;
				if(mtr_contacts[i] & mask)
					(*num_mtr_contacts)++;
				if(symm_diff & mask)
					(*hamming_distance)++;
			}
		}
	}*/

	hash_map<long,float,hash<long>,eqlong> mtr_residue_distances;
	for(vector<Transformation*>::iterator titr = (mtr->transformations).begin(); titr != (mtr->transformations).end(); titr++){
		Transformation* tr = *titr;

		for(hash_map<long,float,hash<long>,eqlong>::iterator rditr = tr->details->contact_energy_by_residue.begin();
				rditr != tr->details->contact_energy_by_residue.end(); rditr++){
			long index = rditr->first;
			if(mtr_residue_distances.count(index) == 0)
				mtr_residue_distances[index] = rditr->second;
			else
				mtr_residue_distances[index] = minimum(mtr_residue_distances[index],rditr->second);
		}
	}
	*num_mtr_contacts = mtr_residue_distances.size();

	//*out << (residue_distances->size()) << " " << mtr_residue_distances.size() << " ";

	for(hash_map<long,float,hash<long>,eqlong>::iterator rditr = mtr_residue_distances.begin();
			rditr != mtr_residue_distances.end(); rditr++){
		if(residue_distances->count(rditr->first) == 0)
			(*hamming_distance)++;
	}

	for(hash_map<long,float,hash<long>,eqlong>::iterator rditr = residue_distances->begin();
			rditr != residue_distances->end(); rditr++){
		if(mtr_residue_distances.count(rditr->first) == 0)
			(*hamming_distance)++;
	}

	//*out << *num_mtr_contacts << " hd " << hamming_distance << endl;
}

void MultiTransformation::score_addition(float *num_new_contacts, float *num_contacts, float *vdw_energy, float *dock_score, Transformation *tr){
	hash_map<long,float,hash<long>,eqlong> residue_energy1;
	hash_map<int, short, hash<int>,eqint> lresidue_vdw_energy1;
	hash_map<int, short, hash<int>,eqint> rresidue_vdw_energy1;
	//*out << "check" << endl;
	compute_details(&residue_energy1, &lresidue_vdw_energy1 , &rresidue_vdw_energy1);
	//*out << "checking " << residue_distances1.size() << " " << lresidue_vdw_energy1.size() << " " << rresidue_vdw_energy1.size() << " " << endl;
	//out->flush();
	*dock_score = 0;
	for(hash_map<long,float,hash<long>,eqlong>::iterator rditr = residue_energy1.begin(); rditr != residue_energy1.end(); rditr++){
		if(tr->details->contact_energy_by_residue.count(rditr->first) == 0)
			(*dock_score) += rditr->second;
		else
			(*dock_score) += minimum(rditr->second,(tr->details->contact_energy_by_residue)[rditr->first]);
	}

	*num_contacts = residue_energy1.size();
	*num_new_contacts = 0;
	for(hash_map<long,float,hash<long>,eqlong>::iterator rditr = tr->details->contact_energy_by_residue.begin();
			rditr != tr->details->contact_energy_by_residue.end(); rditr++){
		if(residue_energy1.count(rditr->first) == 0){
			(*num_new_contacts)++;
			(*dock_score) += rditr->second;
		}
	}
	*num_contacts = *num_contacts + *num_new_contacts;
	tr->details->num_new_contacts = *num_new_contacts;

	//*out << " #contacts " << *num_contacts << " ";

	*vdw_energy = 0;
	for(hash_map<int, short, hash<int>,eqint>::iterator ritr = lresidue_vdw_energy1.begin();
			ritr != lresidue_vdw_energy1.end(); ritr++){
		int residue = ritr->first;
		if(tr->details->lresidue_vdw_energy.count(residue) > 0)
			*vdw_energy += maximum(ritr->second,(tr->details->lresidue_vdw_energy)[residue]);
	}

	for(hash_map<int, short, hash<int>,eqint>::iterator ritr = rresidue_vdw_energy1.begin();
			ritr != rresidue_vdw_energy1.end(); ritr++){
		int residue = ritr->first;
		if(tr->details->rresidue_vdw_energy.count(residue) > 0)
			*vdw_energy += maximum(ritr->second,(tr->details->rresidue_vdw_energy)[residue]);
	}
	//*out << "done";
}

MultiTransformation::MultiTransformation(char* buf){
	stringstream ss (stringstream::in | stringstream::out);
	ss << buf;

	ss >> num_transformations;

	stringstream ss2 (stringstream::in | stringstream::out);
	long trid;
	for(int i = 0 ; i < num_transformations ; i++){
		ss >> trid;
		transformation_ids.push_back(trid);
		if(i > 0)
			ss2 << " ";
		ss2 << trid;
	}
	int idcsize = sizeof(long)*num_transformations + 32; 
	char idc[idcsize];
	ss2.getline(idc,idcsize);
	id = *(new string(idc));

	//*out << num_transformations << " " << id << endl;

	ss >> votes;
	ss >> num_contacts;
	ss >> eVdw;
	ss >> zscore;
	ss >> curvature_score;
	ss >> cscore;
	ss >> vdw_energy;
	ss >> residues_in_surface;
	ss >> residues_in_intermediate;
	ss >> residues_in_interior;
	ss >> num_neighbors;

	/*ss >> n;
	for(int i = 0 ; i < n ; i++){
		long index;
		float d;
		ss >> index;
		ss >> d;
		residue_distances[index] = d;
	}

	ss >> n;
	for(int i = 0 ; i < n ; i++){
		int index;
		float d;
		ss >> index;
		ss >> d;
		lresidue_vdw_energy[index] = d;
	}

	ss >> n;
	for(int i = 0 ; i < n ; i++){
		int index;
		float d;
		ss >> index;
		ss >> d;
		rresidue_vdw_energy[index] = d;
	}*/

	float x,y,z;
	ss >> x; ss >> y; ss >> z;
	delta_r = new Vector(x,y,z);
	ss >> rmsd; ss >> lrmsd; ss >> irmsd; ss >> frac_native_contacts_predicted; 

	contacts = NULL;
}

void MultiTransformation::print_details(ostream *out){
	if(transformations.size() > 0)
		*out << transformations.size() << " " << id << " ";
	else{
		*out << num_transformations << " ";
		for(vector<long>::iterator itr = transformation_ids.begin(); itr != transformation_ids.end(); itr++)
			*out << *itr << " ";
	}
	*out << votes << " " << num_contacts << " " << eResiduepair << " ";
	*out << zscore << " ";
	*out << curvature_score << " ";
	*out << cscore << " ";
	*out << vdw_energy << " " << residues_in_surface << " " << residues_in_intermediate
			<< " " << residues_in_interior << " " ;
	*out << num_neighbors << "\t";

	/*
	hash_map<long,float,hash<long>,eqlong> residue_distances;
	hash_map<int, float, hash<int>,eqint> lresidue_vdw_energy;
	hash_map<int, float, hash<int>,eqint> rresidue_vdw_energy;
	compute_details(&residue_distances, &lresidue_vdw_energy , &rresidue_vdw_energy);

	 *out << residue_distances.size() << " ";
	for(hash_map<long,float,hash<long>,eqlong>::iterator rditr = residue_distances.begin();
		rditr != residue_distances.end(); rditr++)
	 *out << rditr->first << " " << rditr->second << " ";

	 *out << lresidue_vdw_energy.size() << " ";
	for(hash_map<int, float, hash<int>,eqint>::iterator ritr = lresidue_vdw_energy.begin();
		ritr != lresidue_vdw_energy.end(); ritr++)
	 *out << ritr->first << " " << ritr->second << " ";

	 *out << rresidue_vdw_energy.size() << " ";
	for(hash_map<int, float, hash<int>,eqint>::iterator ritr = rresidue_vdw_energy.begin();
		ritr != rresidue_vdw_energy.end(); ritr++)
	 *out << ritr->first << " " << ritr->second << " ";
	 */

	if(delta_r != NULL)
		*out << delta_r->x << " " << delta_r->y << " " << delta_r->z << " ";
	else
		*out << "0 0 0 ";
	*out << rmsd << " " << lrmsd << " " << irmsd << " " << frac_native_contacts_predicted; 

	*out << endl;
}

void MultiTransformation::print_results(ostream *out){
	*out << delta_r->norm() << "\t";
	*out << rmsd << "\t" << lrmsd << "\t";
	*out << irmsd << "\t" << frac_native_contacts_predicted << "\t" << num_contacts << "\t" 
			<< eResiduepair << "\t";
	*out << id << endl;
}

void Transformation::update(Vector ut, Vector uex, Vector uey, int num_matching_points){
	translation->x = ut.x;
	translation->y = ut.y;
	translation->z = ut.z;

	ex->x = uex.x;
	ex->y = uex.y;
	ex->z = uex.z;

	ey->x = uey.x;
	ey->y = uey.y;
	ey->z = uey.z;

	Vector uez = ex->cross(ey);
	ez->x = uez.x;
	ez->y = uez.y;
	ez->z = uez.z;

	num_contacts = num_matching_points;
}


//void Reference_Frame::write_as_pdb(Complex *c, string chains, char chain_start, bool drop_hydrogen, string filename, bool regular){
void Reference_Frame::write_as_pdb(Complex *c, string chains, bool drop_hydrogen, string filename, bool regular){
	fstream pdbout;
	pdbout.open(filename.c_str(), fstream::out);
	pdbout << setiosflags(ios::fixed) << setprecision(3);

	char prevchain=-1;
	for( int i = 0; i < c->num_atoms; i++){
		Atom *a = c->atom[i];
		if(!drop_hydrogen || (drop_hydrogen && (a->name).c_str()[0] != 'H')){
			//cout << c->pdbcode << a->index << "\t" << a->aacindex << " " << a->aaindex << endl;cout.flush();
			Aminoacid *aa = c->aminoacid[a->monocindex];
			if((chains=="-" || chains.find(aa->chain)!=string::npos)){
				//FORMAT('ATOM',I7,2X,A3,1X,A3,I6,4X,3F8.3,2(1X,F5.2))
				float coord[3];
				Vector v;
				if(regular)
					v = transform(*(a->position));
				else
					v = inverse_transform(*(a->position));
				coord[0] = v.x;
				coord[1] = v.y;
				coord[2] = v.z;
				char ci[3][5];
				char cf[3][7];
				for(int i=0;i<3;i++){
					int c = (int) coord[i];
					sprintf(ci[i],"%4d",c);
					float f = coord[i] - c;
					sprintf(cf[i],"%6.3f",f);
					if(c==0 && f < 0)	ci[i][2] = '-';
				}

				char buf[80];
				if(aa->name.size() > 3)
					sprintf(buf,"ATOM  %5d %-4s%-4s %c%4s    %s.%s%s.%s%s.%s %c%6.4f %6.4f",
							(a->cindex)+1,(a->name).c_str(),(aa->name).c_str(),aa->chain,(aa->index).c_str(),
							ci[0],&cf[0][3],ci[1],&cf[1][3],ci[2],&cf[2][3],
							((a->charge < 0)?'-':' '),((a->charge)<0?(0-a->charge):a->charge), a->radius);
				else
					sprintf(buf,"ATOM  %5d  %-4s%-4s%c%4s    %s.%s%s.%s%s.%s %c%6.4f %6.4f",
							(a->cindex)+1,(a->name).c_str(),(aa->name).c_str(),aa->chain,(aa->index).c_str(),
							ci[0],&cf[0][3],ci[1],&cf[1][3],ci[2],&cf[2][3],
							((a->charge < 0)?'-':' '),((a->charge)<0?(0-a->charge):a->charge), a->radius);

				if(prevchain == -1)	prevchain = aa->chain;
				else if( prevchain != aa->chain){
					prevchain = aa->chain;
					pdbout << "TER" << endl;
				}

				pdbout << buf << endl;
			}
		}
	}
	pdbout << "TER" << endl;
	pdbout.close();
}

void Transformation::write_binary(ostream *out,unsigned short type){
	out->write((char*) &frame_number,sizeof(long));

	Vector v = *translation;
	out->write((char*) &(v.x),sizeof(float));
	out->write((char*) &(v.y),sizeof(float));
	out->write((char*) &(v.z),sizeof(float));

	v = *ex;
	out->write((char*) &(v.x),sizeof(float));
	out->write((char*) &(v.y),sizeof(float));
	out->write((char*) &(v.z),sizeof(float));

	v = *ey;
	out->write((char*) &(v.x),sizeof(float));
	out->write((char*) &(v.y),sizeof(float));
	out->write((char*) &(v.z),sizeof(float));

	out->write((char*) &votes,sizeof(float));
	out->write((char*) &curvature_score,sizeof(float));

	out->write((char*) &eVdw,sizeof(float));
	out->write((char*) &num_contacts,sizeof(unsigned short));
	out->write((char*) &eResiduepair,sizeof(float));
	out->write((char*) &eElectrostatic,sizeof(float));
	out->write((char*) &sEvolution_interface,sizeof(float));
	out->write((char*) &eSolvation,sizeof(float));
	out->write((char*) &num_neighbors,sizeof(unsigned short));	

	if(type == TN_VERIFY || type == TN_PROTPROT_DETAILED_SCORE_VERIFY || type == TN_PROTRNA_DETAILED_SCORE_VERIFY){
		out->write((char*) &(vmetrics->rmsd),sizeof(float));
		out->write((char*) &(vmetrics->lrmsd),sizeof(float));
		out->write((char*) &(vmetrics->irmsd),sizeof(float));
		out->write((char*) &(vmetrics->frac_native_contacts_predicted),sizeof(float));
		out->write((char*) &(vmetrics->frac_nonnative_contacts_predicted),sizeof(float));
		out->write((char*) &(vmetrics->delta_r),sizeof(float));
		out->write((char*) &(vmetrics->delta_U),sizeof(float));
	}

	if(type == TN_PROTPROT_DETAILED_SCORE_VERIFY){
		out->write((char*) &num_clashes,sizeof(unsigned short));
		out->write((char*) &num_bbclashes,sizeof(unsigned short));
		out->write((char*) &delta_sasa,sizeof(float));
		out->write((char*) &freeE_rigidmotion,sizeof(float));

		ProtProtDetailedScoringMetrics *detailed_scores = (ProtProtDetailedScoringMetrics*) this->detailed_scores;
		out->write((char*) &(detailed_scores->points_vs_distance[0]),sizeof(unsigned short)*NUM_DTRANSFORM_DIVISIONS);
		out->write((char*) &(detailed_scores->atoms_vs_distance[0]),sizeof(unsigned short)*NUM_DTRANSFORM_DIVISIONS);
		for(int i = 0; i < NUM_DTRANSFORM_DIVISIONS; i++)
			out->write((char*) &(detailed_scores->grid_contacts[i][0]),sizeof(unsigned short)*NUM_DTRANSFORM_DIVISIONS);

		out->write((char*) &(detailed_scores->delta_vdwsa[0]),sizeof(float)*NUM_RESIDUE_TYPES);
		out->write((char*) &(detailed_scores->delta_sasa[0]),sizeof(float)*NUM_RESIDUE_TYPES);
		//out->write((char*) &(detailed_scores->delta_vol[0]),sizeof(float)*NUM_RESIDUE_TYPES);
		//out->write((char*) &(detailed_scores->delta_sevol[0]),sizeof(float)*NUM_RESIDUE_TYPES);
		int num_residue_types = NUM_RESIDUE_TYPES;
		/*{
			for(int i = 0; i < num_residue_types; i++)
				out->write((char*) &(detailed_scores->residue_contacts_core[i][0]),sizeof(unsigned short)*num_residue_types);
			for(int i = 0; i < num_residue_types; i++)
				out->write((char*) &(detailed_scores->residue_contacts_rim[i][0]),sizeof(unsigned short)*num_residue_types);
		}*/
		{
			for(int i = 0; i < num_residue_types; i++)
				out->write((char*) &(detailed_scores->residue_contacts_core[i][0]),sizeof(float)*num_residue_types);
			for(int i = 0; i < num_residue_types; i++)
				out->write((char*) &(detailed_scores->residue_contacts_rim[i][0]),sizeof(float)*num_residue_types);
		}
		for(int i = 0; i < NUM_ENTROPY_DIVISIONS; i++)
			out->write((char*) &(detailed_scores->conserv_contacts[i][0]),sizeof(unsigned short)*NUM_ENTROPY_DIVISIONS);

		out->write((char*) &(detailed_scores->elec[0]),sizeof(float)*NUM_ELEC_FEATURES);

		num_residue_types = NUM_COARSE_RTYPES;
		for(int i = 0; i < num_residue_types; i++)
			for(int j = 0; j < num_residue_types; j++)
				out->write((char*) &(detailed_scores->threebody_contacts[i][j][0]),sizeof(unsigned short)*num_residue_types);

		out->write((char*) &(detailed_scores->delta_vdwsa_atom[0]),sizeof(float)*NUM_ATOM_TYPES);
		out->write((char*) &(detailed_scores->delta_sasa_atom[0]),sizeof(float)*NUM_ATOM_TYPES);
		//out->write((char*) &(detailed_scores->delta_vol_atom[0]),sizeof(float)*NUM_ATOM_TYPES);
		//out->write((char*) &(detailed_scores->delta_sevol_atom[0]),sizeof(float)*NUM_ATOM_TYPES);
		for(int i = 0; i < NUM_ATOM_TYPES; i++)
			if(NUM_ATOM_DISTANCE_DIVISIONS <= 1)
				out->write((char*) &(detailed_scores->atom_contacts[i][0]),sizeof(float)*NUM_ATOM_TYPES);
			else {
				for(int j = 0; j < NUM_ATOM_TYPES; j++)
					out->write((char*) &(detailed_scores->atom_dcontacts[i][j][0]),sizeof(float)*NUM_ATOM_DISTANCE_DIVISIONS);
			}

		out->write((char*) &(detailed_scores->atom_vdw_repulsion[0]),sizeof(unsigned short)*NUM_ATOM_VDWREPULSION_DIVISIONS);
	} else if(type == TN_PROTRNA_DETAILED_SCORE_VERIFY){
		out->write((char*) &num_clashes,sizeof(unsigned short));
		out->write((char*) &num_bbclashes,sizeof(unsigned short));
		out->write((char*) &delta_sasa,sizeof(float));
		out->write((char*) &freeE_rigidmotion,sizeof(float));

		ProtRnaDetailedScoringMetrics *detailed_scores = (ProtRnaDetailedScoringMetrics*) this->detailed_scores;
		for(int i = 0; i < NUM_RESIDUE_TYPES; i++){
			out->write((char*) &(detailed_scores->coarse_contacts[i][0]),sizeof(float)*NUM_RNA_PARTICLE_TYPES);
		}
		/*for(int i = 0; i < NUM_RESIDUE_TYPES; i++)
			for(int j = 0; j < NUM_RNA_PARTICLE_TYPES; j++)
				cout << frame_number << " wbinary " << i << " " << j << " " << detailed_scores->coarse_contacts[i][j] << endl;*/
	}
}

/*
 * Every time this is changed, remember to change the clustering of decoys and creation of inequalities
 */
void Transformation::print_details(ostream *out, unsigned short type){
	if(type == TN_VERIFY || type == TN_PROTPROT_DETAILED_SCORE_VERIFY || type == TN_PROTRNA_DETAILED_SCORE_VERIFY){
		*out << vmetrics->rmsd << " " << vmetrics->lrmsd << " " << vmetrics->irmsd << " ";
		*out << vmetrics->frac_native_contacts_predicted << " ";
		*out << vmetrics->frac_nonnative_contacts_predicted << " ";
		*out << vmetrics->delta_r << " " << vmetrics->delta_U << " ";
	}
	*out << frame_number << " ";

	if(type != TN_VERIFY){
		// translation ex ey scale votes #lresidues lresidues #rresidues rresidues
		Vector v = *translation;
		*out << v.x << " " << v.y << " " << v.z << " ";			
		v = *ex;
		*out << v.x << " " << v.y << " " << v.z << " ";
		v = *ey;
		*out << v.x << " " << v.y << " " << v.z << " ";

		*out << votes << " " ; //<< curvature_score << "\t";

		*out << num_contacts << " " << num_clashes << " " << num_bbclashes << " "
				<< eVdw << " " << eResiduepair << " " << eElectrostatic << " "
				<< eSolvation << " " << delta_sasa << " " << freeE_rigidmotion << " " << freeE_sidechain << " "
				<< sEvolution_interface << " " << num_neighbors << "\t";

		if(type == TN_PROTPROT_DETAILED_SCORE_VERIFY){
			ProtProtDetailedScoringMetrics *detailed_scores = (ProtProtDetailedScoringMetrics*) this->detailed_scores;

			for(unsigned short i = 0; i < NUM_DTRANSFORM_DIVISIONS; i++)
				for(unsigned short j = i; j < NUM_DTRANSFORM_DIVISIONS; j++)
					*out << detailed_scores->grid_contacts[i][j] << " ";

			for(unsigned short i = 0; i < 20/*NUM_RESIDUE_TYPES*/; i++)
				*out << detailed_scores->delta_sasa[i] << " ";

			int num_residue_types = NUM_RESIDUE_TYPES;				
			for(unsigned short i = 0; i < num_residue_types; i++)
				for(unsigned short j = i; j < num_residue_types; j++)
					if(i % NUM_RESIDUE_TYPES < 20 && j % NUM_RESIDUE_TYPES < 20)
						*out << detailed_scores->residue_contacts_core[i][j] << " ";

			for(unsigned short i = 0; i < num_residue_types; i++)
				for(unsigned short j = i; j < num_residue_types; j++)
					if(i % NUM_RESIDUE_TYPES < 20 && j % NUM_RESIDUE_TYPES < 20)
						*out << detailed_scores->residue_contacts_rim[i][j] << " ";

			for(unsigned short i = 0; i < 20; i++)
				*out << detailed_scores->residue_contacts_core[NTER][i] << " ";
			for(unsigned short i = 0; i < 20; i++)
				*out << detailed_scores->residue_contacts_core[CTER][i] << " ";
			*out << detailed_scores->residue_contacts_core[NTER][NTER] << " " << detailed_scores->residue_contacts_core[NTER][CTER]
			                                                                                                                  << " " << detailed_scores->residue_contacts_core[CTER][CTER] << " ";

			for(unsigned short i = 0; i < NUM_ENTROPY_DIVISIONS; i++)
				for(unsigned short j = i; j < NUM_ENTROPY_DIVISIONS; j++)
					*out << detailed_scores->conserv_contacts[i][j] << " ";

			for(unsigned short i = 0; i < NUM_ELEC_FEATURES; i++)
				*out << detailed_scores->elec[i] << " ";

			num_residue_types = NUM_COARSE_RTYPES;
			//num_residue_types = NUM_RESIDUE_TYPES;
			for(int i = 0; i < num_residue_types; i++)
				for(int j = i ; j < num_residue_types; j++)
					for(int k = j ; k < num_residue_types; k++)
						*out << detailed_scores->threebody_contacts[i][j][k] << " ";

			for(unsigned short i = 0; i < 20; i++)
				*out << detailed_scores->delta_vdwsa[i] << " ";
			/*for(unsigned short i = 0; i < 20; i++)
			 *out << detailed_scores->delta_vol[i] << " ";
			for(unsigned short i = 0; i < 20; i++)
			 *out << detailed_scores->delta_sevol[i] << " ";
			/*for(unsigned short i = 0; i < NUM_ATOM_TYPES; i++){
			 *out << detailed_scores->delta_sasa_atom[i] << " ";
				//cout << frame_number << " " << detailed_scores->delta_sasa_atom[i] << endl;
			}*/

			for(unsigned short i = 0; i < NUM_ATOM_TYPES; i++)
				*out << detailed_scores->delta_vdwsa_atom[i] << " ";

			/*	for(unsigned short i = 0; i < NUM_ATOM_TYPES; i++)
			 *out << detailed_scores->delta_vol_atom[i] << " ";
			for(unsigned short i = 0; i < NUM_ATOM_TYPES; i++)
			 *out << detailed_scores->delta_sevol_atom[i] << " ";*/

			eT32S3=0;	
			for(unsigned short i = 0; i < NUM_ATOM_TYPES; i++)
				for(unsigned short j = i; j < NUM_ATOM_TYPES; j++)
					if(NUM_ATOM_DISTANCE_DIVISIONS <=1)
						*out << detailed_scores->atom_contacts[i][j] << " ";
					else
						for(int k = 0; k < NUM_ATOM_DISTANCE_DIVISIONS; k++){
							*out << detailed_scores->atom_dcontacts[i][j][k] << " ";
							//eT32S3 += (detailed_scores->atom_dcontacts[i][j][k]*atom32_dpotential[i][j][k]);
						}
			//if(NUM_ATOM_DISTANCE_DIVISIONS > 1)	*out << eT32S3 << " ";

			for(unsigned short i=0; i < NUM_ATOM_VDWREPULSION_DIVISIONS; i++)
				*out << detailed_scores->atom_vdw_repulsion[i] << " ";
		} else if(type == TN_PROTRNA_DETAILED_SCORE_VERIFY){
			ProtRnaDetailedScoringMetrics *detailed_scores = (ProtRnaDetailedScoringMetrics*) this->detailed_scores;

			for(unsigned short i = 0; i < NUM_RESIDUE_TYPES; i++)
				for(unsigned short j = 0; j <NUM_RNA_PARTICLE_TYPES ; j++)
					if(i % NUM_RESIDUE_TYPES < 20)
						*out << detailed_scores->coarse_contacts[i][j] << " " ;
		}
	} 

	*out << endl;
}

/*
 * Every time this is changed, remember to change the clustering of decoys and creation of inequalities
 */
void Transformation::print_pie_score(unsigned short type,string paramfile){

	// PIE parameters vector
	double pieParams[NUM_PIE_PARAMS], pie_energy = 0.0;
	FILE *fp;
	int count=0;

	/* Get the PIE potential parameters */
	if ((fp = fopen(paramfile.c_str(), "r")) == NULL) {
		printf("Cannot open PIE parameters file.\n");
		return;
	}

	/* Initialize the parameters matrix */
	for (int i = 0; i < NUM_PIE_PARAMS; i++)
	{	fscanf(fp, "%lf", &pieParams[i]); }

	fclose(fp);

	pie_energy += (eVdw - 9*sEvolution_interface)*pieParams[count++] ; // van der Waals energy

	if(type == TN_PROTPROT_DETAILED_SCORE_VERIFY){
		ProtProtDetailedScoringMetrics *detailed_scores = (ProtProtDetailedScoringMetrics*) this->detailed_scores;

		int num_residue_types = NUM_RESIDUE_TYPES;
		for(unsigned short i = 0; i < num_residue_types; i++)
			for(unsigned short j = i; j < num_residue_types; j++)
				if(i % NUM_RESIDUE_TYPES < 20 && j % NUM_RESIDUE_TYPES < 20)
					pie_energy += detailed_scores->residue_contacts_core[i][j]*pieParams[count++];

		for(unsigned short i = 0; i < 20; i++)
			pie_energy += detailed_scores->residue_contacts_core[NTER][i]*pieParams[count++];
		for(unsigned short i = 0; i < 20; i++)
			pie_energy += detailed_scores->residue_contacts_core[CTER][i]*pieParams[count++];;

		pie_energy += detailed_scores->residue_contacts_core[NTER][NTER]*pieParams[count++];
		pie_energy += detailed_scores->residue_contacts_core[NTER][CTER]*pieParams[count++];
		pie_energy += detailed_scores->residue_contacts_core[CTER][CTER]*pieParams[count++];

	}

	printf("%.6lf\n",pie_energy);

}

float Transformation::compute_score(int function){
	float score=0;
	int index;

	ProtProtDetailedScoringMetrics *detailed_scores = (ProtProtDetailedScoringMetrics*) this->detailed_scores;
	if(function == RESIDUE_CONTACT_3B){
		index = 35;
		int num_residue_types = NUM_RESIDUE_TYPES;
		for(unsigned short i = 0; i < num_residue_types; i++)
			for(unsigned short j = i; j < num_residue_types; j++)
				if(i % NUM_RESIDUE_TYPES < 20 && j % NUM_RESIDUE_TYPES < 20)
					score += (detailed_scores->residue_contacts_core[i][j]+detailed_scores->residue_contacts_rim[i][j]) * score_param[index++];

		index = 272;			
		num_residue_types = NUM_COARSE_RTYPES;
		for(int i = 0; i < num_residue_types; i++)
			for(int j = i ; j < num_residue_types; j++)
				for(int k = j ; k < num_residue_types; k++)
					score += detailed_scores->threebody_contacts[i][j][k] * score_param[index++];

		index = 307;		
		for(unsigned short i = 0; i < 20; i++)
			score += detailed_scores->delta_vdwsa[i] * score_param[index++];
	} else if(function == ATOM20_CONTACT){
		index = 20;
		for(int j=0; j< 20; j++)
			for(int k=j; k< 20; k++)
				score += detailed_scores->atom_contacts[j][k] * score_param[index++];
	} else if(function == VDW_RES_BKBN_SPLINE){
		score = vdw_weight * (eVdw - 9.0* eVdw_repulsion);
		for(int j=0; j< 22; j++)
			for(int k=j; k< 22; k++)
				score += detailed_scores->residue_contacts_core[j][k] *residue_bkbn_potential[j][k];
	} /*else {
		for(int j=0; j< 20; j++){
			for(int k=j; k< 20; k++){
				score += detailed_scores->atom_contacts[j][k]*atom20_potential[j][k];
			}
			score += detailed_scores->delta_vdwsa_atom[j]*atom20_potential[20][j];
		}
	}*/

	return score;
}


bool Object::filter3(Object* ligand, Transformation *tr){
	ligand->compute_energy_approx3(this,ligand,tr,(ProtProtDetailedScoringMetrics*)tr_details);
	if(tr->num_contacts >= MIN_RESIDUE_CONTACTS(1)
			// && tr->eResiduepair <= MAX_CONTACT_ENERGY(1)
	)
		return true;
	else
		return false; 
}

ProtProtDetailedScoringMetrics::ProtProtDetailedScoringMetrics(){
	type = PROT_PROT_DETAILS;
	int num_residue_types = NUM_RESIDUE_TYPES;
	/*{
		residue_contacts_core = (unsigned short **) malloc(sizeof(unsigned short*)*num_residue_types);
		residue_contacts_rim = (unsigned short **) malloc(sizeof(unsigned short*)*num_residue_types);
		for(int i = 0; i < num_residue_types; i++){
			residue_contacts_core[i] = (unsigned short *) malloc(sizeof(unsigned short)*num_residue_types);
			residue_contacts_rim[i] = (unsigned short *) malloc(sizeof(unsigned short)*num_residue_types);
		}
	}*/
	{
		residue_contacts_core = (float **) malloc(sizeof(float*)*num_residue_types);
		residue_contacts_rim = (float **) malloc(sizeof(float*)*num_residue_types);
		for(int i = 0; i < num_residue_types; i++){
			residue_contacts_core[i] = (float *) malloc(sizeof(float)*num_residue_types);
			residue_contacts_rim[i] = (float *) malloc(sizeof(float)*num_residue_types);
		}
	}
	conserv_contacts = (unsigned short **) malloc(sizeof(unsigned short*)*NUM_ENTROPY_DIVISIONS);
	for(int i = 0; i < NUM_ENTROPY_DIVISIONS; i++)
		conserv_contacts[i] = (unsigned short *) malloc(sizeof(unsigned short)*NUM_ENTROPY_DIVISIONS);

	atom_contacts = NULL; atom_dcontacts = NULL;
	if(NUM_ATOM_DISTANCE_DIVISIONS <= 1){
		atom_contacts = (float **) malloc(sizeof(float*)*NUM_ATOM_TYPES);
		for(int i = 0; i < NUM_ATOM_TYPES; i++)
			atom_contacts[i] = (float *) malloc(sizeof(float)*NUM_ATOM_TYPES);
	} else {
		atom_dcontacts = (float ***) malloc(sizeof(float**)*NUM_ATOM_TYPES);
		for(int i = 0; i < NUM_ATOM_TYPES; i++){
			atom_dcontacts[i] = (float **) malloc(sizeof(float*)*NUM_ATOM_TYPES);
			for(int j = 0; j < NUM_ATOM_TYPES; j++)
				atom_dcontacts[i][j] = (float *) malloc(sizeof(float*)*NUM_ATOM_DISTANCE_DIVISIONS);
		}
	}
}

ProtProtDetailedScoringMetrics::~ProtProtDetailedScoringMetrics(){
	int num_residue_types = NUM_RESIDUE_TYPES;
	for(int i = 0; i < num_residue_types; i++){
		if(residue_contacts_core[i] != NULL)	free(residue_contacts_core[i]);
		if(residue_contacts_rim[i] != NULL)	free(residue_contacts_rim[i]);
	}
	if(residue_contacts_core != NULL)	free(residue_contacts_core);
	if(residue_contacts_rim != NULL)	free(residue_contacts_rim);

	for(int i = 0; i < NUM_ENTROPY_DIVISIONS; i++)
		if(conserv_contacts[i] != NULL)	free(conserv_contacts[i]);
	if(conserv_contacts != NULL)	free(conserv_contacts);

	if(NUM_ATOM_DISTANCE_DIVISIONS <= 1){
		for(int i = 0; i < NUM_ATOM_TYPES; i++)
			if(atom_contacts[i] != NULL)	free(atom_contacts[i]);
		if(atom_contacts != NULL)	free(atom_contacts);
	} else {
		for(int i = 0; i < NUM_ATOM_TYPES; i++){
			for(int j = 0; j < NUM_ATOM_TYPES; j++)
				if(atom_dcontacts[i][j] != NULL)	free(atom_dcontacts[i][j]);
			if(atom_dcontacts[i] != NULL)	free(atom_dcontacts[i]);
		}
		if(atom_dcontacts != NULL)	free(atom_dcontacts);
	}
}



/*
 * Receptor atoms are written as is, ligand atoms are transformed
 * chains are started at 'A' and increased alphabetically
 * 	receptor first chain is 'A'
 * 
 * WARNIING - using hash_map iterator aminoacids not sorted by index
 */
void write_as_pdb(Complex* receptor, Complex* ligand, Transformation* tr, string filename){
	string alpha_carbon = "CA";

	fstream pdbout;
	pdbout.open(filename.c_str(), fstream::out);
	pdbout << setiosflags(ios::fixed) << setprecision(3);

	int num_atoms = 1;
	char chain = 'A';
	for(hash_map<char, Molecule*, hash<char>, eqint>::iterator mitr = receptor->molecules.begin(); mitr != receptor->molecules.end(); mitr++){
		Molecule *m = mitr->second;
		int natoms = (m->atom).size();
		Atom *m_atoms[natoms + 1];
		int count = 0;
		for(hash_map<unsigned short, Atom*, hash<unsigned short>,eqint>::iterator aitr = (m->atom).begin(); aitr != (m->atom).end(); aitr++){
			m_atoms[count++] = aitr->second;
		}
		// sort atoms by index
		for(int i=0;i<natoms;i++)
			for(int j=i+1;j<natoms;j++)
				if(m_atoms[i]->cindex > m_atoms[j]->cindex){
					Atom *a = m_atoms[i];
					m_atoms[i] = m_atoms[j];
					m_atoms[j] = a;
				}
		for(int i=0;i<natoms;i++){
			Atom *a = m_atoms[i];
			Aminoacid *aa = receptor->aminoacid[a->monocindex];
			if((a->name).size() <= 3){
				//FORMAT('ATOM',I7,2X,A3,1X,A3,I6,4X,3F8.3,2(1X,F5.2))
				float coord[3];
				coord[0] = a->position->x;
				coord[1] = a->position->y;
				coord[2] = a->position->z;
				char ci[3][5];
				char cf[3][7];
				for(int i=0;i<3;i++){
					int c = (int) coord[i];
					sprintf(ci[i],"%4d",c);
					float f = coord[i] - c;
					sprintf(cf[i],"%6.3f",f);
					if(c==0 && f < 0)	ci[i][2] = '-';
				}

				char buf[256], c = aa->index.c_str()[aa->index.length()-1];
				if(aa->name.size() > 3)
					if(c>='A' && c <= 'Z')
						sprintf(buf,"ATOM  %5d %-4s%-4s %c%5s   %s.%s%s.%s%s.%s %c%6.4f %6.4f",
								num_atoms++,(a->name).c_str(),(aa->name).c_str(),chain,(aa->index).c_str(),
								ci[0],&cf[0][3],ci[1],&cf[1][3],ci[2],&cf[2][3],
								((a->charge < 0)?'-':' '),((a->charge)<0?(0-a->charge):a->charge), a->radius);
					else
						sprintf(buf,"ATOM  %5d %-4s%-4s %c%4s    %s.%s%s.%s%s.%s %c%6.4f %6.4f",
								num_atoms++,(a->name).c_str(),(aa->name).c_str(),chain,(aa->index).c_str(),
								ci[0],&cf[0][3],ci[1],&cf[1][3],ci[2],&cf[2][3],
								((a->charge < 0)?'-':' '),((a->charge)<0?(0-a->charge):a->charge), a->radius);
				else
					if(c>='A' && c <= 'Z')
						sprintf(buf,"ATOM  %5d  %-4s%-4s%c%5s   %s.%s%s.%s%s.%s %c%6.4f %6.4f",
								num_atoms++,(a->name).c_str(),(aa->name).c_str(),chain,(aa->index).c_str(),
								ci[0],&cf[0][3],ci[1],&cf[1][3],ci[2],&cf[2][3],
								((a->charge < 0)?'-':' '),((a->charge)<0?(0-a->charge):a->charge), a->radius);
					else
						sprintf(buf,"ATOM  %5d  %-4s%-4s%c%4s    %s.%s%s.%s%s.%s %c%6.4f %6.4f",
								num_atoms++,(a->name).c_str(),(aa->name).c_str(),chain,(aa->index).c_str(),
								ci[0],&cf[0][3],ci[1],&cf[1][3],ci[2],&cf[2][3],
								((a->charge < 0)?'-':' '),((a->charge)<0?(0-a->charge):a->charge), a->radius);
				pdbout << buf << endl;
			}
		}
		chain++;
	}

	for(hash_map<char, Molecule*, hash<char>, eqint>::iterator mitr = ligand->molecules.begin(); mitr != ligand->molecules.end(); mitr++){
		Molecule *m = mitr->second;
		int natoms = (m->atom).size();
		Atom *m_atoms[natoms + 1];
		int count = 0;
		for(hash_map<unsigned short, Atom*, hash<unsigned short>,eqint>::iterator aitr = (m->atom).begin(); aitr != (m->atom).end(); aitr++){
			m_atoms[count++] = aitr->second;
		}
		// sort atoms by index
		for(int i=0;i<natoms;i++)
			for(int j=i+1;j<natoms;j++)
				if(m_atoms[i]->cindex > m_atoms[j]->cindex){
					Atom *a = m_atoms[i];
					m_atoms[i] = m_atoms[j];
					m_atoms[j] = a;
				}
		for(int i=0;i<natoms;i++){
			Atom *a = m_atoms[i];
			Aminoacid *aa = ligand->aminoacid[a->monocindex];		
			if((a->name).size() <= 3){
				Vector v = ((Reference_Frame*) tr)->inverse_transform(*(a->position));
				float coord[3];
				coord[0] = v.x;
				coord[1] = v.y;
				coord[2] = v.z;
				char ci[3][5];
				char cf[3][7];
				for(int i=0;i<3;i++){
					int c = (int) coord[i];
					sprintf(ci[i],"%4d",c);
					float f = coord[i] - c;
					sprintf(cf[i],"%6.3f",f);
					if(c==0 && f < 0)	ci[i][2] = '-';
				}

				char buf[256], c = aa->index.c_str()[aa->index.length()-1];
				if(aa->name.size() > 3)
					if(c>='A' && c <= 'Z')
						sprintf(buf,"ATOM  %5d %-4s%-4s %c%5s   %s.%s%s.%s%s.%s %c%6.4f %6.4f",
								num_atoms++,(a->name).c_str(),(aa->name).c_str(),chain,(aa->index).c_str(),
								ci[0],&cf[0][3],ci[1],&cf[1][3],ci[2],&cf[2][3],
								((a->charge < 0)?'-':' '),((a->charge)<0?(0-a->charge):a->charge), a->radius);
					else
						sprintf(buf,"ATOM  %5d %-4s%-4s %c%4s    %s.%s%s.%s%s.%s %c%6.4f %6.4f",
								num_atoms++,(a->name).c_str(),(aa->name).c_str(),chain,(aa->index).c_str(),
								ci[0],&cf[0][3],ci[1],&cf[1][3],ci[2],&cf[2][3],
								((a->charge < 0)?'-':' '),((a->charge)<0?(0-a->charge):a->charge), a->radius);
				else
					if(c>='A' && c <= 'Z')
						sprintf(buf,"ATOM  %5d  %-4s%-4s%c%5s   %s.%s%s.%s%s.%s %c%6.4f %6.4f",
								num_atoms++,(a->name).c_str(),(aa->name).c_str(),chain,(aa->index).c_str(),
								ci[0],&cf[0][3],ci[1],&cf[1][3],ci[2],&cf[2][3],
								((a->charge < 0)?'-':' '),((a->charge)<0?(0-a->charge):a->charge), a->radius);
					else
						sprintf(buf,"ATOM  %5d  %-4s%-4s%c%4s    %s.%s%s.%s%s.%s %c%6.4f %6.4f",
								num_atoms++,(a->name).c_str(),(aa->name).c_str(),chain,(aa->index).c_str(),
								ci[0],&cf[0][3],ci[1],&cf[1][3],ci[2],&cf[2][3],
								((a->charge < 0)?'-':' '),((a->charge)<0?(0-a->charge):a->charge), a->radius);
				pdbout << buf << endl;
			}
		}
		chain++;
	}
	pdbout.close();
}


void optimize_sidechains(Complex* receptor, Complex* ligand, Transformation* tr){
	char scwrl_in_file[128];
	sprintf(scwrl_in_file,"%s/%s%ld",scratch_dir,SCWRL_INPUT_FILE,tr->frame_number); 

	write_as_pdb(receptor, ligand, tr, string(scwrl_in_file));

	char scwrl_out_file[128];
	sprintf(scwrl_out_file,"models/%s%ld.pdb",SCWRL_OUTPUT_FILE,tr->frame_number); 
	char scwrl_log_file[128];
	sprintf(scwrl_log_file,"%s/%s%ld",scratch_dir,SCWRL_LOG_FILE,tr->frame_number); 

	/*pthread_t scwrl_thread;
	int iret = pthread_create( &scwrl_thread, NULL, run_scwrl, &(tr->frame_number));*/
	char scwrl_command[512];
	sprintf(scwrl_command, "%s/%s -i %s -o %s > %s",piedock_home.c_str(),SCWRL_EXECUTABLE,scwrl_in_file,
			scwrl_out_file,	scwrl_log_file);

	int iret = system(scwrl_command);
	*out << tr->frame_number << " scwrlstatus " << iret << endl;
	//pthread_join(scwrl_thread, NULL);
	//*out << "check " << endl;
}

void run_scwrl(void *fnumber){
	long frame_number = *((long *) fnumber);
	char scwrl_in_file[128];
	sprintf(scwrl_in_file,"%s%ld",SCWRL_INPUT_FILE,frame_number); 
	char scwrl_out_file[128];
	sprintf(scwrl_out_file,"%s%ld",SCWRL_OUTPUT_FILE,frame_number); 
	char scwrl_log_file[128];
	sprintf(scwrl_log_file,"%s%ld",SCWRL_LOG_FILE,frame_number); 

	int execlpstatus = execlp(string(SCWRL_EXECUTABLE).c_str(),string(" -i ").c_str(),
			scwrl_in_file,string(" -o ").c_str(),
			scwrl_out_file,string(" > ").c_str(),
			scwrl_log_file,(char*) 0);
	if(execlpstatus == -1)
		*out << "scwrl status " << errno << endl;

	return;
}

void Object::score(Object *ligand, MultiTransformation *mtr){
	stringstream ss (stringstream::in | stringstream::out);
	stringstream ss2 (stringstream::in | stringstream::out);
	ss << string((mtr->id).c_str());
	for(int i = 0; i < (mtr->transformations).size(); i++){
		if(i > 0)
			ss2 << ".";
		long tid;
		ss >> tid;
		ss2 << tid;
	}
	string mid;
	ss2 >> mid;
	mid = string(mid.c_str());
	ss2.clear();
	ss2 << "models/" << mid;
	string s;
	ss2 >> s;
	string chains = "", rchains = "", lchains = "";
	int num_rchains = (c->molecules).size(), num_lchains = (ligand->c->molecules).size();
	int num_chains = num_rchains + num_lchains;
	char chain = 'A';
	for(int i = 0 ; i < num_chains; i++){
		chains.append(1,chain);
		if(i < num_rchains)
			rchains.append(1,chain);
		else
			lchains.append(1,chain);
		chain++;
	}
	Complex *model = new Complex(s,chains, PDB);

	float model_contact_energy=0;
	hash_set<long,hash<long>,eqlong> model_interface_contacts;
	hash_map<int,int,hash<int>, eqint> rresidue_type, lresidue_type;
	char rchain = 'A';
	for(int i = 0 ; i < num_rchains; i++){
		Molecule *mr = model->molecules[rchain];
		if(mr->type == PROTEIN){
			Protein *pr = (Protein*) mr;
			for(hash_map<const char*, Aminoacid*, hash<const char*>,eqstr>::iterator raitr = pr->aminoacid.begin(); raitr != pr->aminoacid.end(); raitr++){
				Aminoacid *ra = raitr->second;
				char lchain = 'A' + num_rchains;
				for(int i = 0 ; i < num_lchains; i++){
					Molecule *ml = model->molecules[lchain];
					if(ml->type == PROTEIN){
						Protein *pl = (Protein*) ml;
						for(hash_map<const char*, Aminoacid*, hash<const char*>,eqstr>::iterator laitr = pl->aminoacid.begin(); laitr != pl->aminoacid.end(); laitr++){
							Aminoacid *la = laitr->second;
							if(Vector::distance(ra->alpha_carbon->position, la->alpha_carbon->position) < SS_CUTOFF || Vector::distance(ra->centroid, la->centroid) < SS_CUTOFF){
								if(rresidue_type.count(ra->cindex) == 0)
									rresidue_type[ra->cindex] = ra->type;
								if(lresidue_type.count(la->cindex) == 0)
									lresidue_type[la->cindex] = la->type;
								long index = (la->cindex)*MAX_ATOMS + ra->cindex;
								if(model_interface_contacts.count(index) == 0){
									model_interface_contacts.insert(index);
									model_contact_energy += residue_potential[ra->type][la->type];
								}
							}
						}
					}
					lchain++;
				}
			}
		}
		rchain++;
	}

	//compute z-score
	float sum = 0, sum_squares = 0;
	int nr = rresidue_type.size();
	int rresidues[nr+1];
	int count = 0;
	hash_map<int,int,hash<int>, eqint> rresidue_vs_array;
	for(hash_map<int,int,hash<int>, eqint>::iterator itr = rresidue_type.begin(); itr != rresidue_type.end(); itr++){
		rresidues[count] = itr->first;
		rresidue_vs_array[itr->first] = count;
		count++;
	}
	int nl = lresidue_type.size();
	int lresidues[nr+1];
	count = 0;
	hash_map<int,int,hash<int>, eqint> lresidue_vs_array;
	for(hash_map<int,int,hash<int>, eqint>::iterator itr = lresidue_type.begin(); itr != lresidue_type.end(); itr++){
		lresidues[count] = itr->first;
		lresidue_vs_array[itr->first] = count;
		count++;
	}

	for(int i = 0 ; i < ZSCORE_N; i++){
		// permute receptor indices
		int pr[nr+1];
		permute(nr,pr);
		/*out << "nr: " << nr << " ";
		for(int j = 0 ; j < nr; j++)
		 *out << pr[j] << " ";
		 *out << endl;*/
		int pl[nl+1];
		permute(nl,pl);

		float e = 0;
		for(hash_set<long,hash<long>,eqlong>::iterator citr = model_interface_contacts.begin(); citr != model_interface_contacts.end(); citr++){
			long index = *citr;
			int rindex = index % MAX_ATOMS;
			int lindex = index % MAX_ATOMS;
			int prindex = rresidues[pr[rresidue_vs_array[rindex]]];
			int plindex = lresidues[pl[lresidue_vs_array[lindex]]];
			e += residue_potential[rresidue_type[prindex]][lresidue_type[plindex]];
		}
		sum += e;
		sum_squares += (e*e);
	}
	float avg = sum/ZSCORE_N;
	float sigma = sqrt(sum_squares/ZSCORE_N  - avg*avg);

	mtr->zscore = (model_contact_energy - avg)/sigma;

	mtr->num_contacts = model_interface_contacts.size();
	mtr->eResiduepair = model_contact_energy;
	*out << model_interface_contacts.size() << " " << model_contact_energy << " " << (mtr->zscore) << endl;
	//out->flush();
	//delete model;
}

float Object::irmsd(Object *ligand, Transformation *tr1, Transformation *tr2){
	hash_set<long,hash<long>,eqlong> reference_interface_contacts, reference_interface_contacts_capri;
	for(short ci = 0; ci < c->num_aminoacids; ci++){
		Aminoacid *ra = c->aminoacid[ci];
		for(short cj = 0; cj < ligand->c->num_aminoacids; cj++){
			Aminoacid *la = ligand->c->aminoacid[cj];
			float d_ca, d_cd;
			long cindex = (ra->cindex)*MAX_ATOMS + la->cindex;

			if((ra->alpha_carbon != NULL && la->alpha_carbon != NULL && (d_ca = Vector::distance(ra->alpha_carbon->position, tr1->inverse_transform(la->alpha_carbon->position))) < SS_CUTOFF) ||
					(ra->centroid != NULL && la->centroid != NULL && (d_cd = Vector::distance(ra->centroid,  tr1->inverse_transform(la->centroid))) < SS_CUTOFF)){
				if(reference_interface_contacts.count(cindex) == 0){
					reference_interface_contacts.insert(cindex);
				}

				// defining the interface used for computing irmsd
				for(hash_map<const char*, Atom*, hash<const char*>, eqstr>::iterator raitr = ra->atom.begin(); raitr != ra->atom.end(); raitr++){
					Atom *ratom = (Atom *) raitr->second;
					if(ratom->name.c_str()[0] != 'H'){
						for(hash_map<const char*, Atom*, hash<const char*>, eqstr>::iterator laitr = la->atom.begin(); laitr != la->atom.end(); laitr++){
							Atom *latom = (Atom *) laitr->second;
							if(latom->name.c_str()[0] != 'H'){
								if(Vector::distance(ratom->position,  tr1->inverse_transform(latom->position)) <= CAPRI_INTERFACE_CUTOFF){
									reference_interface_contacts_capri.insert(cindex);
								}
							}
						}
					}
				}
			}
		}
	}

	hash_set<int,hash<int>,eqint> added_ligand_residue, added_receptor_residue;
	Vector tr1_interface_points[c->num_aminoacids + ligand->c->num_aminoacids + 1],tr2_interface_points[c->num_aminoacids + ligand->c->num_aminoacids + 1];
	int ca_point_index = 0;

	for(hash_set<long,hash<long>,eqlong>::iterator itr = reference_interface_contacts_capri.begin(); itr != reference_interface_contacts_capri.end(); itr++){
		long index = *itr;
		int rindex = index/MAX_ATOMS;

		if(added_receptor_residue.count(rindex) == 0){
			added_receptor_residue.insert(rindex);
			Aminoacid *aa = c->aminoacid[rindex];

			if(aa->atom.count("CA")> 0){
				Atom *ratom = (Atom *) aa->atom["CA"];
				tr1_interface_points[ca_point_index] = Vector(ratom->position);
				tr2_interface_points[ca_point_index] = Vector(ratom->position);
				ca_point_index++;
			}
		}
	}

	for(hash_set<long,hash<long>,eqlong>::iterator itr = reference_interface_contacts_capri.begin(); itr != reference_interface_contacts_capri.end(); itr++){
		long index = *itr;
		int lindex = index % MAX_ATOMS;

		if(added_ligand_residue.count(lindex) == 0){
			added_ligand_residue.insert(lindex);
			Aminoacid *aa = ligand->c->aminoacid[lindex];

			if(aa->atom.count("CA")> 0){
				Atom *ratom = (Atom *) aa->atom["CA"];
				tr1_interface_points[ca_point_index] = tr1->inverse_transform(ratom->position);
				tr2_interface_points[ca_point_index] = tr2->inverse_transform(ratom->position);
				ca_point_index++;
			} 
		}
	}
	int num_tr1_interface_points_ca = ca_point_index;

	// get the rmsd between the for the interface
	float irmsd = compute_rmsd(num_tr1_interface_points_ca, &tr1_interface_points[0], &tr2_interface_points[0]);
	*out << tr1->frame_number << " " << tr2->frame_number << " " << irmsd << endl;
	return irmsd;
}

/*
 * Score the reference structure - compute votes, cscore, ...
 */
Transformation** Object::score_reference(Object* ligand, Complex* reference, Complex* referenceH, string *rchains, string *lchains, hash_map<long, long, hash<long>, eqlong> *receptor_vs_reference,
		hash_map<long, long, hash<long>, eqlong> *reference_vs_receptor, hash_map<long, long, hash<long>, eqlong> *ligand_vs_reference,
		hash_map<long, long, hash<long>, eqlong> *reference_vs_ligand, int procid){
	int max_num_atoms = c->num_atoms + ligand->c->num_atoms;
	Vector receptor_points[max_num_atoms+1], ligand_points[max_num_atoms+1];
	Vector reference_points[max_num_atoms+1];

	/*Vector relative_r_reference = *(reference->molecules[ligand_chain]->center_of_mass) - *(reference->molecules[receptor_chain]->center_of_mass); 
	cout << relative_r_reference.norm() << "\t";*/ 
	if(procid == 0)	cout << receptor_vs_reference->size() << "|" << reference_vs_receptor->size() << " " << ligand_vs_reference->size() << "|" << reference_vs_ligand->size() << endl;
	//for(hash_map<long, long, hash<long>, eqlong>::iterator itr = reference_vs_receptor->begin(); itr != reference_vs_receptor->end(); itr++)
	//	cout << (long) itr->first << "|" << (long) itr->second << endl;

	int point_index = 0;
	for(int i = 0; i < c->num_aminoacids; i++){
		Aminoacid *aa = c->aminoacid[i];
		if(receptor_vs_reference->count(aa->cindex) > 0){
			int ref_residue_index = (*receptor_vs_reference)[aa->cindex];
			Aminoacid *aa_ref = reference->aminoacid[ref_residue_index];

			if(aa->atom.count("N") > 0 && aa_ref->atom.count("N") > 0){
				reference_points[point_index] = Vector(aa_ref->atom["N"]->position);
				receptor_points[point_index] = Vector(aa->atom["N"]->position);
				point_index++;
			} else if(procid==0)
				cout << point_index << " " << aa_ref << " " << aa << " N missing " << aa->name << " " << aa->index << endl;

			if(aa->atom.count("CA") > 0 && aa_ref->atom.count("CA") > 0){
				reference_points[point_index] = Vector(aa_ref->atom["CA"]->position);
				receptor_points[point_index] = Vector(aa->atom["CA"]->position);
				point_index++;
			} else if(procid==0)
				cout << point_index << " " << aa_ref << " " << aa << " CA missing " << aa->name << " " << aa->index << endl;

			if(aa->atom.count("C") > 0 && aa_ref->atom.count("C") > 0){
				reference_points[point_index] = Vector(aa_ref->atom["C"]->position);
				receptor_points[point_index] = Vector(aa->atom["C"]->position);
				point_index++;
			} else if(procid==0)
				cout << point_index << " " << aa_ref << " " << aa << " C missing " << aa->name << " " << aa->index << endl;

			if(aa->atom.count("O") > 0 && aa_ref->atom.count("O") > 0){
				reference_points[point_index] = Vector(aa_ref->atom["O"]->position);
				receptor_points[point_index] = Vector(aa->atom["O"]->position);
				point_index++;
			} else if(procid==0)
				cout << point_index << " " << aa_ref << " " << aa << " O missing " << aa->name << " " << aa->index << endl;

			/*if((aa_ref->alpha_carbon == NULL || aa->alpha_carbon == NULL ) && procid == 0)
				cout << point_index << " " << aa_ref << " " << aa << " " << aa_ref->alpha_carbon << " " << aa->alpha_carbon << " " << aa->name << " " << aa->index << endl;
			reference_points[point_index] = Vector(aa_ref->alpha_carbon->position);
			receptor_points[point_index] = Vector(aa->alpha_carbon->position);
			point_index++;*/
		}
	}
	int ligand_start = point_index;

	for(int i = 0; i < ligand->c->num_aminoacids; i++){
		Aminoacid *aa = ligand->c->aminoacid[i];
		if(ligand_vs_reference->count(aa->cindex) > 0){
			int ref_residue_index = (*ligand_vs_reference)[aa->cindex];
			Aminoacid *aa_ref = reference->aminoacid[ref_residue_index];

			if(aa->atom.count("N") > 0 && aa_ref->atom.count("N") > 0){
				reference_points[point_index] = Vector(aa_ref->atom["N"]->position);
				ligand_points[point_index - ligand_start] = Vector(aa->atom["N"]->position);
				point_index++;
			} else if(procid==0)
				cout << point_index << " " << aa_ref << " " << aa << " N missing " << aa->name << " " << aa->index << endl;

			if(aa->atom.count("CA") > 0 && aa_ref->atom.count("CA") > 0){
				reference_points[point_index] = Vector(aa_ref->atom["CA"]->position);
				ligand_points[point_index - ligand_start] = Vector(aa->atom["CA"]->position);
				point_index++;
			} else if(procid==0)
				cout << point_index << " " << aa_ref << " " << aa << " CA missing " << aa->name << " " << aa->index << endl;

			if(aa->atom.count("C") > 0 && aa_ref->atom.count("C") > 0){
				reference_points[point_index] = Vector(aa_ref->atom["C"]->position);
				ligand_points[point_index - ligand_start] = Vector(aa->atom["C"]->position);
				point_index++;
			} else if(procid==0)
				cout << point_index << " " << aa_ref << " " << aa << " C missing " << aa->name << " " << aa->index << endl;

			if(aa->atom.count("O") > 0 && aa_ref->atom.count("O") > 0){
				reference_points[point_index] = Vector(aa_ref->atom["O"]->position);
				ligand_points[point_index - ligand_start] = Vector(aa->atom["O"]->position);
				point_index++;
			} else if(procid==0)
				cout << point_index << " " << aa_ref << " " << aa << " O missing " << aa->name << " " << aa->index << endl;

			/*if((aa_ref->alpha_carbon == NULL || aa->alpha_carbon == NULL ) && procid == 0)
				cout << point_index << " " << aa_ref << " " << aa << " " << aa_ref->alpha_carbon << " " << aa->alpha_carbon << " " << aa->name << " " << aa->index << endl;
			reference_points[point_index] = Vector(aa_ref->alpha_carbon->position);
			ligand_points[point_index - ligand_start] = Vector(aa->alpha_carbon->position);
			point_index++;*/
		}
	}
	*out << "score reference - computing lrmsd " << ligand_start << " " << point_index << endl; out->flush();

	float lrmsd_ref_vs_model, lrmsd_ref_vs_orig, lrmsd_orig_vs_model;
	float rrmsd_ref_vs_model, rrmsd_ref_vs_orig, rrmsd_orig_vs_model;
	Transformation *receptor_orig_to_ref = new Transformation(new Vector(0,0,0), new Vector(1,0,0), new Vector(0,1,0), 1.0, 0,0);
	Transformation *ligand_orig_to_ref = new Transformation(new Vector(0,0,0), new Vector(1,0,0), new Vector(0,1,0), 1.0, 0, 0);
	rrmsd_ref_vs_orig = compute_rmsd(ligand_start, &reference_points[0], &receptor_points[0], receptor_orig_to_ref);
	int num_ligand_points = point_index - ligand_start;
	lrmsd_ref_vs_orig = compute_rmsd(num_ligand_points, &reference_points[ligand_start], &ligand_points[0], ligand_orig_to_ref);
	Transformation *tr = generate_transformation(this,receptor_orig_to_ref,ligand,ligand_orig_to_ref,0);
	float combined_rmsd = sqrt((rrmsd_ref_vs_orig*rrmsd_ref_vs_orig*ligand_start+lrmsd_ref_vs_orig*lrmsd_ref_vs_orig*num_ligand_points)/(num_ligand_points+ligand_start));
	if(procid == 0)	cout << "rmsd rec " << rrmsd_ref_vs_orig << "\tlig\t" << lrmsd_ref_vs_orig << "\tcombined\t" << combined_rmsd << "\t";

	hash_set<long,hash<long>,eqlong> reference_interface_contacts, reference_interface_contacts_capri;
	for(short ci = 0; ci < rchains->length(); ci++)
		for(int i = 0; i < reference->num_aminoacids; i++){
			Aminoacid *ra = reference->aminoacid[i];
			if(ra->chain == rchains->c_str()[ci]){
				for(short cj = 0; cj < lchains->length(); cj++)
					for(int j = 0; j < reference->num_aminoacids; j++){
						Aminoacid *la = reference->aminoacid[j];
						if(la->chain == lchains->c_str()[cj]){
							long cindex = (ra->cindex)*MAX_ATOMS + la->cindex;
							float d_ca, d_cd;

							if((ra->alpha_carbon != NULL && la->alpha_carbon != NULL && (d_ca = Vector::distance(ra->alpha_carbon->position, la->alpha_carbon->position)) < SS_CUTOFF) ||
									(ra->centroid != NULL && la->centroid != NULL && (d_cd = Vector::distance(ra->centroid, la->centroid)) < SS_CUTOFF)){
								if(reference_interface_contacts.count(cindex) == 0){
									reference_interface_contacts.insert(cindex);
									if(procid == 0){
										cout <<endl<<la->chain<<":"<<la->get_symbol()<<":"<<la->index << ":" << la->cindex << " " <<ra->chain<<":"<<ra->get_symbol()<<":"<<ra->index << ":" << ra->cindex
												<< " " << d_ca << " " << d_cd;
									}
								}
							}

							// defining the interface used for computing irmsd
							for(hash_map<const char*, Atom*, hash<const char*>, eqstr>::iterator raitr = ra->atom.begin(); raitr != ra->atom.end(); raitr++){
								Atom *ratom = (Atom *) raitr->second;
								if(ratom->name.c_str()[0] != 'H'){
									for(hash_map<const char*, Atom*, hash<const char*>, eqstr>::iterator laitr = la->atom.begin(); laitr != la->atom.end(); laitr++){
										Atom *latom = (Atom *) laitr->second;
										if(latom->name.c_str()[0] != 'H'){
											if(Vector::distance(ratom->position, latom->position) <= CAPRI_INTERFACE_CUTOFF){
												reference_interface_contacts_capri.insert(cindex);
											}
										}
									}
								}
							}
						}
					}
			}
		}

	if(procid == 0){
		float reference_contact_energy_residue=0,reference_contact_energy_coarse=0;
		//compute the residues in receptor and ligand in the interface
		unsigned short contacts_numgaps = 0;
		hash_set<int,hash<int>,eqint> reference_interface_ligand_residues,reference_interface_receptor_residues, ref_vs_lig_interface_missing, ref_vs_rec_interface_missing;
		for(hash_set<long,hash<long>,eqlong>::iterator itr = reference_interface_contacts_capri.begin(); itr != reference_interface_contacts_capri.end(); itr++){
			long index = *itr;
			int refrindex = index/MAX_ATOMS;
			int reflindex = index % MAX_ATOMS;
			if(reference_vs_ligand->count(reflindex) > 0){
				int laindex = (*reference_vs_ligand)[reflindex];
				reference_interface_ligand_residues.insert(laindex);
				if(reference_vs_receptor->count(refrindex) > 0){
					int raindex = (*reference_vs_receptor)[refrindex];
					Aminoacid *ra = c->aminoacid[raindex];
					Aminoacid *la = ligand->c->aminoacid[laindex];
					//cout << refrindex << " " << reflindex << " " << raindex << " " << laindex << "\t"; cout.flush();
					//cout << ra->name << " " << ra->type << " " << ra->index << "\t" << la->name << " " << la->type << " " << la->index << endl; cout.flush();  
					if(ra->type >= 0 && la->type >= 0){
						reference_contact_energy_residue += residue_potential[ra->type][la->type];
						int ss1 = SSTYPE(ra->sstructure), ss2 = SSTYPE(la->sstructure);
						int crs1 = RTYPE(ra->type) *NUM_COARSE_SSTYPES + ss1;
						int crs2 = RTYPE(la->type) *NUM_COARSE_SSTYPES + ss2;
						if(ss1 != DSSP_U && ss2 != DSSP_U){
							//reference_contact_energy_coarse += coarse_potential[crs1][crs2];
						}
					}
				} else {
					contacts_numgaps++;
				}
			} else {
				ref_vs_lig_interface_missing.insert(reflindex);
				contacts_numgaps++;
			}
			if(reference_vs_receptor->count(refrindex) > 0){
				int raindex = (*reference_vs_receptor)[refrindex];
				reference_interface_receptor_residues.insert(raindex);
				//cout << refrindex << ":" << raindex << endl;
			} else {
				ref_vs_rec_interface_missing.insert(refrindex);
			}
		}

		// changes in stability?
		/*c->compute_stability();
		ligand->c->compute_stability();
		reference->compute_contacts();
		reference->compute_stability();*/

		//cout << endl << receptor_vs_reference->size() << "|" << reference_vs_receptor->size() << " " << ligand_vs_reference->size() << "|" << reference_vs_ligand->size() << endl;
		float eInterface = 0;
		float interface_seq_stability = 0, interface_seq_stability_docked=0;
		cout << endl << "receptor residues " << endl;
		bool rinterface[c->num_aminoacids+1], linterface[ligand->c->num_aminoacids+1];
		for(int i = 0; i < c->num_aminoacids; i++)	rinterface[i] = false;
		for(int i = 0; i < ligand->c->num_aminoacids; i++)	linterface[i] = false;
		for(hash_set<int,hash<int>,eqint>::iterator itr = reference_interface_receptor_residues.begin(); itr != reference_interface_receptor_residues.end(); itr++){
			int index = *itr;
			cout << c->aminoacid[index]->chain << "|" << c->aminoacid[index]->cindex << "|" << (*receptor_vs_reference)[index] << ":" << c->aminoacid[index]->get_symbol() << ":" << reference->aminoacid[(*receptor_vs_reference)[index]]->get_symbol() 
						<< " " << c->aminoacid[index]->sequence_stability << ":" << reference->aminoacid[(*receptor_vs_reference)[index]]->sequence_stability << endl;
			eInterface += c->aminoacid[index]->eInterface;
			interface_seq_stability += c->aminoacid[index]->sequence_stability;
			interface_seq_stability_docked += reference->aminoacid[(*receptor_vs_reference)[index]]->sequence_stability;
		}
		cout << endl << "ligand residues " << endl;
		for(hash_set<int,hash<int>,eqint>::iterator itr = reference_interface_ligand_residues.begin(); itr != reference_interface_ligand_residues.end(); itr++){
			int index = *itr;
			cout << ligand->c->aminoacid[index]->chain << "|" << ligand->c->aminoacid[index]->cindex << "|" << (*ligand_vs_reference)[index] << ":" << ligand->c->aminoacid[index]->get_symbol() << ":" << reference->aminoacid[(*ligand_vs_reference)[index]]->get_symbol() 
						<< " " << ligand->c->aminoacid[index]->sequence_stability << ":" << reference->aminoacid[(*ligand_vs_reference)[index]]->sequence_stability << endl;
			eInterface += ligand->c->aminoacid[index]->eInterface;
			interface_seq_stability += ligand->c->aminoacid[index]->sequence_stability;
			interface_seq_stability_docked += reference->aminoacid[(*ligand_vs_reference)[index]]->sequence_stability;
		}
		cout << endl;

		for(hash_set<long,hash<long>,eqlong>::iterator itr = reference_interface_contacts.begin(); itr != reference_interface_contacts.end(); itr++){
			long index = *itr;
			int refrindex = index/MAX_ATOMS;
			int reflindex = index % MAX_ATOMS;
			if(reference_vs_receptor->count(refrindex) > 0){
				int index = (*reference_vs_receptor)[refrindex];
				rinterface[index] = true;
			}
			if(reference_vs_ligand->count(reflindex) > 0){
				int index = (*reference_vs_ligand)[reflindex];
				linterface[index] = true;
			}
		}
		cout << "Receptor:\t";
		for(int j = 0; j < c->num_aminoacids; j++)	cout << c->aminoacid[j]->get_symbol();
		cout << endl << "Receptor:\t";
		for(int j = 0; j < c->num_aminoacids; j++)	
			if(rinterface[j])	
				if(c->aminoacid[j]->rsa >=0.1)	cout << "|";
				else	cout << "B";	
			else cout << "-";
		cout << endl;
		cout << "Ligand:\t";
		for(int j = 0; j < ligand->c->num_aminoacids; j++)	cout << ligand->c->aminoacid[j]->get_symbol();
		cout << endl << "Ligand:\t";
		for(int j = 0; j < ligand->c->num_aminoacids; j++)	
			if(linterface[j])	
				if(ligand->c->aminoacid[j]->rsa >=0.1)	cout << "|";
				else	cout << "B";	
			else cout << "-";
		cout << endl;

		cout << "sequence stability " << interface_seq_stability << " change " << interface_seq_stability_docked-interface_seq_stability << endl;
		cout << "contactsmissing " << contacts_numgaps << " interface#gaps rec " <<  ref_vs_rec_interface_missing.size() << " lig " <<  ref_vs_lig_interface_missing.size() << endl;
		cout << reference->pdbcode << " reference interface " << reference_interface_contacts.size() << "\t" << reference_contact_energy_residue << "\t" << reference_contact_energy_coarse << endl;

		// compute the change of interface
		hash_set<int,hash<int>,eqint> added_ligand_residue, added_receptor_residue;
		Vector ref_rec_interface_points_bb[reference->num_atoms + 1],ref_lig_interface_points_bb[reference->num_atoms + 1],
		rec_interface_points_bb[reference->num_atoms + 1],lig_interface_points_bb[reference->num_atoms + 1];
		Vector ref_interface_points_ca[reference->num_atoms + 1],trunbound_interface_points_ca[reference->num_atoms + 1],
		rec_interface_points_ca[reference->num_atoms + 1],lig_interface_points_ca[reference->num_atoms + 1];
		Vector ref_rec_interface_points_scheavy[reference->num_atoms + 1],ref_lig_interface_points_scheavy[reference->num_atoms + 1],
		rec_interface_points_scheavy[reference->num_atoms + 1],lig_interface_points_scheavy[reference->num_atoms + 1];
		int bb_point_index = 0, scheavy_point_index = 0, ca_point_index = 0;

		for(hash_set<long,hash<long>,eqlong>::iterator itr = reference_interface_contacts_capri.begin(); itr != reference_interface_contacts_capri.end(); itr++){
			long index = *itr;
			int refrindex = index/MAX_ATOMS;

			if(reference_vs_receptor->count(refrindex) > 0)	
				if(added_receptor_residue.count(refrindex) == 0){
					added_receptor_residue.insert(refrindex);
					int rindex = (*reference_vs_receptor)[refrindex];
					Aminoacid *aa_ref = reference->aminoacid[refrindex];
					Aminoacid *aa = c->aminoacid[rindex];

					for(hash_map<const char*, Atom*, hash<const char*>, eqstr>::iterator raitr = aa->atom.begin(); raitr != aa->atom.end(); raitr++){
						Atom *ratom = (Atom *) raitr->second;
						if(ratom->name.c_str()[0] != 'H'){
							if(aa_ref->atom.count(ratom->name.c_str()) > 0){
								if(ratom->name == "N" || ratom->name == "CA" || ratom->name == "C" || ratom->name == "O"){
									ref_rec_interface_points_bb[bb_point_index] = Vector(aa_ref->atom[ratom->name.c_str()]->position);
									rec_interface_points_bb[bb_point_index] = Vector(ratom->position);
									bb_point_index++;
								} else {
									ref_rec_interface_points_scheavy[scheavy_point_index] = Vector(aa_ref->atom[ratom->name.c_str()]->position);
									rec_interface_points_scheavy[scheavy_point_index] = Vector(ratom->position);
									scheavy_point_index++;
								}
								if(ratom->name == "CA"){
									ref_interface_points_ca[ca_point_index] = Vector(aa_ref->atom[ratom->name.c_str()]->position);
									rec_interface_points_ca[ca_point_index] = Vector(ratom->position);
									trunbound_interface_points_ca[ca_point_index] = receptor_orig_to_ref->transform(*(ratom->position));
									ca_point_index++;
								}
							}
						}
					}
				}
		}
		int num_rec_interface_points_bb = bb_point_index,num_rec_interface_points_scheavy = scheavy_point_index, num_rec_interface_points_ca = ca_point_index;
		bb_point_index=0; scheavy_point_index=0, ca_point_index=0; 
		for(hash_set<long,hash<long>,eqlong>::iterator itr = reference_interface_contacts_capri.begin(); itr != reference_interface_contacts_capri.end(); itr++){
			long index = *itr;
			int reflindex = index % MAX_ATOMS;

			if(reference_vs_ligand->count(reflindex) > 0)	
				if(added_ligand_residue.count(reflindex) == 0){
					added_ligand_residue.insert(reflindex);
					int lindex = (*reference_vs_ligand)[reflindex];
					Aminoacid *aa_ref = reference->aminoacid[reflindex];
					Aminoacid *aa = ligand->c->aminoacid[lindex];

					for(hash_map<const char*, Atom*, hash<const char*>, eqstr>::iterator raitr = aa->atom.begin(); raitr != aa->atom.end(); raitr++){
						Atom *ratom = (Atom *) raitr->second;
						if(ratom->name.c_str()[0] != 'H'){
							if(aa_ref->atom.count(ratom->name.c_str()) > 0){
								if(ratom->name == "N" || ratom->name == "CA" || ratom->name == "C" || ratom->name == "O"){
									ref_lig_interface_points_bb[bb_point_index] = Vector(aa_ref->atom[ratom->name.c_str()]->position);
									lig_interface_points_bb[bb_point_index] = Vector(ratom->position);
									bb_point_index++;
								} else {
									ref_lig_interface_points_scheavy[scheavy_point_index] = Vector(aa_ref->atom[ratom->name.c_str()]->position);
									lig_interface_points_scheavy[scheavy_point_index] = Vector(ratom->position);
									scheavy_point_index++;
								}
								if(ratom->name == "CA"){
									ref_interface_points_ca[ca_point_index+num_rec_interface_points_ca] = Vector(aa_ref->atom[ratom->name.c_str()]->position);
									lig_interface_points_ca[ca_point_index] = Vector(ratom->position);
									trunbound_interface_points_ca[ca_point_index+num_rec_interface_points_ca] = ligand_orig_to_ref->transform(*(ratom->position));
									ca_point_index++;
								}
							}
						}
					}
				}
		}
		int num_lig_interface_points_bb = bb_point_index,num_lig_interface_points_scheavy = scheavy_point_index,num_lig_interface_points_ca = ca_point_index;

		// get the rmsd between the for the interface
		float irmsd_bb_rec = compute_rmsd(num_rec_interface_points_bb, &ref_rec_interface_points_bb[0], &rec_interface_points_bb[0]);
		float irmsd_scheavy_rec = compute_rmsd(num_rec_interface_points_scheavy, &ref_rec_interface_points_scheavy[0], &rec_interface_points_scheavy[0]);
		float irmsd_bb_lig = compute_rmsd(num_lig_interface_points_bb, &ref_lig_interface_points_bb[0], &lig_interface_points_bb[0]);
		float irmsd_scheavy_lig = compute_rmsd(num_lig_interface_points_scheavy, &ref_lig_interface_points_scheavy[0], &lig_interface_points_scheavy[0]);

		Transformation *rec_orig_to_ref_interfaceca = new Transformation(new Vector(0,0,0), new Vector(1,0,0), new Vector(0,1,0), 1.0, 0,0);
		Transformation *lig_orig_to_ref_interfaceca = new Transformation(new Vector(0,0,0), new Vector(1,0,0), new Vector(0,1,0), 1.0, 0,0);
		float irmsd_ca_rec = compute_rmsd(num_rec_interface_points_ca, &ref_interface_points_ca[0], &rec_interface_points_ca[0],rec_orig_to_ref_interfaceca);
		float irmsd_ca_lig = compute_rmsd(num_lig_interface_points_ca, &ref_interface_points_ca[num_rec_interface_points_ca], &lig_interface_points_ca[0],lig_orig_to_ref_interfaceca);

		//fill up the scheavy arrarys
		for(int i = 0; i < num_rec_interface_points_bb; i++){
			ref_rec_interface_points_scheavy[num_rec_interface_points_scheavy + i] = ref_rec_interface_points_bb[i];
			rec_interface_points_scheavy[num_rec_interface_points_scheavy + i] = rec_interface_points_bb[i];
		}
		for(int i = 0; i < num_lig_interface_points_bb; i++){
			ref_lig_interface_points_scheavy[num_lig_interface_points_scheavy + i] = ref_lig_interface_points_bb[i];
			lig_interface_points_scheavy[num_lig_interface_points_scheavy + i] = lig_interface_points_bb[i];
		}
		float irmsd_allheavy_rec = compute_rmsd(num_rec_interface_points_scheavy + num_rec_interface_points_bb, &ref_rec_interface_points_scheavy[0], &rec_interface_points_scheavy[0]);
		float irmsd_allheavy_lig = compute_rmsd(num_lig_interface_points_scheavy + num_lig_interface_points_bb, &ref_lig_interface_points_scheavy[0], &lig_interface_points_scheavy[0]);

		combined_rmsd = sqrt((irmsd_bb_rec*irmsd_bb_rec*num_rec_interface_points_bb+irmsd_bb_lig*irmsd_bb_lig*num_lig_interface_points_bb)/(num_rec_interface_points_bb+num_lig_interface_points_bb));
		cout << "irmsd bb rec " << irmsd_bb_rec << " lig " << irmsd_bb_lig << " both " << combined_rmsd << endl;
		combined_rmsd = sqrt((irmsd_ca_rec*irmsd_ca_rec*num_rec_interface_points_ca+irmsd_ca_lig*irmsd_ca_lig*num_lig_interface_points_ca)/(num_rec_interface_points_ca+num_lig_interface_points_ca));
		cout << "irmsd ca rec " << irmsd_ca_rec << " lig " << irmsd_ca_lig << " both " << combined_rmsd << endl;
		combined_rmsd = sqrt((irmsd_scheavy_rec*irmsd_scheavy_rec*num_rec_interface_points_scheavy+irmsd_scheavy_lig*irmsd_scheavy_lig*num_lig_interface_points_scheavy)/(num_rec_interface_points_scheavy+num_lig_interface_points_scheavy));
		cout << "irmsd scheavy rec " << irmsd_scheavy_rec << " lig " << irmsd_scheavy_lig << " both " << combined_rmsd << endl;
		combined_rmsd = sqrt((irmsd_allheavy_rec*irmsd_allheavy_rec*(num_rec_interface_points_scheavy + num_rec_interface_points_bb)
				+irmsd_allheavy_lig*irmsd_allheavy_lig*(num_lig_interface_points_scheavy + num_lig_interface_points_bb))/(num_rec_interface_points_scheavy + num_rec_interface_points_bb+num_lig_interface_points_scheavy + num_lig_interface_points_bb));
		cout << "irmsd allheavy rec " << irmsd_allheavy_rec << " lig " << irmsd_allheavy_lig << " both " << combined_rmsd << endl;

		// output the fraction of native contacts and nonnative contacts
		// overlay receptor and ligand based on the interface calpha and compute the contacts
		hash_set<long,hash<long>,eqlong> unbound_interface_contacts;
		for(int i = 0; i < c->num_aminoacids; i++){
			Aminoacid *ra = c->aminoacid[i];
			for(int j = 0; j < ligand->c->num_aminoacids; j++){
				Aminoacid *la = ligand->c->aminoacid[j];
				long cindex = (ra->cindex)*MAX_ATOMS + la->cindex;

				// defining the interface used for computing irmsd
				for(hash_map<const char*, Atom*, hash<const char*>, eqstr>::iterator raitr = ra->atom.begin(); raitr != ra->atom.end(); raitr++){
					Atom *ratom = (Atom *) raitr->second;
					if(ratom->name.c_str()[0] != 'H'){
						for(hash_map<const char*, Atom*, hash<const char*>, eqstr>::iterator laitr = la->atom.begin(); laitr != la->atom.end(); laitr++){
							Atom *latom = (Atom *) laitr->second;
							if(latom->name.c_str()[0] != 'H'){
								Vector vr = receptor_orig_to_ref->transform(*(ratom->position));
								Vector vl = ligand_orig_to_ref->transform(*(latom->position));
								if(Vector::distance(vr, vl) <= 10.0){
									unbound_interface_contacts.insert(cindex);
									//cout << ra->cindex << "|" << la->cindex << " " << cindex << endl;
								}
							}
						}
					}
				}
			}
		}

		unsigned short unbound_num_correct_contacts = 0;
		for(hash_set<long,hash<long>,eqlong>::iterator itr = reference_interface_contacts_capri.begin(); itr != reference_interface_contacts_capri.end(); itr++){
			long index = *itr;
			int reflindex = index % MAX_ATOMS;
			int refrindex = index/MAX_ATOMS;

			if(reference_vs_receptor->count(refrindex) > 0 && reference_vs_ligand->count(reflindex) > 0){
				int raindex = (*reference_vs_receptor)[refrindex];
				int laindex = (*reference_vs_ligand)[reflindex];
				long uindex = raindex*MAX_ATOMS + laindex;
				if(unbound_interface_contacts.count(uindex) > 0)
					unbound_num_correct_contacts++;
			}
		}
		float trunbound_irmsd_ca = compute_rmsd(num_rec_interface_points_ca+num_lig_interface_points_ca, &ref_interface_points_ca[0], &trunbound_interface_points_ca[0]);
		cout << unbound_num_correct_contacts << " " << reference_interface_contacts_capri.size() << " " << unbound_interface_contacts.size() << endl;
		cout << trunbound_irmsd_ca << " ";
		cout << "frac nat " << ((float)unbound_num_correct_contacts)/reference_interface_contacts_capri.size() << " frac nonnat " << (unbound_interface_contacts.size()-unbound_num_correct_contacts)/((float)(unbound_interface_contacts.size())) << endl; 

		/*float eElectrostatic = 0;
		for(int i = 0; i < rchains->length(); i++){
			Molecule *rchain = referenceH->molecules[rchains->at(i)];
			//cout << rchain->chain << "\t" << rchain->atom.size() << endl;
			for(int j = 0; j < lchains->length(); j++){
				Molecule *lchain = referenceH->molecules[lchains->at(j)];
				//cout << lchain->chain << "\t" << lchain->atom.size() << endl;
	 			for(hash_map<unsigned short, Atom*, hash<unsigned short>,eqint>::iterator ritr = rchain->atom.begin(); ritr != rchain->atom.end(); ritr++){
	 				Atom *ar = ritr->second;
					if(ar->charge != 0){
						for(hash_map<unsigned short, Atom*, hash<unsigned short>,eqint>::iterator litr = lchain->atom.begin(); litr != lchain->atom.end(); litr++){
			 				Atom *al = litr->second;		
							if(al->charge != 0){
								float d = Vector::distance(*(ar->position), ((Reference_Frame*) tr)->inverse_transform(*(al->position)));
								if(d < 1.0*(ar->radius + al->radius)){
									//cout << ar->cindex << " " << al->cindex << " " << d << "\t";
									d = 1.0*(ar->radius + al->radius);
								}	

								eElectrostatic += (ar->charge * al->charge)/d;
								//cout << ar->cindex << " " << al->cindex << " " << eElectrostatic << endl;
							}
						}
					}
				}
			}
		}
		cout << eElectrostatic << "\t";*/
		//cout << eInterface << "\t";
		cout << endl;

		tr->vmetrics = new VerificationMetrics();
		// geometric hashing based docking of just the interface regions
		/*{
			// compute votes
			int num_votes = 0;
			float d[ligand->num_faces][num_faces];
			float max_d;
			vector<long> matches;
			for(int i = 0 ; i < ligand->num_faces; i++){
				Vector vl = tr->inverse_transform(*(ligand->face[i]->point));
				//Vector vl = ligand_orig_to_ref->transform(*(ligand->face[i]->point));
				for(int j = 0; j < num_faces; j++){
					Vector vr = *(face[j]->point);
					//Vector vr = receptor_orig_to_ref->transform(*(face[j]->point));
					float dis = d[i][j] = Vector::distance(vl,vr);
					if(i == 0 && j == 0){
						max_d = dis;
					} else {
						if(dis > max_d)	max_d = dis;
					}
					if(dis < INCREMENT){
						matches.push_back(i*MAX_ATOMS + j);
						num_votes++;
					}
				}
			}
			tr->votes = num_votes;
			tr->vmetrics->rmsd = tr->vmetrics->lrmsd = tr->vmetrics->irmsd = tr->vmetrics->delta_r = tr->vmetrics->delta_U = 0;
			tr->print_details(&cout,TN_BASIC);

			cout << "#rpoints " << num_faces << " #lpoints " << ligand->num_faces << " grid_spacing " << INCREMENT << " #votes " << num_votes << endl;	cout.flush();
			//cout << "max d " << max_d << endl;
			cout << "vote details" << endl;
			unsigned short aligned_normals=0;
			for(vector<long>::iterator itr = matches.begin();itr != matches.end();itr++){
				long l = (long) *itr;
				int pl = l/MAX_ATOMS;
				int pr = l%MAX_ATOMS;
				Vector nl = Vector(ligand->face[pl]->normal->dot(tr->ex),ligand->face[pl]->normal->dot(tr->ey),ligand->face[pl]->normal->dot(tr->ez));
				Vector cnl = Vector(ligand->face[pl]->coarser_normal->dot(tr->ex),ligand->face[pl]->coarser_normal->dot(tr->ey),ligand->face[pl]->coarser_normal->dot(tr->ez));	
				cout << pl << " " << pr << " " << d[pl][pr] << " " << 0 - nl.dot(face[pr]->normal) << " " << 0 - cnl.dot(face[pr]->coarser_normal) << endl;
				if( nl.dot(face[pr]->normal) <= -0.7 || cnl.dot(face[pr]->coarser_normal) <= -0.7)
					aligned_normals++;
			}

			float box_size = INCREMENT;
			int d_divisions = (int) (max_d/box_size) + 1;
			int d_distribution[d_divisions];
			for(int i = 0 ; i < d_divisions; i++)	d_distribution[i] = 0;

			for(int i = 0 ; i < ligand->num_faces; i++)
				for(int j = 0; j < num_faces; j++)
					d_distribution[(int) (d[i][j]/box_size)]++;

			cout << "points distance distribution: " << "\t";
			for(int i = 0 ; i*INCREMENT < 5.0; i++)	cout << d_distribution[i] << " ";
			cout << endl;

			for(int i = 0 ; i < d_divisions; i++)	d_distribution[i] = 0;
			int votepindex[matches.size() + 1], count=0;
			max_d=0;
			for(vector<long>::iterator itr = matches.begin();itr != matches.end();itr++){
				long l = (long) *itr;
				int pr = l%MAX_ATOMS;
				votepindex[count++] = pr;
			}
			for(int i = 0 ; i < count; i++)
				for(int j = i+1; j < count; j++){
					float d = Vector::distance(*(face[votepindex[i]]->point), *(face[votepindex[j]]->point));
					if(d > max_d)	max_d = d;
					d_distribution[(int) (d/box_size)]++;
				}

			cout << "vote spread " << "\t";
			for(int i = 0 ; i*INCREMENT < max_d; i++)	cout << d_distribution[i] << " ";
			cout << endl;

			cout << "signal " << num_votes << " " << aligned_normals << endl;
		}*/
	}

	Transformation **result = (Transformation **) malloc(sizeof(Transformation*)*3);
	result[0] = receptor_orig_to_ref;
	result[1] = ligand_orig_to_ref;
	result[2] = tr;
	return result;
}
/*
 * Compute the lrmsd between the orignal, model and reference
 *  the number of contacts correctly identified
 *  the rmsd between contacts identified between original, model and reference
 */
void Object::verify_combinations(Object* ligand, Complex* reference, hash_map<long, long, hash<long>, eqlong> *receptor_vs_reference,
		hash_map<long, long, hash<long>, eqlong> *reference_vs_receptor, hash_map<long, long, hash<long>, eqlong> *ligand_vs_reference,
		hash_map<long, long, hash<long>, eqlong> *reference_vs_ligand, vector<MultiTransformation*> *combinations){
	int max_num_atoms = c->num_atoms + ligand->c->num_atoms;

	//Vector relative_r_reference = *(reference->molecules[ligand_chain]->center_of_mass) - *(reference->molecules[receptor_chain]->center_of_mass);

	Vector reference_points[max_num_atoms+1], model_points[max_num_atoms+1];
	int point_index = 0;
	for(int i = 0; i < c->num_aminoacids; i++){
		Aminoacid *aa = c->aminoacid[i];
		if(receptor_vs_reference->count(aa->cindex) > 0){
			int ref_residue_index = (*receptor_vs_reference)[aa->cindex];
			Aminoacid *aa_ref = reference->aminoacid[ref_residue_index];	
			reference_points[point_index] = Vector(aa_ref->alpha_carbon->position);
			point_index++;
		}
	}
	int ligand_start = point_index;

	for(int i = 0; i < ligand->c->num_aminoacids; i++){
		Aminoacid *aa = ligand->c->aminoacid[i];
		if(ligand_vs_reference->count(aa->cindex) > 0){
			int ref_residue_index = (*ligand_vs_reference)[aa->cindex];
			Aminoacid *aa_ref = reference->aminoacid[ref_residue_index];
			reference_points[point_index] = Vector(aa_ref->alpha_carbon->position);
			point_index++;
		}
	}

	hash_set<int,hash<int>,eqint> reference_interface_ligand_residues,reference_interface_receptor_residues;
	hash_set<long,hash<long>,eqlong> reference_interface_contacts;
	for(int i = 0; i < reference->num_aminoacids; i++){
		Aminoacid *ra = reference->aminoacid[i];
		if(reference_vs_receptor->count(ra->cindex) > 0){
			int raindex = (*reference_vs_receptor)[ra->cindex];
			for(int j = 0; j < reference->num_aminoacids; j++){
				Aminoacid *la = reference->aminoacid[j];
				if(reference_vs_ligand->count(la->cindex) > 0){
					int laindex = (*reference_vs_ligand)[la->cindex];
					if(Vector::distance(ra->alpha_carbon->position, la->alpha_carbon->position) < SS_CUTOFF){
						reference_interface_ligand_residues.insert(laindex);
						reference_interface_receptor_residues.insert(raindex);
						long cindex = (ra->cindex)*MAX_ATOMS + la->cindex;
						if(reference_interface_contacts.count(cindex) == 0){
							reference_interface_contacts.insert(cindex);
						}
					}
				}
			}
		}
	}

	Vector reference_interface_points[max_num_atoms+1];
	hash_set<int,hash<int>,eqint> added_ligand_residue, added_receptor_residue;
	point_index = 0;
	for(hash_set<long,hash<long>,eqlong>::iterator itr = reference_interface_contacts.begin(); itr != reference_interface_contacts.end(); itr++){
		//*out << *itr << "-";
		long index = *itr;
		int rindex = index/MAX_ATOMS;
		int lindex = index % MAX_ATOMS;

		if(added_receptor_residue.count(rindex) == 0){
			added_receptor_residue.insert(rindex);
			int ref_residue_index = (*receptor_vs_reference)[rindex];
			Aminoacid *aa_ref = reference->aminoacid[ref_residue_index];
			reference_interface_points[point_index] = Vector(aa_ref->alpha_carbon->position);
			point_index++;
		}
		if(added_ligand_residue.count(lindex) == 0){
			added_ligand_residue.insert(lindex);
			int ref_residue_index = (*ligand_vs_reference)[lindex];
			Aminoacid *aa_ref = reference->aminoacid[ref_residue_index];
			reference_interface_points[point_index] = Vector(aa_ref->alpha_carbon->position);
			point_index++;
		}
	}

	for(vector<MultiTransformation*>::iterator mitr = combinations->begin(); mitr != combinations->end(); mitr++){
		MultiTransformation *mtr = *mitr;

		stringstream ss (stringstream::in | stringstream::out);
		stringstream ss2 (stringstream::in | stringstream::out);
		ss << string((mtr->id).c_str());
		for(int i = 0; i < (mtr->transformations).size(); i++){
			if(i > 0)
				ss2 << ".";
			long tid;
			ss >> tid;
			ss2 << tid;
		}
		string mid;
		ss2 >> mid;
		mid = string(mid.c_str());
		ss2.clear();
		ss2 << "models/" << mid;
		string id;
		ss2 >> id;
		Complex *model = new Complex(id,"AB", PDB);

		// check the relative displacement
		/*model->molecules['A']->compute_motions();
		model->molecules['B']->compute_motions();
		Vector relative_r_model = *(model->molecules['B']->center_of_mass) - *(model->molecules['A']->center_of_mass);
		mtr->delta_r = new Vector(relative_r_model - relative_r_reference);*/

		point_index = 0;
		for(int i = 0; i < c->num_aminoacids; i++){
			Aminoacid *aa = c->aminoacid[i];
			if(receptor_vs_reference->count(aa->cindex) > 0){
				Aminoacid *aa_model = model->aminoacid[aa->cindex];
				model_points[point_index] = Vector(aa_model->alpha_carbon->position);
				point_index++;
			}
		}

		for(int i = 0; i < ligand->c->num_aminoacids; i++){
			Aminoacid *aa = ligand->c->aminoacid[i];
			if(ligand_vs_reference->count(aa->cindex) > 0){
				Aminoacid *aa_model = model->aminoacid[c->num_aminoacids + aa->cindex];
				// output receptor chains before ligand chains in the model file
				model_points[point_index] = Vector(aa_model->alpha_carbon->position);
				point_index++;
			}
		}
		*out << "lrmsd # " << point_index << endl;
		//out->flush();
		mtr->lrmsd = compute_lrmsd(point_index, ligand_start, &reference_points[0],&model_points[0], &(mtr->rmsd));

		// how good are your contacts
		float model_contact_energy=0;
		hash_set<int,hash<int>,eqint> model_interface_ligand_residues,model_interface_receptor_residues;
		hash_set<long,hash<long>,eqlong> model_interface_contacts;
		for(int i = 0; i < c->num_aminoacids; i++){
			Aminoacid *ra = model->aminoacid[i];
			for(int j = 0; j < ligand->c->num_aminoacids; j++){
				Aminoacid *la = model->aminoacid[c->num_aminoacids + j];
				//cout << Vector::distance(ra->alpha_carbon->position, la->alpha_carbon->position) << "-";
				if(Vector::distance(ra->alpha_carbon->position, la->alpha_carbon->position) < SS_CUTOFF){
					model_interface_ligand_residues.insert(la->cindex);
					model_interface_receptor_residues.insert(ra->cindex);
					long cindex = (ra->cindex)*MAX_ATOMS + la->cindex;
					if(model_interface_contacts.count(cindex) == 0){
						model_interface_contacts.insert(cindex);
						model_contact_energy += residue_potential[ra->type][la->type];
					}
				}
			}
		}

		int num_native_contacts_predicted=0;
		Vector model_interface_points[max_num_atoms+1];
		hash_set<int,hash<int>,eqint> added_ligand_residue, added_receptor_residue;
		point_index = 0;
		for(hash_set<long,hash<long>,eqlong>::iterator itr = reference_interface_contacts.begin(); itr != reference_interface_contacts.end(); itr++){
			//*out << *itr << "-";
			long index = *itr;
			int rindex = index/MAX_ATOMS;
			int lindex = index % MAX_ATOMS;

			if(added_receptor_residue.count(rindex) == 0){
				added_receptor_residue.insert(rindex);
				Aminoacid *aa_model = model->aminoacid[rindex];
				model_interface_points[point_index] = Vector(aa_model->alpha_carbon->position);
				point_index++;
			}
			if(added_ligand_residue.count(lindex) == 0){
				added_ligand_residue.insert(lindex);
				Aminoacid *aa_model = model->aminoacid[lindex];
				model_interface_points[point_index] = Vector(aa_model->alpha_carbon->position);
				point_index++;
			}

			if(model_interface_contacts.count(index) > 0){
				num_native_contacts_predicted++;
			}
		}
		*out << "irmsd #" << point_index << endl;
		mtr->irmsd = compute_rmsd(point_index, &reference_interface_points[0], &model_interface_points[0]);
		mtr->frac_native_contacts_predicted = num_native_contacts_predicted;
		mtr->num_contacts = model_interface_contacts.size();
		mtr->eResiduepair = model_contact_energy;	
	}
}

void refine_sidechains(Object* receptor, Object *ligand, vector<Transformation*> matching_frames){
	*out << "refining side chains (transform ligand)\t scwrl time limit " << SCWRL_TIME_LIMIT << endl;

	//set a timelimit
	char command[512];
	sprintf(command, "ulimit -t %d",SCWRL_TIME_LIMIT);
	int iret = system(command);

	for(vector<Transformation*>::iterator itr = matching_frames.begin(); itr != matching_frames.end(); itr++){
		Transformation *tr = *itr;
		optimize_sidechains(receptor->c,ligand->c,tr);
	}
}

void verify(Object* receptor, Object *ligand, Complex* reference, vector<Transformation*> matching_frames, int atoms_included){
	//*out << "Verifying ..." << endl;

	int max_num_atoms = reference->num_atoms;
	Vector reference_points[max_num_atoms+1], receptor_points_refined[max_num_atoms+1], ligand_points_refined[max_num_atoms+1];
	Vector receptor_points[max_num_atoms+1], ligand_points[max_num_atoms+1];
	int atom_index, receptor_atom_index, ligand_atom_index;

	*out << "rmsd_r lrmsd_r ";
	*out << "rmsd lrmsd frame_no residue_contacts dock_score atoms_interior atoms_intermediate atoms_surface points_interior points_intermediate " 
			<< "points_surface cscore "
			<< "negneg negneu negpos neuneg neuneu neupos posneg posneu pospos "
			<< "cp_contacts votes " << endl;

	for(vector<Transformation*>::iterator itr = matching_frames.begin(); itr != matching_frames.end(); itr++){
		Transformation *tr = *itr;
		atom_index = 0; receptor_atom_index = 0; ligand_atom_index = 0;
		/*char scwrl_out_file[128];
		fstream scwrlout;
		sprintf(scwrl_out_file,"%s%ld",SCWRL_OUTPUT_FILE,tr->frame_number); 
		scwrlout.open(scwrl_out_file, fstream::in);
		//*out << "check " << tr->frame_number << " " << scwrl_out_file << endl;
		if(!scwrlout.is_open()){
		 *out << "Error opening file " << scwrl_out_file << " " << errno << "\n";
		} else {
			hash_map<char,char,hash<char>,eqint> rm_vs_m, lm_vs_m;
			char chain = 'A';
			for(hash_map<char, Molecule*, hash<char>, eqint>::iterator mitr = c->molecules.begin(); mitr != molecules.end(); mitr++){
				Molecule *m = mitr->second;

			}
			char buf[8192];
			do{
				scwrlout.getline(buf,8192);
				//if(scwrlout.gcount() == 0)
				//	scwrlout.getline(buf,8192,'\r');
				//*out << "a? " << buf << " " << scwrlout.gcount() << endl;
			}while((string(buf)).find("ATOM") == string::npos);

			while (!scwrlout.eof()){
				//*out << buf << endl;
				stringstream* line = new stringstream(buf,stringstream::in);
				string format;
		 *line >> format;
			    int aindex;
		 *line >> aindex;
			    string atom_name;
		 *line >> atom_name;
			    string aminoacid_name;
		 *line >> aminoacid_name;
			    char chain;
		 *line >> chain;
				int aaindex;
		 *line >> aaindex;

				float x,y,z;
		 *line >> x;
		 *line >> y;
		 *line >> z;

				Aminoacid *aa;
				// assuming atoms of molecule A appear before atoms of molecule B
				if((chain - 'A') < c->molecules.size()){
					aa = (o->m->aminoacid)[aaindex];
					Atom *a = aa->atoms[atom_name.c_str()];
					if((atoms_included = CALPHA && a->name == ALPHA_CARBON) || atoms_included == ALLATOMS){
						reference_points[atom_index] = Vector(a->position);
						receptor_points[atom_index] = Vector(a->position);
						receptor_points_refined[atom_index] = Vector(x,y,z);
						atom_index++;
						receptor_atom_index++;
					}
				} else {
					aa = (m->aminoacid)[aaindex];
					Atom *a = aa->atoms[atom_name.c_str()];
					if((atoms_included = CALPHA && a->name == ALPHA_CARBON) || atoms_included == ALLATOMS){
						reference_points[atom_index] = Vector(a->position);
						ligand_points[ligand_atom_index] = Vector(a->position);
						ligand_points_refined[ligand_atom_index] = Vector(x,y,z);
						atom_index++;
						ligand_atom_index++;
					}
				}

				if(!scwrlout.eof())
					scwrlout.getline(buf,8192);
			}
			//*out << "eof " << scwrlout.eof() << " open " << scwrlout.is_open() << endl;
			scwrlout.close();
			//*out << atom_index << " " << receptor_atom_index << " " << ligand_atom_index << endl;

			compute_rmsd(atom_index, receptor_atom_index, &reference_points[0], 
				&receptor_points_refined[0], &ligand_points_refined[0], tr);
			compute_rmsd(atom_index, receptor_atom_index, &reference_points[0],
				&receptor_points[0], &ligand_points[0], tr);
			tr->print_results(out);
		}*/
	}
}

float tr_irmsd_distance(Complex* receptor, Object* ligand_obj, Transformation* tr1, Transformation *tr2){
	Complex *ligand = ligand_obj->c;
	Transformation *tr[2];
	tr[0] = tr1;
	tr[1] = tr2;


	hash_set<int,hash<int>,eqint> interface_receptor_residues, interface_ligand_residues;
	for(int ai = 0 ; ai < receptor->num_aminoacids; ai++){
		Aminoacid *ar = receptor->aminoacid[ai];

		if(ar->centroid != NULL){
			bool contact[2];
			for(short ti = 0; ti < 2; ti++){
				contact[ti] = false;

				Vector vr = tr[ti]->transform(*(ar->centroid));
				Vector vv = (vr - *(ligand_obj->grid_origin)) - Vector(SS_CUTOFF/AA_GRID_SPACING,SS_CUTOFF/AA_GRID_SPACING,SS_CUTOFF/AA_GRID_SPACING);
				int x,y,z;
				x = (int) (vv.x/AA_GRID_SPACING);
				y = (int) (vv.y/AA_GRID_SPACING);
				z = (int) (vv.z/AA_GRID_SPACING);
				if(x >= 0 && x < ligand_obj->aagrid_num_xdivisions && y >= 0 && y < ligand_obj->aagrid_num_ydivisions && z >= 0 && z < ligand_obj->aagrid_num_zdivisions){
					unsigned int index = (x*ligand_obj->aagrid_num_ydivisions + y)*ligand_obj->aagrid_num_zdivisions + z;
					if(ligand_obj->aagrid[index] != NULL){
						unsigned short *aagcontact = ligand_obj->aagrid[index];
						for(int agi = 0; agi < ligand_obj->aagrid_size[index]; agi++){
							int laindex = aagcontact[agi];
							Aminoacid *la = ligand->aminoacid[laindex];
							if(la->centroid != NULL && Vector::distance_squared(vr,*(la->centroid)) < SS_CUTOFF*SS_CUTOFF){
								contact[ti]=true;
								interface_ligand_residues.insert(la->cindex);
							}
						}
					}
				}
			}

			if(contact[0] || contact[1]) {
				interface_receptor_residues.insert(ar->cindex);
			}
		}
	}
	// cout << "d call " << tr1->frame_number << " " << tr2->frame_number << endl; 


	if(interface_receptor_residues.size() == 0 || interface_ligand_residues.size() == 0)	return 100;
	else{
		short nr = interface_receptor_residues.size();
		short nl = interface_ligand_residues.size();
		Vector tr1_interface_points[nr+nl], tr2_interface_points[nr+nl];

		int npt = 0;
		for(hash_set<int,hash<int>,eqint>::iterator itr = interface_receptor_residues.begin(); itr != interface_receptor_residues.end(); itr++){
			int aaindex = *itr;
			Aminoacid *aa = receptor->aminoacid[aaindex];

			if(aa->centroid != NULL){
				tr1_interface_points[npt] = tr2_interface_points[npt] = Vector(aa->centroid);
				npt++;
			}
		}

		for(hash_set<int,hash<int>,eqint>::iterator itr = interface_ligand_residues.begin(); itr != interface_ligand_residues.end(); itr++){
			int aaindex = *itr;
			Aminoacid *aa = ligand->aminoacid[aaindex];

			cout << aa->centroid << endl;
			if(aa->centroid != NULL){
				tr1_interface_points[npt] = tr1->inverse_transform(Vector(aa->centroid)); 


				tr2_interface_points[npt] = tr2->inverse_transform(Vector(aa->centroid));;
				npt++;
			}
		}


		float irmsd = compute_rmsd(npt, &tr1_interface_points[0], &tr2_interface_points[0]);
		// *out << tr1->frame_number << " " << tr2->frame_number << " " << irmsd << endl;

		return irmsd;
	}
}

float tr_irmsd_distance(Complex* receptor, Complex* ligand, Transformation* tr1, Transformation *tr2){
	Transformation *tr[2];
	tr[0] = tr1;
	tr[1] = tr2;

	hash_set<int,hash<int>,eqint> interface_receptor_residues, interface_ligand_residues;
	for(int i = 0; i < receptor->num_aminoacids; i++){
		Aminoacid *ra = receptor->aminoacid[i];
		for(int j = 0; j < ligand->num_aminoacids; j++){
			Aminoacid *la = ligand->aminoacid[j];
			bool contact[2];
			for(short ti = 0; ti < 2; ti++){
				contact[ti] = false;
				if(ra->centroid != NULL && la->centroid != NULL){
					Vector v = Vector(*(la->centroid));
					v = tr[ti]->inverse_transform(v);
					if(Vector::distance(ra->centroid,v) < SS_CUTOFF)	contact[ti]=true;
				} else if(ra->alpha_carbon != NULL && la->alpha_carbon != NULL){
					Vector v = Vector(*(la->alpha_carbon->position));
					v = tr[ti]->inverse_transform(v);
					if(Vector::distance(ra->alpha_carbon->position,v) < SS_CUTOFF)	contact[ti]=true;
				}
			}

			if(contact[0] || contact[1]) {
				interface_ligand_residues.insert(la->cindex);
				interface_receptor_residues.insert(ra->cindex);
			}
		}
	}
	//cout << "d call " << tr1->frame_number << " " << tr2->frame_number << endl; 

	if(interface_receptor_residues.size() == 0 || interface_ligand_residues.size() == 0)	return 100;
	else{
		short nr = interface_receptor_residues.size();
		short nl = interface_ligand_residues.size();
		Vector tr1_interface_points[nr+nl], tr2_interface_points[nr+nl];

		int npt = 0;
		for(hash_set<int,hash<int>,eqint>::iterator itr = interface_receptor_residues.begin(); itr != interface_receptor_residues.end(); itr++){
			int aaindex = *itr;
			Aminoacid *aa = receptor->aminoacid[aaindex];

			if(aa->centroid != NULL){
				tr1_interface_points[npt] = tr2_interface_points[npt] = Vector(aa->centroid);
				npt++;
			}
		}

		for(hash_set<int,hash<int>,eqint>::iterator itr = interface_ligand_residues.begin(); itr != interface_ligand_residues.end(); itr++){
			int aaindex = *itr;
			Aminoacid *aa = ligand->aminoacid[aaindex];

			if(aa->centroid != NULL){
				tr1_interface_points[npt] = tr1->inverse_transform(Vector(aa->centroid)); 
				tr2_interface_points[npt] = tr2->inverse_transform(Vector(aa->centroid));;
				npt++;
			}
		}

		float irmsd = compute_rmsd(npt, &tr1_interface_points[0], &tr2_interface_points[0]);
		//*out << tr1->frame_number << " " << tr2->frame_number << " " << irmsd << endl;

		return irmsd;
	}
}

float tr_frac_shared_contacts(Complex* receptor, Complex* ligand, Transformation* tr1, Transformation *tr2){
	Transformation *tr[2];
	tr[0] = tr1;
	tr[1] = tr2;

	int num_contacts1=0, num_contacts2=0, num_contacts_both=0;
	for(int i = 0; i < receptor->num_aminoacids; i++){
		Aminoacid *ra = receptor->aminoacid[i];
		for(int j = 0; j < ligand->num_aminoacids; j++){
			Aminoacid *la = ligand->aminoacid[j];
			bool contact[2];
			for(short ti = 0; ti < 2; ti++){
				contact[ti] = false;
				if(ra->centroid != NULL && la->centroid != NULL){
					Vector v = Vector(*(la->centroid));
					v = tr[ti]->inverse_transform(v);
					if(Vector::distance(ra->centroid,v) < SS_CUTOFF)	contact[ti]=true;
				} else if(ra->alpha_carbon != NULL && la->alpha_carbon != NULL){
					Vector v = Vector(*(la->alpha_carbon->position));
					v = tr[ti]->inverse_transform(v);
					if(Vector::distance(ra->alpha_carbon->position,v) < SS_CUTOFF)	contact[ti]=true;
				}

				if(contact[ti]){
					if(ti == 0)
						num_contacts1++;
					else
						num_contacts2++;
				}

				if(contact[0] && contact[1])
					num_contacts_both++;
			}
		}
	}

	if(num_contacts1 == 0 || num_contacts2 == 0)
		return 0;
	else {
		int max_contacts = maximum(num_contacts1,num_contacts2);
		float frac = ((float) num_contacts_both)/max_contacts;
		return frac;
	}
}

void Object::verify_transformations(Object* ligand, Complex* reference, string *refrchains, string *reflchains, hash_map<long, long, hash<long>, eqlong> *receptor_vs_reference,
		hash_map<long, long, hash<long>, eqlong> *reference_vs_receptor, hash_map<long, long, hash<long>, eqlong> *ligand_vs_reference,
		hash_map<long, long, hash<long>, eqlong> *reference_vs_ligand, vector<Transformation*> *ligand_vs_receptor_symmetry, vector<Transformation*> *receptor_symmetry, 
		vector<Transformation*> *ligand_symmetry, int atoms_included, vector<Transformation*> *matching_frames, Transformation** trs){
	int max_num_atoms = c->num_atoms + ligand->c->num_atoms;

	unsigned short num_receptor_symmetries = 1 + receptor_symmetry->size();
	unsigned short num_ligand_symmetries = 1 + ligand_symmetry->size();
	unsigned short num_ligand_vs_receptor_symmetries = 1 + ligand_vs_receptor_symmetry->size();
	unsigned int num_symmetries = num_receptor_symmetries*num_ligand_symmetries*num_ligand_vs_receptor_symmetries;
	bool interface_capridefn=true, use_all_bbatom_for_rmsd=false;

	Vector reference_points[num_symmetries][max_num_atoms+1], model_points[max_num_atoms+1], ligand_points[max_num_atoms+1];
	int point_index = 0;
	for(int i = 0; i < c->num_aminoacids; i++){
		Aminoacid *aa = c->aminoacid[i];
		if(receptor_vs_reference->count(aa->cindex) > 0){
			int ref_residue_index = (*receptor_vs_reference)[aa->cindex];
			Aminoacid *aa_ref = reference->aminoacid[ref_residue_index];

			if(aa->atom.count("N") > 0 && aa_ref->atom.count("N") > 0 && use_all_bbatom_for_rmsd){
				reference_points[0][point_index] = Vector(aa_ref->atom["N"]->position);
				model_points[point_index] = Vector(aa->atom["N"]->position);
				point_index++;
			} 
			if(aa->atom.count("CA") > 0 && aa_ref->atom.count("CA") > 0){
				reference_points[0][point_index] = Vector(aa_ref->atom["CA"]->position);
				model_points[point_index] = Vector(aa->atom["CA"]->position);
				point_index++;
			} 
			if(aa->atom.count("C") > 0 && aa_ref->atom.count("C") > 0 && use_all_bbatom_for_rmsd){
				reference_points[0][point_index] = Vector(aa_ref->atom["C"]->position);
				model_points[point_index] = Vector(aa->atom["C"]->position);
				point_index++;
			} 
			if(aa->atom.count("O") > 0 && aa_ref->atom.count("O") > 0 && use_all_bbatom_for_rmsd){
				reference_points[0][point_index] = Vector(aa_ref->atom["O"]->position);
				model_points[point_index] = Vector(aa->atom["O"]->position);
				point_index++;
			}
		}
	}
	int ligand_start = point_index;

	for(int i = 0; i < ligand->c->num_aminoacids; i++){
		Aminoacid *aa = ligand->c->aminoacid[i];
		if(ligand_vs_reference->count(aa->cindex) > 0){
			int ref_residue_index = (*ligand_vs_reference)[aa->cindex];
			Aminoacid *aa_ref = reference->aminoacid[ref_residue_index];

			if(aa->atom.count("N") > 0 && aa_ref->atom.count("N") > 0 && use_all_bbatom_for_rmsd){
				reference_points[0][point_index] = Vector(aa_ref->atom["N"]->position);
				ligand_points[point_index - ligand_start] = Vector(aa->atom["N"]->position);
				point_index++;
			} 
			if(aa->atom.count("CA") > 0 && aa_ref->atom.count("CA") > 0){
				reference_points[0][point_index] = Vector(aa_ref->atom["CA"]->position);
				ligand_points[point_index - ligand_start] = Vector(aa->atom["CA"]->position);
				point_index++;
			} 
			if(aa->atom.count("C") > 0 && aa_ref->atom.count("C") > 0 && use_all_bbatom_for_rmsd){
				reference_points[0][point_index] = Vector(aa_ref->atom["C"]->position);
				ligand_points[point_index - ligand_start] = Vector(aa->atom["C"]->position);
				point_index++;
			} 
			if(aa->atom.count("O") > 0 && aa_ref->atom.count("O") > 0 && use_all_bbatom_for_rmsd){
				reference_points[0][point_index] = Vector(aa_ref->atom["O"]->position);
				ligand_points[point_index - ligand_start] = Vector(aa->atom["O"]->position);
				point_index++;
			}

			//cout << "reference " << i << " " << point_index << "\t" << v.x << " " << v.y << " " << v.z << endl;
		}
	}
	int num_ligand_points = point_index - ligand_start;

	int count=0;
	Transformation *ligand_vs_receptor_symmetries[num_ligand_vs_receptor_symmetries];
	for(vector<Transformation*>::iterator titr = ligand_vs_receptor_symmetry->begin(); titr != ligand_vs_receptor_symmetry->end(); titr++)
		ligand_vs_receptor_symmetries[count++] = ((Transformation *) *titr); 

	count=0;
	Transformation *receptor_symmetries[num_receptor_symmetries];
	for(vector<Transformation*>::iterator titr = receptor_symmetry->begin(); titr != receptor_symmetry->end(); titr++)
		receptor_symmetries[count++] = ((Transformation *) *titr);

	count=0;
	Transformation *ligand_symmetries[num_ligand_symmetries];
	for(vector<Transformation*>::iterator titr = ligand_symmetry->begin(); titr != ligand_symmetry->end(); titr++)
		ligand_symmetries[count++] = ((Transformation *) *titr);

	hash_set<long,hash<long>,eqlong> reference_interface_contacts[num_symmetries];
	hash_set<long,hash<long>,eqlong> reference_interface_contacts_capri[num_symmetries];
	Vector ref_interface_points[num_symmetries][max_num_atoms+1], reflig_interface_points[num_symmetries][max_num_atoms+1];;
	int num_refrec_interface_points[num_symmetries], num_reflig_interface_points[num_symmetries];

	// applying symmetry transformations on the reference
	for(unsigned int orientation=0; orientation < num_symmetries; orientation++){
		unsigned short homodimer_orientation = (orientation/(num_receptor_symmetries*num_ligand_symmetries));
		unsigned short receptor_orientation = (orientation/num_ligand_symmetries) % num_receptor_symmetries;
		unsigned short ligand_orientation = orientation % num_ligand_symmetries;

		if(orientation > 0)
			for(int i = 0; i < num_ligand_points+ligand_start; i++){
				reference_points[orientation][i] = reference_points[0][i]; 
			}

		for(int i = 0; i < num_ligand_points; i++){
			point_index = ligand_start + i;	

			// ligand symmetry
			if(ligand_orientation > 0)
				reference_points[orientation][point_index] = ligand_symmetries[ligand_orientation-1]->transform(reference_points[orientation][point_index]);

			// ligand vs receptor symmetry
			if(homodimer_orientation > 0){
				Transformation *htr = ligand_vs_receptor_symmetries[homodimer_orientation-1];
				reference_points[orientation][point_index] = htr->transform(reference_points[orientation][point_index]);
				reference_points[orientation][point_index] = htr->transform(reference_points[orientation][point_index]);
			} 

			// receptor symmetry
			if(receptor_orientation > 0)
				reference_points[orientation][point_index] = receptor_symmetries[receptor_orientation-1]->inverse_transform(reference_points[orientation][point_index]);

			if(orientation > 0)	
				*out << point_index << " - (" << reference_points[0][point_index].x << "," << reference_points[0][point_index].x << "," <<  reference_points[0][point_index].x << ")\t(" <<
				reference_points[1][point_index].x << "," << reference_points[1][point_index].x << "," <<  reference_points[1][point_index].x << ")" << endl;
		}

		for(int i = 0; i < reference->num_aminoacids; i++){
			Aminoacid *ra = reference->aminoacid[i];
			if(reference_vs_receptor->count(ra->cindex) > 0){
				for(int j = 0; j < reference->num_aminoacids; j++){
					Aminoacid *la = reference->aminoacid[j];
					if(reference_vs_ligand->count(la->cindex) > 0){
						long cindex = (ra->cindex)*MAX_ATOMS + la->cindex;

						Vector vca, vcd;
						if(la->centroid != NULL){
							vcd = Vector(la->centroid);
							if(ligand_orientation > 0)
								vcd = ligand_symmetries[ligand_orientation-1]->transform(vcd);
							if(homodimer_orientation > 0){
								vcd = ligand_vs_receptor_symmetries[homodimer_orientation-1]->transform(vcd);
								vcd = ligand_vs_receptor_symmetries[homodimer_orientation-1]->transform(vcd);
							}
							if(receptor_orientation > 0)
								vcd = receptor_symmetries[receptor_orientation-1]->inverse_transform(vcd);
						}
						if(la->alpha_carbon != NULL){
							vca = Vector(la->alpha_carbon->position);
							if(ligand_orientation > 0)
								vca = ligand_symmetries[ligand_orientation-1]->transform(vca);
							if(homodimer_orientation > 0){
								vca = ligand_vs_receptor_symmetries[homodimer_orientation-1]->transform(vca);
								vca = ligand_vs_receptor_symmetries[homodimer_orientation-1]->transform(vca);
							}
							if(receptor_orientation > 0)
								vca = receptor_symmetries[receptor_orientation-1]->inverse_transform(vca);
						}

						if((ra->centroid != NULL && la->centroid != NULL && (Vector::distance(ra->centroid,vcd) < SS_CUTOFF)) ||
								(ra->alpha_carbon != NULL && la->alpha_carbon != NULL && (Vector::distance(ra->alpha_carbon->position,vca) < SS_CUTOFF))){
							reference_interface_contacts[orientation].insert(cindex);
						}

						// defining the interface used for computing irmsd
						for(hash_map<const char*, Atom*, hash<const char*>, eqstr>::iterator raitr = ra->atom.begin(); raitr != ra->atom.end(); raitr++){
							Atom *ratom = (Atom *) raitr->second;
							if(ratom->name.c_str()[0] != 'H'){
								for(hash_map<const char*, Atom*, hash<const char*>, eqstr>::iterator laitr = la->atom.begin(); laitr != la->atom.end(); laitr++){
									Atom *latom = (Atom *) laitr->second;
									if(latom->name.c_str()[0] != 'H'){
										Vector v = Vector(*(latom->position));
										if(ligand_orientation > 0)
											v = ligand_symmetries[ligand_orientation-1]->transform(v);
										if(homodimer_orientation > 0){
											v = ligand_vs_receptor_symmetries[homodimer_orientation-1]->transform(v);
											v = ligand_vs_receptor_symmetries[homodimer_orientation-1]->transform(v);
										}
										if(receptor_orientation > 0)
											v = receptor_symmetries[receptor_orientation-1]->inverse_transform(v);

										if(Vector::distance(*(ratom->position), v) <= CAPRI_INTERFACE_CUTOFF){
											reference_interface_contacts_capri[orientation].insert(cindex);
										}
									}
								}
							}
						}
					}
				}
			}
		}

		hash_set<long,hash<long>,eqlong> *ref_interface_contacts;
		if(interface_capridefn) 	
			ref_interface_contacts = &(reference_interface_contacts_capri[orientation]);
		else
			ref_interface_contacts = &(reference_interface_contacts[orientation]);

		hash_set<int,hash<int>,eqint> added_ligand_residue, added_receptor_residue;
		num_refrec_interface_points[orientation] =  num_reflig_interface_points[orientation] = 0;
		for(hash_set<long,hash<long>,eqlong>::iterator itr = ref_interface_contacts->begin(); itr != ref_interface_contacts->end(); itr++){
			long index = *itr;
			int refrindex = index/MAX_ATOMS;
			int reflindex = index % MAX_ATOMS;

			if(added_receptor_residue.count(refrindex) == 0){
				added_receptor_residue.insert(refrindex);
				int rindex = (*reference_vs_receptor)[refrindex];
				*out << refrindex << ":" << rindex << endl; out->flush();
				Aminoacid *aa_ref = reference->aminoacid[refrindex];
				Aminoacid *aa = c->aminoacid[rindex];

				if(aa->atom.count("N") > 0 && aa_ref->atom.count("N") > 0 && use_all_bbatom_for_rmsd){
					ref_interface_points[orientation][num_refrec_interface_points[orientation]++] = Vector(aa_ref->atom["N"]->position);
				} 
				if(aa->atom.count("CA") > 0 && aa_ref->atom.count("CA") > 0){
					ref_interface_points[orientation][num_refrec_interface_points[orientation]++] = Vector(aa_ref->atom["CA"]->position);
				} 
				if(aa->atom.count("C") > 0 && aa_ref->atom.count("C") > 0 && use_all_bbatom_for_rmsd){
					ref_interface_points[orientation][num_refrec_interface_points[orientation]++] = Vector(aa_ref->atom["C"]->position);
				} 
				if(aa->atom.count("O") > 0 && aa_ref->atom.count("O") > 0 && use_all_bbatom_for_rmsd){
					ref_interface_points[orientation][num_refrec_interface_points[orientation]++] = Vector(aa_ref->atom["O"]->position);
				}
			}
			if(added_ligand_residue.count(reflindex) == 0){
				added_ligand_residue.insert(reflindex);
				int lindex = (*reference_vs_ligand)[reflindex];
				Aminoacid *aa_ref = reference->aminoacid[reflindex];
				Aminoacid *aa = ligand->c->aminoacid[lindex];

				if(aa->atom.count("N") > 0 && aa_ref->atom.count("N") > 0 && use_all_bbatom_for_rmsd){
					reflig_interface_points[orientation][num_reflig_interface_points[orientation]++] = Vector(aa_ref->atom["N"]->position);
				} 
				if(aa->atom.count("CA") > 0 && aa_ref->atom.count("CA") > 0){
					reflig_interface_points[orientation][num_reflig_interface_points[orientation]++] = Vector(aa_ref->atom["CA"]->position);
				} 
				if(aa->atom.count("C") > 0 && aa_ref->atom.count("C") > 0 && use_all_bbatom_for_rmsd){
					reflig_interface_points[orientation][num_reflig_interface_points[orientation]++] = Vector(aa_ref->atom["C"]->position);
				} 
				if(aa->atom.count("O") > 0 && aa_ref->atom.count("O") > 0 && use_all_bbatom_for_rmsd){
					reflig_interface_points[orientation][num_reflig_interface_points[orientation]++] = Vector(aa_ref->atom["O"]->position);
				} 
			}
		}

		// add ligand interface points at the end of ref_interface_points
		for(int i = 0;i < num_reflig_interface_points[orientation]; i++){
			if(orientation > 0){
				if(ligand_orientation > 0)
					reflig_interface_points[orientation][i] = ligand_symmetries[ligand_orientation-1]->transform(reflig_interface_points[orientation][i]);
				if(homodimer_orientation > 0){
					reflig_interface_points[orientation][i] = ligand_vs_receptor_symmetries[homodimer_orientation-1]->transform(reflig_interface_points[orientation][i]);
					reflig_interface_points[orientation][i] = ligand_vs_receptor_symmetries[homodimer_orientation-1]->transform(reflig_interface_points[orientation][i]);
				}
				if(receptor_orientation > 0)
					reflig_interface_points[orientation][i] = receptor_symmetries[receptor_orientation-1]->inverse_transform(reflig_interface_points[orientation][i]);
			}

			ref_interface_points[orientation][num_refrec_interface_points[orientation]+i] = reflig_interface_points[orientation][i];
		}
	}

	//Transformation *receptor_orig_to_ref = new Transformation(new Vector(0,0,0), new Vector(1,0,0), new Vector(0,1,0), 1.0, 0,0);
	//float rrmsd_ref_vs_orig = compute_rmsd(ligand_start, &reference_points[0][0], &model_points[0], receptor_orig_to_ref);

	// collect the points corresponding to the interface
	Vector receptor_interface_points[num_symmetries][max_num_atoms+1];
	Vector ligand_interface_points[num_symmetries][max_num_atoms+1];
	int num_receptor_interface_points[num_symmetries],num_ligand_interface_points[num_symmetries];
	for(unsigned int orientation=0; orientation < num_symmetries; orientation++){
		unsigned short homodimer_orientation = (orientation/(num_receptor_symmetries*num_ligand_symmetries));
		unsigned short receptor_orientation = (orientation/num_ligand_symmetries) % num_receptor_symmetries;
		unsigned short ligand_orientation = orientation % num_ligand_symmetries;
		//*out << orientation << " " << receptor_orientation << " " << ligand_orientation << " " << homodimer_orientation << endl; out->flush();

		hash_set<long,hash<long>,eqlong> *ref_interface_contacts;
		if(interface_capridefn) 	
			ref_interface_contacts = &(reference_interface_contacts_capri[orientation]);
		else
			ref_interface_contacts = &(reference_interface_contacts[orientation]);

		num_receptor_interface_points[orientation] = 0;
		num_ligand_interface_points[orientation] = 0;

		hash_set<int,hash<int>,eqint> added_ligand_residue, added_receptor_residue;
		// assuming that the hash_set iterator lists elements in the same order every time
		for(hash_set<long,hash<long>,eqlong>::iterator itr = ref_interface_contacts->begin(); itr != ref_interface_contacts->end(); itr++){
			long index = *itr;
			int refrindex = index/MAX_ATOMS;
			int reflindex = index % MAX_ATOMS;

			if(added_receptor_residue.count(refrindex) == 0){
				added_receptor_residue.insert(refrindex);
				int rindex = (*reference_vs_receptor)[refrindex];
				Aminoacid *aa_ref = reference->aminoacid[refrindex];
				Aminoacid *aa = c->aminoacid[rindex];

				if(aa->atom.count("N") > 0 && aa_ref->atom.count("N") > 0 && use_all_bbatom_for_rmsd){
					receptor_interface_points[orientation][num_receptor_interface_points[orientation]++] = Vector(aa->atom["N"]->position);
				} 
				if(aa->atom.count("CA") > 0 && aa_ref->atom.count("CA") > 0){
					receptor_interface_points[orientation][num_receptor_interface_points[orientation]++] = Vector(aa->atom["CA"]->position);
				} 
				if(aa->atom.count("C") > 0 && aa_ref->atom.count("C") > 0 && use_all_bbatom_for_rmsd){
					receptor_interface_points[orientation][num_receptor_interface_points[orientation]++] = Vector(aa->atom["C"]->position);
				} 
				if(aa->atom.count("O") > 0 && aa_ref->atom.count("O") > 0 && use_all_bbatom_for_rmsd){
					receptor_interface_points[orientation][num_receptor_interface_points[orientation]++] = Vector(aa->atom["O"]->position);
				}
			}
			if(added_ligand_residue.count(reflindex) == 0){
				added_ligand_residue.insert(reflindex);
				int lindex = (*reference_vs_ligand)[reflindex];
				Aminoacid *aa_ref = reference->aminoacid[reflindex];
				Aminoacid *aa = ligand->c->aminoacid[lindex];

				if(aa->atom.count("N") > 0 && aa_ref->atom.count("N") > 0 && use_all_bbatom_for_rmsd){
					ligand_interface_points[orientation][num_ligand_interface_points[orientation]++] = Vector(aa->atom["N"]->position);
				} 
				if(aa->atom.count("CA") > 0 && aa_ref->atom.count("CA") > 0){
					ligand_interface_points[orientation][num_ligand_interface_points[orientation]++] = Vector(aa->atom["CA"]->position);
				} 
				if(aa->atom.count("C") > 0 && aa_ref->atom.count("C") > 0 && use_all_bbatom_for_rmsd){
					ligand_interface_points[orientation][num_ligand_interface_points[orientation]++] = Vector(aa->atom["C"]->position);
				} 
				if(aa->atom.count("O") > 0 && aa_ref->atom.count("O") > 0 && use_all_bbatom_for_rmsd){
					ligand_interface_points[orientation][num_ligand_interface_points[orientation]++] = Vector(aa->atom["O"]->position);
				}
			}
		}
	}

	short ref_num_contacts[num_symmetries];
	for(unsigned int orientation=0; orientation < num_symmetries; orientation++){
		/*if(interface_capridefn)	ref_num_contacts[orientation] = (reference_interface_contacts_capri[orientation].size()>0?reference_interface_contacts_capri[orientation].size():1);
		else*/ ref_num_contacts[orientation] = (reference_interface_contacts[orientation].size()>0?reference_interface_contacts[orientation].size():1);
	}

	Transformation *reference_tr[num_symmetries];
	Transformation *ligand_orig_to_ref = trs[1];
	Transformation *receptor_orig_to_ref = trs[0];	
	Reference_Frame *receptor_ref_to_orig = Reference_Frame::invert(receptor_orig_to_ref);
	Vector reftr_ligandcm[num_symmetries];
	for(unsigned int orientation=0; orientation < num_symmetries; orientation++){
		unsigned short homodimer_orientation = (orientation/(num_receptor_symmetries*num_ligand_symmetries));
		unsigned short receptor_orientation = (orientation/num_ligand_symmetries) % num_receptor_symmetries;
		unsigned short ligand_orientation = orientation % num_ligand_symmetries;

		Reference_Frame *rf = new Reference_Frame((Reference_Frame*) ligand_orig_to_ref);
		if(ligand_orientation > 0)
			rf = Reference_Frame::compose(ligand_symmetries[ligand_orientation-1],rf);
		if(homodimer_orientation > 0){
			rf = Reference_Frame::compose(ligand_vs_receptor_symmetries[homodimer_orientation-1],rf);
			rf = Reference_Frame::compose(ligand_vs_receptor_symmetries[homodimer_orientation-1],rf);
		}
		if(receptor_orientation > 0){
			Reference_Frame *inverse = Reference_Frame::invert(receptor_symmetries[receptor_orientation-1]);
			rf = Reference_Frame::compose(inverse,rf);
		}
		rf = Reference_Frame::compose(receptor_ref_to_orig,rf);

		rf = Reference_Frame::invert(rf);
		reference_tr[orientation] = new Transformation(rf->translation, rf->ex, rf->ey, rf->scale, 0, 0-orientation);
		reftr_ligandcm[orientation] = reference_tr[orientation]->inverse_transform(*(ligand->c->center_of_mass)); 

		/*if(orientation == 0){
			trs[2]->print_details(&cout,TN_BASIC);
			Transformation *tr = generate_transformation(this,receptor_orig_to_ref, ligand, ligand_orig_to_ref,0);
			tr->print_details(&cout,TN_BASIC);
			rf = new Reference_Frame(receptor_orig_to_ref);
			rf = Reference_Frame::compose(Reference_Frame::invert(ligand_orig_to_ref),rf);
			tr = new Transformation(rf->translation, rf->ex, rf->ey, rf->scale, 0, 0);
			tr->print_details(&cout,TN_BASIC);
			reference_tr[0]->print_details(&cout,TN_BASIC);
		}*/
	}

	bool fast_rmsd = false;
	double R_receptor_interface[num_symmetries][3][3], R_ligand_interface[num_symmetries][3][3];
	Vector sum_receptor_interface[num_symmetries],sum_ligand_interface[num_symmetries],sum_ligaligned_interface[num_symmetries]; 
	Vector avg_aligned_interface[num_symmetries];
	float avg_aligned_interface_array[num_symmetries][3], sum_ligaligned_interface_array[num_symmetries][3];
	float reference_tr_U[num_symmetries][3][3], reference_tr_t[num_symmetries][3];
	double sum_square_interface_pre[num_symmetries];

	if(fast_rmsd){
		for(unsigned int orientation=0; orientation < num_symmetries; orientation++){
			// compute the R matrix for receptor points
			for(int i = 0; i <3; i++)
				for(int j = 0; i <3; i++)
					R_receptor_interface[orientation][i][j] = 0;
			sum_receptor_interface[orientation] = Vector(0,0,0);

			for(int i = 0 ; i < num_receptor_interface_points[orientation]; i++){
				Vector v = receptor_interface_points[orientation][i];
				R_receptor_interface[orientation][0][0] += v.x*v.x;
				R_receptor_interface[orientation][0][1] += v.x*v.y;
				R_receptor_interface[orientation][0][2] += v.x*v.z;
				R_receptor_interface[orientation][1][1] += v.y*v.y;
				R_receptor_interface[orientation][1][2] += v.y*v.z;
				R_receptor_interface[orientation][2][2] += v.z*v.z;
				sum_receptor_interface[orientation] = sum_receptor_interface[orientation] + v; 
			}
			R_receptor_interface[orientation][1][0] = R_receptor_interface[orientation][0][1];
			R_receptor_interface[orientation][2][0] = R_receptor_interface[orientation][0][2];
			R_receptor_interface[orientation][2][1] = R_receptor_interface[orientation][1][2];

			// compute the ligand points for the "correct" solution
			for(int i = 0; i <3; i++)
				for(int j = 0; i <3; i++)
					R_ligand_interface[orientation][i][j] = 0;

			sum_ligaligned_interface[orientation] = Vector(0,0,0);
			sum_ligand_interface[orientation] = Vector(0,0,0);

			unsigned short homodimer_orientation = (orientation/(num_receptor_symmetries*num_ligand_symmetries));
			unsigned short receptor_orientation = (orientation/num_ligand_symmetries) % num_receptor_symmetries;
			unsigned short ligand_orientation = orientation % num_ligand_symmetries;

			for(int i = 0 ; i < num_ligand_interface_points[orientation]; i++){
				Vector v = ligand_interface_points[orientation][i];
				R_ligand_interface[orientation][0][0] += v.x*v.x;
				R_ligand_interface[orientation][0][1] += v.x*v.y;
				R_ligand_interface[orientation][0][2] += v.x*v.z;
				R_ligand_interface[orientation][1][1] += v.y*v.y;
				R_ligand_interface[orientation][1][2] += v.y*v.z;
				R_ligand_interface[orientation][2][2] += v.z*v.z;
				sum_ligand_interface[orientation] = sum_ligand_interface[orientation] + v;

				v = ligand_orig_to_ref->transform(v);
				if(ligand_orientation > 0)
					v = ligand_symmetries[ligand_orientation-1]->transform(v);
				if(homodimer_orientation > 0){
					v = ligand_vs_receptor_symmetries[homodimer_orientation-1]->transform(v);
					v = ligand_vs_receptor_symmetries[homodimer_orientation-1]->transform(v);
				}
				if(receptor_orientation > 0)
					reflig_interface_points[orientation][i] = receptor_symmetries[receptor_orientation-1]->inverse_transform(reflig_interface_points[orientation][i]);

				v = receptor_orig_to_ref->inverse_transform(v);

				sum_ligaligned_interface[orientation] = sum_ligaligned_interface[orientation] + v;
			}
			R_ligand_interface[orientation][1][0] = R_ligand_interface[orientation][0][1];
			R_ligand_interface[orientation][2][0] = R_ligand_interface[orientation][0][2];
			R_ligand_interface[orientation][2][1] = R_ligand_interface[orientation][1][2];

			avg_aligned_interface[orientation] = (sum_receptor_interface[orientation]+sum_ligaligned_interface[orientation]) * (1.0/(num_receptor_interface_points[orientation]+num_ligand_interface_points[orientation]));

			// convert the Vector representation to matrices for more readable code
			reference_tr[orientation]->translation->toarray(reference_tr_t[orientation]);
			reference_tr[orientation]->ex->toarray(reference_tr_U[num_symmetries][0]);
			reference_tr[orientation]->ey->toarray(reference_tr_U[num_symmetries][1]);
			reference_tr[orientation]->ez->toarray(reference_tr_U[num_symmetries][2]); 

			sum_ligaligned_interface[orientation].toarray(sum_ligaligned_interface_array[orientation]); 
			avg_aligned_interface[orientation].toarray(avg_aligned_interface_array[orientation]); 

			sum_square_interface_pre[orientation] = 2*(R_receptor_interface[orientation][0][0]+R_receptor_interface[orientation][1][1]+R_receptor_interface[orientation][2][2]);
			double quadratic_term = 0;
			for(int i = 0; i < 3; i++)
				for(int j = 0; j < 3; j++)
					for(int i1 = 0; i1 < 3; i1++)
						for(int j1 = 0; j1 < 3; j1++)
							quadratic_term += reference_tr_U[orientation][i1][i] * reference_tr_U[orientation][j1][j] * R_ligand_interface[orientation][i1][j1];

			double linear_term = Vector(reference_tr_U[orientation][0][0], reference_tr_U[orientation][1][0], reference_tr_U[orientation][2][0]).dot(sum_ligand_interface[orientation]) * reference_tr[orientation]->translation->x 
					+ Vector(reference_tr_U[orientation][0][1], reference_tr_U[orientation][1][1], reference_tr_U[orientation][2][1]).dot(sum_ligand_interface[orientation]) * reference_tr[orientation]->translation->y
					+ Vector(reference_tr_U[orientation][0][2], reference_tr_U[orientation][1][2], reference_tr_U[orientation][2][2]).dot(sum_ligand_interface[orientation]) * reference_tr[orientation]->translation->z;

			sum_square_interface_pre[orientation] += quadratic_term + linear_term + reference_tr[orientation]->translation->norm_squared()*(num_ligand_interface_points[orientation]);
		}
	}

	// process transformations

	bool compute_contact_metrics = true; //false;
	float rmsd[num_symmetries], lrmsd[num_symmetries],irmsd[num_symmetries], model_contact_energy_residue, model_contact_energy_coarse;
	short num_native_contacts_predicted[num_symmetries], num_contacts;
	for(vector<Transformation*>::iterator titr = matching_frames->begin(); titr != matching_frames->end(); titr++){
		Transformation *tr = *titr;
		tr->vmetrics = new VerificationMetrics();
		*out << tr->frame_number << " ";

		// convert vector representation to matrix representation
		float U[3][3], translation[3];
		tr->ex->toarray(U[0]);
		tr->ey->toarray(U[1]);
		tr->ez->toarray(U[2]);
		tr->translation->toarray(translation);

		// find irmsd taking symmetry into account
		for(unsigned int orientation=0; orientation < num_symmetries; orientation++){
			if(fast_rmsd){
				// compute R using precomputed values
				Vector v = sum_ligand_interface[orientation];
				v = tr->inverse_transform(v);
				Vector sum_model = sum_receptor_interface[orientation] + v;

				float v_array[3];
				v.toarray(v_array);
				float sum_model_array[3];
				sum_model.toarray(sum_model_array);

				double R[3][3];
				for(int i = 0; i < 3; i++)
					for(int j = 0; j < 3; j++){
						double quadratic_term = 0;
						for(int i1 = 0; i1 < 3; i1++)
							for(int j1 = 0; j1 < 3; j1++)
								quadratic_term += reference_tr_U[orientation][i1][i] * U[j1][j] * R_ligand_interface[orientation][i1][j1];

						R[i][j] = R_receptor_interface[orientation][i][j] + quadratic_term + sum_ligaligned_interface_array[orientation][i]*translation[j] 
						                                                                                                                                + v_array[j]*reference_tr_t[orientation][i] - avg_aligned_interface_array[orientation][i]*sum_model_array[j];
					}

				double sum_square_interface=0;
				{
					double quadratic_term = 0;
					for(int i = 0; i < 3; i++)
						for(int j = 0; j < 3; j++)
							for(int i1 = 0; i1 < 3; i1++)
								for(int j1 = 0; j1 < 3; j1++)
									quadratic_term += U[i1][i] * U[j1][j] * R_ligand_interface[orientation][i1][j1];

					double linear_term = Vector(U[0][0], U[1][0], U[2][0]).dot(sum_ligand_interface[orientation]) * tr->translation->x 
							+ Vector(U[0][1], U[1][1], U[2][1]).dot(sum_ligand_interface[orientation]) * tr->translation->y
							+ Vector(U[0][2], U[1][2], U[2][2]).dot(sum_ligand_interface[orientation]) * tr->translation->z;

					double translation_term = tr->translation->norm_squared();
					sum_square_interface += quadratic_term + linear_term + translation_term*(num_ligand_interface_points[orientation]);
				}

				cubicsolve_rmsd_fromRE(R, sum_square_interface_pre[orientation]+sum_square_interface, num_receptor_interface_points[orientation]+num_ligand_interface_points[orientation], &irmsd[orientation]);
				*out << "fast " << irmsd[orientation];
			}

			Vector model_interface_points[max_num_atoms+1];
			for(int i = 0 ; i < num_receptor_interface_points[orientation]; i++)
				model_interface_points[i] = receptor_interface_points[orientation][i];
			for(int i = 0 ; i < num_ligand_interface_points[orientation]; i++)
				model_interface_points[i+num_receptor_interface_points[orientation]] = tr->inverse_transform(ligand_interface_points[orientation][i]);

			irmsd[orientation] = compute_rmsd(num_receptor_interface_points[orientation]+num_ligand_interface_points[orientation], ref_interface_points[orientation], &model_interface_points[0]);
			if(fast_rmsd)	*out << "\tregular " << irmsd[orientation] << endl;
		}

		short orientationselected=0;
		float minirmsd=irmsd[0];
		for(short orientation=0; orientation < num_symmetries; orientation++){
			if(irmsd[orientation] < minirmsd){
				orientationselected=orientation;
				minirmsd=irmsd[orientation];
			}
			*out << "|" << irmsd[orientation] << "| ";
		}
		*out << orientationselected << endl;

		for(int i = 0; i < num_ligand_points; i++){
			point_index = ligand_start + i;
			model_points[point_index] = ligand_points[i];

			// dock transformation
			model_points[point_index] = tr->inverse_transform(model_points[point_index]);

			// computing L_rms transform ligand such that receptor is aligned
			//model_points[point_index] = receptor_orig_to_ref->transform(model_points[point_index]);
		}

		if(minirmsd < 10.0){
			unsigned int orientation=orientationselected;

			//cout << ligand_vs_reference->size() << " ";
			//*out << "lrmsd # " << ligand_start << " " << point_index << " ";//out->flush();
			lrmsd[orientation] = compute_lrmsd(ligand_start + num_ligand_points, ligand_start, reference_points[orientation],&model_points[0], &(rmsd[orientation]));
			// no minimization for L_rms
			/*double crmsd = 0;
			for(int i = ligand_start; i < ligand_start + num_ligand_points; i++){
				Vector v = reference_points[orientation][i];
				Vector w = model_points[i];
				crmsd += Vector::distance_squared(v,w);
			}
			rmsd[orientation] = lrmsd[orientation] = sqrt(crmsd/(point_index - ligand_start));*/
			//*out << tr->vmetrics->rmsd << " " << tr->vmetrics->lrmsd << endl;
		} else {
			rmsd[orientationselected] = lrmsd[orientationselected] = irmsd[orientationselected];
		}

		tr->vmetrics->rmsd = rmsd[orientationselected];
		tr->vmetrics->lrmsd = lrmsd[orientationselected];
		tr->vmetrics->irmsd = irmsd[orientationselected];

		// computing residue contacts is currently the most expensive step so do it on demand
		model_contact_energy_residue=0,model_contact_energy_coarse=0;
		hash_set<long,hash<long>,eqlong> model_interface_contacts;
		if(compute_contact_metrics){
			hash_set<int,hash<int>,eqint> model_interface_receptor_residues, model_interface_ligand_residues;
			for(int i = 0; i < c->num_aminoacids; i++){
				Aminoacid *ra = c->aminoacid[i];
				if(receptor_vs_reference->count(ra->cindex) > 0){
					for(int j = 0; j < ligand->c->num_aminoacids; j++){
						Aminoacid *la = ligand->c->aminoacid[j];
						if(ligand_vs_reference->count(la->cindex) > 0){
							bool contact = false;
							/*if(interface_capridefn){
					 			for(hash_map<const char*, Atom*, hash<const char*>, eqstr>::iterator raitr = ra->atom.begin(); (raitr != ra->atom.end()) && !contact ; raitr++){
					 		  		Atom *ratom = (Atom *) raitr->second;
					 		   		for(hash_map<const char*, Atom*, hash<const char*>, eqstr>::iterator laitr = la->atom.begin(); (laitr != la->atom.end()) && !contact; laitr++){
					 					Atom *latom = (Atom *) laitr->second;
					 					Vector v = Vector(*(latom->position));
										v = tr->inverse_transform(v);
					 					if(Vector::distance(*(ratom->position), v) <= 10.0)	contact = true;
		 					   		}
						 		}
				 			} else*/ {
								if(ra->centroid != NULL && la->centroid != NULL){
									Vector v = Vector(*(la->centroid));
									v = tr->inverse_transform(v);
									if(Vector::distance(ra->centroid,v) < SS_CUTOFF)	contact=true;
								}
								if(ra->alpha_carbon != NULL && la->alpha_carbon != NULL){
									Vector v = Vector(*(la->alpha_carbon->position));
									v = tr->inverse_transform(v);
									if(Vector::distance(ra->alpha_carbon->position,v) < SS_CUTOFF)	contact=true;
								}
				 			}
				 			if(contact) {
				 				model_interface_ligand_residues.insert(la->cindex);
				 				model_interface_receptor_residues.insert(ra->cindex);
				 				long cindex = (*receptor_vs_reference)[ra->cindex]*MAX_ATOMS + (*ligand_vs_reference)[la->cindex];
				 				if(model_interface_contacts.count(cindex) == 0){
				 					model_interface_contacts.insert(cindex);
				 					if(ra->type >= 0 && la->type>= 0)
				 						model_contact_energy_residue += residue_potential[ra->type][la->type];
				 					int ss1 = SSTYPE(ra->sstructure), ss2 = SSTYPE(la->sstructure);
				 					int crs1 = RTYPE(ra->type) *NUM_COARSE_SSTYPES + ss1;
				 					int crs2 = RTYPE(la->type) *NUM_COARSE_SSTYPES + ss2;
				 					if(ss1 != DSSP_U && ss2 != DSSP_U){
				 						//model_contact_energy_coarse += coarse_potential[crs1][crs2];
				 					}
				 				}
				 			}
						}
					}
				}
			}
			num_contacts = model_interface_contacts.size();

			num_native_contacts_predicted[orientationselected]=0;
			for(hash_set<long,hash<long>,eqlong>::iterator itr = reference_interface_contacts[orientationselected].begin(); itr != reference_interface_contacts[orientationselected].end(); itr++){
				long index = *itr;	
				if(model_interface_contacts.count(index) > 0)	num_native_contacts_predicted[orientationselected]++;
			}

			*out << " pairp " << model_contact_energy_coarse << " "; //out->flush();
			*out << num_native_contacts_predicted[orientationselected] << " " << ref_num_contacts[orientationselected] << endl;

			tr->vmetrics->frac_native_contacts_predicted = ((float)num_native_contacts_predicted[orientationselected])/ref_num_contacts[orientationselected];
			tr->vmetrics->frac_nonnative_contacts_predicted = ((float)(num_contacts - num_native_contacts_predicted[orientationselected]))/num_contacts;
			tr->num_contacts = num_contacts;
			tr->eResiduepair = model_contact_energy_residue;
		}

		Vector tr_ligandcm = tr->inverse_transform(*(ligand->c->center_of_mass));
		tr->distance(&(tr->vmetrics->delta_U),reference_tr[orientationselected]);
		tr->vmetrics->delta_r = Vector::distance(reftr_ligandcm[orientationselected],tr_ligandcm);

		if(count++ % 100000 == 0)	*out << "verified " << count << endl;
	}
}

/*
 * U applied to reference_points
 * Tr applied to model points
 */
float compute_rmsd(int num_points, Vector* reference_points, Vector* model_points, Transformation *tr){
	float ref_xlist[num_points+1][3], mov_xlist[num_points+1][3];
	float mov_com[3], mov_to_ref[3], U[3][3];
	float rmsd;

	for(int i = 0; i < num_points; i++){
		Vector v = reference_points[i];
		ref_xlist[i][0] = v.x;
		ref_xlist[i][1] = v.y;
		ref_xlist[i][2] = v.z;

		v = model_points[i];
		mov_xlist[i][0] = v.x;
		mov_xlist[i][1] = v.y;
		mov_xlist[i][2] = v.z;
	}

	// function writes on the arrays passed
	int rank;
	calculate_rotation_rmsd(ref_xlist, mov_xlist, num_points, mov_com, mov_to_ref,&rank,U,&rmsd);

	Vector *ex = new Vector(U[0][0],U[0][1],U[0][2]); 
	Vector *ey = new Vector(U[1][0],U[1][1],U[1][2]); 
	ex->normalize();ey->normalize();
	Vector *ez = new Vector(ex->cross(ey));
	tr->ex = ex; tr->ey =ey; tr->ez = ez;

	float ref_com[3];
	for (int i=0; i<3; i++)	ref_com[i] = mov_to_ref[i]+ mov_com[i]; 

	Vector tmov = Vector(mov_com[0],mov_com[1],mov_com[2]);
	Vector tref = Vector(ref_com[0],ref_com[1],ref_com[2]);
	float x = tref.dot(Vector(ex->x,ey->x,ez->x));
	float y = tref.dot(Vector(ex->y,ey->y,ez->y));
	float z = tref.dot(Vector(ex->z,ey->z,ez->z));
	Vector Ui_tref = Vector(x,y,z);

	tr->translation = new Vector(tmov - Ui_tref);

	/*{
	 * cout << mov_to_ref[0] << " " << mov_to_ref[1] << " " << mov_to_ref[2] << endl; cout.flush();
		cout << ex->x << " " << ex->y << " " << ex->z << endl; cout.flush();
		cout << ey->x << " " << ey->y << " " << ey->z << endl; cout.flush();

		double crmsd = 0;
		for(int i = 0; i < num_points; i++){
			Vector v = reference_points[i];
			Vector w = tr->transform(model_points[i]);
			float d = Vector::distance(v,w);  
			crmsd += d*d;
		}
		cout << rmsd << " " << sqrt(crmsd/num_points) << endl;
	}*/

	return rmsd;
}

void compute_rmsd(int total_num_atoms, int ligand_start,Vector* reference_points, Vector* receptor_points,
		Vector* ligand_points, Transformation *tr){
	float ref_xlist[total_num_atoms+1][3], mov_xlist[total_num_atoms+1][3];
	float mov_com[3], mov_to_ref[3], U[3][3];
	float rmsd;
	int ligand_total_atoms = total_num_atoms - ligand_start;	

	//*out << total_num_atoms << " " << ligand_start << endl;
	for(int i = 0; i < total_num_atoms; i++){
		Vector v = reference_points[i];
		ref_xlist[i][0] = v.x;
		ref_xlist[i][1] = v.y;
		ref_xlist[i][2] = v.z;
	}
	for(int i = 0; i < ligand_start; i++){
		Vector v = receptor_points[i];
		mov_xlist[i][0] = v.x;
		mov_xlist[i][1] = v.y;
		mov_xlist[i][2] = v.z;
	}
	for(int i = 0 ; i < ligand_total_atoms; i++){
		Vector v = ((Reference_Frame*) tr)->inverse_transform(ligand_points[i]);
		int atom_index = i + ligand_start;
		mov_xlist[atom_index][0] = v.x;
		mov_xlist[atom_index][1] = v.y;
		mov_xlist[atom_index][2] = v.z;
		//*out << v.x << " " << v.y << " " << v.z << endl;
	}

	// function writes on the arrays passed
	int rank;
	calculate_rotation_rmsd(ref_xlist, mov_xlist, total_num_atoms, mov_com, mov_to_ref,&rank,U,&rmsd);

	// compute the rmsd between ligand points
	float lrmsd = 0;
	if(rank == 3){
		Vector ux = Vector(U[0][0],U[0][1],U[0][2]);
		ux.normalize();
		Vector uy = Vector(U[1][0],U[1][1],U[1][2]);
		uy.normalize();
		Vector uz = Vector(U[2][0],U[2][1],U[2][2]);
		uz.normalize();
		for(int i = 0 ; i < ligand_total_atoms; i++){
			Vector v = ((Reference_Frame*) tr)->inverse_transform(ligand_points[i]);
			v = v + (Vector(-mov_com[0],-mov_com[1],-mov_com[2]));
			float tvx,tvy,tvz;
			tvx = v.dot(ux);
			tvy = v.dot(uy);
			tvz = v.dot(uz);
			v = reference_points[i + ligand_start];
			v = v + (Vector(-(mov_to_ref[0]+mov_com[0]),
					-(mov_to_ref[1]+mov_com[1]),
					-(mov_to_ref[2]+mov_com[2])));
			lrmsd += ((tvx - v.x)*(tvx - v.x) + (tvy - v.y)*(tvy - v.y) + (tvz - v.z)*(tvz - v.z));
		}
		lrmsd = lrmsd/total_num_atoms;
		lrmsd = sqrt(lrmsd);

		//*out << rank << " " << rmsd << " " << lrmsd << endl;
		tr->vmetrics->rmsd = rmsd;
		tr->vmetrics->lrmsd = lrmsd;
	}
}

/*
transform receptor to ligand's reference frame ...
to compute 	t2i s2i r2i r1 s1 t1
		= t2i r2i r1 s2i s1 t1
		= (t2i (r2i r1 s2i s1 t1)) (r2i r1) (s2i s1)
		= t r s
		= r s translation	
try local minimization over the matching points?	
 */
Transformation* Object::generate_transformation(Object *receptor, Reference_Frame *rf1, Object *ligand, Reference_Frame *rf2, int votes){
	//, vector<int> lresidue_indices,vector<int> rresidue_indices){
	float scale = 1; /*(rf2->scale / rf1->scale);
	if(scale > SCALE_MAX || scale < SCALE_MIN)
		return NULL;*/

	Vector ux = Vector((rf1->ex)->x,(rf1->ey)->x,(rf1->ez)->x);
	Vector uy = Vector((rf1->ex)->y,(rf1->ey)->y,(rf1->ez)->y);
	Vector uz = Vector((rf1->ex)->z,(rf1->ey)->z,(rf1->ez)->z);

	Vector vx = Vector((rf2->ex)->x,(rf2->ey)->x,(rf2->ez)->x);
	Vector vy = Vector((rf2->ex)->y,(rf2->ey)->y,(rf2->ez)->y);
	Vector vz = Vector((rf2->ex)->z,(rf2->ey)->z,(rf2->ez)->z);

	Vector ex = Vector( vx.dot(ux), vx.dot(uy), vx.dot(uz));
	Vector ey = Vector( vy.dot(ux), vy.dot(uy), vy.dot(uz));
	Vector ez = Vector( vz.dot(ux), vz.dot(uy), vz.dot(uz)); // = ex.cross(ey))

	Vector t = Vector(*(rf2->translation));
	t = t * scale;
	float x = t.dot(Vector(ex.x,ey.x,ez.x));
	float y = t.dot(Vector(ex.y,ey.y,ez.y));
	float z = t.dot(Vector(ex.z,ey.z,ez.z));
	t = Vector(rf1->translation) - Vector(x,y,z) ;
	Vector *translation = new Vector(t);

	Transformation *tr = new Transformation(translation, new Vector(ex), new Vector(ey), scale, votes, transformation_index++);	
	tr->compute_coordinates((void *)receptor, (void *)ligand);
	return tr;
}

/*
 * Compute the reduced representation of the rotation matrix
 * represent ex and ey as points on the unit sphrere using theta and phi
 * Also compute compute cmr in spherical coordinates
 */
void Transformation::compute_coordinates(void *receptor, void *ligand){	
	cmr = new Vector(*(((Object*) ligand)->c->center_of_mass) - transform(*(((Object*) receptor)->c->center_of_mass))); 
	Vector v = Vector(cmr);
	v.normalize();

	cmr_phi = acos(v.z);
	if(v.x != 0){
		cmr_theta = atan(v.y/v.x);
	} else
		cmr_theta = v.y > 0 ? PI/2.0 : -PI/2.0;

	if(cmr_theta > 0){
		if(v.x < 0)
			cmr_theta = PI - cmr_theta;
	} else {
		if(v.x < 0)
			cmr_theta =  PI - cmr_theta;
		else
			cmr_theta = 2*PI + cmr_theta;
	}
}

/*
 * Compute residue in contact, residue contact energy, residue eInterface
 */
void Object::compute_energy_approx3(Object *receptor, Object *ligand, Transformation *tr, ProtProtDetailedScoringMetrics *details){
	int num_residue_types = NUM_RESIDUE_TYPES;
	for(int i = 0; i < num_residue_types; i++)
		for(int j = 0 ; j < num_residue_types; j++)
			details->residue_contacts_core[i][j] = details->residue_contacts_rim[i][j] = 0;
	for(int i = 0; i < NUM_ENTROPY_DIVISIONS; i++)
		for(int j = 0 ; j < NUM_ENTROPY_DIVISIONS; j++)
			details->conserv_contacts[i][j] = 0;
	num_residue_types = NUM_COARSE_RTYPES;
	for(int i = 0; i < num_residue_types; i++)
		for(int j = 0 ; j < num_residue_types; j++)
			for(int k = 0 ; k < num_residue_types; k++)
				details->threebody_contacts[i][j][k] = 0;

	vector<long> aminoacid_contacts;
	tr->sEvolution_interface = 0;

	// variations of residue contact potential used in the fourier transform
	{
		for(int laindex = 0; laindex < ligand->c->num_aminoacids; laindex++){
			Aminoacid *al = ligand->c->aminoacid[laindex];

			if(al->type >= 0 && al->centroid != NULL){
				Vector vl = tr->inverse_transform(*(al->centroid));
				for(int raindex = 0; raindex < receptor->c->num_aminoacids; raindex++){
					Aminoacid *ar = receptor->c->aminoacid[raindex];
					float d, d2, factor=0,nfactor;
					bool contact=false;

					// residue contacts (centroid centroid)
					if(ar->type >= 0 && ar->centroid != NULL){
						d2 = Vector::distance_squared(vl,*(ar->centroid));
#ifdef 	STEP_POTENTIAL 			
						if(d2 < AA_CUTOFF_SQUARED)
#endif
#ifdef LINEAR_SPLINE_POTENTIAL										 			
							if(d2 < AA_SMTHP_CUTOFF_SQUARED)
#endif												
							{
								contact=true;
#ifdef 	STEP_POTENTIAL 	
								factor = 1.0;
#endif
#ifdef LINEAR_SPLINE_POTENTIAL	 									
								if(d2 < AA_CUTOFF_SQUARED)	factor = 1.0;
								else{ float d=sqrt(d2);	factor = AA_SMTHP_FACTOR(d); }
								//cout << ar->index << " " << al->index << " " << d2 << " factor " << factor << endl;
#endif						
							}
					}//*/

					/*/ residue contacts (heavy atom centroid)
#ifdef 	STEP_POTENTIAL					
					for(hash_map<const char*, Atom*, hash<const char*>, eqstr>::iterator aitr = ar->atom.begin(); !contact && (aitr != ar->atom.end()); aitr++)
#endif
#ifdef LINEAR_SPLINE_POTENTIAL
					for(hash_map<const char*, Atom*, hash<const char*>, eqstr>::iterator aitr = ar->atom.begin(); aitr != ar->atom.end(); aitr++)
#endif
					{
						Atom *a = aitr->second;
						d2 = Vector::distance_squared(vl,*(a->position));
#ifdef 	STEP_POTENTIAL 			
						if(d2 < AA_CUTOFF_SQUARED)
#endif
#ifdef LINEAR_SPLINE_POTENTIAL										 			
						if(d2 < AA_SMTHP_CUTOFF_SQUARED)
#endif												
						{	
							contact=true;								
#ifdef 	STEP_POTENTIAL 	
							factor = 1.0;
#endif
#ifdef LINEAR_SPLINE_POTENTIAL	 									
				 			if(d2 < AA_CUTOFF_SQUARED)	nfactor = 1.0;
							else{ float d=sqrt(d2);	nfactor = AA_SMTHP_FACTOR(d); }
							if(factor < nfactor)	factor = nfactor;
#endif									
						}
					}//*/

					if(contact){
						short ratype=ar->type, latype=al->type;
						if(ratype >= 0 && latype>= 0)
							if(ratype <= latype)
								details->residue_contacts_core[ratype][latype] += factor;
							else
								details->residue_contacts_core[latype][ratype] += factor;
						aacontact_core[laindex][raindex] = true;
						aminoacid_contacts.push_back(laindex*MAX_ATOMS + raindex);
						if(!aarcontact[raindex]){
							aarcontact[raindex] = true;
						}
						if(!aalcontact[laindex]){
							aalcontact[laindex] = true;
						}
					}
				}
			}

			// backbone centroid and backbone backbone contacts
			{
				Vector *vaa[3];
				vaa[0] = al->centroid;
				vaa[1] = (al->amide_nitrogen == NULL) ? NULL : al->amide_nitrogen->position;
				vaa[2] = (al->carbonyl_oxygen == NULL) ? NULL : al->carbonyl_oxygen->position;
				for(int aapi = 0; aapi < 3; aapi++) 
					if(vaa[aapi]  != NULL){
						Vector vl = tr->inverse_transform(*(vaa[aapi]));
						for(int raindex = 0; raindex < receptor->c->num_aminoacids; raindex++){
							Aminoacid *ra = receptor->c->aminoacid[raindex];
							if(ra->type >= 0)
								for(int aapj = 0; aapj < 3; aapj++){
									float d, d2, factor=0;
									if((aapi == 0 && aapj != 0 && al->type >= 0) ||
											(aapi != 0 && aapj == 0 && ra->centroid != NULL)) {
										bool computefactor=false;
										if(aapi == 0 && aapj == 1 && ra->amide_nitrogen != NULL){
											computefactor=true;
											d2=Vector::distance_squared(vl,*(ra->amide_nitrogen->position));
										}
										if(aapi == 0 && aapj == 2 && ra->carbonyl_oxygen != NULL){
											computefactor=true;
											d2=Vector::distance_squared(vl,*(ra->carbonyl_oxygen->position));
										}
										if(aapi != 0){
											computefactor=true;
											d2=Vector::distance_squared(vl,*(ra->centroid));
										}

										if(computefactor)
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
													if(d2 < BS_CUTOFF_SQUARED)	factor = 1.0;
													else{ float d=sqrt(d2);	factor = BS_SMTHP_FACTOR(d); }
#endif						
												}
										if(aapi == 0 && aapj != 0 && al->type >= 0){
											details->residue_contacts_core[al->type][19+aapj] += factor;
											details->residue_contacts_core[19+aapj][al->type] += factor;
										} else {
											details->residue_contacts_core[ra->type][19+aapi] += factor;
											details->residue_contacts_core[19+aapi][ra->type] += factor;
										}
									}

									if(aapi > 0 && aapj > 0){
										bool computefactor=false;
										if(aapj == 1 && ra->amide_nitrogen != NULL){
											computefactor=true;
											d2=Vector::distance_squared(vl,*(ra->amide_nitrogen->position));
										}
										if(aapj == 2 && ra->carbonyl_oxygen != NULL){
											computefactor=true;
											d2=Vector::distance_squared(vl,*(ra->carbonyl_oxygen->position));
										}

										if(computefactor)
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
													if(d2 < BB_CUTOFF_SQUARED)	factor = 1.0;
													else{ float d=sqrt(d2);	factor = BB_SMTHP_FACTOR(d); }
#endif		
												}
										details->residue_contacts_core[19+aapi][19+aapj] += factor;
										if(aapi != aapj)	details->residue_contacts_core[19+aapj][19+aapi] += factor;
									}
								} // aapj
						}
					}
			}
		}
		/*{
			float sum = (tr->eVdw - 9.0 * tr->eVdw_repulsion)*vdw_weight;
			for(unsigned short i = 0; i < 22; i++)
				for(unsigned short j = i; j < 22; j++)
					sum += details->residue_contacts_rim[i][j] * residue_bkbn_potential[i][j];

			if(tr->frame_number == 0){
				cout << "pot " << vdw_weight << endl;
				for(unsigned short i = 0; i < 20; i++)
					for(unsigned short j = i; j < 20; j++)	cout << "pot " << residue_bkbn_potential[i][j] << endl;
				for(unsigned short i = 0; i < 20; i++)	cout << "pot " << residue_bkbn_potential[20][i] << endl;
				for(unsigned short i = 0; i < 20; i++)	cout << "pot " << residue_bkbn_potential[21][i] << endl;
				cout << "pot " << residue_bkbn_potential[20][20] << "\t" << residue_bkbn_potential[20][21] << "\t"
					<< residue_bkbn_potential[21][21] << endl;
			}

			cout << "lcscore " << tr->frame_number << " " << sum << endl;
		}*/
		return;
	}

	float sppider_score_rec = 0, sppider_score_lig = 0; 
	vector<int> aalcontactv, aarcontactv;
	for(unsigned short i = 0; i < ligand->c->num_aminoacids; i++)
		if(aalcontact[i]){
			aalcontactv.push_back(i);
			sppider_score_lig += ligand->c->aminoacid[i]->pInterface;
		}
	for(unsigned short i = 0; i < receptor->c->num_aminoacids; i++)
		if(aarcontact[i]){
			aarcontactv.push_back(i);
			sppider_score_rec += receptor->c->aminoacid[i]->pInterface;
			for(vector<int>::iterator itr = aalcontactv.begin(); itr != aalcontactv.end(); itr++){
				int lindex = *itr;
				if(aacontact_core[lindex][i] || aacontact_rim[lindex][i])	aminoacid_contacts.push_back(lindex*MAX_ATOMS + i);
				//if(aacontact_core[lindex][i] && aacontact_rim[lindex][i])
				//	*out << "double count " << lindex << ":" << i << endl; 
			}
		}
	*out << "interface #residues receptor: " << aarcontactv.size() << " ligand: " << aalcontactv.size() << endl;
	tr->sEvolution_interface = sppider_score_rec/aarcontactv.size() + sppider_score_lig/aalcontactv.size(); 

	tr->num_contacts = aminoacid_contacts.size();
	for(vector<long>::iterator itr = aminoacid_contacts.begin(); itr != aminoacid_contacts.end(); itr++){
		long index = *itr;
		int laindex = index / MAX_ATOMS;
		int raindex = index % MAX_ATOMS;
		Aminoacid *ra = receptor->c->aminoacid[raindex], *la = ligand->c->aminoacid[laindex];

		short ratype, latype;
		// 3 body term
		ratype = RTYPE(ra->type); latype = RTYPE(la->type);
		for(vector<int>::iterator itr = aalcontactv.begin(); itr != aalcontactv.end(); itr++){
			int aa3index = (int) *itr;
			if(aa3index > laindex && ligand->c->aacontact[laindex][aa3index] && (aacontact_core[aa3index][raindex] || aacontact_rim[aa3index][raindex])){
				short aa3type = ligand->c->aminoacid[aa3index]->type;
				aa3type = RTYPE(aa3type);
				if(ratype >= 0 && latype >= 0 && aa3type >= 0){
					short type[3];
					type[0] = ratype; type[1] = latype; type[2] = aa3type;
					for(short i = 0; i < 3; i++)
						for(short j = i+1; j < 3; j++)
							if(type[j] < type[i]){
								short t = type[i];
								type[i] = type[j];
								type[j] = t;
							}
					details->threebody_contacts[type[0]][type[1]][type[2]]++;
					//cout << "ligand " << la->cindex << " " << ra->cindex << " " << ligand->c->aminoacid[aa3index]->cindex << "\t";
					//cout << type[0] << " " << type[1] << " " << type[2] << endl;
				}
			}
		}
		for(vector<int>::iterator itr = aarcontactv.begin(); itr != aarcontactv.end(); itr++){
			int aa3index = (int) *itr;
			if(aa3index > raindex && receptor->c->aacontact[raindex][aa3index] && (aacontact_core[laindex][aa3index] || aacontact_rim[laindex][aa3index])){
				short aa3type = receptor->c->aminoacid[aa3index]->type;
				aa3type = RTYPE(aa3type);
				if(ratype >= 0 && latype >= 0 && aa3type >= 0){
					short type[3];
					type[0] = ratype; type[1] = latype; type[2] = aa3type;
					for(short i = 0; i < 3; i++)
						for(short j = i+1; j < 3; j++)
							if(type[j] < type[i]){
								short t = type[i];
								type[i] = type[j];
								type[j] = t;
							}
					details->threebody_contacts[type[0]][type[1]][type[2]]++;
					//cout << "receptor " << la->cindex << " " << ra->cindex << " " << receptor->c->aminoacid[aa3index]->cindex << "\t";
					//cout << type[0] << " " << type[1] << " " << type[2] << endl;
				}
			}
		}

		ratype = ENTROPY(ra->entropy_sequence_percentile);
		latype = ENTROPY(la->entropy_sequence_percentile);
		if(ratype <= latype)
			details->conserv_contacts[ratype][latype]++;
		else
			details->conserv_contacts[latype][ratype]++;

		/*int ss1 = SSTYPE(ra->sstructure), ss2 = SSTYPE(la->sstructure);
		int crs1 = RTYPE(ra->type) *NUM_COARSE_SSTYPES + ss1;
		int crs2 = RTYPE(la->type) *NUM_COARSE_SSTYPES + ss2;
		if(ss1 != DSSP_U && ss2 != DSSP_U)
			tr->eResiduepair += coarse_potential[crs1][crs2];*/

		aacontact_core[laindex][raindex] = aacontact_rim[laindex][raindex] = false;
		aarcontact[raindex] = false; aalcontact[laindex] = false;
	}
}


bool Object::consistent_with_info(Object *receptor, Object *ligand, ModelInfo *modelinfo, Transformation *tr){
	bool consistent = true;
	unsigned short *aagcontact;
	bool interface_capridefn = true;
	bool residue_ininterface;

	for(hash_set<unsigned short, hash<short>, eqint>::iterator itr = modelinfo->rinterface.begin(); itr != modelinfo->rinterface.end(); itr++){
		Aminoacid *ar = receptor->c->aminoacid[*itr];
		residue_ininterface = false;

		if(interface_capridefn){
			for(hash_map<const char*, Atom*, hash<const char*>, eqstr>::iterator raitr = ar->atom.begin(); raitr != ar->atom.end(); raitr++){
				Atom *ratom = (Atom *) raitr->second;
				if(ratom->name.c_str()[0] != 'H'){
					Vector vr = tr->transform(*(ratom->position));
					for(int i=0; i< ligand->c->num_atoms; i++){
						Atom *latom = ligand->c->atom[i];
						if(latom->name.c_str()[0] != 'H'){
							if(Vector::distance_squared(vr, *(latom->position)) <= 5.0*5.0){//CAPRI_INTERFACE_CUTOFF*CAPRI_INTERFACE_CUTOFF){
								residue_ininterface = true;
							}
						}
					}
				}
			}
		} else if(ar->centroid != NULL){
			Vector vr = tr->transform(*(ar->centroid));
			Vector vv = (vr - *(ligand->grid_origin)) - Vector(SS_CUTOFF/AA_GRID_SPACING,SS_CUTOFF/AA_GRID_SPACING,SS_CUTOFF/AA_GRID_SPACING);
			int x,y,z;
			x = (int) (vv.x/AA_GRID_SPACING);
			y = (int) (vv.y/AA_GRID_SPACING);
			z = (int) (vv.z/AA_GRID_SPACING);
			if(x >= 0 && x < ligand->aagrid_num_xdivisions && y >= 0 && y < ligand->aagrid_num_ydivisions && z >= 0 && z < ligand->aagrid_num_zdivisions){
				unsigned int index = (x*ligand->aagrid_num_ydivisions + y)*ligand->aagrid_num_zdivisions + z;
				//*out << laindex << " rgrid " << index << " " << ligand->aagrid[index] << endl;
				if(ligand->aagrid[index] != NULL){
					aagcontact = ligand->aagrid[index];
					for(int agi = 0; agi < ligand->aagrid_size[index]; agi++){
						int laindex = aagcontact[agi];
						Aminoacid *la = ligand->c->aminoacid[laindex];
						if(la->centroid != NULL && Vector::distance_squared(vr,*(la->centroid)) < SS_CUTOFF*SS_CUTOFF){
							residue_ininterface = true; 
						}
					}
				}
			}
		}
		consistent &= residue_ininterface;
	}

	for(hash_set<unsigned short, hash<short>, eqint>::iterator itr = modelinfo->linterface.begin(); itr != modelinfo->linterface.end(); itr++){	
		Aminoacid *al = ligand->c->aminoacid[*itr];
		residue_ininterface = false;

		if(interface_capridefn){
			for(hash_map<const char*, Atom*, hash<const char*>, eqstr>::iterator laitr = al->atom.begin(); laitr != al->atom.end(); laitr++){
				Atom *latom = (Atom *) laitr->second;
				if(latom->name.c_str()[0] != 'H'){
					Vector vl = tr->inverse_transform(*(latom->position));
					for(int i=0; i< receptor->c->num_atoms; i++){
						Atom *ratom = receptor->c->atom[i];
						if(ratom->name.c_str()[0] != 'H'){
							if(Vector::distance_squared(vl, *(ratom->position)) <= 5.0*5.0){//CAPRI_INTERFACE_CUTOFF*CAPRI_INTERFACE_CUTOFF){
								residue_ininterface = true;
							}
						}
					}
				}
			}
		} else if(al->centroid != NULL){
			Vector vl = tr->inverse_transform(*(al->centroid));
			Vector vv = (vl - *(receptor->grid_origin)) - Vector(SS_CUTOFF/AA_GRID_SPACING,SS_CUTOFF/AA_GRID_SPACING,SS_CUTOFF/AA_GRID_SPACING);
			int x,y,z;
			x = (int) (vv.x/AA_GRID_SPACING);
			y = (int) (vv.y/AA_GRID_SPACING);
			z = (int) (vv.z/AA_GRID_SPACING);
			//*out << "check " << laindex << " (" << x << "," << y << "," << z << ")\n";
			if(x >= 0 && x < receptor->aagrid_num_xdivisions && y >= 0 && y < receptor->aagrid_num_ydivisions && z >= 0 && z < receptor->aagrid_num_zdivisions){
				unsigned int index = (x*receptor->aagrid_num_ydivisions + y)*receptor->aagrid_num_zdivisions + z;
				//*out << laindex << " rgrid " << index << " " << receptor->aagrid[index] << endl;
				if(receptor->aagrid[index] != NULL){
					aagcontact = receptor->aagrid[index];
					for(int agi = 0; agi < receptor->aagrid_size[index]; agi++){
						int raindex = aagcontact[agi];
						Aminoacid *ra = receptor->c->aminoacid[raindex];
						if(ra->centroid != NULL && Vector::distance_squared(vl,*(ra->centroid)) < SS_CUTOFF*SS_CUTOFF){
							residue_ininterface = true; 
						}
					}
				}
			}
		}
		consistent &= residue_ininterface;
	}

	return consistent;
}

bool Object::consistent_with_info_T50(Object *receptor, Object *ligand, ModelInfo *modelinfo, Transformation *tr){
	int num_rresidues_consistent=0, num_lresidues_consistent=0;
	unsigned short *aagcontact;
	bool interface_capridefn = true;
	bool residue_ininterface;

	for(hash_set<unsigned short, hash<short>, eqint>::iterator itr = modelinfo->rinterface.begin(); itr != modelinfo->rinterface.end(); itr++){
		Aminoacid *ar = receptor->c->aminoacid[*itr];
		residue_ininterface = false;

		if(interface_capridefn){
			for(hash_map<const char*, Atom*, hash<const char*>, eqstr>::iterator raitr = ar->atom.begin(); raitr != ar->atom.end(); raitr++){
				Atom *ratom = (Atom *) raitr->second;
				if(ratom->name.c_str()[0] != 'H'){
					Vector vr = tr->transform(*(ratom->position));
					for(int i=0; i< ligand->c->num_atoms; i++){
						Atom *latom = ligand->c->atom[i];
						if(latom->name.c_str()[0] != 'H'){
							if(Vector::distance_squared(vr, *(latom->position)) <= 5.0*5.0){//CAPRI_INTERFACE_CUTOFF*CAPRI_INTERFACE_CUTOFF){
								residue_ininterface = true;
							}
						}
					}
				}
			}
		} else if(ar->centroid != NULL){
			Vector vr = tr->transform(*(ar->centroid));
			Vector vv = (vr - *(ligand->grid_origin)) - Vector(SS_CUTOFF/AA_GRID_SPACING,SS_CUTOFF/AA_GRID_SPACING,SS_CUTOFF/AA_GRID_SPACING);
			int x,y,z;
			x = (int) (vv.x/AA_GRID_SPACING);
			y = (int) (vv.y/AA_GRID_SPACING);
			z = (int) (vv.z/AA_GRID_SPACING);
			if(x >= 0 && x < ligand->aagrid_num_xdivisions && y >= 0 && y < ligand->aagrid_num_ydivisions && z >= 0 && z < ligand->aagrid_num_zdivisions){
				unsigned int index = (x*ligand->aagrid_num_ydivisions + y)*ligand->aagrid_num_zdivisions + z;
				//*out << laindex << " rgrid " << index << " " << ligand->aagrid[index] << endl;
				if(ligand->aagrid[index] != NULL){
					aagcontact = ligand->aagrid[index];
					for(int agi = 0; agi < ligand->aagrid_size[index]; agi++){
						int laindex = aagcontact[agi];
						Aminoacid *la = ligand->c->aminoacid[laindex];
						if(la->centroid != NULL && Vector::distance_squared(vr,*(la->centroid)) < SS_CUTOFF*SS_CUTOFF){
							residue_ininterface = true;
						}
					}
				}
			}
		}
		if(residue_ininterface) num_rresidues_consistent++;
	}

	for(hash_set<unsigned short, hash<short>, eqint>::iterator itr = modelinfo->linterface.begin(); itr != modelinfo->linterface.end(); itr++){
		Aminoacid *al = ligand->c->aminoacid[*itr];
		residue_ininterface = false;

		if(interface_capridefn){
			for(hash_map<const char*, Atom*, hash<const char*>, eqstr>::iterator laitr = al->atom.begin(); laitr != al->atom.end(); laitr++){
				Atom *latom = (Atom *) laitr->second;
				if(latom->name.c_str()[0] != 'H'){
					Vector vl = tr->inverse_transform(*(latom->position));
					for(int i=0; i< receptor->c->num_atoms; i++){
						Atom *ratom = receptor->c->atom[i];
						if(ratom->name.c_str()[0] != 'H'){
							if(Vector::distance_squared(vl, *(ratom->position)) <= 5.0*5.0){//CAPRI_INTERFACE_CUTOFF*CAPRI_INTERFACE_CUTOFF){
								residue_ininterface = true;
							}
						}
					}
				}
			}
		} else if(al->centroid != NULL){
			Vector vl = tr->inverse_transform(*(al->centroid));
			Vector vv = (vl - *(receptor->grid_origin)) - Vector(SS_CUTOFF/AA_GRID_SPACING,SS_CUTOFF/AA_GRID_SPACING,SS_CUTOFF/AA_GRID_SPACING);
			int x,y,z;
			x = (int) (vv.x/AA_GRID_SPACING);
			y = (int) (vv.y/AA_GRID_SPACING);
			z = (int) (vv.z/AA_GRID_SPACING);
			//*out << "check " << laindex << " (" << x << "," << y << "," << z << ")\n";
			if(x >= 0 && x < receptor->aagrid_num_xdivisions && y >= 0 && y < receptor->aagrid_num_ydivisions && z >= 0 && z < receptor->aagrid_num_zdivisions){
				unsigned int index = (x*receptor->aagrid_num_ydivisions + y)*receptor->aagrid_num_zdivisions + z;
				//*out << laindex << " rgrid " << index << " " << receptor->aagrid[index] << endl;
				if(receptor->aagrid[index] != NULL){
					aagcontact = receptor->aagrid[index];
					for(int agi = 0; agi < receptor->aagrid_size[index]; agi++){
						int raindex = aagcontact[agi];
						Aminoacid *ra = receptor->c->aminoacid[raindex];
						if(ra->centroid != NULL && Vector::distance_squared(vl,*(ra->centroid)) < SS_CUTOFF*SS_CUTOFF){
							residue_ininterface = true;
						}
					}
				}
			}
		}
		if(residue_ininterface)	num_lresidues_consistent++;
	}

	float fr = (num_rresidues_consistent + 1.0) / ( modelinfo->rinterface.size() + 1.0);
	bool consistent = false;
	if(fr>=0.6){
		float fl = (num_lresidues_consistent + 1.0) / ( modelinfo->linterface.size() + 1.0);
		if(fl>=0.6){
			consistent = true;
			*out << "consistency " << tr->frame_number << " " << fr << " " << fl << endl;
		}
	}
	return consistent;
}


void Object::compute_tobibahar_energy(Object *receptor, Object *ligand, Transformation *tr, ProtProtDetailedScoringMetrics *details){
	int num_residue_types = NUM_RESIDUE_TYPES;
	for(int i = 0; i < num_residue_types; i++)
		for(int j = 0 ; j < num_residue_types; j++)
			details->residue_contacts_core[i][j] = details->residue_contacts_rim[i][j] = 0;

	unsigned short *aagcontact;
	for(int laindex = 0; laindex < ligand->c->num_aminoacids; laindex++){
		Aminoacid *al = ligand->c->aminoacid[laindex];

		// the aminoacid grid was computed on the basis of the side chain center of mass, can be used to compute side chain contacts 
		// and contacts between CO and NH of ligand and the side chain of receptor 
		Vector *vaa[3];
		vaa[0] = al->centroid;
		vaa[1] = (al->amide_nitrogen == NULL) ? NULL : al->amide_nitrogen->position;
		vaa[2] = (al->carbonyl_oxygen == NULL) ? NULL : al->carbonyl_oxygen->position;
		for(int aapi = 0; aapi < 3; aapi++) 
			if(vaa[aapi]  != NULL){
				Vector vl = tr->inverse_transform(*(vaa[aapi]));
				for(int raindex = 0; raindex < receptor->c->num_aminoacids; raindex++){
					Aminoacid *ra = receptor->c->aminoacid[raindex]; 
					if(aapi == 0) {
						if(ra->centroid != NULL && Vector::distance_squared(vl,*(ra->centroid)) < 6.8*6.8){
							//*out << "gridc\t" << laindex << " " << raindex << " " << Vector::distance(vl,*(receptor->c->aminoacid[raindex]->centroid)) << endl;
							short ratype=ra->type, latype=al->type;
							if(ratype >= 0 && latype>= 0)
								if(ratype <= latype)
									details->residue_contacts_core[ratype][latype]++;
								else
									details->residue_contacts_core[latype][ratype]++;
						}
						if(ra->amide_nitrogen != NULL && Vector::distance_squared(vl,*(ra->amide_nitrogen->position)) < 5.6*5.6)
							details->residue_contacts_core[al->type][NTER]++;
						if(ra->carbonyl_oxygen != NULL && Vector::distance_squared(vl,*(ra->carbonyl_oxygen->position)) < 5.6*5.6)
							details->residue_contacts_core[al->type][CTER]++;
					} else if(aapi == 1) {
						if(ra->centroid != NULL && Vector::distance_squared(vl,*(ra->centroid)) < 5.6*5.6)
							details->residue_contacts_core[ra->type][NTER]++;
						if(ra->amide_nitrogen != NULL && Vector::distance_squared(vl,*(ra->amide_nitrogen->position)) < 16.0)
							details->residue_contacts_core[NTER][NTER]++;
						if(ra->carbonyl_oxygen != NULL && Vector::distance_squared(vl,*(ra->carbonyl_oxygen->position)) < 16.0){
							details->residue_contacts_core[CTER][NTER]++;
							details->residue_contacts_core[NTER][CTER]++;
						}
					} else if(aapi == 2) {
						if(ra->centroid != NULL && Vector::distance_squared(vl,*(ra->centroid)) < 5.6*5.6)
							details->residue_contacts_core[ra->type][CTER]++;
						if(ra->amide_nitrogen != NULL && Vector::distance_squared(vl,*(ra->amide_nitrogen->position)) < 16.0){
							details->residue_contacts_core[CTER][NTER]++;
							details->residue_contacts_core[NTER][CTER]++;
						} if(ra->carbonyl_oxygen != NULL && Vector::distance_squared(vl,*(ra->carbonyl_oxygen->position)) < 16.0)
							details->residue_contacts_core[CTER][CTER]++;
					}
				}
			}
	}
}

void Object::compute_luskolnick_energy(Object *receptor, Object *ligand, Transformation *tr, ProtProtDetailedScoringMetrics *details){
	int num_residue_types = NUM_RESIDUE_TYPES;
	for(int i = 0; i < num_residue_types; i++)
		for(int j = 0 ; j < num_residue_types; j++)
			details->residue_contacts_core[i][j] = details->residue_contacts_rim[i][j] = 0;

	bool **aacontact = aacontact_core;
	for(int laindex = 0; laindex < ligand->c->num_aminoacids; laindex++)
		for(int raindex = 0; raindex < receptor->c->num_aminoacids; raindex++)
			aacontact[laindex][raindex] = false;

	for(int laindex = 0; laindex < ligand->c->num_aminoacids; laindex++){
		Aminoacid *la = ligand->c->aminoacid[laindex];
		for(hash_map<const char*, Atom*, hash<const char*>, eqstr>::iterator laitr = la->atom.begin(); laitr != la->atom.end(); laitr++){
			Atom *latom = (Atom *) laitr->second;
			if(latom->name.c_str()[0] != 'H'){
				Vector v = tr->inverse_transform(*(latom->position)) - *(receptor->grid_origin);
				int cellx = (int) (v.x/GRID_SPACING);
				int celly = (int) (v.y/GRID_SPACING);
				int cellz = (int) (v.z/GRID_SPACING);

				if(cellx >= 0 && cellx < receptor->grid_num_xdivisions && celly >= 0 && celly < receptor->grid_num_ydivisions && cellz >= 0 && cellz < receptor->grid_num_zdivisions){
					unsigned int index = (ATOM_NEIGHBOR_GS_BY_GS*(cellx/ATOM_NEIGHBOR_GS_BY_GS)*receptor->grid_num_ydivisions + ATOM_NEIGHBOR_GS_BY_GS*(celly/ATOM_NEIGHBOR_GS_BY_GS))*receptor->grid_num_zdivisions + ATOM_NEIGHBOR_GS_BY_GS*(cellz/ATOM_NEIGHBOR_GS_BY_GS);

					if(receptor->grid[index] != NULL){
						GridCell *coarse_cell = receptor->grid[index]; 
						unsigned short* atom_neighbors[2];
						atom_neighbors[0] = coarse_cell->atom_overlap_neighbors;
						atom_neighbors[1] = coarse_cell->atom_contact_neighbors;
						int num_atom_neighbors[2];
						num_atom_neighbors[0] = coarse_cell->num_atom_overlap_neighbors;
						num_atom_neighbors[1] = coarse_cell->num_atom_contact_neighbors;
						for(unsigned short anli = 0; anli < 2; anli++)
							if(atom_neighbors[anli] != NULL)
								for(int ani = 0; ani < num_atom_neighbors[anli]; ani++){
									int ratom_cindex = atom_neighbors[anli][ani];
									Atom *ratom = receptor->c->atom[ratom_cindex];
									if(ratom->name.c_str()[0] != 'H' && !aacontact[laindex][ratom->monocindex]){
										double d2 = Vector::distance_squared(tr->inverse_transform(*(latom->position)),ratom->position);
										if(d2 <= 4.50*4.50){
											Aminoacid *ra = receptor->c->aminoacid[ratom->monocindex]; 
											aacontact[laindex][ra->cindex]=true;
											short ratype=ra->type, latype=la->type;
											if(ratype >= 0 && latype>= 0)
												if(ratype <= latype)
													details->residue_contacts_core[ratype][latype]++;
												else
													details->residue_contacts_core[latype][ratype]++;			
										}
									}
								}	 		
					}
				}
			}
		}
	}
}

void Object::compute_dpp_energy(Object *receptor, Object *ligand, Transformation *tr){
	hash_map<long,unsigned short,hash<long>,eqlong> pair_d[4];
	compute_dpp_contacts(receptor,ligand,tr,pair_d);

	tr->eResiduepair = 0;
	for(int i=0; i<4; i++){
		for(hash_map<long,unsigned short,hash<long>,eqlong>::iterator itr = pair_d[i].begin(); itr != pair_d[i].end(); itr++){
			long index = itr->first;
			unsigned short d = itr->second;
			int laindex = index / MAX_ATOMS;
			int raindex = index % MAX_ATOMS;
			Aminoacid *ra = receptor->c->aminoacid[raindex], *la = ligand->c->aminoacid[laindex];
			short ratype=ra->type,latype=la->type;
			if(ratype >= 0 && latype >= 0 && ratype < 20 && latype < 20){
				float delta;
				switch(i){
				case 0:
					if(ratype <= latype)	delta = cacar_potential[ratype][latype][d];
					else	delta = cacar_potential[latype][ratype][d];
					break;
				case 3:
					if(ratype <= latype)	delta = cmcmr_potential[ratype][latype][d];
					else	delta = cmcmr_potential[latype][ratype][d];
					break;
				case 1:
					delta = cmcar_potential[ratype][latype][d];
					break;
				case 2:
					delta = cmcar_potential[latype][ratype][d];
					break;
				}
				//*out << latype << " " << ratype << " " << d << " " << delta << endl;
				tr->eResiduepair += delta;
			}
		}
	}
}

void Object::compute_dpp_contacts(Object *receptor, Object *ligand, Transformation *tr, hash_map<long,unsigned short,hash<long>,eqlong> *pair_d){
	unsigned short *aagcontact;
	for(int laindex = 0; laindex < ligand->c->num_aminoacids; laindex++){
		Aminoacid *al = ligand->c->aminoacid[laindex];
		if(al->type != -1){
			Vector v[2];
			bool vnotnull[2];
			for(int i = 0; i <2; i++)	vnotnull[i] = false;
			float mass = 0;
			Vector sc_cm(0,0,0);
			for(hash_map<const char*, Atom*, hash<const char*>, eqstr>::iterator itr = al->atom.begin(); itr != al->atom.end(); itr++){
				Atom *a = itr->second;
				if(a->name == "CA")	{
					v[0] = tr->inverse_transform(*(a->position));
					vnotnull[0] = true;
				} else if(a->name != "N" && a->name != "C" && a->name != "O" && a->name.c_str()[0] != 'H') {  
					sc_cm = sc_cm + *(a->position) * a->mass;
					mass += a->mass;
				}
			}
			if(mass != 0){
				v[1] = tr->inverse_transform((sc_cm*(1.0/mass)));
				vnotnull[1] = true;
			}

			for(int i = 0 ; i < 2 ; i++){
				bool examinedaa[receptor->c->num_aminoacids];
				for(int j=0; j< receptor->c->num_aminoacids; j++)	examinedaa[j]=false;
				if(vnotnull[i]){
					Vector vv = (v[i] - *(receptor->grid_origin)) - Vector(SS_CUTOFF/AA_GRID_SPACING,SS_CUTOFF/AA_GRID_SPACING,SS_CUTOFF/AA_GRID_SPACING);
					int vvx,vvy,vvz;
					vvx = (int) (vv.x/AA_GRID_SPACING);
					vvy = (int) (vv.y/AA_GRID_SPACING);
					vvz = (int) (vv.z/AA_GRID_SPACING);
					for(int xi = -1; xi < 1; xi++)
						for(int yi = -1; yi < 1; yi++)
							for(int zi = -1; zi < 1; zi++){
								int x,y,z;
								x = vvx+xi; y = vvy+yi; z = vvz+zi;
								if(x >= 0 && x < receptor->aagrid_num_xdivisions && y >= 0 && y < receptor->aagrid_num_ydivisions && z >= 0 && z < receptor->aagrid_num_zdivisions){
									unsigned int index = (x*receptor->aagrid_num_ydivisions + y)*receptor->aagrid_num_zdivisions + z;
									if(receptor->aagrid[index] != NULL){
										aagcontact = receptor->aagrid[index];
										for(int agi = 0; agi < receptor->aagrid_size[index]; agi++){
											int raindex = aagcontact[agi];
											if(!examinedaa[raindex]){
												examinedaa[raindex] = true;
												Aminoacid *ar = receptor->c->aminoacid[raindex];
												if(ar->type != -1){
													Vector w[2];
													bool wnotnull[2];
													for(int j=0; j<2; j++)	wnotnull[j] = false;
													float mass = 0;
													Vector rsc_cm(0,0,0);
													for(hash_map<const char*, Atom*, hash<const char*>, eqstr>::iterator itr = ar->atom.begin(); itr != ar->atom.end(); itr++){
														Atom *a = itr->second;
														if(a->name == "CA")	{
															w[0] = *(a->position);
															wnotnull[0] = true;
														} else if(a->name != "N" && a->name != "C" && a->name != "O" && a->name.c_str()[0] != 'H') {  
															rsc_cm = rsc_cm + *(a->position) * a->mass;
															mass += a->mass;
														}
													}
													if(mass != 0){
														w[1] = (rsc_cm*(1.0/mass));
														wnotnull[1] = true;
													}
													for(int j=0; j<2; j++)
														if(wnotnull[j]){
															float d = Vector::distance(v[i],w[j]);
															int discrete_d = 10*(d-0.05);
															if(discrete_d < 0)	discrete_d = 0;
															long index = laindex*MAX_ATOMS + raindex;
															if(discrete_d < 135){
																if(i==0 &&  j==0)	pair_d[0][index] = (unsigned short) discrete_d;
																else if(i==1 &&  j==1)	pair_d[3][index] = (unsigned short) discrete_d;
																else if(i==1 &&  j==0)	pair_d[2][index] = (unsigned short) discrete_d;
																else if(i==0 &&  j==1)	pair_d[1][index] = (unsigned short) discrete_d;
															}
														}
												}
											}
										}
									}
								}
							}
				}
			}
		}
	}
}


#ifndef DSSP_EXECUTABLE
#define DSSP_EXECUTABLE "external/dsspcmbi"
#endif
bool sasa_dssp_computed_unbound=false;
hash_map<const char*, float, hash<const char*>,eqstr> rsasa, lsasa;
/*
 * Use dssp to compute the change in surface area
 * Assume all chains are proteins no nucleic acids
 */
void Object::compute_deltasas_dssp(Object *receptor, Object *ligand, Transformation *tr, ProtProtDetailedScoringMetrics *details){
	stringstream ss (stringstream::in | stringstream::out);
	string filename;
	fstream dsspin;
	char buf[8192],command[1024];
	int ret;

	if(!sasa_dssp_computed_unbound){
		ss << getenv(string("SCRATCH").c_str()) << "/receptor" << tr->frame_number << ".pdb";
		ss >> filename;
		receptor->c->write_as_pdb(filename);
		ss.clear();
		ss << getenv(string("SCRATCH").c_str()) << "/receptor" << tr->frame_number << ".dssp";
		string rdsspoutfile; ss >> rdsspoutfile;
		sprintf(command,"%s/%s %s %s",piedock_home.c_str(),string(DSSP_EXECUTABLE).c_str(), filename.c_str(), rdsspoutfile.c_str());
		*out << command << " "; out->flush();
		ret = system(command);
		*out << ret << endl; out->flush();

		ss.clear();
		ss << getenv(string("SCRATCH").c_str()) << "/ligand" << tr->frame_number << ".pdb";
		ss >> filename;
		ligand->c->write_as_pdb(filename);
		ss.clear();
		ss << getenv(string("SCRATCH").c_str()) << "/ligand" << tr->frame_number << ".dssp";
		string ldsspoutfile; ss >> ldsspoutfile;
		sprintf(command,"%s/%s %s %s",piedock_home.c_str(),string(DSSP_EXECUTABLE).c_str(), filename.c_str(), ldsspoutfile.c_str());
		*out << command << " "; out->flush();
		ret = system(command);
		*out << ret << endl; out->flush();

		dsspin.open(rdsspoutfile.c_str(), fstream::in);
		do {
			dsspin.getline(buf,8192);
		} while (((string(buf)).find("STRUCTURE") == string::npos) && !dsspin.eof());
		while(!dsspin.eof()){
			dsspin.getline(buf,8192);
			char chain = buf[11];
			if(receptor->c->chains == "-")	chain = '-';
			if((receptor->c->chains.find(chain) != string::npos) && (buf[13] != '!')){
				string aacid_index;
				{
					stringstream line(buf+5,stringstream::in);
					line >> aacid_index;
				}
				if(buf[10] != ' '){
					aacid_index = aacid_index.substr(0,aacid_index.find(buf[10])+1);
				}
				float sasa;
				{
					stringstream line(buf+34,stringstream::in);
					line >> sasa;
				}
				ss.clear();
				ss << chain << aacid_index;
				string caaindex; ss >> caaindex;
				caaindex = *(new string(caaindex.c_str()));
				rsasa[caaindex.c_str()] = sasa;
				*out << caaindex << " " << sasa << endl;
			}
		}
		dsspin.close();
		dsspin.clear();

		dsspin.open(ldsspoutfile.c_str(), fstream::in);
		do {
			dsspin.getline(buf,8192);
		} while (((string(buf)).find("STRUCTURE") == string::npos) && !dsspin.eof());
		while(!dsspin.eof()){
			dsspin.getline(buf,8192);
			char chain = buf[11];
			if(ligand->c->chains == "-")	chain = '-';
			if((ligand->c->chains.find(chain) != string::npos) && (buf[13] != '!')){
				string aacid_index;
				{
					stringstream line(buf+5,stringstream::in);
					line >> aacid_index;
				}
				if(buf[10] != ' '){
					aacid_index = aacid_index.substr(0,aacid_index.find(buf[10])+1);
				}
				float sasa;
				{
					stringstream line(buf+34,stringstream::in);
					line >> sasa;
				}
				ss.clear();
				ss << chain << aacid_index;
				string caaindex; ss >> caaindex;
				caaindex = *(new string(caaindex.c_str()));
				lsasa[caaindex.c_str()] = sasa;
				*out << caaindex << " " << sasa << endl;
			}
		}
		dsspin.close();
		dsspin.clear();

		*out << rsasa.size() << " " << lsasa.size() << endl;
		sasa_dssp_computed_unbound=true;
	}

	ss.clear();
	ss << getenv(string("SCRATCH").c_str()) << "/" << "complex" << tr->frame_number << ".pdb";
	ss >> filename;
	write_as_pdb(receptor->c,ligand->c,tr,filename);
	ss.clear();
	ss << getenv(string("SCRATCH").c_str()) << "/" << "complex" << tr->frame_number << ".dssp";
	string cdsspoutfile; ss >> cdsspoutfile;
	sprintf(command,"%s/%s %s %s",piedock_home.c_str(),string(DSSP_EXECUTABLE).c_str(), filename.c_str(), cdsspoutfile.c_str());
	*out << command << " "; out->flush();
	ret = system(command);
	*out << ret << endl; out->flush();

	for(int i = 0 ; i < NUM_RESIDUE_TYPES; i++)	details->delta_sasa[i] = 0;
	hash_map<const char*, float, hash<const char*>,eqstr> csasa;
	dsspin.open(cdsspoutfile.c_str(), fstream::in);
	do {
		dsspin.getline(buf,8192);
	} while (((string(buf)).find("STRUCTURE") == string::npos) && !dsspin.eof());
	while(!dsspin.eof()){
		dsspin.getline(buf,8192);
		char fchain = buf[11];
		if(!dsspin.eof() && buf[13] != '!'){
			char chain;
			string aacid_index;
			{
				stringstream line(buf+5,stringstream::in);
				line >> aacid_index;
			}
			if(buf[10] != ' '){
				aacid_index = aacid_index.substr(0,aacid_index.find(buf[10])+1);
			}
			float sasa;
			{
				stringstream line(buf+34,stringstream::in);
				line >> sasa;
			}
			ss.clear();

			Complex *c;
			char cchainstart;
			if(fchain < 'A'+receptor->c->chains.size()){
				cchainstart='A';
				c = receptor->c;
			}else{
				cchainstart='A'+receptor->c->chains.size();
				c = ligand->c;
			}

			Molecule *m;
			if(c->chains == "-"){
				chain = '-';
				m = (Molecule *) (c->molecules.begin())->second;
			} else {
				unsigned short moleculeno=0;
				for(hash_map<char, Molecule*, hash<char>, eqint>::iterator mitr = c->molecules.begin(); mitr != c->molecules.end(); mitr++){
					Molecule *mc = mitr->second;
					if(fchain == cchainstart + moleculeno){
						m = mc;
						chain=m->chain;
					}
					moleculeno++;
				}
			}

			ss << chain << aacid_index;
			string caaindex; ss >> caaindex;
			float delta_sasa;
			if(fchain < 'A'+receptor->c->chains.size()){
				delta_sasa = rsasa[caaindex.c_str()] - sasa;
				// sometimes delta_sasa is negative
				//*out << caaindex << " " << rsasa[caaindex.c_str()] << " " << sasa << " " << delta_sasa << endl;
			} else {
				delta_sasa = lsasa[caaindex.c_str()] - sasa;
				//*out << caaindex << " " << lsasa[caaindex.c_str()] << " " << sasa << " " << delta_sasa << endl;
			}
			//*out << aacid_index << " " << m->chain << " " << m->aminoacid.count(aacid_index.c_str()) << endl;
			if(((Protein*) m)->aminoacid.count(aacid_index.c_str()) > 0){
				Aminoacid *aa = ((Protein*) m)->aminoacid[aacid_index.c_str()];
				if(aa->type >= 0)
					details->delta_sasa[aa->type] += delta_sasa;
			}
		}
	}
	dsspin.close();

	tr->delta_sasa = 0;
	for(int i = 0 ; i < NUM_RESIDUE_TYPES; i++){
		tr->delta_sasa += details->delta_sasa[i];
	}
	*out << tr->frame_number << " sasa " << tr->delta_sasa << endl;
}

void read_usphere_solution(){    
	string filename = piedock_home + "/" + string(DOCK_DIR) + "usphere.10";
	fstream fusphere(filename.c_str());

	int n = 2*UNIT_GRID_PHI_DIVISIONS*UNIT_GRID_PHI_DIVISIONS;
	char buf[8192];
	while(fusphere.good()){
		int index;
		fusphere.getline(buf,8192);
		if(fusphere.gcount() > 0){
			stringstream ss (stringstream::in | stringstream::out);
			ss << buf;
			ss >> index;
			for(int i = 0 ;i < n; i++){
				ss >> sorted_cell[index][i];
			}
		}
		fusphere.getline(buf,8192);
		if(fusphere.gcount() > 0){
			stringstream ss (stringstream::in | stringstream::out);
			ss << buf;
			for(int i = 0 ;i < n; i++){
				ss >> distance_bound[index][i];
			}
		}
	}
	fusphere.close();
}
