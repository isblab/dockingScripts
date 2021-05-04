#ifndef INCLUDED_MOLECULE
	#define INCLUDED_MOLECULE
	#include "molecule.hpp"
	#include "evolution.hpp"
#endif

#ifndef INCLUDED_ERRNO
#define INCLUDED_ERRNO
#include <errno.h>
#endif

#include <pthread.h>

#define NUM_PIE_PARAMS 254

#define MAX_BINS 4096
#define MIN_COORDINATE -256
#define MAX_COORDINATE 256
#define NUM_ADJACENT_CELLS_EXAMINED 2
#define INCREMENT 1.2
#define MAX_MATCH_DISTANCE INCREMENT
#define MIN_NORMAL_DOT 0//0.3
/* the bins are distributed from -32 to 31.875*/
#define MAX_POINTS 65536

#define ATOM_NEIGHBOR_GS_BY_GS 2
#define AA_GRID_SPACING (3*GRID_SPACING)

#define INTERIOR_ATOM_SPREAD 1
//#define CLOSEST_POINT_CUTOFF 3.0
#define CONTACT_MIN_ANGLE 0.5
#define CONTACT_CUTOFF 0

#define MIN_VOTES( n1 , n2 ) (min(6, (int) (sqrt(sqrt(n1*n2))*0.4)))

#define RIGID_MAX_ATOMS_INTERIOR 2
#define RIGID_MAX_ATOMS_INTERMEDIATE 6
#define RIGID_MAX_AREA_INTERIOR (3*RIGID_MAX_ATOMS_INTERIOR)
#define RIGID_MAX_AREA_INTERMEDIATE (3*RIGID_MAX_ATOMS_INTERMEDIATE)

// allow upto 3 misplaced side chains 
#define FLEX_MAX_ATOMS_INTERIOR 6
#define FLEX_MAX_ATOMS_INTERMEDIATE 24
#define FLEX_MAX_AREA_INTERIOR (3*FLEX_MAX_ATOMS_INTERIOR)
#define FLEX_MAX_AREA_INTERMEDIATE (3*FLEX_MAX_ATOMS_INTERMEDIATE)

// (-inf,-3.5),[-3.5,-2.0),[-2.0,-0.5),[-0.5,0.5) [0.5,inf) 
// d is the minimum distance to a surface cell, it is negated for cells inside the protein
// surface is the SAS, no probe radius addition
#define NUM_DTRANSFORM_DIVISIONS 5
#define DTRANSFORM_INDEX(d) ((d > 0.5) ? 0 : ((d > -0.5) ? 1 : ((d > -2.0) ? 2 : ((d > -3.5) ? 3 : 4))))
#define MIN_VDW_FILTER1 30 //20
#define MIN_VDW_FILTER2 50//30
#define MIN_DELTA_SASA 500.0
// vdwrepulsion 10 means that an atom of a residue is entirely buried in the other molecule
#define MAX_VDW_REPULSION 40 

#define UNIT_GRID_DIVISIONS 10
#define UNIT_GRID_PHI_DIVISIONS UNIT_GRID_DIVISIONS
#define UNIT_GRID_THETA_DIVISIONS (2*UNIT_GRID_DIVISIONS)
#define UNIT_GRID_THETA_SPACING (PI/UNIT_GRID_DIVISIONS)
#define UNIT_GRID_PHI_SPACING (PI/UNIT_GRID_DIVISIONS)

#define MIN_RESIDUE_CONTACTS(n) 15 //(20 * n)
#define MAX_CONTACT_ENERGY(n) (20.0  - 2.0 * n)
#define MIN_NEW_CONTACTS 20

#define ATOM_DISTANCE_DIVISION(d2) ((d2 > 42.25) ? -1 : ((d2 > 25) ? 2 : ((d2 > 12.25) ? 1 : ((d2>4.0) ? 0 : -1))))

//#define SS_CUTOFF_COMPARE_NATIVE 8.0
#define CAPRI_INTERFACE_CUTOFF 10.0

#define SCALE_MIN 0.99
#define SCALE_MAX 1.01

#define NUM_CHARGE_TYPES 3
#define NEGATIVE 0
#define NEUTRAL 1
#define POSITIVE 2
#define MAX_NEG_CHARGE -0.3
#define MIN_POS_CHARGE 0.3

#define RIGID_DOCK 0
#define FLEX_SIDECHAINS 1
#define FLEX_BACKBONE 2

#define RUN_EWALD "scripts/run_ewald.sh"

#define SCWRL_TIME_LIMIT 300
#define SCWRL_EXECUTABLE "external/scwrl/scwrl3"//Scwrl4"
#define SCWRL_INPUT_FILE "scwrlinp"
#define SCWRL_OUTPUT_FILE "scwrlout"
#define SCWRL_LOG_FILE "scwrllog"

#define TR_NUM_LOCAL_NEIGHBORS 729

#ifndef CALPHA
#define CALPHA 0
#define BACKBONE 1
#define ALLATOMS 2
#endif

#define INTERMEDIATE_IMAG 3.0
#define INTERIOR_IMAG 8.0

#define NUM_ELEC_FEATURES 20
#define NUM_ATOM_VDWREPULSION_DIVISIONS 20

#define DETAILS 1
#define RESULTS 2

#define ZSCORE_N 100

#define RESIDUE_CONTACT_3B 1
#define ATOM20_CONTACT 2 
#define VDW_RES_BKBN_SPLINE 3

//#define STEP_POTENTIAL
#define LINEAR_SPLINE_POTENTIAL

extern float *vdw_energy,*residue_solvation_energy,**residue_potential, **residue_bkbn_potential, **typeconservation_potential;
extern float ***cacar_potential, ***cmcar_potential, ***cmcmr_potential, *score_param;

class Reference_Frame{
	public:
		Vector *translation;
		Vector *ex, *ey, *ez;
		//float theta, phi, psi; // euler angles
		float scale;
		long frame_number;
		
	Reference_Frame();
	Reference_Frame(Vector *t, Vector *ex, Vector *ey, float s, long frame_number);
	Reference_Frame(Reference_Frame*);
	~Reference_Frame();

	Vector rotate(Vector p);
	Vector inverse_rotate(Vector p);
	Vector transform(Vector p);
	Vector inverse_transform(Vector p);
	
	void distance(float*,float*,Reference_Frame *tr, Vector);
	void distance(float*,float*,Reference_Frame *tr);
	void distance(float*,Reference_Frame *tr);
	
	void eulerangles(float*,float*,float*);
	static Reference_Frame* compose(Reference_Frame*,Reference_Frame*);
	static Reference_Frame* invert(Reference_Frame*);
	
	void write_as_pdb(Complex *c, string chains, bool, string filename, bool);
	//void write_as_pdb(Complex *c, string chains, char, bool, string filename, bool);
};

/*
 * Pattern has sufficient information to define a reference frame
 * the frame_number of the reference frame is also the identifier of the pattern
 *
class Pattern{
	public:
		short type;
		
	Reference_Frame get_reference_frame();
	void locally_maximize();
};

class GeometricPattern : public Pattern{
};*/

class VerificationMetrics {
	public:
		float rmsd,lrmsd,irmsd,delta_r,delta_U;
		// delta_r is the difference in the center of mass of the ligand from correct tr, delta_U is difference in rotations about the center of mass of the ligand
		// the rotation matrices are the same when we move the origin to the center of mass of the ligand 
		float frac_native_contacts_predicted, frac_nonnative_contacts_predicted;
};

#define PROT_PROT_DETAILS 0
#define PROT_RNA_DETAILS 1

class DetailedScoringMetrics {
	protected:
		short type;
		
	public:
		float elec[NUM_ELEC_FEATURES];
	
		//DetailedScoringMetrics();
		virtual ~DetailedScoringMetrics() {};
};

class ProtProtDetailedScoringMetrics : public DetailedScoringMetrics {
	public:
		unsigned short points_vs_distance[NUM_DTRANSFORM_DIVISIONS],atoms_vs_distance[NUM_DTRANSFORM_DIVISIONS],
			grid_contacts[NUM_DTRANSFORM_DIVISIONS][NUM_DTRANSFORM_DIVISIONS],
			threebody_contacts[NUM_COARSE_RTYPES][NUM_COARSE_RTYPES][NUM_COARSE_RTYPES],
			atom_vdw_repulsion[NUM_ATOM_VDWREPULSION_DIVISIONS+1]
			;
		
		float delta_vdwsa[NUM_RESIDUE_TYPES], delta_sasa[NUM_RESIDUE_TYPES], //delta_vol[NUM_RESIDUE_TYPES],delta_sevol[NUM_RESIDUE_TYPES], 
			delta_vdwsa_atom[NUM_ATOM_TYPES], delta_sasa_atom[NUM_ATOM_TYPES], //delta_vol_atom[NUM_ATOM_TYPES],delta_sevol_atom[NUM_ATOM_TYPES], 
			**atom_contacts, ***atom_dcontacts;
			;
		unsigned short /**residue_contacts_core, **residue_contacts_rim,*/ **conserv_contacts;
		float **residue_contacts_core, **residue_contacts_rim;
		
		ProtProtDetailedScoringMetrics();
		~ProtProtDetailedScoringMetrics();
};

class ProtRnaDetailedScoringMetrics : public DetailedScoringMetrics {
	public:
		float coarse_contacts[NUM_RESIDUE_TYPES][NUM_RNA_PARTICLE_TYPES];
	
		ProtRnaDetailedScoringMetrics() {};
		~ProtRnaDetailedScoringMetrics() {};
};

class ScoringMetrics {
	public:
		float votes, curvature_score;
		
		float eVdw, eVdw_repulsion, eElectrostatic, sEvolution_interface, eSolvation, eResiduepair, delta_sasa, delta_vdwsa, delta_vdwsapair;
		float eT32S3;
		unsigned short num_contacts, num_neighbors, num_clashes, num_bbclashes;
		
		float freeE_rigidmotion, freeE_sidechain;
		DetailedScoringMetrics *detailed_scores;
};

#define TN_BASIC 0
#define TN_VERIFY 1
#define TN_PROTPROT_DETAILED_SCORE_VERIFY 2
#define TN_PROTRNA_DETAILED_SCORE_VERIFY 3

class Transformation : public Reference_Frame , public ScoringMetrics {
	public:
		class TransformationDetails{
			public:
				// convert to float or ...
				hash_map<int,float, hash<int>,eqint> lresidue_vdw_energy, rresidue_vdw_energy;
				hash_map<long,float,hash<long>,eqlong> contact_energy_by_residue;
					
				vector<Transformation*> nearest_neighbors;
				float num_new_contacts; // used while generating combinations
		};
		
		Vector *cmr; // use cmr for clustering since translation can be very large
		float cmr_theta, cmr_phi;
		
		TransformationDetails *details;
		VerificationMetrics *vmetrics;
		
		static const int basic_byte_size = 76,
				verify_byte_size = basic_byte_size + 28,
				protprot_detailed_scores_verify_byte_size = verify_byte_size + 12 + sizeof(unsigned short)*NUM_DTRANSFORM_DIVISIONS*(2+NUM_DTRANSFORM_DIVISIONS) + sizeof(float)*NUM_RESIDUE_TYPES*2 
					//+ sizeof(unsigned short)*NUM_RESIDUE_TYPES*NUM_RESIDUE_TYPES*2
					+ sizeof(float)*NUM_RESIDUE_TYPES*NUM_RESIDUE_TYPES*2
					+ sizeof(unsigned short)*NUM_ENTROPY_DIVISIONS*NUM_ENTROPY_DIVISIONS
					+ sizeof(float)*NUM_ELEC_FEATURES
					+ sizeof(unsigned short)*NUM_COARSE_RTYPES*NUM_COARSE_RTYPES*NUM_COARSE_RTYPES
					+ sizeof(float)*NUM_ATOM_TYPES*NUM_ATOM_TYPES*NUM_ATOM_DISTANCE_DIVISIONS
					+ sizeof(float)*NUM_ATOM_TYPES*2
					+ sizeof(unsigned short)*NUM_ATOM_VDWREPULSION_DIVISIONS
					,
				protrna_detailed_scores_verify_byte_size = verify_byte_size + 12 + sizeof(float)*NUM_RESIDUE_TYPES*NUM_RNA_PARTICLE_TYPES
					;

	Transformation(Vector *t, Vector *ex, Vector *ey, float s, float votes, long);
	Transformation(float,float,float,long);
	Transformation(char* , unsigned short );
	
	~Transformation();
	
	void compute_coordinates(void *receptor, void *ligand); // cmr from ligand to receptor
	
	void update(Vector ut, Vector uex, Vector uey, int matching_points);
	
	float compute_score(int);
	
	void write_binary(ostream*,unsigned short);
	void print_details(ostream*, unsigned short);
	void print_pie_score(unsigned short,string);
};

class MultiTransformation : public ScoringMetrics, public VerificationMetrics{
	public:
		string id;
		vector<Transformation*> transformations;
		int num_transformations;
		vector<long> transformation_ids;
		
		unsigned char *contacts;
		hash_map<long,float,hash<long>,eqlong> *residue_distances;
		hash_map<int, short, hash<int>,eqint> *lresidue_vdw_energy, *rresidue_vdw_energy;
		float vdw_energy, cscore, curvature_score, zscore, votes;
		int residues_in_interior, residues_in_intermediate, residues_in_surface;
		short cluster_id; // is the index of the cluster of which this combination is the median
		
		Vector* delta_r;
		
	MultiTransformation(string,Transformation*);
	MultiTransformation(string id, MultiTransformation *cn, Transformation *tr);
	MultiTransformation(char* buf);
	
	void compute_details(hash_map<long,float,hash<long>,eqlong> *residue_energy,
		hash_map<int, short, hash<int>,eqint> *lresidue_vdw_energy,
		hash_map<int, short, hash<int>,eqint> *rresidue_vdw_energy);
			
	void score();
	void score_addition(float*, float*, float*,float*, Transformation*);
	
	void distance(float *d, float *ud, MultiTransformation *mtr);
	void hamming_distance(int*, int*, MultiTransformation*,int,int);
			
	void print_details(ostream*);
	void print_results(ostream*);
};

class Reference_Frame_Point_Tuple{
	public:
		Reference_Frame* rf;
		int point_index;
		
	Reference_Frame_Point_Tuple(Reference_Frame*,int);	
};

class Reference_Frame_Point_Tuple_List{	
	public:
		list<Reference_Frame_Point_Tuple*> rflist;
	
	Reference_Frame_Point_Tuple_List();
};

class Hash_Table{
	public:
		float extent;
		int maxindex; // determines the array size
		Reference_Frame_Point_Tuple_List **lists;

	Hash_Table(float);
	
	void insert(Vector *p, int pindex, Reference_Frame* rf);
	void retrieve(Face *,Reference_Frame *rf, Face **, unsigned short *);
	void retrieve(Face *,int, Reference_Frame *rf, Face **, unsigned int **, unsigned int **);
};

/* using a critical point notation of the object - each object is a set of critical points and the associated normals */
class Object{
	public:
		Complex *c, *cH; // cH complex with hydrogen atoms - using becos surface generation not working with hydrogens

		int num_faces;
		Cluster** face;
		float maxpcmd, minpcmd;
		
		// for each face the type of the residue+sstructure
		int *crstype;
		float *curvature;

		Hash_Table *hash_table;
		hash_map<unsigned int,Reference_Frame*,hash<int>,eqint> reference_frame;
		int num_frames;

		int grid_num_xdivisions,grid_num_ydivisions,grid_num_zdivisions;
		int aagrid_num_xdivisions,aagrid_num_ydivisions,aagrid_num_zdivisions;
		int contactgrid_num_xdivisions,contactgrid_num_ydivisions,contactgrid_num_zdivisions;
		
		Vector *grid_origin;
		GridCell **grid;
		unsigned int grid_size, *grid_non_null_index;
		float max_grid_point_distance; 		// only among the surface grid cells
		
		// list of points and atoms in grid cell, used for computing neighbor lists
		hash_map<unsigned int,vector<int>*,hash<int>,eqint> grid_points_l;
		hash_map<unsigned int,vector<int>*,hash<int>,eqint> grid_atom_contact_neighbors_l,grid_atom_overlap_neighbors_l;
		
		// partition the surface grid cells on the basis of the closest point
		unsigned int **point_vs_grid;
		unsigned int *point_vs_grid_size;
		bool *contactgrid;	// true for cells within ss-cutoff from an aminoacid
		//hash_map<long,float,hash<long>,eqlong> smooth_grid;
		//hash_map<long,,hash<long>,eqlong> grid_points;
		
		unsigned short **aagrid, *aagrid_size;
		unsigned short **atom_unit_grid, *atom_unit_grid_size;
		unsigned int **point_unit_grid, *point_unit_grid_size;
		
		long transformation_index;
		unsigned short *match_size;
		unsigned int **match_range, **match_domain;
		DetailedScoringMetrics *tr_details;
		unsigned int *atoms_in_contact; unsigned int num_atoms_in_contact;
		GridCell **partner_cell_of_atom;
		Vector *transformed_atom_position;
		
	Object(Complex *c, Complex *cH);

	void build_hash_table();

	void match_surfaces(Object* o, int i_start, int i_end, vector<Transformation*> *, bool, long*);
	
	void match_surfaces_at_interface(Object* o, vector<Transformation*> *);
	
	bool can_build_reference_frame(int, int);

	vector<Reference_Frame*> build_reference_frames(int, int, bool);
	
	bool can_process_point(int, int, int);
	
	GridCell* access_point_in_grid(Vector *point, Vector *normal);
	
	int find_close_points(Vector *point, int*);

	Transformation* generate_transformation(Object *receptor, Reference_Frame *receptor_rf, Object *ligand, Reference_Frame *ligand_rf, 
						int votes);// vector<int>,vector<int>);
		
	void compute_energy_approx3(Object *receptor, Object *ligand, Transformation *tr, ProtProtDetailedScoringMetrics *details);
	
	void compute_tobibahar_energy(Object *receptor, Object *ligand, Transformation *tr, ProtProtDetailedScoringMetrics *details);
	
	void compute_luskolnick_energy(Object *receptor, Object *ligand, Transformation *tr, ProtProtDetailedScoringMetrics *details);
	
	void compute_electrostatic_energy_sphere(Object *receptor, Object *ligand, Transformation *tr, DetailedScoringMetrics *details);
	
	void compute_electrostatic_energy_pb(Object *receptor, Object *ligand, Transformation *tr, DetailedScoringMetrics *details);
	
	void compute_electrostatic_energy_coloumb(Object *receptor, Object *ligand, Transformation *tr);
	
	void compute_electrostatic_energy_ewald(Object *receptor, Object *ligand, Transformation *tr);
	
	//void compute_deltavolsa_approx1(Object *receptor, Object *ligand, Vector *, unsigned int*, unsigned int, GridCell **, Transformation *, ProtProtDetailedScoringMetrics *details);
	
	void compute_deltavolsa_approx2(Object *receptor, Object *ligand, Vector *, unsigned int*, unsigned int, GridCell **, Transformation *, ProtProtDetailedScoringMetrics *details);

	void compute_deltasas_dssp(Object *receptor, Object *ligand, Transformation *, ProtProtDetailedScoringMetrics *details);
	
	void compute_deltasa_byresidue(Object *receptor, Object *ligand, Vector *, unsigned int*, unsigned int, GridCell **, Transformation *, float *, float *,
		float *, float *);
  
	void compute_dpp_energy(Object *receptor, Object *ligand, Transformation *tr);
	
	void compute_dpp_contacts(Object *receptor, Object *ligand, Transformation *tr, hash_map<long,unsigned short,hash<long>,eqlong> *pair_d);
	
	void compute_seqavg_dpp_energy(Object *, Object *, Transformation *,hash_map<long, long, hash<long>, eqlong>*, hash_map<long, long, hash<long>, eqlong>*, 
		hash_map<const char*, Sequence*, hash<const char*>, eqstr> *,hash_map<const char*, Sequence*, hash<const char*>, eqstr> *,int,int);
		
	void compute_seqavg_energy(Object *, Object *, Transformation *,float**,float**,unsigned short**,unsigned short**,float *,float*,float *, float *,
		hash_map<long, long, hash<long>, eqlong>*, hash_map<long, long, hash<long>, eqlong>*, Species**, int, ProtProtDetailedScoringMetrics***);
	
	void compute_seqavg_spline_energy(Object *, Object *, Transformation *,float *,float*,float *, float *,
		hash_map<long, long, hash<long>, eqlong>*, hash_map<long, long, hash<long>, eqlong>*, Species**, int, ProtProtDetailedScoringMetrics***);
	
	void compute_detailed_seqavg_energy(Object *, Object *, Transformation *, hash_map<long, long, hash<long>, eqlong>*, hash_map<long, long, hash<long>, eqlong>*,	
		Species**, int, ProtProtDetailedScoringMetrics ***);
	
	bool filter3(Object* ligand, Transformation *tr);
	
	float irmsd(Object*, Transformation*, Transformation*);
	
	bool consistent_with_info(Object *receptor, Object *ligand, ModelInfo *, Transformation *tr);
	
	bool consistent_with_info_T50(Object *receptor, Object *ligand, ModelInfo *, Transformation *tr);

	void compute_contacts(Object *receptor, Object *ligand, Transformation *tr );
	
	void compute_detailed_scores(Object* ligand, Transformation *tr);
	
	vector<Transformation*> score_transformations(char*transformationsfile,Object* receptor);
	
	void score(Object *o, MultiTransformation*);
	
	Transformation **score_reference(Object* o, Complex* reference, Complex *referenceH, string *rchain, string *lchain, hash_map<long, long, hash<long>, eqlong>*,
			hash_map<long, long, hash<long>, eqlong>*, hash_map<long, long, hash<long>, eqlong>*,
			hash_map<long, long, hash<long>, eqlong>*, int);
			
	void verify_transformations(Object* o, Complex* reference, string *rchain, string *lchain, hash_map<long, long, hash<long>, eqlong>*,
			hash_map<long, long, hash<long>, eqlong>*, hash_map<long, long, hash<long>, eqlong>*,
			hash_map<long, long, hash<long>, eqlong>*, vector<Transformation*>*, vector<Transformation*>*, vector<Transformation*>*,
			int, vector<Transformation*>*, Transformation **);
	
	void verify_combinations(Object* o, Complex* reference, hash_map<long, long, hash<long>, eqlong>*,
			hash_map<long, long, hash<long>, eqlong>*, hash_map<long, long, hash<long>, eqlong>*,
			hash_map<long, long, hash<long>, eqlong>*, vector<MultiTransformation*>*);
	
	void print_grid(char* filename, short format);
	
	
	
};

class ProtRnaObject : public Object{
public:
	ProtRnaObject(Complex *c, Complex *cH);

	void compute_detailed_scores_protrna(ProtRnaObject* ligand, Transformation *tr);

	//virtual void compute_detailed_scores(Object* ligand, Transformation *tr);

	void compute_prot_rna_energy_approx1(ProtRnaObject *receptor, ProtRnaObject *ligand, Transformation *tr, ProtRnaDetailedScoringMetrics *details);


};

float compute_rmsd(int num_points, Vector* reference_points, Vector* model_points, Transformation *tr);
	
void compute_rmsd(int total_num_atoms, int ligand_start,Vector* reference_points, Vector* receptor_points, Vector* ligand_points, Transformation *tr);

float tr_irmsd_distance(Complex* , Complex* , Transformation* , Transformation*);

float tr_irmsd_distance(Complex* , Object* , Transformation* , Transformation*);

float tr_frac_shared_contacts(Complex* , Complex* , Transformation* , Transformation*);

//float tr_frac_shared_contacts(Complex* , Object* , Transformation* , Transformation*);

void write_as_pdb(Complex* receptor, Complex* ligand, Transformation* tr, string filename);

void write_as_crd(Complex* receptor, Complex* ligand, Transformation* tr);

void optimize_sidechains(Complex* receptor, Complex* ligand, Transformation* tr);

Transformation *optmize_transformation_bruteforce(Object *, Object *, Transformation *tr);
	
void refine_sidechains(Object*,Object*, vector<Transformation*> matching_frames);
			
void verify(Object*,Object*, Complex* reference, vector<Transformation *>, int);
	
void run_scwrl(void *);

void read_usphere_solution();

extern ostream *out;

