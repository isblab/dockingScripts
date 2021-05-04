/* PIEDOCK                                                            *
 * Authors: D. V. S. Ravikant                                         *
 * (C) 2011 D. V. S. Ravikant and Ron Elber                           */

#include "utilities.hpp"
#include "rmsd.h"

#include <fstream>
#include <sstream>
#include <string>
#include <cstdlib>
#include <algorithm>
#include <cctype>
#include <iomanip>

#include <ext/hash_set>

#ifndef INCLUDED_ERRNO
#define INCLUDED_ERRNO
#include <errno.h>
#endif

#define MAX_ATOMS 20000
#define MAX_ATOM_RADIUS 5.0

//#define probe_sigma 0.8 
#define PROBE_RADIUS 1.3 //1.5
//#define MAX_RADIUS 2048
#define triangle_area 0.125
#define triangle_side_length 0.5
#define max_saddle_breadth triangle_side_length
#define max_saddle_length triangle_side_length
#define convex_face_inc triangle_side_length

#define STRAIGHT_EDGE 1
#define CONVEX_EDGE 2
#define CIRCLE_EDGE 3

#define TRIANGLE_FACE 11
#define SADDLE_FACE 12
#define TORUS_FACE 13 // the torus is free
#define ATOM_FACE 14
#define CLUSTER_FACE 15

#define CONVEX 20
#define CONCAVE 21
#define PLANE 22

//without hydrogen atoms
#define NUM_CLUSTERS(concave,convex) (0.2*concave + 0.7*convex)
//with hydrogens added
//#define NUM_CLUSTERS(concave,convex) (0.2*concave + 0.6*convex)
#define NEIGHBOR_CUTOFF 5.0

// focus on pruning combinations using a different filter, normals become unstable when convex and concave regions are combined
//#define NUM_CLUSTERS(concave,convex) concave + convex
#define KMEDIANS_TRIALS 10
#define MAX_D_USED_FOR_NORMAL 5.0

#define DISTANCE_NORMAL_WEIGHT 1 //3
#define DISTANCE_CURVATURE_WEIGHT 10
#define DENSITY_WEIGHT 1
#define CONCAVE_FACE_MIN_AREA_RATIO 0.01
#define CONCAVE_FACE_MAX_AREA_RATIO 0.1
#define CONVEX_FACE_MIN_AREA_RATIO 0.05
#define CONCAVE_FACE_MAX_CLUSTER_RADIUS 1.5

#define PROTEIN 0
#define NA 1
#define DNA 2
#define RNA 3

#define ALA 0
#define ARG 1
#define ASN 2
#define ASP 3
#define CYS 4
#define GLN 5
#define GLU 6
#define GLY 7
#define HIS 8
#define ILE 9
#define LEU 10
#define LYS 11
#define MET 12
#define PHE 13
#define PRO 14
#define SER 15
#define THR 16
#define TRP 17
#define TYR 18
#define VAL 19
#define NTER 20
#define CTER 21

#define IS_AMINOACID(t) ((t>=0 && t<20) ? 1 : 0)

#define GAP -1

#define DA 100
#define DT 101
#define DC 102
#define DG 103
#define RA 104
#define RU 105
#define C 106
#define G 107

#define IS_NUCLEOTIDE(t) ((t>=100 && t<110) ? 1 : 0)

#define AMIDE_N 0
#define CARBONYL_O 1

#ifndef NUM_RESIDUE_TYPES
#define NUM_COARSE_ATOM_TYPES 2
#define NUM_ATOM_TYPES 20 //32 //85 
#define NUM_ATOM_DISTANCE_DIVISIONS 1 //2
#define NUM_RESIDUE_TYPES 22
#define NUM_COARSE_RTYPES 5
#define NUM_COARSE_SSTYPES 4
#define NUM_ENTROPY_DIVISIONS 3
// blast nr profile entropy 0 - 0.51 - 0.91 - 4
#define ENTROPY_INT1_MAX 0.51
#define ENTROPY_INT2_MAX 0.91
//#define ENTROPY(e) ((e <= ENTROPY_INT1_MAX)? 0 : ((e <= ENTROPY_INT2_MAX)? 1 :2))
#define ENTROPY(e) ((e <= 1)? 0 : ((e <= 2)? 1 :2))
#define NUM_RESIDUE_TYPES_CONSERV 20
#endif

#define POL 0 //51 GLY PRO ASN GLN SER THR CYS
#define BAS 1 //52 LYS ARG
#define ACI 2 //53 ASP GLU
#define ALI 3 //54 ALA VAL LEU ILE MET
#define ARO 4 //55 PHE TYR TRP HIS
#define RTYPE(r) ((r<0)? r : ((r == LYS || r == ARG) ? BAS : ((r == ASP || r == GLU) ? ACI : ((r == PHE || r == TYR || r == TRP || r == HIS) ? ARO : ((r == ALA || r == VAL || r == ILE || r == LEU || r == MET) ? ALI : POL)))))

#define DSSP_B 0
#define DSSP_H 1
#define DSSP_S 2
#define DSSP_U 3
#define SSTYPE(s) ((s == 'T' || s == 'S') ? DSSP_S : ((s == 'B' || s == 'E') ? DSSP_B : ((s == 'H' || s == 'G' || s == 'I') ? DSSP_H : DSSP_U)))

#define CRD 0
#define PDB 1
#define PROCESSED 2
#define PQR 3

#define ALPHA_CARBON "CA"

#define CURVATURE_SPHERE_RADIUS 5

#ifndef KMEANS
#define KMEANS 1
#define DENSITY_CLUSTER 2
#endif

#define PDB_DIR "pdb"

#define MSMS_EXECUTABLE "external/msms"
#define msms_density 2.0
#define MSMS_EXPANSION_FACTOR 1.30

#define DISTANCE_DIVISIONS_FOR_VOLUME 10

#define ATOM_CLASH_CUTOFF_SQUARED 9.0
#define BBATOM_CLASH_CUTOFF_SQUARED 6.25 //4.0

#define VDW_ATTR(s) ((s>1.8) ? 0 : ((s > 1.32623) ? 0.15*(1.8-s)/0.47377 : ((s>1.1) ? 0.15 : ((s>0.8)? 0.5*(s-0.8) : 0))))
#define VDW_REPUL(s) ((s<0.8) ? (0.8-s) : 0)

//#define AA_CUTOFF 6.5
//#define AA_CUTOFF_SQUARED 42.25

#define AA_CUTOFF 6.0
#define AA_CUTOFF_SQUARED 36.0
#define AA_SMTHP_CUTOFF 8.0
#define AA_SMTHP_CUTOFF_SQUARED 64.0
#define AA_SMTHP_FACTOR(d) ((d < AA_CUTOFF) ? 1.0 : ((d < AA_SMTHP_CUTOFF) ? (AA_SMTHP_CUTOFF-d)/2.0 : 0 ))
// the computation of aminoacid contacts in Transformation::compute_details() is correct only if the atom_cutoff is >= aa_cutoff

#define BB_CUTOFF 4.0 //4.5
#define BB_CUTOFF_SQUARED 16.0
#define BB_SMTHP_CUTOFF 6.0
#define BB_SMTHP_CUTOFF_SQUARED 36.0
#define BB_SMTHP_FACTOR(d) ((d < BB_CUTOFF) ? 1.0 : ((d < BB_SMTHP_CUTOFF) ? (BB_SMTHP_CUTOFF-d)/2.0 : 0 ))

#define BS_CUTOFF 5.0 //5.6
#define BS_CUTOFF_SQUARED 25.0
#define BS_SMTHP_CUTOFF 7.0
#define BS_SMTHP_CUTOFF_SQUARED 49.0
#define BS_SMTHP_FACTOR(d) ((d < BS_CUTOFF) ? 1.0 : ((d < BS_SMTHP_CUTOFF) ? (BS_SMTHP_CUTOFF-d)/2.0 : 0 ))

#define SS_CUTOFF 6.8 //6.4
#define RES_ATOM_CUTOFF 3

#define ATOM_STEPP_DCUTOFF 6.0
#define ATOM_STEPP_DCUTOFF_SQUARED 36.0
#define ATOM_STEPP_DMIN_SQUARED 4.0 // 6.25
#define ATOM_STEPP_DMAX_SQUARED 42.25 // 30.25

#define ATOM_SMTHP_DCUTOFF 8.0 //11.0
#define ATOM_SMTHP_DCUTOFF_SQUARED 64.0 //121.0
// linear interpolate, adjust factor 5 if cutoffs are adjusted
//#define ATOM_SMTHP_FACTOR(d) ((d < ATOM_STEPP_DCUTOFF) ? 1.0 : ((d < ATOM_SMTHP_DCUTOFF) ? (ATOM_SMTHP_DCUTOFF-d)/5.0 : 0 ))
#define ATOM_SMTHP_FACTOR(d) ((d < ATOM_STEPP_DCUTOFF) ? 1.0 : ((d < ATOM_SMTHP_DCUTOFF) ? (ATOM_SMTHP_DCUTOFF-d)/2.0 : 0 ))

#define NUM_RNA_PARTICLE_TYPES 6
#define RA_CUTOFF 8.0
#define RA_CUTOFF_SQUARED 64.0
#define RA_SMTHP_CUTOFF 10.0
#define RA_SMTHP_CUTOFF_SQUARED 100.0
#define RA_SMTHP_FACTOR(d) ((d < RA_CUTOFF) ? 1.0 : ((d < RA_SMTHP_CUTOFF) ? (RA_SMTHP_CUTOFF-d)/2.0 : 0 )) 

#define PA_CUTOFF 9.0
#define PA_CUTOFF_SQUARED 81.0
#define PA_SMTHP_CUTOFF 13.0
#define PA_SMTHP_CUTOFF_SQUARED 169.0
#define PA_SMTHP_FACTOR(d) ((d < PA_CUTOFF) ? 1.0 : ((d < PA_SMTHP_CUTOFF) ? (PA_SMTHP_CUTOFF-d)/4.0 : 0 ))

#define SUGA_CUTOFF 9.0
#define SUGA_CUTOFF_SQUARED 81.0
#define SUGA_SMTHP_CUTOFF 11.0
#define SUGA_SMTHP_CUTOFF_SQUARED 121.0
#define SUGA_SMTHP_FACTOR(d) ((d < PA_CUTOFF) ? 1.0 : ((d < SUGA_SMTHP_CUTOFF) ? (SUGA_SMTHP_CUTOFF-d)/2.0 : 0 ))

#define EPS_PROTEIN 2.0
#define EPS_WATER 80.0

class Face{
	public:
		unsigned int face_number;
		short type, orientation;
		// normal always points outward
		Vector *point, *normal;
		float radius,area,curvature;
		unsigned short num_neighbors;
		
		Face* root_face;

		bool is_buried;
		
		// keep track of the residues that form each face
		vector<unsigned short> face_monomer_indices;
		vector<unsigned short> face_atom_indices;
		
	void print_details();
};

class Cluster : public Face{
	public:
		Vector *coarser_normal;
		vector<Face *> faces;
		Face *median;
		float normalstability;
	
	Cluster();
	Cluster(Face*);
	Cluster(char*);
	void compute_area_point_normal_curvature_atoms();
	void compute_coarsenormal_radius(Face **, int);
	float get_spread();
	//bool can_accomodate(Face *f);
	//void add(Cluster*);
	void print_details(ostream*);
};

// each edge is either a straight line or an arc
class Edge{
	public:
		Vector *start, *end;
		short type;
		Vector *center, *axis;
		float radius;
		
		// breakup of the edges of a given edge
		vector<Edge*> refinement;
		Face *left_face, *right_face;

	Edge(Vector* s, Vector* e, Vector* c, Vector* axis, float r);
	
	Edge(Vector* c, Vector* axis, float r);

	//~Edge();
};

class Cycle{
	public:
		list<Edge *> edges;
		unsigned short atom;
		Vector* center;
		float radius;

		vector<Vector*> vertices_on_other_cycles;
		
	Cycle(Vector* c, float r);

	//~Cycle();

	Vector inverse_stereographic_project(Vector *pole, Vector *point);

	bool contains(Cycle* c);

	void print_details();
};

class Triangle : public Face{
	public:
		Vector *v1,*v2,*v3;
		Edge *e12,*e23,*e31;
		Vector *center;

	Triangle();

	Triangle(Vector *v1, Vector *v2, Vector *v3, Vector* center, float radius,short orientation);

	void triangulate(vector<Triangle*> *);
		
	float compute_area();
	
	void print_details();
};

class Atom : public Cluster{
	public:
		unsigned short index, cindex;
		string name;
		const char* monoindex;
		unsigned short monocindex;
		short atom_type;
		Vector* position;
		float mass, charge, sigma, eps, sigma_cubed, sqrt_eps; //, radius - there is a radius in Face;
		float vdwsa, sasa, volume, overlapvolume, seoverlapvol, sasa_complex, vdwsa_complex;
		bool isbbatom;

		unsigned int	mask[8], pmask[8]; 	/* SASA masks */
		unsigned int	vmask[DISTANCE_DIVISIONS_FOR_VOLUME][8],sevmask[DISTANCE_DIVISIONS_FOR_VOLUME][8];
		
		vector<Atom*> neighbors;

		vector<Triangle *> triangles;
		
		vector<Face*> toruses;
		vector<Face*> convexfaces;
		//vector<Face *> trianglefaces;		

	Atom();

	Atom(unsigned short index, string name, const char* aaindex, float, float,float, float sigma, float eps, float charge, float mass, short);

	~Atom();

	static float get_max_radius();
	
	static bool get_atom_details(string atom, string amino_acid, float*, float*,float*,float*,short*);

	static vector<Atom *> get_common_neighbors(Atom* a1, Atom* a2);

	// requires atoms in common neighbors to be sorted
	vector<Atom *> get_neighbors_in(vector<Atom *> common_neighbors);

	// uses the fact that the neighbors are sorted by index
	static vector<Atom *> get_intersection(vector<Atom*> v1, vector<Atom*> v2);

	void build_convex_faces();

	void triangulate();
	
	void compute_point();
	
	static int get_max_cycles_on_convexface(vector<Face*> cfaces);	

	void compute_sas_burial(Atom *a, float d, float *b, float *bprime, float *fraclost);
	
	void print_details(ostream*);
};

// Looks anti-clockwise when viewed from outside the molecule
class Triangle_Face : public Triangle{
	public:
		Atom *a1, *a2, *a3;
		
	Triangle_Face(Atom *a1, Atom *a2, Atom *a3, Vector *v1, Vector *v2, Vector *v3, Vector *center);

	~Triangle_Face();
	
	void compute_point();

	void print_details();
};

// the Face::radius is not defined
// Looks clockwise when viewed from outside the molecule
class Saddle : public Face{
	public:
		Vector *v1,*v2,*v3,*v4;
		Edge *e12,*e23,*e34, *e41;
		vector<Triangle*> triangles;

	Saddle(Vector *v1, Vector *v2, Vector *v3, Vector* v4, Edge *e12, Edge *e23, Edge *e34, Edge *e41, vector<unsigned short>* );

	~Saddle();
	
	void compute_point(Vector* axis, Vector* center, float radius);

	void print_details();
};

class Torus : public Face{
	public:
		Atom *a1,*a2;
		float distance;
		Vector* center;
		//float radius;
		// axis from a1 to a2
		Vector* axis;
	
		bool is_free;
		vector<Vector*> vertices_on_c1, vertices_on_c2;	
		vector<Triangle_Face *> trianglefaces;

		Vector *ccc1,*ccc2;
		float ccr1,ccr2;
		vector<Saddle*> saddles;

		vector<Triangle *> triangles;

	Torus(Atom* a1, Atom* a2, float distance);

	~Torus();

	void build_saddle_faces();

	void triangulate();

	void triangulate_saddle(Saddle* s, float angle, Vector v1, Vector v4);

	void print_details();
};

// Cycles look anti-clockwise when viewed from outside the molecule
class Convex_Face : public Face{
	public:
		vector<Cycle *> cycles;
		Atom* atom;
		
		Convex_Face *container;
		vector<Convex_Face*> refinement;

	Convex_Face(Atom* a);

	//~Convex_Face();

	bool contains(Convex_Face *cf);
	
	void refine_edges(vector<Edge*> *edges);
	
	void triangulate(vector<Triangle*> *triangles);
		
	void print_details();
};

class Monomer{
	public:
		hash_map<const char*, Atom*, hash<const char*>, eqstr> atom;
		char chain, sstructure;
		string name, index;
		unsigned short cindex;
		short type;
		
		float bfactor, pInterface, eInterface, entropy_sequence, blast_pseudocount;
		
		Vector *centroid;
		
		//short num_neighbors;
		float sa, rsa; // rsa - relative solvent accessibility
		
		vector<Atom*> pseudoatoms;
		short num_pseudoatoms;
		
	Monomer(string index, char chain, string name);
	
	string get_name();
	static char get_symbol(short);
	char get_symbol();
		
	float distance(Monomer*);
	void compute_reference_points();
	
	void compute_pseudoatoms();
};

class Aminoacid : public Monomer{
	public:
		Atom *amide_nitrogen, *carbonyl_oxygen, *alpha_carbon;
		
		short profile_typecountsp[21], profile_typecountsn[21];
		double profile_probabilities_p[21], profile_probabilities_n[21], interaction_evolutionary_pressure; 
		// computed over residues on the surface of the complex
		float entropy_sequence_percentile, sequence_stability;
		
	Aminoacid(string index, char chain, string name);
	
	void compute_reference_points();
	
	void compute_pseudoatoms();
	
	void compute_evolutionary_pressure(float);
	
	void print_details(ostream*);
};

class Nucleotide : public Monomer{
	public:
		Atom *phosphate;
		Vector *sugar_centroid;
	
	Nucleotide(string index, char chain, string name);
	
	void compute_reference_points();
};

class Molecule{
	public:
		string pdbcode; // assume is in lower case
		char chain;
		
		short type;
		
		hash_map<unsigned short, Atom*, hash<unsigned short>,eqint> atom;
		hash_map<unsigned short, Atom*, hash<unsigned short>,eqint> pseudoatom;
		unsigned short num_atoms;
		float mass, charge, sa;
		
		// array of size 2 - average, variance
		hash_map<long, float*, hash<long>, eqlong> distance_constraints;
	
	~Molecule();
	
	//void write_as_pdb(ostream *);
	virtual void print_details(ostream*){};
};

class Protein : public Molecule{
	public:
		string compute_aasequence();
	
		hash_map<const char*, Aminoacid*, hash<const char*>,eqstr> aminoacid;
		unsigned short num_aminoacids;
	
	Protein(fstream*, string pdbid ,char chainid,short filetype);
	void print_details(ostream*);
};

class NucleicAcid : public Molecule{
	public:
		hash_map<const char*, Nucleotide*, hash<const char*>,eqstr> nucleotide;
		unsigned short num_nucleotides;
	
	NucleicAcid(fstream*, string pdbid ,char chainid,short filetype);
};

#define GRID_SPACING 1.0 //1.0 
#define ABSENT -1000
#define INTERIOR 5
#define INTERIOR_BB 15
#define INTERMEDIATE 3
#define SURFACE 1
#define SAS_NOT_SES 6
#define OUTER_SHELL 99
#define EXTERIOR_ATOM_NEIGHBOR 100
#define POINT_ATOM_TYPE 3
#define CENTROID_TYPE 4

class GridCell{
	public:
		unsigned int index;
		short type;
		// distance from surface negative for cells inside the molecule and positive for cells outside the molecule
		float distance_from_surface;
		unsigned int *points;
		unsigned short num_points;
		unsigned short *atom_overlap_neighbors, num_atom_overlap_neighbors;
		// contact neighbors are neighbors in contact and not overlapping - avoid double saving for overlap neighbors
		unsigned short *atom_contact_neighbors, num_atom_contact_neighbors;
		unsigned short closest_point;
		float curvature;
		
	GridCell(unsigned int);
	GridCell(unsigned int,short);
	void build_atom_neighbor_array(vector<int>*,vector<int>*);
	void build_point_array(vector<int>*);
};

/*
 * Assuming within a pdb file the atom index is unique, do not require chainid to uniquely identify an atom given index
 */
class Complex{
	public:
		string pdbcode, chains;
		short filetype;
		hash_map<char, Molecule*, hash<char>, eqint> molecules;
		int num_atoms, num_monomers;
		
		Atom **atom;
		Monomer **monomer;
		hash_map<char, unsigned short, hash<char>, eqint> mmonostart,mmonoend;
		hash_map<const long, Torus *, hash<long>,eqlong> toruses;
		hash_map<const long, Triangle_Face *, hash<long>,eqlong> trianglefaces;
		
		vector<Triangle *> triangles;
		vector<Vector *> triangle_centroids;
		
		vector<Cluster*> clusters;
		Vector *center_of_mass, *center, *dipole_moment;
		float mass, charge, sa;
		float diameter /*atom center to atom center*/ , variance, rgyration, max_distance_from_cm /*cm to surface*/, min_distance_from_cm;
		
		Vector *min_sphere_center;
		float min_sphere_radius;
		Atom **charged_atom;
		unsigned short num_charged_atoms;
		
		short type;
		
		int num_aminoacids;
		Aminoacid **aminoacid;
		bool **aacontact_core, **aacontact_rim, **aacontact;
		
			
	Complex(string id, string chains, short filetype);
	
	Complex(Atom**, float*, float*,int,int);
	
	~Complex();
	
	void compute_motions();
	
	void compute_surface();

	void eliminate_faces_in_holes();
	
	void triangulate_surface();

	void eliminate_triangles_in_holes();
	
	void compute_convex_surface_of_atom(int);
	
	void cleanup_except_atom(int atom_index);
	
	void triangulate_surface(int atom_index);
		
	void compute_curvature(unsigned int, Face **);
	
	void compute_sphere_based_curvature(Vector *grid_origin, int grid_num_xdivisions, int grid_num_ydivisions, int grid_num_zdivisions, GridCell **grid, 
		unsigned int, Face **);
	
	void compute_points(int);
	
	void compute_points_msms(Vector *grid_origin, int grid_num_xdivisions, int grid_num_ydivisions, int grid_num_zdivisions, GridCell **grid, int);
	
	void cluster_points_KMedians(Face **, int , int);
	
	void cluster_points_agglomerative(Face**, int , int);
	
	//void compute_sas_burial();
	
//	void compute_volume_sasa(bool);
	
	void compute_volume_sasa(bool,bool);
	
	void compute_electrostatics();
	
	void write_as_pdb(string filename);
	
	void print_details(ostream*);	
	
	void compute_aacontacts();
	
	void compute_stability();
	
	void read_sppider_predictions(const char *);	
};
	

class ModelInfo {
	public: 
		hash_set<unsigned short, hash<short>, eqint> rinterface, linterface;
	
		ModelInfo();
};

namespace std{
	template <> class less<Monomer*> {
		public:
			bool operator()(Monomer* const& aa1, Monomer* const& aa2);
	};
	template <> class less<Aminoacid*> {
		public:
			bool operator()(Aminoacid* const& aa1, Aminoacid* const& aa2);
	};
	template <> class less<Atom*> {
		public:
			bool operator()(Atom* const& a1, Atom* const& a2){
				return a1->cindex < a2->cindex;
			};
	};
}

#ifndef DOCK_DIR
#define DOCK_DIR "moil.dock/"
#endif

#ifndef MOP_DIR
#define MOP_DIR "moil.mop/"
#endif

#ifndef EQSTR
#define EQSTR
struct eqstr
{
  bool operator()(const char* s1, const char* s2) const
  {
    return strcmp(s1,s2) == 0;
  }
};
#endif

extern float probe_radius;

void read_molecule_config(short);

void read_molecule_config();

string get_aaname(char);

void update_cluster(Cluster*,hash_map<int,int, hash<int>,eqint>);

float cluster_density(Cluster* c1, Cluster* c2);

void get_edges(vector<Face*> *convexfaces,vector<Edge*>*);

Vector* compute_gravitational_center(vector<Triangle*> triangles, float *area);

#define OBJ 12
#define VRML1 13
#define VRML2 14

void print_triangles(vector<Triangle*>, char*, short);

void print_points(vector<Cluster*>*, char*, short);

void print_spheres(float **, unsigned int, char*, short);

short get_aatype(char);

float compute_rmsd(int, Vector*, Vector*);

float compute_lrmsd(int, int, Vector*, Vector*, float*);

extern float **atom18_potential,**atom20_potential;

extern ostream *out;

