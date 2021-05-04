/* PIEDOCK                                                            *
 * Authors: D. V. S. Ravikant                                         *
 * (C) 2011 D. V. S. Ravikant and Ron Elber                           */

#include "molecule.hpp"

/*
 * KNOWN ERRORS
 * edge directions are flipped on saddles - Atom::build_convex_faces
 */

hash_map<const char*, float, hash<const char*>, eqstr> particle_sigma, particle_charge, particle_mass, particle_eps;
hash_map<const char*, hash_map<const char*, const char*, hash<const char*>, eqstr>, hash<const char*>, eqstr> monomer_atom_vs_particle;
hash_map<const char*, short, hash<const char*>, eqstr> particle_type;
hash_map<const char*, short, hash<const char*>, eqstr> monomer_name_vs_type;

float probe_radius = PROBE_RADIUS;

float *vdw_energy,*residue_solvation_energy,**residue_potential, **residue_bkbn_potential, **typeconservation_potential;//, **coarse_potential ; 
float *score_param, vdw_weight, res_bkbn_vs_atomp_weight;
float ***cacar_potential, ***cmcar_potential, ***cmcmr_potential;
float **atom18_potential, **atom20_potential, ***atom32_dpotential;

char* tmp_dir;
string piedock_home;

void coarsen_atomtypesto18();
void coarsen_atomtypesto20();
void coarsen_atomtypesto32();
void coarsen_atomtypesto46();


bool Atom::get_atom_details(string atom, string aminoacid, float *sigma, float *eps, float *charge, float* mass, short *atype){
	*atype = -1;
	if(monomer_atom_vs_particle.count(aminoacid.c_str()) > 0){
		hash_map<const char*, const char*, hash<const char*>, eqstr> atom_vs_particle = monomer_atom_vs_particle[aminoacid.c_str()];
		const char *aname;
		aname = atom.c_str();
		//if(atom == "OXT")
		//	aname = (string("O")).c_str();
		if(atom_vs_particle.count(aname) > 0){
			const char* particle_name = atom_vs_particle[aname];
			*sigma = particle_sigma[particle_name];
			*charge = particle_charge[particle_name];
			*mass = particle_mass[particle_name];
			*eps = particle_eps[particle_name];
			*atype = particle_type[particle_name];
			return true;
		}
		 /* else {
			 *out << "WARNING: did not find the atom ";
		} */
	} else {
		// HETATM records
		const char* particle_name = atom.c_str(); 
		if(particle_sigma.count(particle_name) > 0){
			*sigma = particle_sigma[particle_name];
			*charge = particle_charge[particle_name];
			*mass = particle_mass[particle_name];
			*eps = particle_eps[particle_name];
			*atype = particle_type[particle_name];
			return true;	
		}
		/* else
			 *out << " WARNING: did not find the aminoacid ";
		 */

	}
	// *out << "in " << aminoacid << " " << atom << endl;
	//exit(0);
	return false;
}

Atom::Atom(unsigned short index, string name, const char* aaindex, float x, float y, float z, float sigma, float eps, float charge, float mass, short atom_type){
	this->index = index;
	this->name = name; //*(new string(name));
	this->monoindex = aaindex;
	position = new Vector(x,y,z);
	this->radius = sigma/2.0;
	//if(radius == 0.15)	radius = 1.20;
		
	this->charge = charge;
	this->mass = mass;
	this->sigma = sigma;
	this->eps = eps;
	sigma_cubed = sigma*sigma*sigma;
	sqrt_eps = sqrt(eps);
	this->atom_type = atom_type;
	//*out << name << " " << radius << endl;
	
	type = ATOM_FACE;
	orientation = CONVEX;
	
	((Face*) this)->face_atom_indices.push_back(cindex);
	is_buried = false;
	sasa_complex = vdwsa_complex = 0;
}

vector<Atom *> Atom::get_common_neighbors(Atom* a1, Atom* a2){
	return get_intersection(a1->neighbors, a2->neighbors);
}

// requires atoms in common neighbors to be sorted
vector<Atom *> Atom::get_neighbors_in(vector<Atom *> common_neighbors){
	return get_intersection(neighbors, common_neighbors);
}

// uses the fact that the neighbors are sorted by index
vector<Atom *> Atom::get_intersection(vector<Atom*> v1, vector<Atom*> v2){
	vector<Atom *> result;

	vector<Atom *>::iterator index1 = v1.begin();
	vector<Atom *>::iterator index2 = v2.begin();
     
	while(index1 != v1.end() && index2 != v2.end()){
	    int p1 = ((Atom*) *index1)->cindex;
		int p2 = ((Atom*) *index2)->cindex;
		if(p1 < p2)
			index1++;
		else if(p1 > p2)
			index2++;
		else { // p1 == p2{
			result.push_back(*index1);	
		
			index1++;
			index2++;
		}
    }

	return result;
}

/* build convex faces on the atom
 * each edge links to the neighboring faces
 */
void Atom::build_convex_faces(){
	//*out << "build convex faces " << cindex << " " << radius << endl;
	vector<Edge *> edges;
	// build the list of edges from the toruses

	//*out << "#toruses " << toruses.size() << endl;
	for(vector<Face*>::iterator itr = toruses.begin(); itr != toruses.end(); itr++){
		Torus *ts = (Torus*) *itr;
		bool is_a1_of_ts = (cindex == (ts->a1)->cindex);
	
		if(ts->is_free){
			Edge *c;
			if(is_a1_of_ts)
				c = new Edge(ts->ccc1, ts->axis, ts->ccr1); 
			else
				c = new Edge(ts->ccc2, new Vector(Vector(0,0,0)-*(ts->axis)), ts->ccr2);
			//*out << "torus face " << ts->ccr1 << " " << ts->ccr2 << " " << c->radius << endl;
			c->right_face = ts;
			edges.push_back(c);
		} else if(!ts->is_buried) {
			//*out << ts->a2->cindex << " " << ts->is_buried << " #saddles " << (ts->saddles).size() << " #triangles " << ts->trianglefaces.size() << endl;
			for(vector<Saddle*>::iterator sitr = (ts->saddles).begin(); sitr != (ts->saddles).end(); sitr++){
				Saddle *s = *sitr;
				Edge *e;
				
				if(is_a1_of_ts)
					e = s->e12;
				else
					e = s->e34;
				
				Edge *en = new Edge(new Vector(e->start), new Vector(e->end), new Vector(e->center), new Vector(e->axis), e->radius);
				en->right_face = s;
				edges.push_back(en);
			}
		}
	}

	// find cycles of edges
	int num_edges = edges.size();
	Edge *edge[num_edges+1];
	int i = 0;
	for(vector<Edge*>::iterator eitr = edges.begin(); eitr != edges.end(); eitr++){
		edge[i++] = *eitr;
	}
	
	int next_edge[num_edges+1], previous_edge[num_edges+1];
	bool examined_edge[num_edges+1];
	for(int i = 0 ; i < num_edges ; i++){
		next_edge[i] = previous_edge[i] = -1;
		examined_edge[i] = false;
	}
	
	int num_cycles, num_circles = 0, num_non_circles = 0;
	Cycle* cycle[num_edges + 1];
	int current_cycle = 0;
	bool edge_flipped = false;
	for(int i = 0 ; i < num_edges ; i++){
		Edge *e = edge[i];
		
		if(e->type == CIRCLE_EDGE){
			Cycle* c = new Cycle(position,radius);
			c->edges.push_back(e);
			cycle[current_cycle++] = c;
			num_circles++;
			examined_edge[i] = true;
		} else {
			if(!examined_edge[i]){
				Cycle* c = new Cycle(position,radius);
				
				int current_edge = i;
				while(!examined_edge[current_edge]){
					examined_edge[current_edge] = true;
					e = edge[current_edge];
					c->edges.push_back(e);

					// find the next edge
					//*out << "finding the next edge of " << current_edge << endl;
					bool found_next = false;
					for(int j = 0 ; !found_next && j < num_edges ; j++){
						if( edge[j]->type != CIRCLE_EDGE && j != current_edge && *(e->end) == *(edge[j]->start)){
							if(previous_edge[j] != current_edge || *(e->start) == *(edge[j]->end)){
								//*out << "next edge " << current_edge << " " << j << "\n";
								found_next = true;
								next_edge[current_edge] = j; previous_edge[j] = current_edge;
								current_edge = j;
							}
						}
					}
					if(!found_next){
						for(int j = 0 ; !found_next && j < num_edges ; j++){
							if( edge[j]->type != CIRCLE_EDGE && j != current_edge && *(e->end) == *(edge[j]->end)){
								if(next_edge[j] != current_edge || *(e->start) == *(edge[j]->start)){
									edge_flipped = true;
									//*out << "next edge (flipped) " << current_edge << " " << j << "\n";
									found_next = true;
									Vector* v = edge[j]->start;
									edge[j]->start = edge[j]->end;
									edge[j]->end = v;
									next_edge[current_edge] = j; previous_edge[j] = current_edge;
									current_edge = j;
								}
							}
						}
					}
				}
				if(c->edges.size() == 1){
					*out << "ERROR: break in cycle on atom " << cindex << endl;
					out->flush(); exit(-1);
				}
				cycle[current_cycle++] = c;
				num_non_circles++;
			}
		}
	}
	if(edge_flipped){
		*out << "WARNING: edge flipped in Atom::build_convex_faces" << endl;
		*out << "Atom " << cindex << " #toruses " << toruses.size() << " #edges " << edges.size() << " details " << endl;
		for(int i = 0 ; i < num_edges ; i++){
			Edge *e = edge[i];
			if(e->type != CIRCLE_EDGE)
				*out << "(" << e->start->x << "," << e->start->y << "," << e->start->z << ") -> ("
					<< e->end->x << "," << e->end->y << "," << e->end->z << ")\t";
			else
				*out << "circle\t";
		}
		*out << endl;
	}
		
	num_cycles = num_circles + num_non_circles;
	//*out << " num cycles " << num_cycles << " " << current_cycle << endl;

	bool interior[num_cycles+1][num_cycles+1];
	for(int i = 0 ; i < num_cycles ; i++)
		for(int j = 0 ; j < num_cycles ; j++)
			if(i != j){
				interior[i][j] =  cycle[j]->contains(cycle[i]);
				//*out << i << "(" << j << ")-" << interior[i][j] << endl;
			}

	// find complete cliques in interior and construct the convex faces
	bool examined_cycle[num_cycles+1];
	for(int i = 0 ; i < num_cycles; i++)
		examined_cycle[i] = false;

	int cycles_in_face[num_cycles+1];
	int num_cycles_in_face = 0;
	for(int i = 0 ; i < num_cycles ; i++){
		if(!examined_cycle[i]){
			examined_cycle[i] = true;

			Convex_Face *cf = new Convex_Face(this);
			convexfaces.push_back(cf);
			cf->cycles.push_back(cycle[i]);

			for(int j = 0 ; j < num_cycles ; j++){
				cycles_in_face[j] = -1;
			}
			num_cycles_in_face = 0;
			cycles_in_face[0] = i;
			num_cycles_in_face++;

			// using the fact that one cycle cannot be part of multiple faces
			for(int j = 0 ; j < num_cycles ; j++){
				if(!examined_cycle[j]){
					bool adjacent_to_all_cycles = true;
					for(int k = 0 ; k < num_cycles_in_face; k++){
						adjacent_to_all_cycles &= ( interior[j][cycles_in_face[k]] && interior[cycles_in_face[k]][j] );
					}
					if(adjacent_to_all_cycles){
						cf->cycles.push_back(cycle[j]);
						cycles_in_face[num_cycles_in_face++] = j;
						examined_cycle[j] = true;
					}
				}
			}
		}
	}

	//*out << "Atom " << index << " " << cindex << " #convex faces " << convexfaces.size() << endl;
	// " details" << endl;
	// convex face is incident on all edges
	for(vector<Face*>::iterator cfitr = convexfaces.begin(); cfitr != convexfaces.end(); cfitr++){
		Convex_Face *cf = (Convex_Face*) *cfitr;
		//cf->print_details();
		for(vector<Cycle*>::iterator cycleitr = (cf->cycles).begin(); cycleitr != (cf->cycles).end(); cycleitr++){
			Cycle *c = *cycleitr;
			for(list<Edge*>::iterator eitr = (c->edges).begin(); eitr != (c->edges).end(); eitr++){
				Edge *e = *eitr;
				e->left_face = cf;
			}
		}
	}
}
	
void Atom::print_details(ostream *out){
    *out << index << " " << name << " " << position->x << " " << position->y << " " << position->z << " " << monoindex << " " << is_buried << endl;
    
    //*out << "% #convex faces " << convexfaces.size() << endl;
	/*fout << "Neighbors ";
	vector<Atom *>::iterator itr2 = neighbors.begin();
	while(itr2 != neighbors.end()){
		fout << ((Atom*) *itr2)->index << " " ;
		itr2++;
	}
	fout << endl;*/
	
	/**out << "Convexfaces ";
	vector<Face*>::iterator itr3 = convexfaces.begin();
	while(itr3 != convexfaces.end()){
		((Face*) *itr2)->print_details() ;
		itr3++;
	}
	fout << endl;*/
}


Edge::Edge(Vector* s, Vector* e, Vector* c, Vector* axis, float r){
	start = s;
	end = e;
	center = c;
	this->axis = axis;
	radius = r;

	type = CONVEX_EDGE;
}

// assumes that the circle is counterclockwise when seen down the axis
Edge::Edge(Vector* c, Vector* axis, float r){
	center = c;
	this->axis = axis;
	radius = r;
	
	type = CIRCLE_EDGE;
}


/*
Edge::~Edge(){
	//delete start;
	//delete end;
	//delete center;
	delete axis;
}
*/

Cycle::Cycle(Vector* c, float r){
	center = c;
	radius = r;
}

/*
Cycle::~Cycle(){
		edges.clear();
	}
*/

/*
 * returns a point in the plane tangent to the sphere at the base
 */
Vector Cycle::inverse_stereographic_project(Vector *pole, Vector *point){
	Vector u = *pole - *center;
	u.normalize();
	Vector r = *point - *pole;
	float f = - r.dot(u);
	Vector R = *pole + r * ((2*radius)/f);
	
	return R;
}

/*
 * Determines containment by checking if a point in the interior cycle is inside the outer cycle
 *  -- independent of the orientation of the inner cycle
 */
bool Cycle::contains(Cycle* c){
	Vector pole;
	Edge* e = (Edge* ) (c->edges).front();
	if(e->type == CIRCLE_EDGE){
		Vector v = Vector(1,0,0);
		Vector axis_perpendicular = v - *(e->axis) * (v.dot(e->axis));
		if(axis_perpendicular == Vector(0,0,0)){
			v = Vector(0,1,0);
			axis_perpendicular = v - *(e->axis) * (v.dot(e->axis));
		}
		axis_perpendicular.normalize();

		pole = *(e->center) + axis_perpendicular*(e->radius);
		//*out << "check " << radius << " " << (e->radius) << " " << Vector::distance(*center, *(e->center)) << endl;  
	} else {
		pole = *(e->start);
	}

	Vector u = pole - *center;
	u.normalize();
	
	e = (Edge* ) edges.front();
	int num_vertices;
	if(e->type == CIRCLE_EDGE || edges.size() == 2)
		num_vertices = 4;
	else
		num_vertices = edges.size();
	
	Vector v[num_vertices];
	if(e->type == CIRCLE_EDGE){
		Vector vv = Vector(1,0,0);
		Vector axis_perpendicular = vv - *(e->axis) * (vv.dot(e->axis));
		if(axis_perpendicular == Vector(0,0,0)){
			vv = Vector(0,1,0);
			axis_perpendicular = vv - *(e->axis) * (vv.dot(e->axis));
		}
		axis_perpendicular.normalize();

		Vector u1 = Vector(axis_perpendicular);
		Vector u2 = u1.cross(*(e->axis));

		v[0] = *(e->center) + (u1)*(e->radius);
		v[1] = *(e->center) + (u2)*(e->radius);
		v[2] = *(e->center) - (u1)*(e->radius);
		v[3] = *(e->center) - (u2)*(e->radius);
		
		/*out << "check " << (v[0] - pole).dot(e->axis) << "\t";
		*out << (v[1] - pole).dot(e->axis) << "\t";
		*out << (v[2] - pole).dot(e->axis) << "\t";
		*out << (v[3] - pole).dot(e->axis) << endl;*/
	} else if(edges.size() == 2){
		v[0] = *e->start;
		Vector vv = *e->start * 0.5 + *e->end * 0.5;
		float d = (vv - *e->center).norm();
		v[1] = vv * (e->radius/d) - (*e->center) * ((e->radius - d)/d);

		e = (Edge *) edges.back();
		v[2] = *e->start;
		d = (vv - *e->center).norm();
		v[3] = vv * (e->radius/d) - (*e->center) * ((e->radius - d)/d);
		
		/**out << "distance " << Vector::distance(v[1],v[3]) << endl;
		Edge* en1 = (Edge* ) (c->edges).front();
		Edge* en2 = (Edge* ) (c->edges).back();
		Edge* eo1 = (Edge* ) edges.front();
		Edge* eo2 = (Edge* ) edges.back();
		
		*out << "check edge 1 " << (*eo1->center - *en1->center).dot(*eo1->axis) << " " << (*eo1->axis == *en1->axis) << endl;
		*out << "check edge 2 " << (*eo2->center - *en2->center).dot(*eo2->axis) << " " << (*eo2->axis == *en2->axis) << endl;
		
		*out << "check end points " << (*eo1->start - *en1->start).dot(*eo1->axis) << " " << " ";
		*out << " " << (*eo1->end - *en1->end).dot(*eo1->axis) << " " << endl;
		
		*out << "edges parallel? " << (*eo1->start - *eo1->end).dot(*en1->start - *en1->end) << endl;
		/*for(int i = 0; i < 4; i++)
			*out << "(" << v[i].x << "," << v[i].y << "," << v[i].z << ")  (" 
					<< Vector::distance(v[i],*center) <<")-> ";
		*out << endl;*/
	} else {
		int index = 0;
		for(list<Edge*>::iterator itr = edges.begin(); itr != edges.end(); itr++){
			v[index++] = *(((Edge *) *itr)->end);
		}
	}
	
	Vector pv[num_vertices];
	for(int i = 0 ; i < num_vertices; i++){
		pv[i] = inverse_stereographic_project(&pole,&v[i]);
		//*out << "("<< v[i].x << ","<< v[i].y << ","<< v[i].z << ")->("<< pv[i].x << ","<< pv[i].y << ","<< pv[i].z << ")" << endl;
	}

	return ( (pv[1] - pv[0]).cross(pv[2] - pv[1]).dot(u) > 0);
}

void Cycle::print_details(){
	list<Edge *>::iterator eitr = edges.begin(), eend = edges.end();
	while(eitr != eend){
		Edge *e = *eitr;
		if(e->type != CIRCLE_EDGE)
			*out << "(" << e->start->x << "," << e->start->y << "," << e->start->z << ") -> ("
				<< e->end->x << "," << e->end->y << "," << e->end->z << ") "
				<< "(" << (*e->start - *e->end).norm()
				<<  " c = (" << e->center->x << "," << e->center->y << "," << e->center->z << ") "
				<<  " a = (" << e->axis->x << "," << e->axis->y << "," << e->axis->z << ") "
				<< " r = " << e->radius
				<< ") => ";
		else
			*out << " circle edge, center (" << e->center->x << "," << e->center->y << "," << e->center->z << ") radius " << e->radius << " ";
			eitr++;
		}
		*out << endl;
	}

void Face::print_details(){
	*out << type << " point (" << point->x << "," << point->y << "," << point->z << 
		") normal (" << normal->x << "," << normal->y << "," << normal->z << ") area " << area << endl;
	
	*out << "Atoms #" << face_atom_indices.size() << " ";
	for(vector<unsigned short>::iterator itr = face_atom_indices.begin(); itr != face_atom_indices.end(); itr++){
		*out << *itr << " ";
	}
	*out << endl;
	
	*out << "Residues #" << face_monomer_indices.size() << " ";
	for(vector<unsigned short>::iterator itr = face_monomer_indices.begin(); itr != face_monomer_indices.end(); itr++){
		*out << *itr << " ";
	}
	*out << endl;
}
	

Triangle::Triangle(){
	is_buried = false;
}

Triangle::Triangle(Vector* v1, Vector* v2, Vector* v3, Vector* c, float r, short orientation){
	this->v1 = v1;
	this->v2 = v2;
	this->v3 = v3;
	center = c;
	radius = r;
	this->orientation = orientation;
}

/*
 * Recursively triangulate, return the set of final triangles
 */
void Triangle::triangulate(vector<Triangle*> *result){
	float d12 = Vector::distance(v1,v2), d23 = Vector::distance(v2,v3), d31 = Vector::distance(v3,v1);
	Vector *vm1 = v1, *vm2 = v2, *vm3 = v3; float dmax = d12;
	if(d23 > dmax){
		vm1 = v2; vm2 = v3; vm3 = v1; dmax = d23;
	}
	if(d31 > dmax){
		vm3 = v3; vm2 = v1; vm3 = v2; dmax = d31;
	}
	if(dmax > triangle_side_length){
		Vector *vmid = new Vector( (*vm1 + *vm2) * (1.0/2.0) );

		Triangle *t1 = new Triangle(vm1,vmid,vm3,center,radius,orientation);
		Triangle *t2 = new Triangle(vmid,vm2,vm3,center,radius,orientation);
	
		t1->triangulate(result);
		t2->triangulate(result);
	} else
		result->push_back((Triangle*) this);
}

void Triangle::print_details(){
	//*out << Vector::distance(v1,v2) << " " << Vector::distance(v2,v3) << " " << Vector::distance(v3,v1) << endl;
	*out << "(" << v1->x << "," << v1->y << "," << v1->z << ")\t";
	*out << "(" << v2->x << "," << v2->y << "," << v2->z << ")\t";
	*out << "(" << v3->x << "," << v3->y << "," << v3->z << ")\t";
	*out << compute_area() << "\n";
	/*if(center != NULL){
		Vector vn = (*v1 + *v2 + *v3) * (1.0/3.0) - *center;
		vn.normalize();
		*out << "\t(" << vn.x << "," << vn.y << "," << vn.z << ")";
	}
	*out << endl;*/
}

Triangle_Face::Triangle_Face(Atom *a1, Atom *a2, Atom *a3, Vector *v1, Vector *v2, Vector *v3, Vector *probe_center){
	this->a1 = a1;
	this->a2 = a2;
	this->a3 = a3;
	this->v1 = v1;
	this->v2 = v2;
	this->v3 = v3;

	Vector *a = new Vector((*v2 - *probe_center).cross(*v1 - *probe_center));
	a->normalize();
	e12 = new Edge(v1,v2,probe_center, a, probe_radius);
	
	a = new Vector((*v3 - *probe_center).cross(*v2 - *probe_center));
	a->normalize();
	e23 = new Edge(v2,v3,probe_center, a, probe_radius);
	
	a = new Vector((*v1 - *probe_center).cross(*v3 - *probe_center));
	a->normalize();
	e31 = new Edge(v3,v1,probe_center, a, probe_radius);

	e12->left_face = this;
	e23->left_face = this;
	e31->left_face = this;

	center = probe_center;
	radius = probe_radius;
	orientation = CONCAVE;
	type = TRIANGLE_FACE;
	
	((Face*) this)->face_atom_indices.push_back(a1->cindex);
	((Face*) this)->face_atom_indices.push_back(a2->cindex);
	((Face*) this)->face_atom_indices.push_back(a3->cindex);
}

Triangle_Face::~Triangle_Face(){
	delete e12;
	delete e23;
	delete e31;
}

void Triangle_Face::print_details(){
	*out << "Atoms " << a1->cindex << " " << a2->cindex << " " << a3->cindex << "\t";

	Vector *v = v1;
	*out << "v1 = (" << v->x << "," << v->y << "," << v->z << ")\t";
	v = v2;
	*out << "v2 = (" << v->x << "," << v->y << "," << v->z << ")\t";
	v = v3;
	*out << "v3 = (" << v->x << "," << v->y << "," << v->z << ")\t";

	v = e12->start;
	*out << "e12 = (" << v->x << "," << v->y << "," << v->z << ") -> ";
	v = e12->end;
	*out << "(" << v->x << "," << v->y << "," << v->z << ")\t";
	v = e23->start;
	*out << "e23 = (" << v->x << "," << v->y << "," << v->z << ") -> ";
	v = e23->end;
	*out << "(" << v->x << "," << v->y << "," << v->z << ")\t";
	v = e31->start;
	*out << "e31 = (" << v->x << "," << v->y << "," << v->z << ") -> ";
	v = e31->end;
	*out << "(" << v->x << "," << v->y << "," << v->z << ")\t";

	*out << a1->radius << " " << (*v1 - *(a1->position)).norm() << " ";
	*out << a2->radius << " " << (*v2 - *(a2->position)).norm() << " ";
	*out << a3->radius << " " << (*v3 - *(a3->position)).norm() << " ";
	*out << endl;
}


Saddle::Saddle(Vector *v1, Vector *v2, Vector *v3, Vector* v4, Edge* e12, Edge *e23, Edge* e34, Edge *e41, vector<unsigned short> *atom_indices){
	this->v1 = v1;
	this->v2 = v2;
	this->v3 = v3;
	this->v4 = v4;

	this->e12 = e12;
	this->e23 = e23;
	this->e34 = e34;
	this->e41 = e41;

	e12->right_face = this;
	e23->right_face = this;
	e34->right_face = this;
	e41->right_face = this;
					
	type = SADDLE_FACE;
	
	((Face*) this)->face_atom_indices = *atom_indices;
}

Saddle::~Saddle(){
}

void Saddle::print_details(){
	Vector* v = v1;
	*out << "v1 = (" << v->x << "," << v->y << "," << v->z << ")\t";
	v = v2;
	*out << "v2 = (" << v->x << "," << v->y << "," << v->z << ")\t";
	v = v3;
	*out << "v3 = (" << v->x << "," << v->y << "," << v->z << ")\t";
	v = v4;
	*out << "v4 = (" << v->x << "," << v->y << "," << v->z << ")\t";

	v = e12->start;
	*out << "e12 = (" << v->x << "," << v->y << "," << v->z << ") -> ";
	v = e12->end;
	*out << "(" << v->x << "," << v->y << "," << v->z << ")\t";
	v = e23->start;
	*out << "e23 = (" << v->x << "," << v->y << "," << v->z << ") -> ";
	v = e23->end;
	*out << "(" << v->x << "," << v->y << "," << v->z << ")\t";
	v = e34->start;
	*out << "e34 = (" << v->x << "," << v->y << "," << v->z << ") -> ";
	v = e34->end;
	*out << "(" << v->x << "," << v->y << "," << v->z << ")\t";
	v = e41->start;
	*out << "e41 = (" << v->x << "," << v->y << "," << v->z << ") -> ";
	v = e41->end;
	*out << "(" << v->x << "," << v->y << "," << v->z << ")\t";

	*out << endl;
}


Torus::Torus(Atom* a1, Atom* a2, float distance){
	this->a1 = a1;
	this->a2 = a2;
	this->distance = distance;

	axis = new Vector(*(a2->position) - *(a1->position));
	axis->normalize();

	float t = ((a1->radius+probe_radius)*(a1->radius+probe_radius) - (a2->radius+probe_radius)*(a2->radius+probe_radius))/(distance*distance);
	center = new Vector( *(a1->position)*((1-t)/2.0) + *(a2->position)*((1+t)/2.0) );

	radius = 1/2.0*sqrt((a1->radius+a2->radius + 2*probe_radius)*(a1->radius+a2->radius + 2*probe_radius) - distance*distance)*sqrt(distance*distance - (a1->radius-a2->radius)*(a1->radius-a2->radius))/distance;

	//*out << (distance*distance - (a1->radius-a2->radius)*(a1->radius-a2->radius)) << " " << radius << endl; 
	is_free = true;
	is_buried = false;
	type = TORUS_FACE;

	// contact circle on atom a1
	ccc1 = new Vector(*(a1->position) * ((probe_radius)/(a1->radius + probe_radius)) + *(center) * ((a1->radius)/(a1->radius + probe_radius)));
	ccr1 = radius*(a1->radius)/(a1->radius + probe_radius);

	// contact circle on atom a2
	ccc2 = new Vector(*(a2->position) * ((probe_radius)/(a2->radius + probe_radius)) + *(center) * ((a2->radius)/(a2->radius + probe_radius)));
	ccr2 = radius*(a2->radius)/(a2->radius + probe_radius);

	//*out << ccr1 << " " << ccr2 << endl;
	a1->toruses.push_back(this);	
	a2->toruses.push_back(this);
	
	((Face*) this)->face_atom_indices.push_back(a1->cindex);
	((Face*) this)->face_atom_indices.push_back(a2->cindex);
}

Torus::~Torus(){
	delete axis;
	delete center;
	delete ccc1;
	delete ccc2;
	
	for(vector<Saddle*>::iterator sitr = saddles.begin(); sitr != saddles.end(); sitr++){
		delete *sitr;
	}
}

void Torus::build_saddle_faces(){
	if(!is_free){
		int num_vertices = trianglefaces.size();
		if(num_vertices % 2 == 1){
			*out << "ERROR: torus odd number of trianglefaces " << a1->index << " " << a1->cindex << " " << a2->index << " " << a2->cindex << " " << trianglefaces.size() << endl;
			for(vector<Triangle_Face *>::iterator itr = trianglefaces.begin(); itr != trianglefaces.end(); itr++){
				Triangle_Face *t = *(itr);
				t->print_details();
			}
			exit(-1);
		} else if(num_vertices > 0) {
			/*out << "building saddle faces " << a1->index << " " << a2->index << " " << is_buried << " " << num_vertices << " "
				<< Vector::distance(a1->position, a2->position) << " " << Vector::distance(a1->position, center) << " "
				<< Vector::distance(a2->position, center) << endl;*/
			// form the list of vertices on the torus
			Triangle_Face *tface[num_vertices];
			Vector *a1_vertex[num_vertices], *a2_vertex[num_vertices], tfedge[num_vertices];
			int i = 0;
			for(vector<Triangle_Face *>::iterator itr = trianglefaces.begin(); itr != trianglefaces.end(); itr++){
				Triangle_Face *t = *(itr);			
				tface[i] = t;
				
				short v1_index;
				if(a1->cindex == (t->a1)->cindex){
					a1_vertex[i] = t->v1;
					v1_index = 1;
				} else if(a1->cindex == (t->a2)->cindex){
					a1_vertex[i] = t->v2;
					v1_index = 2;
				} else {
					a1_vertex[i] = t->v3;
					v1_index = 3;
				}

				short v2_index;
				if(a2->cindex == (t->a1)->cindex){
					a2_vertex[i] = t->v1;
					v2_index = 1;
				} else if(a2->cindex == (t->a2)->cindex){
					a2_vertex[i] = t->v2;
					v2_index = 2;
				} else {
					a2_vertex[i] = t->v3;
					v2_index = 3;
				}
				
				if( ( v2_index - (v1_index + 1) ) % 3 == 0 )
					tfedge[i] = *a2_vertex[i] - *a1_vertex[i];
				else
					tfedge[i] = *a1_vertex[i] - *a2_vertex[i];

				i++;
			}
			
			// sort on the basis of angle with the first vertex
			float angle_with_base[num_vertices];
			for(int i = 0 ; i < num_vertices; i++){
				angle_with_base[i] = Vector::angle( *(a1_vertex[0]) - ccc1 , *(a1_vertex[i]) - ccc1 , *axis);
			}
	
			for(int i = 0 ; i < num_vertices; i++)
				for(int j = i+1 ; j < num_vertices; j++)
					if(angle_with_base[j] < angle_with_base[i]){
						float d = angle_with_base[i];
						angle_with_base[i] = angle_with_base[j];
						angle_with_base[j] = d;
	
						Vector* v = a1_vertex[i];
						a1_vertex[i] = a1_vertex[j];
						a1_vertex[j] = v;
				
						v = a2_vertex[i];
						a2_vertex[i] = a2_vertex[j];
						a2_vertex[j] = v;
						
						Triangle_Face *t = tface[i];
						tface[i] = tface[j];
						tface[j] = t;
						
						Vector e = tfedge[i];
						tfedge[i] = tfedge[j];
						tfedge[j] = e;
					}
	
			/*for(int i = 0 ; i < num_vertices; i++)
				*out << angle_with_base[i] << "\t";
			*out << endl;*/
			
			// 0,1 form a saddle, for this, the side of the triangleface 0 has to be opposite the axis used for sorting
			if(tfedge[0].dot(axis)>0){
				Vector *v1 = a1_vertex[0], *v2 = a2_vertex[0], e = tfedge[0];
				Triangle_Face *t = tface[0];
				for(int i = 0 ; i < num_vertices-1; i++){
					a1_vertex[i] = a1_vertex[i+1];
					a2_vertex[i] = a2_vertex[i+1];
					tface[i] = tface[i+1];
					tfedge[i] = tfedge[i+1];
				}
				a1_vertex[num_vertices-1] = v1;
				a2_vertex[num_vertices-1] = v2;
				tface[num_vertices-1] = t;
				tfedge[num_vertices-1] = e;
			}
			
			bool edges_not_alternating = false;
			for(int i = 0 ; i < num_vertices; i = i+2)
				if(tfedge[i].dot(axis) * tfedge[i+1].dot(axis) > 0)
					edges_not_alternating = true;
			if(edges_not_alternating){
				*out << "WARNING: build_saddles edges of trianglefaces not alternating" << endl;
				for(int i = 0 ; i < num_vertices; i++){
					tface[i]->print_details();
					*out << tface[i]->compute_area() << " ";
					Vector *v = a1_vertex[i];
					*out << "(" << v->x << "," << v->y << "," << v->z << ") ";
					v = a2_vertex[i];
					*out << "(" << v->x << "," << v->y << "," << v->z << ") " << endl;
				}
			}
			
			Vector *negaxis = new Vector(Vector(0,0,0) - *axis);;
			// build a saddle face between consecutive triangles
			short direction = 0, forward = 1, backward = -1;
			for(int i = 0 ; i < num_vertices; i = i+2){
				Triangle_Face *t = tface[i];
				Vector* v1 = a1_vertex[i];
				Vector* v2 = a2_vertex[i];
				short v2_index;
				if(a2->cindex == (t->a1)->cindex){
					v2_index = 1;
				} else if(a2->cindex == (t->a2)->cindex){
					v2_index = 2;
				} else {
					v2_index = 3;
				}
				
				int next = i+1;
				Triangle_Face *tnext = tface[next];
				Vector *v0 = a1_vertex[next];
				Vector* v3 = a2_vertex[next];
				short v3_index;
				if(a2->cindex == (tnext->a1)->cindex){
					v3_index = 1;
				} else if(a2->cindex == (tnext->a2)->cindex){
					v3_index = 2;
				} else {
					v3_index = 3;
				}
	
				Saddle *s;
				Edge* e10 = new Edge(v1,v0, ccc1, axis, ccr1);
				if(negaxis == NULL)
					negaxis = new Vector(Vector(0,0,0) - *axis);
				Edge* e32 = new Edge(v3,v2, ccc2, negaxis, ccr2);
	
				Edge *e21, *e03;
				switch(v2_index){
					case 1:
						e21 = t->e12;
						break;
					case 2:
						e21 = t->e23;
						break;
					case 3:
						e21 = t->e31;
						break;
				}
				switch(v3_index){
					case 1:
						e03 = tnext->e31;
						break;
					case 2:
						e03 = tnext->e12;
						break;
					case 3:
						e03 = tnext->e23;
						break;
				}
				/*if(!(*e21->end == *v1) || !(*e03->start == *v0)){
					*out << "ERROR: saddle edges not on triangle face" << endl;
					t->print_details();
					tnext->print_details();
					*out << "(" << v1->x << "," << v1->y << "," << v1->z << ") ";
					*out << "(" << v2->x << "," << v2->y << "," << v2->z << ") " << endl;
				}*/
				
				s = new Saddle(v1,v0,v3,v2, e10, e03, e32, e21, &(((Face*) this)->face_atom_indices));
				saddles.push_back(s);
				//s->print_details();
			}
		}
	}
}

void Torus::triangulate(){
	if(!is_free && !is_buried){
		//*out << "#saddles " << saddles.size() << endl;
		for(vector<Saddle *>::iterator itr = saddles.begin(); itr != saddles.end(); itr++){
			Saddle *s = *(itr);
			// compute the length of the arc on the first circle
			Vector u = *s->v3 - *s->v2;
			u.normalize();
			float angle = Vector::angle( *(s->v1) - *ccc1, *(s->v2) - *ccc1, u);
			//*out << "check axis " << ((*(s->v1) - *ccc1).cross(*(s->v2) - *ccc1)).cross(*axis).norm() << endl;
			triangulate_saddle(s,angle, s->v1, s->v4);
		}
	} else if(!is_buried) {
		Vector v = Vector(1,0,0);
        Vector axis_perpendicular = v - *axis * (v.dot(axis));
        if(axis_perpendicular == Vector(0,0,0)){
                v = Vector(0,1,0);
                axis_perpendicular = v - *axis * (v.dot(axis));
        }
        axis_perpendicular.normalize();
		
		Vector v1 = Vector( *ccc1 + axis_perpendicular*ccr1 );
		Vector v2 = Vector( *ccc2 + axis_perpendicular*ccr2 );
			
		triangulate_saddle(NULL,2*PI, v1, v2);
	}
	
	//*out << "Torus " << a1->cindex << " " << a2->cindex << " free " << is_free << " saddles " << saddles.size() << " triangles " << triangles.size() << endl;
}

/*
Triangles resulting from this triangulation do not lie on a sphere
triangles such that (v2-v1)x(v3-v1) points out?
*/
void Torus::triangulate_saddle(Saddle *s,float angle, Vector v1, Vector v4){
	float d = (*ccc1 - *ccc2).norm();
	int breadth_divisions = (int) (d/max_saddle_breadth);
	if(d/max_saddle_breadth - breadth_divisions < 0.5 && breadth_divisions > 0)
		breadth_divisions--;

	float l1 = ccr1*angle;
	float l2 = ccr2*angle;

	float max_ccr = (ccr1 > ccr2) ? ccr1 : ccr2;
	int length_divisions = (int) (max_ccr*angle/max_saddle_length);
	if(max_ccr*angle/max_saddle_length - length_divisions < 0.5 && length_divisions > 0)
		length_divisions--;

	float theta = angle/(length_divisions + 1.0);
	Vector ur = Vector(v1 - *ccc1);
	ur.normalize();
	Vector ut = ur.cross(axis);
	ut.normalize();

	//*out << "saddle divisions " << length_divisions+1 << " " << breadth_divisions<< endl;
	Vector* vertices[length_divisions+2][breadth_divisions+2];
	for(int i = 0 ; i < length_divisions+2 ; i++){
		Vector vc;
		float r = ccr1;
		float alpha = i*theta;
		for(int j = 0 ; j < breadth_divisions+2 ; j++){
			Vector v;

			vc = Vector( *ccc1 + (*ccc2 - *ccc1) * (j/(breadth_divisions+1.0)) );
			r = ccr1 + ( j/(breadth_divisions+1.0) )*(ccr2 - ccr1);
			alpha = i*theta;
			v = Vector( vc + ur * (r*cos(alpha)) + ut * (r*sin(alpha)) );
			vertices[i][j] = new Vector(v);
		}
	}
	
	for(int i = 0 ; i < length_divisions+1 ; i++){
		if(is_free){
			vertices_on_c1.push_back(vertices[i][0]);
			vertices_on_c2.push_back(vertices[i][breadth_divisions+1]);
		}
		for(int j = 0 ; j < breadth_divisions+1 ; j++){				
			Triangle *t1 = new Triangle(vertices[i][j], vertices[i+1][j], vertices[i+1][j+1],NULL,0,PLANE);
			triangles.push_back(t1);
			if(s != NULL)
				(s->triangles).push_back(t1);

			Triangle* t2 = new Triangle(vertices[i+1][j+1], vertices[i][j+1], vertices[i][j],NULL,0,PLANE);
			triangles.push_back(t2);
			if(s != NULL)
				(s->triangles).push_back(t2);
			
			//*out << i << " " << j << " " << t1->compute_area() << " " << t2->compute_area() << endl;
		}
	}
}
				
void Torus::print_details(){
	*out << "Atoms " << a1->cindex << " " << a2->cindex << " r1 " << a1->radius << " r2 " << a2->radius << " d " << distance << "\tis_free " << is_free << "\tis_buried " << is_buried << "\t num triangles " << trianglefaces.size() 
		<< " ccr1 " << ccr1 << " ccr2 " << ccr2 << " axis (" << axis->x << "," << axis->y << "," << axis->z << ")" << endl;

	*out << "trianglefaces " << endl;
	vector<Triangle_Face*>::iterator titr = trianglefaces.begin(), tend = trianglefaces.end();
	while(titr != tend){
		Triangle_Face *tf = (Triangle_Face*) *titr;
		tf->print_details();
		titr++;
	}

	*out << "saddles " << endl;
	if(!is_free){
		vector<Saddle *>::iterator itr = saddles.begin();
		while(itr != saddles.end()){
			Saddle *s = *(itr);			
			s->print_details();
			itr++;
		}	
		*out << endl;
	}
}


Convex_Face::Convex_Face(Atom *a){
	atom = a;
	//refinement = NULL;
	container = NULL;
	orientation = CONVEX;
	
	((Face*) this)->face_atom_indices.push_back(a->cindex);
}

/*
Convex_Face::~Convex_Face(){
		cycles.clear();
	}
*/

void Convex_Face::print_details(){
	*out << "#cycles " << cycles.size() << endl;
	vector<Cycle *>::iterator citr = cycles.begin(), cend = cycles.end();
	while(citr != cend){
		Cycle* c = *citr;
		c->print_details();
		citr++;
	}
	*out << endl;
}

Monomer::Monomer(string index,char chain, string name){
	this->index = index;
	this->name = name;
	this->chain = chain;
	/*string crindex;
 	crindex.append(1,chain);
 	crindex.append(index);
 	cindex = *(new string(crindex));*/
	if(monomer_name_vs_type.count(name.c_str()) > 0) 
		type = monomer_name_vs_type[name.c_str()];
	else
		type = -1;
	
	sstructure = ' ';
}

Aminoacid::Aminoacid(string index,char chain, string name) : Monomer(index,chain,name){
	entropy_sequence = entropy_sequence_percentile = 0;
}

Nucleotide::Nucleotide(string index,char chain, string name) : Monomer(index,chain,name){
	
}

void Nucleotide::compute_reference_points(){
	phosphate = NULL;	
	Vector centroid = Vector(0,0,0);
	float mass = 0;
	int number_nucleotide_atoms = 0;
	float sugar_mass = 0;
	Vector sugar_centroid = Vector(0,0,0);
	float num_sugar_atoms = 0;
	for(hash_map<const char*, Atom*, hash<const char*>, eqstr>::iterator itr = atom.begin(); itr != atom.end(); itr++){
		Atom *a = itr->second;
		if(a->name == "P"){
			phosphate = a;
		} else if(a->mass > 1.0 && (a->name == "C1'" || a->name == "C2'" || a->name == "C3'" || a->name == "C4'" ||
				a->name == "C5'" || a->name == "O2'" || a->name == "O4'" )) {
			sugar_centroid = sugar_centroid + *(a->position);
			num_sugar_atoms++;

		} else if(a->mass > 1.0 && a->name != "O1P" && a->name != "O2P" && a->name != "O3'" && a->name != "O5'"){
			int namelen = a->name.length();
			char * namec = (char*) a->name.c_str();
			if(namec[namelen-1] != '\''){
				//centroid = centroid + *(a->position) * a->mass;
				centroid = centroid + *(a->position);
				number_nucleotide_atoms++;
			}
		}
		mass += a->mass;
	}
	
	if(mass != 0){
		this->centroid = new Vector(centroid*(1.0/number_nucleotide_atoms));
	} else{
		this->centroid = NULL;
	}

	if(num_sugar_atoms > 0){
		this->sugar_centroid = new Vector(sugar_centroid * (1.0/num_sugar_atoms));
	} else
		this->sugar_centroid = NULL;
}

/*
 * Assumes that centroids of both aminoacids are already computed
 */
float Monomer::distance(Monomer *ma){
	float d;
	d = Vector::distance(ma->centroid, centroid);
	return d;
}

/*
 * centroid is the centroid of the side chain
 */
void Aminoacid::compute_reference_points(){
	amide_nitrogen = carbonyl_oxygen = alpha_carbon = NULL;
	Vector centroid = Vector(0,0,0);
	float mass = 0;
	//*out << "aa " << index << " " << type << " "; out->flush();
	for(hash_map<const char*, Atom*, hash<const char*>, eqstr>::iterator itr = atom.begin(); itr != atom.end(); itr++){
		Atom *a = itr->second;
		if(a->name == "N"){
			amide_nitrogen = a;
			a->type = AMIDE_N;
		}	else if(a->name == "O"){
			carbonyl_oxygen = a;
			a->type = CARBONYL_O;
		} 	else if(a->name == "CA")
			alpha_carbon = a;
		else if(a->name != "C"){ // && a->mass > 1.0){
		//else if(a->mass > 1.0){
			centroid = centroid + *(a->position) * a->mass;
			mass += a->mass;
		}
	}
	if(mass != 0){
		this->centroid = new Vector(centroid*(1.0/mass));
	} else{
		if(alpha_carbon != NULL)
			this->centroid = new Vector(alpha_carbon->position);
		else
			this->centroid = NULL;
	}
	//*out << mass << " " << amide_nitrogen << " " << carbonyl_oxygen << endl;
}

void Aminoacid::print_details(ostream *out){
    *out << index << " " << name << " " << sstructure << " " << pInterface;
    *out << " " << entropy_sequence << " ";// << entropy_sequence_zscore << " ";
    /*out << type << " ";
    Vector v = *centroid;
    *out << v.x << " " << v.y << " " << v.z << " "; 
    *out << amide_nitrogen->index << " " << carbonyl_oxygen->index << " ";*/
    *out << endl;
}

string Monomer::get_name(){
	string name;
	switch(type){
		case ALA:
			name = "ALA";
		break;
		case ARG:
			name = "ARG"; 
		break;
		case ASN:
			name = "ASN"; 
		break; 
		case ASP:
			name = "ASP"; 
		break; 
		case CYS:
			name = "CYS"; 
		break; 
		case GLN:
			name = "GLN"; 
		break; 
		case GLU:
			name = "GLU"; 
		break;
		case GLY:
			name = "GLY"; 
		break; 
		case HIS:	
			name = "HIS"; 
		break; 
		case ILE:
			name = "ILE"; 
		break; 
		case LEU:
			name = "LEU"; 
		break; 
		case LYS:
			name = "LYS"; 
		break; 
		case MET:
			name = "MET"; 
		break; 
		case PHE:
			name = "PHE"; 
		break; 
		case PRO:
			name = "PRO"; 
		break; 
		case SER:
			name = "SER";
		break; 
		case THR:
			name = "THR"; 
		break; 
		case TRP:
			name = "TRP"; 
		break; 
		case TYR:
			name = "TYR"; 
		break; 
		case VAL:
			name = "VAL"; 
		break;
		case CTER:
			name = "CTER";
		break;
		case NTER:
			name = "NTER";
		break;
	}
	return name;
}

string get_aaname(char symbol){
	string name;
	switch(symbol){
		case 'A':
			name = "ALA";
		break;
		case 'R':
			name = "ARG"; 
		break;
		case 'N':
			name = "ASN"; 
		break; 
		case 'D':
			name = "ASP"; 
		break; 
		case 'C':
			name = "CYS"; 
		break; 
		case 'Q':
			name = "GLN"; 
		break; 
		case 'E':
			name = "GLU"; 
		break;
		case 'G':
			name = "GLY"; 
		break; 
		case 'H':	
			name = "HIS"; 
		break; 
		case 'I':
			name = "ILE"; 
		break; 
		case 'L':
			name = "LEU"; 
		break; 
		case 'K':
			name = "LYS"; 
		break; 
		case 'M':
			name = "MET"; 
		break; 
		case 'F':
			name = "PHE"; 
		break; 
		case 'P':
			name = "PRO"; 
		break; 
		case 'S':
			name = "SER";
		break; 
		case 'T':
			name = "THR"; 
		break; 
		case 'W':
			name = "TRP"; 
		break; 
		case 'Y':
			name = "TYR"; 
		break; 
		case 'V':
			name = "VAL"; 
		break;
		case DA : 
			name = "DA"; 
		break;
		case DT :
			name = "DT";
		break;
		case DC : 
			name = "DC";
		break;
		case DG : 
			name = "DG";
		break;
		case RA : 
			name = "A"; 
		break;
		case C : 
			name = "C";
		break;
		case G : 
			name = "G";
		break;
		case RU :
			name = "U";
		break;
	}
	return name;
}

char Monomer::get_symbol(){
	return get_symbol(type);
}

char Monomer::get_symbol(short type){
	char aasymbol='-';
	switch(type){
		case ALA:
			aasymbol = 'A';
		break;
		case ARG:
			aasymbol ='R'; 
		break;
		case ASN:
			aasymbol ='N'; 
		break; 
		case ASP:
			aasymbol ='D'; 
		break; 
		case CYS:
			aasymbol ='C'; 
		break; 
		case GLN:
			aasymbol ='Q'; 
		break; 
		case GLU:
			aasymbol ='E'; 
		break;
		case GLY:
			aasymbol ='G'; 
		break; 
		case HIS:
			aasymbol ='H'; 
		break; 
		case ILE:
			aasymbol ='I'; 
		break; 
		case LEU:
			aasymbol ='L'; 
		break; 
		case LYS:
			aasymbol ='K'; 
		break; 
		case MET:
			aasymbol ='M'; 
		break; 
		case PHE:
			aasymbol ='F'; 
		break; 
		case PRO:
			aasymbol ='P'; 
		break; 
		case SER:
			aasymbol ='S'; 
		break; 
		case THR:
			aasymbol ='T'; 
		break; 
		case TRP:
			aasymbol ='W'; 
		break; 
		case TYR:
			aasymbol ='Y'; 
		break; 
		case VAL:
			aasymbol ='V'; 
		break;
		case DA : case RA:
			aasymbol ='A';
		break;
		case DT :
			aasymbol ='T';
		break;
		case DC : case C:
			aasymbol ='C';
		break;
		case DG : case G:
			aasymbol ='G';
		break;
		case RU :
			aasymbol = 'U';
		break;
	}
	return aasymbol;
}

namespace std{
	bool less<Monomer*>::operator()(Monomer* const& mr1, Monomer* const& mr2){
		bool result;
		
		int a1index, a2index;
		a1index = ((Atom*) mr1->atom.begin()->second)->cindex;
		a2index = ((Atom*) mr2->atom.begin()->second)->cindex;
		result = (a1index < a2index);
		return result;
	};
	bool less<Aminoacid*>::operator()(Aminoacid* const& aa1, Aminoacid* const& aa2){
		bool result;
	  /*{ // the aminoacids have to be sorted not by the index but in the order in which they appear in the file
	  	int i1 = atoi(aa1->index.c_str()), i2 = atoi(aa2->index.c_str());
		if(i1 == i2){
			// 221A < 221B < 221
			if(aa1->index.size() != aa2->index.size())
				result = (aa1->index.size() > aa2->index.size());
			else
				result = (aa1->index < aa2->index);
		} else
			result = i1 < i2;
		//*out << aa1->index << " " << aa2->index << " " << i1 << " " << i2 << " " << result << endl;
	  }*/	
		
		int a1index, a2index;
		a1index = ((Atom*) aa1->atom.begin()->second)->cindex;
		a2index = ((Atom*) aa2->atom.begin()->second)->cindex;
		result = (a1index < a2index);
		return result;
	};
}

/*
 * Missing aminoacids ?
 * Additional aminoacids ?
 */
string Protein::compute_aasequence(){	
	stringstream ss (stringstream::in | stringstream::out);
	Aminoacid* aacid[aminoacid.size()+1];
	int aaindex=0;
	for(hash_map<const char*, Aminoacid*, hash<const char*>,eqstr>::iterator aaitr = aminoacid.begin(); aaitr != aminoacid.end(); aaitr++)
		aacid[aaindex++] = ((Aminoacid*) aaitr->second);
	
	sort(aacid, aacid+aminoacid.size(),less<Aminoacid*>());	
	/*for(int i = 0; i < aaindex; i++)
		for(int j = i+1; j < aaindex; j++)
			if(less<Aminoacid*>()(aacid[j],aacid[i])){
				Aminoacid* aa = aacid[i];
				aacid[i] = aacid[j];
				aacid[j] = aa;
			}*/
	
	for(int i = 0; i < aaindex; i++){
		Aminoacid* aa = aacid[i];
		char aasymbol = aa->get_symbol();
		*out << aa->index << " " << aa->cindex << " " << aa->name << " " << aasymbol << endl;
		//if(i+1 < aaindex)	*out << less<Aminoacid*>()(aa,aacid[i+1]) << endl;
		if(aasymbol != '-')	ss << aasymbol;
	}

	string result;
	ss >> result;
	return result;
}

/*
 * Assumes only 1 chain is present and the chain is a protein
 */
void Complex::read_sppider_predictions(const char *dir){
	stringstream ss (stringstream::in | stringstream::out);
	ss << string(dir) << "/" << pdbcode.substr(3,pdbcode.length()-3) << "_" << chains << ".sppider.pred";
	string filename;
	ss >> filename;
	fstream fin(filename.c_str(), fstream::in);
	*out << filename << " " << fin.is_open() << endl;
	char buf[8192];
	do{
		fin.getline(buf,8192);
	}while((string(buf)).find("SECTION_INT_PROBABILITIES") == string::npos);
	
	//*out << string(buf) << endl;
	int aaindex = 0;
	while (!fin.eof()){
		fin.getline(buf,8192);
		if(fin.gcount() > 0){
			stringstream line(buf,stringstream::in);
			string tag; line >> tag;
			if(tag == ">"){
				//*out << "check " << aastart << " " << aaend << endl;
				
				fin.getline(buf,8192);
				stringstream line(buf,stringstream::in);
				line >> tag;
				int n = tag.length();
				fin.getline(buf,8192);
				stringstream pline(buf,stringstream::in);
				for(int i = 0; i < n; i++){
					int prediction;
					pline >> prediction;
					((Aminoacid*) monomer[aaindex++])->pInterface = prediction;
				}
			}
		}
	}
}

Complex::Complex(string id, string chains, short filetype){
	pdbcode = *(new string(id));
	this->chains = *(new string(chains));
	this->filetype = filetype;
	const char* chains_c = chains.c_str();
	int num_chains = chains.length();
	num_atoms =  num_monomers = 0;
	
	fstream *fin = NULL, pdbfin, localfin, localchainfin;
	if(filetype == PDB || filetype == PQR || filetype == PROCESSED){
		stringstream ss (stringstream::in | stringstream::out);
		ss << pdbcode;
		switch(filetype){
			case PDB:
				ss << ".pdb";
			break;
			case PQR:
				ss << ".pqr";
			break;
			case PROCESSED:
				ss << "_" << chains << ".surface";
			break;
		}
		string filename;
		ss >> filename;
		localfin.open(filename.c_str(), fstream::in);
	
		if(!localfin.is_open()){
			if(filetype == PDB){
				ss.clear();
				ss << pdbcode << "_" << chains << ".pdb";
				ss >> filename;
				localchainfin.open(filename.c_str(), fstream::in);
				if(!localchainfin.is_open()){
					char pdbcodelowercase[pdbcode.length()+1];
					for(short i = 0; i <= pdbcode.length();i++)	pdbcodelowercase[i]=tolower(pdbcode.c_str()[i]);
					ss.clear();
					ss << piedock_home << "/" << PDB_DIR << "/pdb" << ((string) pdbcodelowercase) << ".ent";
					ss >> filename;
					*out << filename << endl;
					pdbfin.open(filename.c_str(), fstream::in);
					fin = &pdbfin;
				} else {
					fin = &localchainfin;
				}
			}
			if(!localchainfin.is_open() && !pdbfin.is_open()) {
				*out << "ERROR: could not open file " << filename << " " << errno << "\n"; out->flush();
				exit(-1);
			}
		} else
			fin = &localfin;
	  	// *out << filename << " " << fin->is_open() << endl;
	}
	
	for(int i = 0; i < num_chains; i++){
		char chain = chains_c[i];
		//*out << chain << endl; out->flush();
		if((filetype == PDB || filetype == PQR) && i > 0){
			fin->clear();
			fin->seekg(ios_base::beg);
			//*out << filename << " " << fin->is_open() << " " << fin->eof() << endl;
		}
		Molecule *m;
		// determine if Protein or Nulceic Acid
		short type = PROTEIN;
		if(filetype == PDB){
			string format;
			char buf[8192];
			while (!fin->eof()){
				fin->getline(buf,8192);
				char pdbchainid = buf[21];
				if(chain == '-' || pdbchainid == chain){
					stringstream ss(buf,stringstream::in);
					string aformat;
					ss >> aformat;
					if(aformat == "ATOM"){
						string atom_name;
				    	string monomer;
				    	{
							stringstream line(buf+13,stringstream::in);
						    line >> atom_name;
						    if(buf[11] == ' ' && buf[12] != ' '){
						    	stringstream line(buf+12,stringstream::in);
							    line >> atom_name;
						      	line >> monomer;
						    } else {
							    if(atom_name.length() > 3){
							    	string s = string(atom_name.c_str());
							    	atom_name = s.substr(0,3);
							    	monomer = s.substr(3,s.length());
							    } else
								   	line >> monomer;
						    }
				    	}
				    	if(monomer.length() > 4){
					    	string s = string(monomer.c_str());
					    	monomer = s.substr(0,4);
					    }
					    if(monomer_name_vs_type.count(monomer.c_str()) > 0) {
							short mtype = monomer_name_vs_type[monomer.c_str()];
							if(IS_AMINOACID(mtype) || IS_NUCLEOTIDE(mtype)){
								if(IS_NUCLEOTIDE(mtype)) type = NA;
								break;
							}
					    }
					}
				}
			}
			fin->clear();
			fin->seekg(ios_base::beg);
		}
		
		if(type==PROTEIN)
		{
				Protein *p = new Protein(fin,id,chain,filetype);
		 		num_monomers += p->num_aminoacids;
		 		m = (Molecule*) p;
		}
		else if(type==NA)
		{
		 		NucleicAcid *n = new NucleicAcid(fin,id,chain,filetype);
		 		num_monomers += n->num_nucleotides;
		 		m = (Molecule*) n;
		 }
		molecules[m->chain] = m;
		num_atoms += m->num_atoms;
	}
	if(filetype == PROCESSED){
		char buf[8192];			
		int num_clusters;
		fin->getline(buf,8192);
		stringstream line(buf,stringstream::in);
		line >> num_clusters;
		*out << "#clusters " << num_clusters << endl;
		clusters = *(new vector<Cluster*>); 
		//*out << "~ point normal cnormal area #atoms atom_indices #residues residue_indices" << endl;
		fin->getline(buf,8192);
		for(int i = 0 ; i < num_clusters; i++){
			fin->getline(buf,8192);
			//*out << buf << endl;
			Cluster* c = new Cluster(buf);
			clusters.push_back(c);
			//c->print_details(out);
		}
		
		int num_triangle_centroids;
		fin->getline(buf,8192);
		{
			stringstream line(buf,stringstream::in);
			line >> num_triangle_centroids;
		}
		triangle_centroids = *(new vector<Vector *>);
		//*out << num_triangle_centroids << endl;
		// ~ v.x v.y v.z
		fin->getline(buf,8192);
		for(int i = 0 ; i < num_triangle_centroids; i++){
			fin->getline(buf,8192);
			//*out << buf << endl;
			stringstream line(buf,stringstream::in);
			float x,y,z;
			line >> x;
			line >> y;
			line >> z;
			triangle_centroids.push_back(new Vector(x,y,z));
		}
		*out << " #clusters " << clusters.size() << " #centroids " << triangle_centroids.size() << endl;
	}
	if(fin != NULL)
		fin->close();
	
	atom = (Atom **) malloc(sizeof(Atom*)*num_atoms);
	int count = 0;
	for(hash_map<char, Molecule*, hash<char>, eqint>::iterator mitr = molecules.begin(); mitr != molecules.end(); mitr++){
		Molecule *m = mitr->second;
		Atom* matom[m->atom.size()+1];
		int maindex=0;
		hash_map<unsigned short, unsigned short, hash<unsigned short>,eqint> originalorder;
		
		for(hash_map<unsigned short, Atom*, hash<unsigned short>,eqint>::iterator aitr = m->atom.begin(); aitr != m->atom.end(); aitr++){
			Atom* a = aitr->second;
			matom[maindex] = a;
			//if(filetype == PROCESSED)
				originalorder[maindex] = aitr->first;
			maindex++;
		}
		
		//if(filetype == PROCESSED)
			for(int i  = 0; i < maindex; i++)
				for(int j = i+1; j < maindex; j++)
					if(originalorder[j] < originalorder[i]){
						Atom *a = matom[i];
						matom[i] = matom[j];
						matom[j] = a;
						unsigned int k = originalorder[j];
						originalorder[j] = originalorder[i];
						originalorder[i] = k;
					}
				
		for(int i = 0; i < maindex; i++){
			Atom* a = matom[i];
			atom[count] = a;
			a->cindex = count;
			count++;
		}
	}
	
	float sumpInterface = 0;
	monomer = (Monomer **) malloc(sizeof(Monomer*)*num_monomers);
	count = 0;
	for(hash_map<char, Molecule*, hash<char>, eqint>::iterator mitr = molecules.begin(); mitr != molecules.end(); mitr++){
		Molecule *m = mitr->second;
		short num_monomers;
		switch(m->type){
			case PROTEIN :
				num_monomers = ((Protein*) m)->num_aminoacids;
			break;
			case NA:
				num_monomers = ((NucleicAcid*) m)->num_nucleotides;
			break;
		}
		Monomer* mmer[num_monomers+1];
		int maaindex=0;
		
		int num_actualmmers;
		
		if (m->type==PROTEIN)
		{
				Protein *p = (Protein*) m;
				for(hash_map<const char*, Aminoacid*, hash<const char*>,eqstr>::iterator aaitr = p->aminoacid.begin(); aaitr != p->aminoacid.end(); aaitr++)
					if(((Aminoacid*) aaitr->second)->atom.size() > 0)
						mmer[maaindex++] = ((Aminoacid*) aaitr->second);
				
				num_actualmmers = maaindex;
		
				for(hash_map<const char*, Aminoacid*, hash<const char*>,eqstr>::iterator aaitr = p->aminoacid.begin(); aaitr != p->aminoacid.end(); aaitr++)
					if(((Aminoacid*) aaitr->second)->atom.size() == 0)
						mmer[maaindex++] = ((Aminoacid*) aaitr->second);
		}
		else if(m->type==NA || m->type==DNA || m->type==RNA)
		{

				NucleicAcid *n = (NucleicAcid*) m;
				for(hash_map<const char*, Nucleotide*, hash<const char*>,eqstr>::iterator aaitr = n->nucleotide.begin(); aaitr != n->nucleotide.end(); aaitr++)
					if(((Nucleotide*) aaitr->second)->atom.size() > 0)
						mmer[maaindex++] = ((Nucleotide*) aaitr->second);
				
				num_actualmmers = maaindex;
				
				for(hash_map<const char*, Nucleotide*, hash<const char*>,eqstr>::iterator aaitr = n->nucleotide.begin(); aaitr != n->nucleotide.end(); aaitr++)
					if(((Nucleotide*) aaitr->second)->atom.size() == 0)
						mmer[maaindex++] = ((Nucleotide*) aaitr->second);
			
		}
		
		sort(mmer, mmer+num_actualmmers,less<Monomer*>());	
		
		mmonostart[m->chain] = count;
		for(int i = 0; i < maaindex; i++){
			Monomer* mr = mmer[i];
			mr->cindex = count;
			monomer[count++] = mr;
			//cout << pdbcode << mr->cindex << " " << mr->index << " " << mr->atom.size() << endl;cout.flush();
			for(hash_map<const char*, Atom*, hash<const char*>, eqstr>::iterator aitr = mr->atom.begin(); aitr != mr->atom.end(); aitr++){
				((Atom*) aitr->second)->monocindex = mr->cindex;
				//cout << pdbcode << ((Atom*) aitr->second)->index << "\t" << mr->cindex << " " << mr->index << endl;cout.flush();
			}
			if(m->type == PROTEIN)
				sumpInterface += ((Aminoacid*) mr)->pInterface;
		}
		mmonoend[m->chain] = count;
	}
	
	int numInterface = 0;
	num_aminoacids = 0;
	for(int i = 0; i < num_monomers; i++)
		if(IS_AMINOACID(monomer[i]->type)){
			num_aminoacids++;
			if(((Aminoacid*) monomer[i])->pInterface > 0)
				numInterface++;
		} /* else {
			*out << "non aminoacid " << monomer[i]->name << " " << monomer[i]->index << " " << monomer[i]->cindex << " " << monomer[i]->type << endl;
		}
		*/
	
	aminoacid = (Aminoacid **) malloc(sizeof(Aminoacid*)*(num_aminoacids+1));
	unsigned short aaindex = 0;
	
	for(int i = 0; i < num_monomers; i++){
		if(IS_AMINOACID(monomer[i]->type)){
			Aminoacid* aa = (Aminoacid*) monomer[i];
			aminoacid[aaindex++] = aa;
			if(aa->pInterface > 0)
				aa->eInterface = 0 - log(aa->pInterface * numInterface/sumpInterface);
			else // buried residues
				aa->eInterface = 0;
		}
	}
	
	for(unsigned int i = 0; i < num_atoms; i++){
		Atom *a = atom[i];
		if(a->name == "N" || a->name == "CA" || a->name == "O" || a->name == "C")
			a->isbbatom = true;
		else	a->isbbatom = false;
	}
	// *out << "cmplx #atoms " << num_atoms << " #monomers " << num_monomers << " #aacids " << num_aminoacids << endl;
}

/*
 * Compute the expected fraction of surface of an atom buried by another atom
 */
void Atom::compute_sas_burial(Atom *a, float d, float *b, float *bprime, float *fraclost){
	float S=4*PI*(radius)*(radius);
	if(d < radius + a->radius && radius + d > a->radius){
		*b = PI*(radius)*(radius + a->radius - d)*(1+(a->radius - radius)/d);
		//*bprime = PI*(radius + probe_radius)*(radius + a->radius - d)*(1+(a->radius - (2*probe_radius + radius))/d);
		//if(*bprime < 0) *bprime = 0;
		//if(*bprime > S) *bprime = S;
		*bprime = 0;
		*fraclost = (*b - *bprime)/S;
	} else {
		if(d > radius + a->radius + 2*probe_radius){
			*b = *bprime = *fraclost = 0;
		}
		if(radius + d < a->radius){
			*bprime = S;
			*b = 0;
			*fraclost = 1;
		}
	}

	if(*fraclost < 0) *fraclost = 0;
	if(*fraclost > 1) *fraclost = 1;
}

/*
 * Compute the burail of each atom
 *
void Complex::compute_sas_burial(){
	sas_b = (float **) malloc(sizeof(float*) * (num_atoms + 1));
	sas_bprime = (float **) malloc(sizeof(float*) * (num_atoms + 1));
	sas_fraclost = (float **) malloc(sizeof(float*) * (num_atoms + 1));
	for(int i = 0 ; i < num_atoms; i++){
		sas_b[i] = (float*) malloc(sizeof(float) * (num_atoms + 1));
		sas_bprime[i] = (float*) malloc(sizeof(float) * (num_atoms + 1));
		sas_fraclost[i] = (float*) malloc(sizeof(float) * (num_atoms + 1));
	}
	
	for(int i = 0 ; i < num_atoms; i++){
		Atom* a1 = atom[i];
		for(int j = 0; j < num_atoms; j++)
			if(i != j){
				Atom* a2 = atom[j];
				float d = Vector::distance((a1->position), (a2->position));
				a1->compute_sas_burial(a2, d , &sas_b[i][j], &sas_bprime[i][j], &sas_fraclost[i][j]);
			}
	}
}*/

Complex::~Complex(){
	for(hash_map<char, Molecule*, hash<char>, eqint>::iterator mitr = molecules.begin(); mitr != molecules.end(); mitr++){
	 	Molecule *m = (Molecule*) mitr->second;
	 	delete m;
	}
	 free(atom);
	 free(monomer);
}

ModelInfo::ModelInfo() {
	rinterface.clear();
	linterface.clear();
};

/*
 * Assumes pdbid in lower case
 * 
 * computes reference points of aminoacids
 */
Protein::Protein(fstream *fin, string pdbid, char chain, short filetype){
	type = PROTEIN;
	
	pdbcode = *(new string(pdbid));
	if(filetype != PROCESSED)
		this->chain = chain;

	bool to_close_fin = false;
	fstream pdbfin, localfin;
	if(fin == NULL){
		to_close_fin = true;
		stringstream ss (stringstream::in | stringstream::out);
		ss << pdbcode;
		if(filetype == PDB)
			ss << ".pdb";
		else if(filetype == CRD)
			ss << "_" << chain << ".crd";
		
		string filename;
		ss >> filename;
		
		localfin.open(filename.c_str(), fstream::in);
		
		if(!localfin.is_open()){
			if(filetype == PDB){
				ss.clear();
				char pdbcodelowercase[pdbcode.length()+1];
				for(short i = 0; i <= pdbcode.length();i++)	pdbcodelowercase[i]=tolower(pdbcode.c_str()[i]);
				ss << piedock_home << "/" << PDB_DIR << "/pdb" << ((string)pdbcodelowercase) << ".ent";
				ss >> filename;
				pdbfin.open(filename.c_str(), fstream::in);
				fin = &pdbfin;
			}
			if(!pdbfin.is_open()) {
				*out << "ERROR: could not open file " << filename << " " << errno << "\n";
				exit(-1);
			}
		} else
			fin = &localfin;
	  	*out << filetype << " " << filename << " " << fin->is_open() << endl;
	}
	
	char buf[8192];
	string format;
	switch(filetype){
		case PDB:
			while (!fin->eof()){
				fin->getline(buf,8192);
				stringstream ss(buf,stringstream::in);
				string aformat;
				ss >> aformat;
				if(aformat == "ATOM" || aformat == "HETATM"){
					//*out << "|" << buf << endl;
				    int position;
					{
						stringstream line(buf,stringstream::in);
						line >> format;
					    line >> position;
					}
					string atom_name;
			    	string amino_acid;
			    	{
						stringstream line(buf+13,stringstream::in);
					    line >> atom_name;
					    if(buf[11] == ' ' && buf[12] != ' '){
					    	stringstream line(buf+12,stringstream::in);
						    line >> atom_name;
					      	line >> amino_acid;
					    } else {
						    if(atom_name.length() > 3){
						    	string s = string(atom_name.c_str());
						    	atom_name = s.substr(0,3);
						    	amino_acid = s.substr(3,s.length());
						    } else
							   	line >> amino_acid;
					    }
			    	}
			    	if(amino_acid.length() > 4){
				    	string s = string(atom_name.c_str());
				    	amino_acid = s.substr(0,4);
				    }
				    
				    // select Aresidues if there are multiple residues
					if(amino_acid.size() == 4 && !(amino_acid == "CTER") && !(amino_acid == "NTER") 
				      && monomer_atom_vs_particle.count(amino_acid.substr(1,3).c_str()) > 0 && amino_acid.c_str()[0] == 'A')
						amino_acid = amino_acid.substr(1,3);
				    
			    	char pdbchainid = buf[21];
				    if((chain == '-' || pdbchainid == chain)
				     && (( aformat == "ATOM" && (monomer_name_vs_type.count(amino_acid.c_str()) > 0 || amino_acid == "CTER" || amino_acid == "NTER"))
				     	|| (aformat == "HETATM" && amino_acid == "MSE")
				     //	|| (aformat == "HETATM" && (amino_acid == "SAM" || amino_acid == "ZN" )) // T46 modification
				     /*	|| ((aformat == "ATOM" || aformat == "HETATM") && (amino_acid == "ADP" || amino_acid == "NHE" || amino_acid == "MTT" || amino_acid == "MN" || amino_acid == "GTP" 
				     		|| amino_acid == "GOL" || amino_acid == "GNP" || amino_acid == "DIO" || amino_acid == "CA" || amino_acid == "APC" || amino_acid == "ANP" 
				     		|| amino_acid == "AGS")) */
				        )
				    ){
//				    	if(amino_acid == "MSE")	amino_acid = "MET";
				    	
						stringstream line(buf+22,stringstream::in);
				    	//*out << "||" << buf << endl;
						string aaindex;
						line >> aaindex;
						
						float x,y,z;
						line >> x;
						line >> y;
						line >> z;
							
						Aminoacid *aa;	
						if(aminoacid.count(aaindex.c_str()) == 0){
							aaindex = *(new string(aaindex));
							amino_acid = *(new string(amino_acid.c_str()));
							aminoacid[aaindex.c_str()] = new Aminoacid(aaindex,chain,amino_acid);
						}
						aa = aminoacid[aaindex.c_str()];
						//*out << aa_index << " " << amino_acid << endl;
						
						if(atom_name != ""){
							float sigma, eps;	short type;
							bool ret = Atom::get_atom_details(atom_name,amino_acid, &sigma, &eps, &charge, &mass, &type);
							if(ret){
								atom_name = *(new string(atom_name.c_str()));
								Atom* a = new Atom(position, atom_name, aa->index.c_str(), x,y,z,sigma, eps,charge,mass,type);
								//a->print_details(out);
								//*out << "ATOM " << atom->radius << " " << atom->charge << endl;
								atom[position] = a;
								aa->atom[(a->name).c_str()] = a;
							}
							/* else
								*out << buf << endl; */
						}
						
						float bfactor;
						line >> bfactor;
						line >> bfactor;
						aa->pInterface = bfactor;
				    }
			    } else if(format == "ENDMDL")
			    	break;
			}
		break;
		case CRD:	
			// first 3 lines do not contain useful information
			fin->getline(buf,8192);
			fin->getline(buf,8192);
			fin->getline(buf,8192);
			while (!fin->eof()){
				fin->getline(buf,8192);
				//*out << buf << endl;
				stringstream line(buf,stringstream::in);
				int position;
			    line >> position;
				string aaindex;
				line >> aaindex;
				string amino_acid;	
				line >> amino_acid;
					
				Aminoacid *aa;
				if(aminoacid.count(aaindex.c_str()) == 0){
					aaindex = *(new string(aaindex));
					amino_acid = *(new string(amino_acid.c_str()));
					aminoacid[aaindex.c_str()] = new Aminoacid(aaindex,chain,amino_acid);
				}
				aa = aminoacid[aaindex.c_str()];
					
				string atom_name;
				line >> atom_name;
				float x,y,z;
				line >> x;
				line >> y;
				line >> z;
			
				if(atom_name != ""){
					float sigma, eps;	short type;
					bool ret = Atom::get_atom_details(atom_name,amino_acid, &sigma, &eps, &charge, &mass, &type);
					if(ret){
						atom_name = *(new string(atom_name.c_str()));
						Atom* a = new Atom(position, atom_name, aa->index.c_str(), x,y,z,sigma, eps,charge,mass,type);
						//*out << "check " << a->radius << " " << a->charge << endl;
						atom[position] = a;
						aa->atom[(a->name).c_str()] = a;
					} /* else
						*out << buf << endl; */
				}
			}
		break;
		case PQR:
			while (!fin->eof()){
				fin->getline(buf,8192);
				stringstream ss(buf,stringstream::in);
				ss >> format;
				if(format == "ATOM"){
					//*out << "|" << buf << endl;
				    int position;
					{
						stringstream line(buf,stringstream::in);
						line >> format;
					    line >> position;
					}
					string atom_name;
			    	string amino_acid;
			    	{
						stringstream line(buf+13,stringstream::in);
					    line >> atom_name;
					    if(buf[11] == ' ' && buf[12] != ' '){
					    	stringstream line(buf+12,stringstream::in);
						    line >> atom_name;
					      	line >> amino_acid;
					    } else {
						    if(buf[12] == ' ' && atom_name.length() > 3){
						    	string s = string(atom_name.c_str());
						    	atom_name = s.substr(0,3);
						    	amino_acid = s.substr(3,s.length());
						    } else
							   	line >> amino_acid;
					    }
			    	}
			    	if(amino_acid.length() > 4){
				    	string s = string(atom_name.c_str());
				    	amino_acid = s.substr(0,4);
				    }
				    
			    	char pdbchainid = buf[21];
				    if(chain == '-' || pdbchainid == chain){
						stringstream line(buf+22,stringstream::in);
				    	//*out << "||" << buf << endl;
						string aaindex;
						line >> aaindex;
						
						string chainstr;
						if(amino_acid == "NTR" || amino_acid == "CTR"){
							char cbuf[2];
							cbuf[0] = pdbchainid; cbuf[1] = '\0';
							chainstr=string(cbuf);
						}
						if(amino_acid == "NTR")
							if(atom_name=="N"||atom_name=="CA")	aaindex=chainstr+"A";
							else	aaindex=chainstr+"B";
						else if(amino_acid == "CTR")
							if(atom_name=="OT2")	aaindex=chainstr+"C";
							else aaindex=chainstr+"D";
						
						float x,y,z;
						line >> x;
						line >> y;
						line >> z;
							
						Aminoacid *aa;	
						if(aminoacid.count(aaindex.c_str()) == 0){
							aaindex = *(new string(aaindex));
							amino_acid = *(new string(amino_acid.c_str()));
							aminoacid[aaindex.c_str()] = new Aminoacid(aaindex,chain,amino_acid);
						}
						aa = aminoacid[aaindex.c_str()];
						
						if(atom_name != ""){
							line >> charge;
							float radius;
							line >> radius;
							atom_name = *(new string(atom_name.c_str()));
							// dont know the mass
							float mass=1; short type=-1;	float b,c, eps;
							bool ret = Atom::get_atom_details(atom_name,amino_acid, &b, &eps, &c, &mass, &type);
							Atom* a = new Atom(position, atom_name, aa->index.c_str(), x,y,z,2*radius,eps,charge,mass,type);
							//a->print_details(out);
							atom[position] = a;
							aa->atom[(a->name).c_str()] = a;
						}
				    }
			    }
			}
		break;
		case PROCESSED:
			fin->getline(buf,8192);
			stringstream line(buf,stringstream::in);
			line >> num_aminoacids;
			fin->getline(buf,8192);
			for(int i = 0 ; i < num_aminoacids; i++){
				fin->getline(buf,8192);
				stringstream line(buf,stringstream::in);
				string aaindex; line >> aaindex;
				aaindex = *(new string(aaindex));
				string name; line >> name;
				name = *(new string(name));
				Aminoacid *aa = new Aminoacid(aaindex,chain,name);
				aminoacid[aaindex.c_str()] = aa;
				char space = line.get();
				aa->sstructure = line.get();
				line >> aa->pInterface;
				line >> aa->entropy_sequence;
				//line >> aa->entropy_sequence_zscore;
				//aa->print_details(out); out->flush();
			}
			fin->getline(buf,8192);
			{
				stringstream line(buf,stringstream::in);
				line >> num_atoms;
			}
			//*out << "~ chain index name position aaindex aaname is_buried" << endl;
			fin->getline(buf,8192);
			for(int i = 0 ; i < num_atoms; i++){
				fin->getline(buf,8192);
				//*out << i << " " << buf << endl; out->flush();
				stringstream line(buf,stringstream::in);
				char chain; line >> chain;
				this->chain = chain;
				int position; line >> position;
				string atom_name; line >> atom_name;
				atom_name = *(new string(atom_name));
				float x,y,z;
				line >> x;
				line >> y;
				line >> z;
				string aaindex; line >> aaindex;
				int is_buried; line >> is_buried;
				
				Aminoacid *aa = aminoacid[aaindex.c_str()];
				if(atom_name != ""){
					float sigma, eps;	short type;
					bool ret = Atom::get_atom_details(atom_name,aa->name, &sigma, &eps, &charge, &mass, &type);
					if(ret){
						Atom* a = new Atom(position, atom_name, aaindex.c_str(), x,y,z,sigma, eps,charge,mass,type);
						a->is_buried = (is_buried!=0);
						//*out << "check " << atom->radius << " " << atom->charge << endl;
						atom[position] = a;
						aa->atom[a->name.c_str()] = a;
						//a->print_details(out); out->flush();
					} /* else
						 *out << buf << endl; */
				}
			}
			//correct the chain assignment for the aminoacids
			for(hash_map<const char*, Aminoacid*, hash<const char*>,eqstr>::iterator aaitr = aminoacid.begin(); aaitr != aminoacid.end(); aaitr++){
				Aminoacid* aa = (aaitr->second);
				aa->chain = this->chain;
			}
			/*if(chain != this->chain){
				cout << "WARNING: chains do not match " << chain << ":" << this->chain << endl;
				this->chain = chain;
			}*/
		break;
	}
	if(to_close_fin){
		fin->close();
		fin = NULL;
	}
	
	num_atoms = atom.size();
	num_aminoacids = aminoacid.size();
	//*out << "#atoms " << num_atoms << " " << "#aminoacids " << num_aminoacids << endl;
	int num_pseudoatoms = 1;
	for(hash_map<const char*, Aminoacid*, hash<const char*>,eqstr>::iterator aaitr = aminoacid.begin(); aaitr != aminoacid.end(); aaitr++){
		Aminoacid* aa = (aaitr->second);
		aa->compute_reference_points();
		/*aa->compute_pseudoatoms();
		vector<Atom*>::iterator aapsitr = aa->pseudoatoms.begin(), 
			aapsend = aa->pseudoatoms.end();
		while(aapsitr != aapsend){
			//*out << "pseudoatom " << num_pseudoatoms << " " << ((Atom*) *aapsitr)->mass << endl;
			pseudoatoms[num_pseudoatoms++] = (Atom*) *aapsitr;
			aapsitr++;
		}*/
	}
	//*out << "initiated protein" << endl;
}

NucleicAcid::NucleicAcid(fstream *fin, string pdbid, char chain, short filetype){
	type = NA;
	
	pdbcode = *(new string(pdbid));
	if(filetype != PROCESSED)
		this->chain = chain;

	bool to_close_fin = false;
	fstream pdbfin, localfin;
	if(fin == NULL){
		to_close_fin = true;
		stringstream ss (stringstream::in | stringstream::out);
		ss << pdbcode;
		if(filetype == PDB)
			ss << ".pdb";
		else if(filetype == CRD)
			ss << "_" << chain << ".crd";
		
		string filename;
		ss >> filename;
		
		localfin.open(filename.c_str(), fstream::in);
		
		if(!localfin.is_open()){
			if(filetype == PDB){
				ss.clear();
				char pdbcodelowercase[pdbcode.length()+1];
				for(short i = 0; i <= pdbcode.length();i++)	pdbcodelowercase[i]=tolower(pdbcode.c_str()[i]);
				ss << piedock_home << "/" << PDB_DIR << "/pdb" << ((string)pdbcodelowercase) << ".ent";
				ss >> filename;
				pdbfin.open(filename.c_str(), fstream::in);
				fin = &pdbfin;
			}
			if(!pdbfin.is_open()) {
				*out << "ERROR: could not open file " << filename << " " << errno << "\n";
				exit(-1);
			}
		} else
			fin = &localfin;
	  	*out << filetype << " " << filename << " " << fin->is_open() << endl;
	}
	
	char buf[8192];
	string format;
	switch(filetype){
		case PDB:
			while (!fin->eof()){
				fin->getline(buf,8192);
				stringstream ss(buf,stringstream::in);
				string aformat;
				ss >> aformat;
				if(aformat == "ATOM" || aformat == "HETATM"){
					//*out << "|" << buf << endl;
				    int position;
					{
						stringstream line(buf,stringstream::in);
						line >> format;
					    line >> position;
					}
					string atom_name;
			    	string nucleotide_name;
			    	{
						stringstream line(buf+13,stringstream::in);
					    line >> atom_name;
					    if(buf[11] == ' ' && buf[12] != ' '){
					    	stringstream line(buf+12,stringstream::in);
						    line >> atom_name;
					      	line >> nucleotide_name;
					    } else {
						    if(atom_name.length() > 3){
						    	string s = string(atom_name.c_str());
						    	atom_name = s.substr(0,3);
						    	nucleotide_name = s.substr(3,s.length());
						    } else
							   	line >> nucleotide_name;
					    }
			    	}
			    	if(nucleotide_name.length() > 4){
				    	string s = string(atom_name.c_str());
				    	nucleotide_name = s.substr(0,4);
				    }
				    
				    // select A if there are multiple coordinated for same monomer
									    
			    	char pdbchainid = buf[21];
				    if((chain == '-' || pdbchainid == chain)
				     //&& (( aformat == "ATOM" && (monomer_name_vs_type.count(nucleotide_name.c_str()) > 0)))
				     && (((aformat == "ATOM"||aformat == "HETATM") && (monomer_name_vs_type.count(nucleotide_name.c_str()) > 0))
				     || ((aformat == "ATOM" || aformat == "HETATM") && (nucleotide_name == "ADP" || nucleotide_name == "NHE" || nucleotide_name == "MTT" || nucleotide_name == "MN" 
				     		|| nucleotide_name == "GTP" || nucleotide_name == "GOL" || nucleotide_name == "GNP" || nucleotide_name == "DIO" || nucleotide_name == "CA" || nucleotide_name == "APC" 
				     		|| nucleotide_name == "ANP" || nucleotide_name == "AGS")))
				    ){				    	
						stringstream line(buf+22,stringstream::in);
				    	//*out << "||" << buf << endl;
						string nindex;
						line >> nindex;
						
						float x,y,z;
						line >> x;
						line >> y;
						line >> z;
							
						Nucleotide *n;	
						if(nucleotide.count(nindex.c_str()) == 0){
							nindex = *(new string(nindex));
							nucleotide_name = *(new string(nucleotide_name.c_str()));
							nucleotide[nindex.c_str()] = new Nucleotide(nindex,chain,nucleotide_name);
						}
						n = nucleotide[nindex.c_str()];
						
						if(atom_name != ""){
							float sigma, eps;	short type;
							bool ret = Atom::get_atom_details(atom_name,nucleotide_name, &sigma, &eps, &charge, &mass, &type);
							if(ret){
								atom_name = *(new string(atom_name.c_str()));
								Atom* a = new Atom(position, atom_name, n->index.c_str(), x,y,z,sigma, eps,charge,mass,type);
								//a->print_details(out);
								//*out << "ATOM " << atom->radius << " " << atom->charge << endl;
								atom[position] = a;
								n->atom[(a->name).c_str()] = a;
							} /* else
								*out << buf << endl; */
						}
						
						float bfactor;
						line >> bfactor;
						line >> bfactor;
						n->bfactor = bfactor;
				    }
			    } else if(format == "ENDMDL")
			    	break;
			}
		break;
	}
	if(to_close_fin){
		fin->close();
		fin = NULL;
	}
	
	num_atoms = atom.size();
	num_nucleotides = nucleotide.size();
	for(hash_map<const char*, Nucleotide*, hash<const char*>,eqstr>::iterator naitr = nucleotide.begin(); naitr != nucleotide.end(); naitr++){
		Nucleotide* na = (naitr->second);
		na->compute_reference_points();
	}
	*out << "#atoms " << num_atoms << " " << "#nucleotides " << num_nucleotides << endl;
	int num_pseudoatoms = 1;
	*out << "initiated nulceic acid" << endl;
}

void Complex::compute_motions(){
	Vector v = Vector(0,0,0);
	mass = charge = 0;
	for(hash_map<char, Molecule*, hash<char>, eqint>::iterator mitr = molecules.begin(); mitr != molecules.end(); mitr++){
		Molecule *m = mitr->second;
		for(hash_map<unsigned short, Atom*, hash<unsigned short>,eqint>::iterator aitr = m->atom.begin(); aitr != m->atom.end(); aitr++){
			Atom* a = aitr->second;
			v = v + *(a->position)*(a->mass);
			mass += a->mass;
			charge += a->charge;
		}
	}
	center_of_mass = new Vector(v * (1.0/mass));
	
	int count = 0;
	for(hash_map<char, Molecule*, hash<char>, eqint>::iterator mitr = molecules.begin(); mitr != molecules.end(); mitr++){
		Molecule *m = mitr->second;
		for(hash_map<unsigned short, Atom*, hash<unsigned short>,eqint>::iterator aitr = m->atom.begin(); aitr != m->atom.end(); aitr++){
			Atom* a = aitr->second;
			Vector v = *(a->position);
			float d = Vector::distance(a->position, center_of_mass) + a->radius;
			if(count++ == 0){
				max_distance_from_cm = min_distance_from_cm = d;
			}
			
			max_distance_from_cm = max(d,max_distance_from_cm);
			min_distance_from_cm = min(d,min_distance_from_cm);
		}
	}
	
	// diameter - surface to surface - computing an upper bound on d
	diameter = 0;
	unsigned short diam_a1, diam_a2;
	for(int i = 0 ; i < num_atoms; i++){
		Atom *a1 = atom[i];
		for(int j = i+1 ; j < num_atoms; j++){
			Atom *a2 = atom[j];
			float d = Vector::distance(a1->position, a2->position) + a1->radius + a2->radius;
			if(diameter < d){
				diameter = d;
				diam_a1 = i;
				diam_a2 = j;
			}
		}
	}
	center = new Vector((*(atom[diam_a1]->position) + *(atom[diam_a2]->position))*0.5);
	*out << "charge " << charge << " diameter " << diameter << " max d from cm " << max_distance_from_cm << " min " << min_distance_from_cm << endl;
}

/*
 * Require solvent accessible areas to label contacts
 */
void Complex::compute_aacontacts(){
	//cout << num_aminoacids << endl; cout.flush();
	aacontact_core = (bool **) malloc(sizeof(bool*)*num_aminoacids);
	aacontact_rim = (bool **) malloc(sizeof(bool*)*num_aminoacids);
	for(int i = 0; i < num_aminoacids; i++){
		aacontact_core[i] = (bool *) malloc(sizeof(bool)*num_aminoacids);
		aacontact_rim[i] = (bool *) malloc(sizeof(bool)*num_aminoacids);
		for(int j = 0; j < num_aminoacids; j++){
			aacontact_core[i][j] = aacontact_rim[i][j] = false;
		}
	}
		
	for(int i = 0; i < num_aminoacids; i++){
		Aminoacid *aa1 = aminoacid[i];
		//if(aa1->centroid != NULL)
			for(int j = i+3; j < num_aminoacids; j++){
				Aminoacid *aa2 = aminoacid[j];
				/*if(aa2->centroid != NULL)
					if(Vector::distance_squared(aa1->centroid,aa2->centroid) < SS_CUTOFF*SS_CUTOFF)*/
				for(hash_map<const char*, Atom*, hash<const char*>, eqstr>::iterator aitr1 = aa1->atom.begin(); aitr1 != aa1->atom.end(); aitr1++){
					Atom *a1 = (Atom*) aitr1->second;
					if((a1->name).c_str()[0] != 'H'){
						for(hash_map<const char*, Atom*, hash<const char*>, eqstr>::iterator aitr2 = aa2->atom.begin(); aitr2 != aa2->atom.end(); aitr2++){
							Atom *a2 = (Atom*) aitr2->second;
							if((a2->name).c_str()[0] != 'H'){
								double d2 = Vector::distance_squared(a1->position,a2->position);
								if(d2 < ATOM_STEPP_DMAX_SQUARED ){ //&& d2 >= ATOM_STEPP_DMIN_SQUARED ){
									if(a1->is_buried && a2->is_buried){
										aacontact_core[i][j] = aacontact_core[j][i] = true;
										aacontact_rim[i][j] = aacontact_rim[j][i] = false;
									} else if(!aacontact_core[i][j])
										aacontact_rim[i][j] = aacontact_rim[j][i] = true;
								}
							}
						}
					}
				}
			}
	}
	
	aacontact = (bool **) malloc(sizeof(bool*)*num_aminoacids);
	for(int i = 0; i < num_aminoacids; i++){
		aacontact[i] = (bool *) malloc(sizeof(bool)*num_aminoacids);
		for(int j = 0; j < num_aminoacids; j++)
			aacontact[i][j] = aacontact_core[i][j] || aacontact_rim[i][j];
	}
}

void Complex::compute_stability(){
	float net_stability = 0;
	for(int i = 0; i < num_aminoacids; i++){
	 	Aminoacid *aa1 = aminoacid[i];
		short type1 = aa1->type;
		if(aa1->centroid != NULL){
			float energy[NUM_RESIDUE_TYPES];
			for(int k = 0; k < NUM_RESIDUE_TYPES; k++)	energy[k] = 0;
			for(int j = 0; j < num_aminoacids; j++){
				if(aacontact[i][j]){
					Aminoacid *aa2 = aminoacid[j];
					short type2 = aa2->type;
					for(int k = 0; k < NUM_RESIDUE_TYPES; k++)	energy[k] += residue_potential[k][type2];
				}
			}
			// is native residue type leading to minimum energy?
			//cout << aa1->cindex << "_" << aa1->get_symbol() << " ";
			float native_energy = energy[type1], min_energy = native_energy, sum_energy=0, sum_squares_energy=0;
			for(int k = 0; k < NUM_RESIDUE_TYPES; k++){
				if(k != type1 && energy[k] < native_energy)
					min_energy = native_energy;
				sum_energy += energy[k];
				sum_squares_energy += energy[k]*energy[k];
				//cout << energy[k] << " ";
			}
			float avg_energy = sum_energy/20;
			float std_energy = sqrt(sum_squares_energy/20 - avg_energy*avg_energy);
			//cout << avg_energy << " " << std_energy << endl;
			if(sum_squares_energy > 0)	// some residues do not form contacts - for our definition of contact
				aa1->sequence_stability = (native_energy - avg_energy) / std_energy;
			else
				aa1->sequence_stability = 0;
		} else
			aa1->sequence_stability = 0;
		net_stability += aa1->sequence_stability;
	}	
	cout << "net stability " << net_stability << endl;
}

void Aminoacid::compute_pseudoatoms(){
	hash_map<const char*, Atom*, hash<const char*>, eqstr>::iterator itr, end;
	/*Vector CA_position, N_position, C_position;
	itr = atoms.begin();
	end = atoms.end();
	while(itr != end){
		Atom *a = (Atom*) *itr;
		if(a->name == "CA"){
			CA_position = Vector(a->position);
		}
		itr++;
	}
	itr = atoms.begin();
	end = atoms.end();
	while(itr != end){
		Atom *a = (Atom*) *itr;
		if(a->name == "N" && (Vector::distance(CA_position, *(a->position))<2)){
			N_position = Vector(a->position);
		} else if(a->name == "C" && (Vector::distance(CA_position, *(a->position))<2)){
			C_position = Vector(a->position);
		}
		itr++;
	}

	// compute the reference frame
	Vector ez = N_position - CA_position;
	ez.normalize();
	*/
	
	Atom* pa1,*pa2,*pa3;
	float mass = 0, charge =0, radius=0;
	int num_atoms=0; 
	Vector position = Vector(0,0,0);
	itr = atom.begin();
	end = atom.end();
	while(itr != end){
		Atom *a = (Atom*) itr->second;
		if((type == -1) || a->name == "CA" || a->name == "N" || a->name == "H" || a->name == "O" || a->name == "CB"){  
			num_atoms++;
			mass += a->mass;
			charge += (a->mass)* (a->charge);
			radius += (a->mass)* (a->radius);
			position = position + *(a->position);
		}
		itr++;
	}
	mass = mass/num_atoms;
	charge = charge/mass;
	radius = radius/mass;
	position = position * (1.0/mass);
	pa1 = new Atom(-1, "",string("").c_str(), position.x, position.y, position.z, 2*radius, 1.0, charge, mass, -1);
	pseudoatoms.push_back(pa1);
	
	if(!(type == -1 || type == ALA || type == GLY) ){
		mass = 0; charge =0; radius=0;
		num_atoms=0; 
		position = Vector(0,0,0);
		float mass3 = 0, charge3 =0, radius3=0;
		int num_atoms3=0; 
		Vector position3 = Vector(0,0,0);
		itr = atom.begin();
		end = atom.end();
		while(itr != end){
			Atom *a = (Atom*) itr->second;
			if((a->name).c_str()[1] == 'G' || (a->name).c_str()[1] == 'D'){  
				num_atoms++;
				mass += a->mass;
				charge += (a->mass)* (a->charge);
				radius += (a->mass)* (a->radius);
				position = position + *(a->position);
			} else {
				num_atoms3++;
				mass3 += a->mass;
				charge3 += (a->mass)* (a->charge);
				radius3 += (a->mass)* (a->radius);
				position3 = position3 + *(a->position);
			}
			itr++;
		}
		mass = mass/num_atoms;
		charge = charge/mass;
		radius = radius/mass;
		position = position * (1.0/mass);
		pa2 = new Atom(-1, "",string("").c_str(), position.x, position.y, position.z, 2*radius,1.0, charge, mass, -1);
		pseudoatoms.push_back(pa2);
		if(num_atoms3 > 0){
			mass3 = mass3/num_atoms3;
			charge3 = charge3/mass3;
			radius3 = radius3/mass3;
			position3 = position3 * (1.0/mass3);
			pa3 = new Atom(-1, "",string("").c_str(), position3.x, position3.y, position3.z, 2*radius3,1.0, charge3, mass3, -1);
			pseudoatoms.push_back(pa3);
		}
	}
	
	//case ASN : case ASP : case CYS : case ILE : case LEU : case PRO : case SER : case THR : case VAL :
	//case ARG : case GLN : case GLU : case HIS : case LYS : case MET : case PHE : case TRP : case TYR :
	
	num_pseudoatoms = pseudoatoms.size();
	//*out << "aminoacid " << index << " #pseudoatoms " << num_pseudoatoms << endl;
	
	/*if(num_pseudoatoms >= 2){
		*out << index << " 1-2 " << type << " " << pa1->radius << " " << pa2->radius << " " << Vector::distance(pa1->position,pa2->position) << " " << pa1->radius + pa2->radius - Vector::distance(pa1->position,pa2->position)  << endl;
	}
	if(num_pseudoatoms >= 3){
		*out << index << " 2-3 " << type << " " << pa2->radius << " " << pa3->radius << " "<< " " << Vector::distance(pa3->position,pa2->position) << " " << pa2->radius + pa3->radius - Vector::distance(pa3->position,pa2->position)  << endl;
	}*/
}


void Complex::compute_surface(){
	float distance[num_atoms][num_atoms];
	for(int i = 0 ; i < num_atoms; i++){
		Atom* a1 = atom[i];
		for(int j = i+1; j < num_atoms; j++){
			Atom* a2 = atom[j];
			distance[i][j] = distance[j][i] = Vector::distance((a1->position), (a2->position));

			a1->is_buried |= (distance[i][j] + a1->radius < a2->radius);
			a2->is_buried |= (distance[i][j] + a2->radius < a1->radius);
		}
	}
	
	for(int i = 0 ; i < num_atoms; i++){
		Atom* a1 = atom[i];
		if(!a1->is_buried && a1->mass > 1){
			for(int j = i+1; j < num_atoms; j++){
				Atom* a2 = atom[j];
				if(!a2->is_buried && a2->mass > 1){
					// pairs of atoms forming tori
					if(distance[i][j] < a1->radius + a2->radius + 2*probe_radius){
						//*out << a1->cindex << " " << a2->cindex << endl;
						Torus* torus = new Torus(a1,a2,distance[i][j]);
						toruses[(a1->cindex)*MAX_ATOMS + (a2->cindex)] = torus;
										
						// update neighbors
						a1->neighbors.push_back(a2);
						a2->neighbors.push_back(a1);
					}
				}
			}
		}
	}
	
	// sort the neighbors of atoms intersection assumes that the neighbors are sorted by index
	for(int i = 0 ; i < num_atoms; i++){
		Atom* a = atom[i];
		sort(a->neighbors.begin(), a->neighbors.end(),less<Atom*>());
		/*out << a->cindex << "- ";
		for(vector<Atom*>::iterator nitr = a->neighbors.begin(); nitr != a->neighbors.end(); nitr++)
			*out << ((Atom*) *nitr)->cindex << " ";
		*out << endl;*/
	}
	
	// for each torus compute the triangles it is involved in
	for(hash_map<const long,Torus*,hash<long>,eqlong>::iterator titr = toruses.begin(); titr != toruses.end(); titr++){
		Torus* t12 = titr->second;
		Atom* a1 = t12->a1;	
		Atom* a2 = t12->a2;

		//*out << "Examining Torus " << a1->cindex << " " << a2->cindex << endl; out->flush();	
		// find mutual neigbhors of the atoms involved
		vector<Atom *> all_common_neighbors = Atom::get_common_neighbors(a1,a2);
		vector<Atom *> common_neighbors = Atom::get_common_neighbors(a1,a2);
		
		/*for(vector<Atom *>::iterator nitr = common_neighbors.begin(); nitr != common_neighbors.end(); nitr++){
			*out << ((Atom*) *nitr)->cindex << " ";
		}
		*out << endl;*/
		
		// ensure that a1->cindex < a2->cindex < a3->cindex
		// uses the fact that the common neighbors are sorted by index
		vector<Atom *>::iterator itr = common_neighbors.begin();
		while( (itr != common_neighbors.end()) && ((*itr)->cindex < a2->cindex) ){
			itr++;
		}
		common_neighbors.erase(common_neighbors.begin(),itr);
		
		for(itr = common_neighbors.begin();itr != common_neighbors.end(); itr++){
			Atom* a3 = *itr;

			Torus* t13 = toruses[a1->cindex * MAX_ATOMS + a3->cindex];
			Torus* t23 = toruses[a2->cindex * MAX_ATOMS + a3->cindex];
			Vector* u12 = t12->axis;
			Vector* u13 = t13->axis;
			float triangle_angle = acos(u12->dot(u13));
			//*out << a1->cindex << " " << a2->cindex << " " << a3->cindex << " " << triangle_angle << "\t";
			if(triangle_angle != 0 && u12->dot(u13) != -1){ // check for colinearity
				Vector base_plane_normal = u12->cross(u13);
				base_plane_normal.normalize();
				Vector base_point = Vector( *(t12->center) + (base_plane_normal.cross(u12)) * ( (*(t13->center) - *(t12->center)).dot(u13)/sin(triangle_angle) ) );
				float probe_height_squared = (a1->radius + probe_radius)*(a1->radius + probe_radius) - (base_point - *(a1->position)).norm() * (base_point - *(a1->position)).norm();
				//*out << " probe height squared " << probe_height_squared << endl;
				
				if(probe_height_squared > 0){
					float probe_height = sqrt(probe_height_squared);
					
					t12->is_free = false;	
					t13->is_free = false;	
					t23->is_free = false;	
					
					// check if the probe position causes collisions
					vector<Atom *> neighbors_a1a2a3 = a3->get_neighbors_in(all_common_neighbors);
					
					bool can_place_probe_above = true, can_place_probe_below = true;
					Vector *probe_position_above = new Vector(base_point + (base_plane_normal * probe_height));
					Vector *probe_position_below = new Vector(base_point - (base_plane_normal * probe_height));	
					/*out << "probe position? " << (Vector::distance(probe_position_above, a1->position) - (a1->radius + probe_radius)) << " " 
						<< (Vector::distance(probe_position_above, a2->position) - (a2->radius + probe_radius)) << " " 
						<< (Vector::distance(probe_position_above, a3->position) - (a3->radius + probe_radius)) << endl;*/
					
					if(neighbors_a1a2a3.size() > 0){
						// check for collision with each atom
						for(vector<Atom *>::iterator itr1 = neighbors_a1a2a3.begin(); itr1 != neighbors_a1a2a3.end(); itr1++){
							Atom* a = *itr1;
							if(can_place_probe_above && Vector::distance(probe_position_above, (a->position)) < a->radius + probe_radius){
								can_place_probe_above = false;
								//*out << "lost above " << a3->cindex << " intersects " << a->cindex << " " << (Vector::distance(probe_position_above, a->position) - (a->radius + probe_radius)) << endl;
								delete probe_position_above;
							}
							if(can_place_probe_below && Vector::distance(probe_position_below, (a->position)) < a->radius + probe_radius){
								can_place_probe_below = false;
								//*out << "lost below " << a3->cindex << " intersects " << a->cindex << " " << (Vector::distance(probe_position_below, a->position) - (a->radius + probe_radius))<< endl;
								delete probe_position_below;
							}
						}
					}
					
					Vector *v1, *v2, *v3;
					if(can_place_probe_above){
						//*out << "obtained triangle above " << a1->cindex << " " << a2->cindex << " " << a3->cindex << endl;
						
						v1 = new Vector(*(a1->position) * ((probe_radius)/(a1->radius + probe_radius)) + *(probe_position_above) * ((a1->radius)/(a1->radius + probe_radius)));
						v2 = new Vector(*(a2->position) * ((probe_radius)/(a2->radius + probe_radius)) + *(probe_position_above) * ((a2->radius)/(a2->radius + probe_radius)));
						v3 = new Vector(*(a3->position) * ((probe_radius)/(a3->radius + probe_radius)) + *(probe_position_above) * ((a3->radius)/(a3->radius + probe_radius)));
						// direct the edges such that they are counter clockwise when seen from the probe position
						float orientation = base_plane_normal.dot( (*v2 - *v1).cross(*v3 - *v1) );
						Triangle_Face* triangle;
						if(orientation > 0)
							triangle = new Triangle_Face(a1,a2,a3,v1,v2,v3, probe_position_above);
						else
							triangle = new Triangle_Face(a1,a3,a2,v1,v3,v2, probe_position_above);
						
						//triangle->print_details();
						
						trianglefaces[2*((a1->cindex * MAX_ATOMS + a2->cindex)*MAX_ATOMS + a3->cindex)] = triangle;
						t12->trianglefaces.push_back(triangle);
						t13->trianglefaces.push_back(triangle);
						t23->trianglefaces.push_back(triangle);
					}

					if(can_place_probe_below){
						//*out << "obtained triangle below " << a1->cindex << " " << a2->cindex << " " << a3->cindex << endl;
						
						v1 = new Vector(*(a1->position) * ((probe_radius)/(a1->radius + probe_radius)) + *(probe_position_below) * ((a1->radius)/(a1->radius + probe_radius)));
                        v2 = new Vector(*(a2->position) * ((probe_radius)/(a2->radius + probe_radius)) + *(probe_position_below) * ((a2->radius)/(a2->radius + probe_radius)));
                        v3 = new Vector(*(a3->position) * ((probe_radius)/(a3->radius + probe_radius)) + *(probe_position_below) * ((a3->radius)/(a3->radius + probe_radius)));
                        // direct the edges such that they are counter clockwise when seen from the probe position
                        float orientation = - base_plane_normal.dot( (*v2 - *v1).cross(*v3 - *v1) );					
						Triangle_Face* triangle;
						if(orientation > 0)
							triangle = new Triangle_Face(a1,a2,a3,v1,v2,v3, probe_position_below);
						else
							triangle = new Triangle_Face(a1,a3,a2,v1,v3,v2, probe_position_below);

						//triangle->print_details();
						
						trianglefaces[2*((a1->cindex * MAX_ATOMS + a2->cindex)*MAX_ATOMS + a3->cindex)+1] = triangle;
						t12->trianglefaces.push_back(triangle);
						t13->trianglefaces.push_back(triangle);
						t23->trianglefaces.push_back(triangle);
					}
				} else {
					// find out if the torii are buried
					Vector u = (*t12->center - *a3->position) - *t12->axis * (*t12->center - *a3->position).dot(*t12->axis);
					u.normalize();
					if((*t12->center + u * t12->radius - *a3->position).norm() <= a3->radius + probe_radius
					 && (*t12->ccc1 + u * t12->ccr1 - *a3->position).norm() <= a3->radius + probe_radius
					 && (*t12->ccc2 + u * t12->ccr2 - *a3->position).norm() <= a3->radius + probe_radius){
						//*out << "Trous buried " << t12->a2->cindex << " " << t12->a2->cindex << endl;
						t12->is_free = false;
						t12->is_buried = true;
					}

					u = (*t23->center - *a1->position) - *t23->axis * (*t23->center - *a1->position).dot(*t23->axis);
					u.normalize();
					if((*t23->center + u * t23->radius - *a1->position).norm() <= a1->radius + probe_radius
					 && (*t23->ccc1 + u * t23->ccr1 - *a1->position).norm() <= a1->radius + probe_radius
					 && (*t23->ccc2 + u * t23->ccr2 - *a1->position).norm() <= a1->radius + probe_radius){
						//*out << "Trous buried " << t23->a1->cindex << " " << t23->a2->cindex << endl;
						t23->is_free = false;
						t23->is_buried = true;
					}

					u = (*t13->center - *a2->position) - *t13->axis * (*t13->center - *a2->position).dot(*t13->axis);
					u.normalize();
					if((*t13->center + u * t13->radius - *a2->position).norm() <= a2->radius + probe_radius
					 && (*t13->ccc1 + u * t13->ccr1 - *a2->position).norm() <= a2->radius + probe_radius
					 && (*t13->ccc2 + u * t13->ccr2 - *a2->position).norm() <= a2->radius + probe_radius){
						//*out << "Trous buried " << t13->a1->cindex << " " << t13->a2->cindex << endl;
						t13->is_free = false;
						t13->is_buried = true;
					}
				} // probe_height_squred
			}
		}
	}
	*out << "#toruses " << toruses.size() << endl;
	
	*out << "#trianglefaces " << trianglefaces.size() << endl;

	//compute saddle faces
	for(hash_map<const long,Torus*,hash<long>,eqlong>::iterator titr = toruses.begin(); titr != toruses.end(); titr++){
		Torus* t = titr->second;
		if(!t->is_buried)//t->trianglefaces.size() > 0)
			t->build_saddle_faces();
		if(!t->is_free)
			t->is_buried = ((t->saddles).size() == 0);
	}

	//compute convex faces
    int num_buried = 0, nonH = 0;
    for(int i = 0 ; i < num_atoms; i++){
   	  	Atom* a = atom[i];
   	  	if(!a->is_buried  && a->mass > 1){
   	  		a->is_buried = true;
	   	  	for(vector<Face*>::iterator titr = a->toruses.begin(); titr != a->toruses.end(); titr++){
    	    	Torus* t = (Torus *) *titr;
        		(a->is_buried) &= (t->is_buried);
   	  		}
   	  	}
   	  	//a->print_details();
   	  	if(!a->is_buried && a->mass > 1)
			a->build_convex_faces();
		else
			num_buried++;
		if(a->mass > 1)
			nonH++;
	}
	*out << "# atoms " << num_atoms << " #nonH " << nonH << " #buried " << num_buried << endl;
}

/*
 * Left_face and right_face information set for all edges
 * 
 * Assumes that the complex is closed - the expsoed surface is connected
 * 
 * Right now triangulating atom is not by triangulating each convex face - so knowing the burial of convex faces
 * is not useful - so need to eliminate triangles in the holes after the triangulation
 */
void Complex::eliminate_faces_in_holes(){
	hash_map<unsigned int,Face*,hash<long>,eqlong> face;
	Graph* face_graph = new Graph();
	hash_map<long,GVertex*,hash<long>,eqlong> face_graph_node;
	unsigned int num_faces = 0;
	GVertex *gv, *gu;
	
	for(hash_map<const long,Torus*,hash<long>,eqlong>::iterator titr = toruses.begin(); titr != toruses.end(); titr++){
		Torus* t = titr->second;
		if(t->is_free){
			{
				t->face_number = num_faces;
				gv = new GVertex(num_faces);
				gv->data = (Face*) t;
				face_graph->insert_vertex(gv);
				face_graph_node[num_faces] = gv;

				num_faces++;
			}
		} else if(!t->is_buried){
			for(vector<Saddle*>::iterator sitr = t->saddles.begin(); sitr != t->saddles.end(); sitr++){
				Saddle *s = *sitr;
				{
					s->face_number = num_faces;
					gv = new GVertex(num_faces);
					gv->data = (Face*) s;
					face_graph->insert_vertex(gv);
					face_graph_node[num_faces] = gv;
	
					num_faces++;
				}
			}
		}
	}
	
	for(hash_map<const long, Triangle_Face *, hash<long>,eqlong>::iterator tfitr = trianglefaces.begin(); tfitr != trianglefaces.end(); tfitr++){
		Triangle_Face* tf = tfitr->second;
		
		{
			tf->face_number = num_faces;
			gv = new GVertex(num_faces);
			gv->data = (Face*) tf;
			face_graph->insert_vertex(gv);
			face_graph_node[num_faces] = gv;
		
			num_faces++;
		}
	}
	
	for(int i = 0 ; i < num_atoms; i++){
   		Atom* a = atom[i];
   		if(!a->is_buried  && a->mass > 1){
   			for(vector<Face*>::iterator cfitr = a->convexfaces.begin(); cfitr != a->convexfaces.end(); cfitr++){
				Convex_Face *cf = (Convex_Face*) *cfitr;
				
				{
					cf->face_number = num_faces;
					gv = new GVertex(num_faces);
					gv->data = (Face*) cf;
					face_graph->insert_vertex(gv);
					face_graph_node[num_faces] = gv;
				
					num_faces++;
				}
			}
   		}
	}
	 
	// add the neighbor information
	for(hash_map<const long, Triangle_Face *, hash<long>,eqlong>::iterator tfitr = trianglefaces.begin(); tfitr != trianglefaces.end(); tfitr++){
		Triangle_Face* tf = tfitr->second;
		
		gu = face_graph_node[tf->face_number];
		
		gv = face_graph_node[tf->e12->right_face->face_number];
		gu->add_neighbor(gv);
		gv->add_neighbor(gu);
		
		gv = face_graph_node[tf->e23->right_face->face_number];
		gu->add_neighbor(gv);
		gv->add_neighbor(gu);
		
		gv = face_graph_node[tf->e31->right_face->face_number];
		gu->add_neighbor(gv);
		gv->add_neighbor(gu);
	}
		
	for(int i = 0 ; i < num_atoms; i++){
   		Atom* a = atom[i];
   		if(!a->is_buried  && a->mass > 1){
   			for(vector<Face*>::iterator cfitr = a->convexfaces.begin(); cfitr != a->convexfaces.end(); cfitr++){
				Convex_Face *cf = (Convex_Face*) *cfitr;
				gu = face_graph_node[cf->face_number];
				for(vector<Cycle*>::iterator cycleitr = (cf->cycles).begin(); cycleitr != (cf->cycles).end(); cycleitr++){
					Cycle *c = *cycleitr;
					for(list<Edge*>::iterator eitr = (c->edges).begin(); eitr != (c->edges).end(); eitr++){
						Edge *e = *eitr;
						gv = face_graph_node[e->right_face->face_number];
						gu->add_neighbor(gv);
						gv->add_neighbor(gu);
					}
				}
   			}
   		}
	}

	// find a face that is definitely exposed using an approximation
	Face *exposed_face = NULL;
	float miny;
	for(hash_map<const long, Triangle_Face *, hash<long>,eqlong>::iterator tfitr = trianglefaces.begin(); tfitr != trianglefaces.end(); tfitr++){
		Triangle_Face* tf = tfitr->second;
		if(exposed_face == NULL || tf->v1->y < miny || tf->v2->y < miny || tf->v3->y < miny){
			exposed_face = (Face*) tf;
			miny = minimum(miny,tf->v1->y);
			miny = minimum(miny,tf->v2->y);
			miny = minimum(miny,tf->v3->y);
		}
	}
		
	// compute the list of faces accessible from outside	
	list<long> connected_component;
	face_graph->compute_connected_component(&connected_component,exposed_face->face_number);
	*out << "#faces " << num_faces << " " << face_graph->num_vertices << " " << connected_component.size() << endl;
	
	int num_convex_faces=0, num_torus_faces=0, num_saddle_faces=0, num_triangle_faces=0;
	for(hash_map<long,GVertex*,hash<long>,eqlong>::iterator fvitr = face_graph_node.begin(); fvitr != face_graph_node.end(); fvitr++){
		Face* f = ((Face*) fvitr->second->data);
		f->is_buried = true;
		switch(f->type){
			case TORUS_FACE:
				num_torus_faces++;
				break;
			case SADDLE_FACE:
				num_saddle_faces++;
				break;
			case TRIANGLE_FACE:
				num_triangle_faces++;
				break;
			default:
				num_convex_faces++;
				break;
		}
	}
	
	for(list<long>::iterator itr = connected_component.begin(); itr != connected_component.end(); itr++){
		Face* f = ((Face*) face_graph_node[*itr]->data);
		f->is_buried = false;
		switch(f->type){
			case TORUS_FACE:
				num_torus_faces--;
				break;
			case SADDLE_FACE:
				num_saddle_faces--;
				break;
			case TRIANGLE_FACE:
				num_triangle_faces--;
				break;
			default:
				num_convex_faces--;
				break;
		}
	}
	*out << "Buried by type torus:" << num_torus_faces << " saddle:" << num_saddle_faces << " triangle:" << num_triangle_faces << " convex:" << num_convex_faces << endl;

	for(hash_map<const long,Torus*,hash<long>,eqlong>::iterator titr = toruses.begin(); titr != toruses.end(); titr++){
		Torus* t = titr->second;
		if(!t->is_buried){
			t->is_buried = true;
			for(vector<Saddle*>::iterator sitr = t->saddles.begin(); sitr != t->saddles.end(); sitr++){
				Saddle *s = *sitr;
				t->is_buried &= s->is_buried;
			}
		}
	}
	
	for(int i = 0 ; i < num_atoms; i++){
   	  	Atom* a = atom[i];
   	  	if(!a->is_buried  && a->mass > 1){
   	  		a->is_buried = true;
	   	  	for(vector<Face*>::iterator titr = a->toruses.begin(); titr != a->toruses.end(); titr++){
    	    	Torus* t = (Torus *) *titr;
        		(a->is_buried) &= (t->is_buried);
   	  		}
   	  	}
   	  	if(!a->is_buried){
   	  		a->is_buried = true;
	   	  	for(vector<Face*>::iterator cfitr = a->convexfaces.begin(); cfitr != a->convexfaces.end(); cfitr++){
				Convex_Face *cf = (Convex_Face*) *cfitr;
				a->is_buried &= cf->is_buried;
	   	  	}
   	  	}
   	  	if(!a->is_buried){
   	  		*out << a->index << "-" << a->convexfaces.size() << "(";
   	  		for(vector<Face*>::iterator cfitr = a->convexfaces.begin(); cfitr != a->convexfaces.end(); cfitr++){
				Convex_Face *cf = (Convex_Face*) *cfitr;
				if(!cf->is_buried) *out << (cf->cycles).size() <<",";
   	  		}
   	  		*out << ") ";
   	  	}
	}
	*out << endl;
}

/*
 * Job is to eliminate triangles that are not exposed to the solvent
 * need to remove them triangles and atom triangles
 * 
 * Does not work, the edges in the graph are incomplete
 */
void Complex::eliminate_triangles_in_holes(){
	hash_map<long,Vector*,hash<long>,eqlong> vertex;
	hash_map<const char*,long,hash<const char*>,eqstr> coord_vs_vertex;
	Graph* triangle_graph = new Graph();
	hash_map<long,GVertex*,hash<long>,eqlong> triangle_graph_node;
	hash_map<long,hash_set<long,hash<long>,eqlong>,hash<long>,eqlong> vertex_vs_triangle;
	
	long num_vertices = 0;
	long num_triangles = 0;
	
	for(vector<Triangle *>::iterator itr = triangles.begin(); itr != triangles.end(); itr++){
		Triangle *t = *itr;
		long vi1,vi2,vi3;
		vi1 = vi2 = vi3 = -1;
		string vdetails;
		stringstream line1,line2,line3;
		Vector v = *(t->v1);
		line1 << "(" << v.x << "," << v.y << "," << v.z << ")";
		vdetails = line1.str();
		if(coord_vs_vertex.count(vdetails.c_str()) > 0){
			vi1 = coord_vs_vertex[vdetails.c_str()];
		} else {
			vi1 = num_vertices;
			vertex[num_vertices] = (t->v1);
			coord_vs_vertex[(new string(vdetails))->c_str()] = vi1;
			vertex_vs_triangle[vi1] = *(new hash_set<long,hash<long>,eqlong>);
			
			num_vertices++;
		}

		v = *(t->v2);
		line2 << "(" << v.x << "," << v.y << "," << v.z << ")";
		vdetails = line2.str();
		if(coord_vs_vertex.count(vdetails.c_str()) > 0){
			vi2 = coord_vs_vertex[vdetails.c_str()];
		} else {
			vi2 = num_vertices;
			vertex[num_vertices] = (t->v2);
			coord_vs_vertex[(new string(vdetails))->c_str()] = vi2;
			vertex_vs_triangle[vi2] = *(new hash_set<long,hash<long>,eqlong>);
			
			num_vertices++;
		}

		v = *(t->v3);
		line3 << "(" << v.x << "," << v.y << "," << v.z << ")";
		vdetails = line3.str();
		if(coord_vs_vertex.count(vdetails.c_str()) > 0){
			vi3 = coord_vs_vertex[vdetails.c_str()];
		} else {
			vi3 = num_vertices;
			vertex[num_vertices] = (t->v3);
			coord_vs_vertex[(new string(vdetails))->c_str()] = vi3;
			vertex_vs_triangle[vi3] = *(new hash_set<long,hash<long>,eqlong>);
			
			num_vertices++;
		}
		
		{
			GVertex *gv = new GVertex(num_triangles);
			gv->data = t;
			triangle_graph->insert_vertex(gv);
			triangle_graph_node[num_triangles] = gv;
		
			vertex_vs_triangle[vi1].insert(num_triangles);
			vertex_vs_triangle[vi2].insert(num_triangles);
			vertex_vs_triangle[vi3].insert(num_triangles);
		
			num_triangles++;
		}
	}	
	
	// compute neighbors of each triangle
	for(hash_map<long,GVertex*,hash<long>,eqlong>::iterator itr = triangle_graph_node.begin(); itr != triangle_graph_node.end(); itr++){
		long tindex = itr->first;
		Triangle *t = ((Triangle*) ((GVertex*) itr->second)->data);
		long vi1,vi2,vi3;
		string vdetails;
		stringstream line1,line2,line3;
		
		Vector v = *(t->v1);
		line1 << "(" << v.x << "," << v.y << "," << v.z << ")";
		vdetails = line1.str();
		vi1 = coord_vs_vertex[vdetails.c_str()];

		v = *(t->v2);
		line2 << "(" << v.x << "," << v.y << "," << v.z << ")";
		vdetails = line2.str();
		vi2 = coord_vs_vertex[vdetails.c_str()];

		v = *(t->v3);
		line3 << "(" << v.x << "," << v.y << "," << v.z << ")";
		vdetails = line3.str();
		vi3 = coord_vs_vertex[vdetails.c_str()];
		
		//find common neighbors of vi1 and vi2
		hash_set<long,hash<long>,eqlong> ts1,ts2,ts3;
		ts1 = vertex_vs_triangle[vi1];
		ts2 = vertex_vs_triangle[vi2];
		ts3 = vertex_vs_triangle[vi3];
		for(hash_set<long,hash<long>,eqlong>::iterator titr = ts1.begin(); titr != ts1.end(); titr++){
			long ntindex = *titr;
			if(ntindex != tindex && (ts2.count(ntindex) >0 || ts3.count(ntindex) > 0))
				triangle_graph_node[tindex]->add_neighbor(triangle_graph_node[ntindex]);
		}
		for(hash_set<long,hash<long>,eqlong>::iterator titr = ts2.begin(); titr != ts2.end(); titr++){
			long ntindex = *titr;
			if(ntindex != tindex && ts3.count(ntindex) > 0)
				triangle_graph_node[tindex]->add_neighbor(triangle_graph_node[ntindex]);
		}
	}
	
	// find a triangle that is definitely exposed
	Triangle *exposed_triangle = NULL;
	long exposed_triangle_index;
	/*for(hash_map<const long, Triangle_Face *, hash<long>,eqlong>::iterator tfitr = trianglefaces.begin(); tfitr != trianglefaces.end(); tfitr++){
		Triangle_Face* tf = tfitr->second;
		if((exposed_triangle == NULL || triangle_graph_node[exposed_triangle_index]->get_degree() == 0) && !tf->is_buried){
			exposed_triangle = (Triangle*) tf;
			{
				long vi1,vi2,vi3;
				string vdetails;
				stringstream line1,line2,line3;
				
				Vector v = *(exposed_triangle->v1);
				line1 << "(" << v.x << "," << v.y << "," << v.z << ")";
				vdetails = line1.str();
				vi1 = coord_vs_vertex[vdetails.c_str()];
		
				v = *(exposed_triangle->v2);
				line2 << "(" << v.x << "," << v.y << "," << v.z << ")";
				vdetails = line2.str();
				vi2 = coord_vs_vertex[vdetails.c_str()];
		
				v = *(exposed_triangle->v3);
				line3 << "(" << v.x << "," << v.y << "," << v.z << ")";
				vdetails = line3.str();
				vi3 = coord_vs_vertex[vdetails.c_str()];
				
				hash_set<long,hash<long>,eqlong> ts1,ts2,ts3;
				ts1 = vertex_vs_triangle[vi1];
				ts2 = vertex_vs_triangle[vi2];
				ts3 = vertex_vs_triangle[vi3];
				for(hash_set<long,hash<long>,eqlong>::iterator titr = ts1.begin(); titr != ts1.end(); titr++){
					long ntindex = *titr;
					if(ts2.count(ntindex) > 0 && ts3.count(ntindex) > 0){
						exposed_triangle_index = ntindex;
						break;
					}
				}
			}
		}
	}*/
	Vector min = ((Vector*) vertex.begin()->second);
	for(hash_map<long,Vector*,hash<long>,eqlong>::iterator itr = vertex.begin(); itr != vertex.end(); itr++){
		Vector *v = itr->second;
		if(v->x < min.x)	min = *v;
	}
	stringstream line;
	line << "(" << min.x << "," << min.y << "," << min.z << ")";
	string vdetails = line.str();
	long vminindex = coord_vs_vertex[vdetails.c_str()];
	exposed_triangle_index = *((vertex_vs_triangle[vminindex]).begin());
	
	*out << "exposed triangle " << exposed_triangle_index << " degree " << triangle_graph_node[exposed_triangle_index]->get_degree() << endl;
	
	list<long> connected_component;
	triangle_graph->compute_connected_component(&connected_component,exposed_triangle_index);
	*out << "#triangles " << triangles.size() << " " << triangle_graph->num_vertices << " " << connected_component.size() << endl;
	
	bool connected[num_triangles];
	for(int i = 0; i < num_triangles; i++){
		connected[i] = false;
		((Triangle*) triangle_graph_node[i]->data)->is_buried = true;
	}
	for(list<long>::iterator itr = connected_component.begin(); itr != connected_component.end(); itr++){
		connected[*itr] = true;
		((Triangle*) triangle_graph_node[*itr]->data)->is_buried = false;
	}
}


Complex::Complex(Atom **atom_details, float* radius, float* eps, int count_atoms, int max_atoms){
	//*out << "new complex " << count_atoms << " ";
	num_atoms = count_atoms;
	
	atom = (Atom **) malloc(sizeof(Atom*)*(max_atoms+1));
	for(int i = 0 ; i < max_atoms; i++)
		atom[i] = NULL;
		
	for(int i = 0 ; i < count_atoms; i++){
		//*out << atom_details[i]->cindex << " ";
		Atom *a = atom_details[i];
		Atom* an = new Atom(a->index, a->name, a->monoindex, a->position->x,a->position->y,a->position->z, 2*radius[i], eps[i], a->charge, a->mass, a->atom_type);
		atom[a->cindex] = an;
		an->cindex = a->cindex;
	}
	//*out << "done" << endl; out->flush();
}

Molecule::~Molecule(){
}

Atom::Atom(){
}

Atom::~Atom(){
	delete position;
}

/*
 * Check if this face contains cf
 */
bool Convex_Face::contains(Convex_Face *cf){
	bool result = true;
	for(vector<Cycle*>::iterator ocitr = cycles.begin(); ocitr != cycles.end(); ocitr++){
		for(vector<Cycle*>::iterator icitr = cf->cycles.begin(); icitr != cf->cycles.end(); icitr++){
			result &= ((Cycle*) *ocitr)->contains((Cycle*) *icitr);
		}
	}
	return result;
}

// if circle edges disappear, then we should not look for the nearest in the inner vertices -- use the center of the cycle instead
// cannot add the center arbitrarily, the center might be closer than the vertices on the inner cycle
void Convex_Face::triangulate(vector<Triangle*> *triangles){
	vector<Edge*> edges;
	refine_edges(&edges);
	
	vector<Vector*> in_vertices;
	
	Cycle *cycle[cycles.size()];
	Vector *cyclecenter[cycles.size()];
	float cyclesize[cycles.size()];
	if(refinement.size() == 0){
		int cycle_index = 0;

		for(vector<Cycle *>::iterator citr = cycles.begin(); citr != cycles.end(); citr++){
			Cycle *c = *citr;
			int num_cycle_vertices = 0;
			Vector vc = Vector(0,0,0);
			for(list<Edge *>::iterator eitr = (c->edges).begin(); eitr != (c->edges).end(); eitr++){
				Edge *e = *eitr;
				for(vector<Edge*>::iterator eritr = e->refinement.begin(); eritr != e->refinement.end(); eritr++){
					Edge *er = *eritr;
					vc = vc + *er->start;
					num_cycle_vertices++;
				}
			}
			vc = vc * (1.0/num_cycle_vertices);
			vc = *(atom->position) + (vc - *(atom->position))*(atom->radius/(vc - *(atom->position)).norm());
			cyclecenter[cycle_index] = new Vector(vc);
			cyclesize[cycle_index] = num_cycle_vertices;
			cycle[cycle_index] = c;
			cycle_index++;
		}
		
		if(cycles.size() == 1)
			in_vertices.push_back(cyclecenter[0]);
	} else {	
		for(vector<Convex_Face*>::iterator icfitr = refinement.begin(); icfitr != refinement.end(); icfitr++){
			Convex_Face *icf = (Convex_Face*) *icfitr;
			icf->triangulate(triangles);
			for(vector<Cycle *>::iterator citr = icf->cycles.begin(); citr != icf->cycles.end(); citr++){
				Cycle *c = *citr;
				for(list<Edge *>::iterator eitr = (c->edges).begin(); eitr != (c->edges).end(); eitr++){
					Edge *e = *eitr;
					for(vector<Edge*>::iterator eritr = e->refinement.begin(); eritr != e->refinement.end(); eritr++){
						Edge *er = *eritr;
						in_vertices.push_back(er->start);
					}
				}
			}
		}
	}
	
	float max_distance = -1, min_distance = 2*radius;
	if(refinement.size() > 0 || cycles.size() == 1){
		for(vector<Cycle *>::iterator citr = cycles.begin(); citr != cycles.end(); citr++){
			Cycle *c = *citr;
			for(list<Edge *>::iterator eitr = (c->edges).begin(); eitr != (c->edges).end(); eitr++){
				Edge *e = *eitr;
				for(vector<Edge*>::iterator eritr = e->refinement.begin(); eritr != e->refinement.end(); eritr++){
					Edge *er = *eritr;
					Vector *v1 = Vector::find_closest(&in_vertices,er->start);
					Vector *v2 = Vector::find_closest(&in_vertices,er->end);
					if(v1 != NULL){
						float distance = Vector::distance(v1,er->start);
						if(distance > max_distance)
							max_distance = distance;
						else if(distance < min_distance)
							min_distance = distance;									
						triangles->push_back(new Triangle(v1,er->start,er->end,atom->position,atom->radius,CONVEX));
						if(v2 != NULL && !(*v1 == *v2))
							triangles->push_back(new Triangle(v1,er->end,v2,atom->position,atom->radius,CONVEX));
					}
				}
			}
		}
	} else {
		// more than one cycle and no refinement
		Cycle *largest_cycle;
		float max_size = 0.0;
		for(int i = 0; i < cycles.size(); i++){
			if(cyclesize[i] > max_size){
				max_size = cyclesize[i];
				largest_cycle = cycle[i];
			}
		}
		
		for(vector<Cycle *>::iterator citr = cycles.begin(); citr != cycles.end(); citr++){
			Cycle *c = *citr;
			if(c != largest_cycle){
				for(list<Edge *>::iterator eitr = (c->edges).begin(); eitr != (c->edges).end(); eitr++){
					Edge *e = *eitr;
					for(vector<Edge*>::iterator eritr = e->refinement.begin(); eritr != e->refinement.end(); eritr++){
						Edge *er = *eritr;
						in_vertices.push_back(er->start);
					}
				}
			}
		}
		
		for(list<Edge *>::iterator eitr = (largest_cycle->edges).begin(); eitr != (largest_cycle->edges).end(); eitr++){
			Edge *e = *eitr;
			for(vector<Edge*>::iterator eritr = e->refinement.begin(); eritr != e->refinement.end(); eritr++){
				Edge *er = *eritr;
				Vector *v1 = Vector::find_closest(&in_vertices,er->start);
				Vector *v2 = Vector::find_closest(&in_vertices,er->end);
				if(v1 != NULL){
					float distance = Vector::distance(v1,er->start);
					if(distance > max_distance)
						max_distance = distance;
					else if(distance < min_distance)
						min_distance = distance;									
					triangles->push_back(new Triangle(v1,er->start,er->end,atom->position,atom->radius,CONVEX));
					if(v2 != NULL && !(*v1 == *v2))
						triangles->push_back(new Triangle(v1,er->end,v2,atom->position,atom->radius,CONVEX));
				}
			}
		}
	}
	//*out << refinement.size() << " " << cycles.size() << " " << triangles->size() << endl;
}
	
void Atom::triangulate(){
	hash_map<int, float, hash<int>,eqint> angle;	
	int num_iterations = 0;
	bool done = (convexfaces.size() == 0);

	vector<Face*> in_convexfaces, out_convexfaces = convexfaces;
	bool refined_outermostfaces = false;
	
	int num_neighbors = toruses.size();
	Atom **atom_details = (Atom**) malloc((num_neighbors+1)*sizeof(Atom*));
	float *atom_radii = (float*) malloc((num_neighbors+1)*sizeof(float));
	float *atom_eps = (float*) malloc((num_neighbors+1)*sizeof(float));
	hash_map<int, Torus*, hash<int>,eqint> torii;
	
	int maxcindex = cindex;
	for(vector<Atom *>::iterator nitr = neighbors.begin(); nitr != neighbors.end(); nitr++){
		Atom *a = *nitr;
		if(a->cindex > maxcindex)
			maxcindex = a->cindex;
	}
	//*out << "triangulate atom " <<  " #neighbors " << toruses.size() << " " << neighbors.size() << " " << maxcindex << endl;
	while(!done){
		*out << "Atom " << index << " " << cindex << " #iterations " << num_iterations << endl;
	
		// increase the radii of neighbors, keep the atoms sorted by index
		// initialize
		if(num_iterations == 0){
			for(vector<Face *>::iterator titr = toruses.begin(); titr != toruses.end(); titr++){
				Torus *t = (Torus*) *titr;
				if(t->a1->cindex == cindex)
					torii[t->a2->cindex] = t;	
				else
					torii[t->a1->cindex] = t;
			}
		}

		int ind = 0;
		bool inserted_this = false;
		for(vector<Atom *>::iterator nitr = neighbors.begin(); nitr != neighbors.end(); nitr++){
			Atom *a = *nitr;
			Torus *t = torii[a->cindex];
			//*out << a->cindex << " " << t << endl;
			float rc, dc;
			if(cindex < a->cindex){
				rc = t->ccr2;
				dc = (t->axis)->dot(*(a->position) - *(t->ccc2));
				
				if(!inserted_this){
					atom_details[ind] = this;
					atom_radii[ind] = radius;
					ind++;
					inserted_this = true;
				}
			} else {
				rc = t->ccr1;
				dc = (t->axis)->dot(*(t->ccc1) - *(a->position));
			}
		
			atom_details[ind] = a;
			if(num_iterations == 0){
				//*out << "atan rc " << rc << " dc " << dc << endl;
				angle[a->cindex] = atan(rc/dc);
				if(angle[a->cindex] < 0)
					angle[a->cindex] = PI + angle[a->cindex];
			} 
			angle[a->cindex] += convex_face_inc/radius;
			done = done || (angle[a->cindex] >= PI );
			
			float d1 = t->distance - (radius + probe_radius)*cos(angle[a->cindex]);
			float d2 = (radius + probe_radius)*sin(angle[a->cindex]);
			atom_radii[ind] = sqrt(d1*d1 + d2*d2) - probe_radius;
			ind++;
		}
		if(!inserted_this){
			atom_details[ind] = this;
			atom_radii[ind] = radius;
			atom_eps[ind] = eps;
			ind++;
		}
		
		/*for(int i = 0 ; i <= num_neighbors; i++){
			*out << " " << atom_details[i]->cindex << " " << atom_radii[i] << endl;
		}*/

		if(!done){
			// ensure that the atoms are sorted by index
			Complex *c = new Complex(atom_details,atom_radii,atom_eps,num_neighbors+1, maxcindex+1);
			c->compute_convex_surface_of_atom(cindex);
			in_convexfaces = c->atom[cindex]->convexfaces;
			
			//*out << "new cfaces " << in_convexfaces.size() << " old " << convexfaces.size() << endl;
			// check if the atom is buried
			done = done || (in_convexfaces.size() == 0);
			
			if(in_convexfaces.size() > 0){
				// check if can link cycles on the basis of containment
				bool contained_in[in_convexfaces.size()][out_convexfaces.size()];
				int ocfindex = 0;
				for(vector<Face*>::iterator cfitr = out_convexfaces.begin(); cfitr != out_convexfaces.end(); cfitr++){
					Convex_Face *ocf = (Convex_Face*) *cfitr;
					int icfindex = 0;
					for(vector<Face*>::iterator icfitr = in_convexfaces.begin(); icfitr != in_convexfaces.end(); icfitr++){
						Convex_Face *icf = (Convex_Face*) *icfitr;
						contained_in[icfindex][ocfindex] = ocf->contains(icf);
						if(contained_in[icfindex][ocfindex]){
							icf->container = ocf;
							ocf->refinement.push_back(icf);
						}
						icfindex++;
					}
					ocfindex++;
				}
				if(out_convexfaces.size() == 1 && in_convexfaces.size() == 1){
					Convex_Face *ocf = (Convex_Face*) *(out_convexfaces.begin());
					Convex_Face *icf = (Convex_Face*) *(in_convexfaces.begin());
					if(!contained_in[0][0]){
						*out << "WARNING: atom-"<< index << " convexface refinement not contained" << endl;
						icf->container = ocf;
						ocf->refinement.push_back(icf);
					}
				}
				for(int i = 0; i < in_convexfaces.size(); i++){
					int outdegree=0;
					for(int j = 0; j < out_convexfaces.size(); j++)
						if(contained_in[i][j]) outdegree++;
					if(outdegree > 1)
						*out << "ERROR: inner convexface in the interior of more than one outer convexfaces" << endl;
				}
				
				out_convexfaces = in_convexfaces;
			}
			
			c->cleanup_except_atom(cindex);
			delete c;
		}
		num_iterations++;
	}

	delete atom_details;
	delete atom_radii;
	
	*out << "triangulating final faces " << out_convexfaces.size() << endl; out->flush();
	// triangulate final convexfaces
	if(convexfaces.size() > 0){
		int cf_index = 0;
		for(vector<Face*>::iterator cfitr = convexfaces.begin(); cfitr != convexfaces.end(); cfitr++){
			Convex_Face *cf = (Convex_Face*) *cfitr;
			if(!cf->is_buried){
				vector<Triangle*> cf_triangles;
				cf->triangulate(&cf_triangles);
				
				*out << " cface " << cf_index << " #triangles " << cf_triangles.size() << endl;
				//" max distance " << max_distance << " min distance " << min_distance << endl;
				for(vector<Triangle*>::iterator titr = cf_triangles.begin(); titr != cf_triangles.end(); titr++){
					triangles.push_back(*titr);
				}
			}
			cf_index++;
		}
	}
	*out << cindex << " " << is_buried << " " << triangles.size() << endl;
}

void Convex_Face::refine_edges(vector<Edge*> *edges){
	//for(vector<Face*>::iterator cfitr = convexfaces->begin(); cfitr != convexfaces->end(); cfitr++){
		Convex_Face *cf = this; //(Convex_Face*) *cfitr;
		for(vector<Cycle*>::iterator citr = cf->cycles.begin(); citr != cf->cycles.end(); citr++){
			Cycle *c = *citr;
			for(list<Edge *>::iterator eitr = (c->edges).begin(); eitr != (c->edges).end(); eitr++){
				Edge *e = *eitr;
				//*out << "radius " << e->radius << " " << e->type << endl;
				float angle;
				Vector start, end;
				if(e->type == CIRCLE_EDGE){
					Vector v = Vector(1,0,0);
					Vector axis_perpendicular = v - *e->axis * (v.dot(e->axis));
					if(axis_perpendicular == Vector(0,0,0)){
						v = Vector(0,1,0);
						axis_perpendicular = v - *e->axis * (v.dot(e->axis));
					}
				    axis_perpendicular.normalize();
					
					angle = 2*PI;
					start = Vector( *e->center + axis_perpendicular*(e->radius) );
					end = Vector( *e->center + axis_perpendicular*(e->radius) );
				} else {
					angle = Vector::angle( *(e->start) - *(e->center), *(e->end) - *(e->center), *e->axis);
					start = *e->start;
					end = *e->end;
				}
				
				int length_divisions = (int) (e->radius*angle/max_saddle_length);
				if(e->radius*angle/max_saddle_length - length_divisions < 0.5 && length_divisions > 0)
					length_divisions--;
				if(e->type == CIRCLE_EDGE && length_divisions < 2)
					length_divisions = 2;
						
				float theta = angle/(length_divisions + 1.0);
				Vector ur = Vector(start - *e->center);
				ur.normalize();
				Vector ut = ur.cross(e->axis);
				ut.normalize();

				//*out << "edge divisions " << length_divisions+1 << " " << angle << " " << e->radius << endl;
				Vector* vertices[length_divisions+2];
				for(int i = 0 ; i < length_divisions+2 ; i++){
					Vector vc = *e->center;
					float r = e->radius;
					float alpha = i*theta;
					Vector v = Vector( vc + ur * (r*cos(alpha)) + ut * (r*sin(alpha)) );
					vertices[i] = new Vector(v);
				}
		
				for(int i = 0 ; i < length_divisions+1 ; i++){
					Edge *edge = new Edge(vertices[i],vertices[i+1],e->center,e->axis,e->radius);
					edges->push_back(edge);
					e->refinement.push_back(edge);
				}
			}
		}
	//}
	//*out << "getedges #edges " << edges->size() << endl;
}
		
int Atom::get_max_cycles_on_convexface(vector<Face*> cfaces){
	int max = 0;
	vector<Face*>::iterator cfitr = cfaces.begin(), cfend = cfaces.end();
	while(cfitr != cfend){
		Convex_Face *cf = (Convex_Face*) *cfitr;

		if(cf->cycles.size() > max)
			max = cf->cycles.size();
		cfitr++;
	}
	return max;
}
	
void Complex::compute_convex_surface_of_atom(int atom_index){
	//*out << "compute convex suface of " << atom_index << endl; out->flush();

	// the atoms indices are not increasing
	Atom *sorted_atoms[num_atoms];
	int count = 0;
	for(int i = 0; i < num_atoms; i++){
		while(atom[count] == NULL) count++;
	   	sorted_atoms[i] = atom[count++];
	   	//*out << count << " " << 	sorted_atoms[i]->cindex << endl; out->flush();
	}

	// hash_map iterator does not preserve order, sort them again
	for(int i = 0 ; i < num_atoms; i++){
		for(int j = i+1 ; j < num_atoms; j++){
			if(sorted_atoms[i]->cindex > sorted_atoms[j]->cindex){
				Atom *a = sorted_atoms[i];
				sorted_atoms[i] = sorted_atoms[j];
				sorted_atoms[j] = a;
			}
		}
	}

	/*out << "sorted indices in m \n";
	for(int i = 0 ; i < num_atoms; i++)
		*out << " " << sorted_atoms[i]->cindex;
	*out << endl; out->flush();*/

	float distance[num_atoms][num_atoms];
	for(int i = 0 ; i < num_atoms; i++){
		Atom* a1 = sorted_atoms[i];
		for(int j = i+1; j < num_atoms; j++){
			Atom* a2 = sorted_atoms[j];
			distance[i][j] = distance[j][i] = Vector::distance((a1->position), (a2->position));

			a1->is_buried |= (distance[i][j] + a1->radius < a2->radius);
			a2->is_buried |= (distance[i][j] + a2->radius < a1->radius);
		}
	}

	vector<Torus*> vtoruses;
	for(int i = 0 ; i < num_atoms; i++){
		Atom* a1 = sorted_atoms[i];
		if(!a1->is_buried && a1->mass > 1){
			for(int j = i+1; j < num_atoms; j++){
				Atom* a2 = sorted_atoms[j];
				if(!a2->is_buried  && a2->mass > 1){
					// pairs of atoms forming tori
					if(distance[i][j] < a1->radius + a2->radius + 2*probe_radius){
						//*out << a1->cindex << " " << a2->cindex << " " << a1->radius << " " << a2->radius << " " << distance[i][j] << endl;	out->flush();
						// create torus we have a1->cindex < a2->cindex
						Torus* torus = new Torus(a1,a2,distance[i][j]);
						toruses[(a1->cindex)*MAX_ATOMS + (a2->cindex)] = torus;
						vtoruses.push_back(torus);
						
						// update neighbors
						a1->neighbors.push_back(a2);
						a2->neighbors.push_back(a1);
					}
				}
			}
		}
	}
	
	// sort the neighbors of atoms intersection assumes that the neighbors are sorted by index
	for(int i = 0 ; i < num_atoms; i++){
		Atom* a = sorted_atoms[i];
		sort(a->neighbors.begin(), a->neighbors.end(),less<Atom*>());
	}
	
	// for each torus compute the triangles it is involved in
	// interested in torii that involve target atom
	count = 0;
	//for(hash_map<const long,Torus*,hash<long>,eqlong>::iterator titr = toruses.begin(); titr != toruses.end(); titr++){
	for(vector<Torus*>::iterator titr = vtoruses.begin(); titr != vtoruses.end(); titr++){
		count++;
		Torus* t12 = (Torus*) (*titr);
		Atom* a1 = t12->a1;
		Atom* a2 = t12->a2;
		bool target_found = (a1->cindex == atom_index ) || (a2->cindex == atom_index);

		// find mutual neigbhors of the atoms involved
		vector<Atom *> all_common_neighbors = Atom::get_common_neighbors(a1,a2);
		vector<Atom *> common_neighbors = Atom::get_common_neighbors(a1,a2);
		
		// ensure that a1->cindex < a2->cindex < a3->cindex
		// uses the fact that the common neighbors are sorted by index
		vector<Atom *>::iterator itr = common_neighbors.begin();
		while( (itr != common_neighbors.end()) && ((*itr)->cindex < a2->cindex) ){
			itr++;
		}
		common_neighbors.erase(common_neighbors.begin(),itr);
		//*out << a1->cindex << " " << a2->cindex << " #cn " << common_neighbors.size() << endl;
		
		for(itr = common_neighbors.begin();itr != common_neighbors.end(); itr++){
			Atom* a3 = *itr;

			if(target_found || (a3->cindex == atom_index)){
				Torus* t13 = toruses[a1->cindex * MAX_ATOMS + a3->cindex];	
				Torus* t23 = toruses[a2->cindex * MAX_ATOMS + a3->cindex];	

				Vector* u12 = t12->axis;
				Vector* u13 = t13->axis;
		
				float triangle_angle = acos(u12->dot(u13));
				//*out << a1->cindex << " " << a2->cindex << " " << a3->cindex << " " << triangle_angle << "\t";
				if(triangle_angle != 0 && u12->dot(u13) != -1){ // check for colinearity
					Vector base_plane_normal = u12->cross(u13);
					base_plane_normal.normalize();
			
					Vector base_point = Vector( *(t12->center) + (base_plane_normal.cross(u12)) * ( (*(t13->center) - *(t12->center)).dot(u13)/sin(triangle_angle) ) );
					float probe_height_squared = (a1->radius + probe_radius)*(a1->radius + probe_radius) - (base_point - *(a1->position)).norm() * (base_point - *(a1->position)).norm();
					//*out << " probe height squared " << probe_height_squared << endl;
					
					if(probe_height_squared > 0){
						float probe_height = sqrt(probe_height_squared);
						t12->is_free = false;
						t13->is_free = false;
						t23->is_free = false;

						// check if the probe position causes collisions
						vector<Atom *> neighbors_a1a2a3 = a3->get_neighbors_in(all_common_neighbors);
						bool can_place_probe_above = true, can_place_probe_below = true;

						Vector* probe_position_above = new Vector(base_point + (base_plane_normal * probe_height));
						Vector* probe_position_below = new Vector(base_point - (base_plane_normal * probe_height));
						if(neighbors_a1a2a3.size() > 0){
							for(vector<Atom *>::iterator itr1 = neighbors_a1a2a3.begin(); itr1 != neighbors_a1a2a3.end(); itr1++){
								// check for collision with each atom
								Atom* a = *itr1;
								if(can_place_probe_above && Vector::distance(probe_position_above, a->position) < a->radius + probe_radius){
									can_place_probe_above = false;
									//*out << "lost above " << a3->cindex << " intersects " << a->cindex << " " << (Vector::distance(probe_position_above, a->position) - (a->radius + probe_radius)) << endl;
									delete probe_position_above;
								}
								if(can_place_probe_below && Vector::distance(probe_position_below, a->position) < a->radius + probe_radius){
									can_place_probe_below = false;
									//*out << "lost below " << a3->cindex << " intersects " << a->cindex << " " << (Vector::distance(probe_position_below, a->position) - (a->radius + probe_radius))<< endl;
									delete probe_position_below;
								}
							}
						}
						
						Vector *v1, *v2, *v3;
						if(can_place_probe_above){
							//*out << "obtained triangle above " << a1->cindex << " " << a2->cindex << " " << a3->cindex << endl;

							v1 = new Vector(*(a1->position) * ((probe_radius)/(a1->radius + probe_radius)) + *(probe_position_above) * ((a1->radius)/(a1->radius + probe_radius)));
							v2 = new Vector(*(a2->position) * ((probe_radius)/(a2->radius + probe_radius)) + *(probe_position_above) * ((a2->radius)/(a2->radius + probe_radius)));
							v3 = new Vector(*(a3->position) * ((probe_radius)/(a3->radius + probe_radius)) + *(probe_position_above) * ((a3->radius)/(a3->radius + probe_radius)));
							// direct the edges such that they are counter clockwise when seen from the probe position
							float orientation = base_plane_normal.dot( (*v2 - *v1).cross(*v3 - *v1) );

							Triangle_Face* triangle;
							if(orientation > 0){
								triangle = new Triangle_Face(a1,a2,a3,v1,v2,v3, probe_position_above);
								trianglefaces[(a1->cindex * MAX_ATOMS + a2->cindex)*MAX_ATOMS + a3->cindex] = triangle;
							} else {
								triangle = new Triangle_Face(a1,a3,a2,v1,v3,v2, probe_position_above);
								trianglefaces[(a1->cindex * MAX_ATOMS + a3->cindex)*MAX_ATOMS + a2->cindex] = triangle;
							}
							
							//*out << orientation << " "; triangle->print_details();
							t12->trianglefaces.push_back(triangle);
							t13->trianglefaces.push_back(triangle);
							t23->trianglefaces.push_back(triangle);
						}

						if(can_place_probe_below){
							//*out << "obtained triangle below " << a1->cindex << " " << a2->cindex << " " << a3->cindex << endl;
							v1 = new Vector(*(a1->position) * ((probe_radius)/(a1->radius + probe_radius)) + *(probe_position_below) * ((a1->radius)/(a1->radius + probe_radius)));
                            v2 = new Vector(*(a2->position) * ((probe_radius)/(a2->radius + probe_radius)) + *(probe_position_below) * ((a2->radius)/(a2->radius + probe_radius)));
	                        v3 = new Vector(*(a3->position) * ((probe_radius)/(a3->radius + probe_radius)) + *(probe_position_below) * ((a3->radius)/(a3->radius + probe_radius)));
        	                // direct the edges such that they are counter clockwise when seen from the probe position
                	        float orientation = - base_plane_normal.dot( (*v2 - *v1).cross(*v3 - *v1) );
							
							Triangle_Face* triangle;
							if(orientation > 0){
								triangle = new Triangle_Face(a1,a2,a3,v1,v2,v3, probe_position_above);
								trianglefaces[(a1->cindex * MAX_ATOMS + a2->cindex)*MAX_ATOMS + a3->cindex] = triangle;
							} else {
								triangle = new Triangle_Face(a1,a3,a2,v1,v3,v2, probe_position_above);
								trianglefaces[(a1->cindex * MAX_ATOMS + a3->cindex)*MAX_ATOMS + a2->cindex] = triangle;
							}
							
							//*out << orientation << " "; triangle->print_details();
							t12->trianglefaces.push_back(triangle);
							t13->trianglefaces.push_back(triangle);
							t23->trianglefaces.push_back(triangle);
						}
					} else {
						Vector u = (*t12->center - *a3->position) - *t12->axis * (*t12->center - *a3->position).dot(*t12->axis);
						u.normalize();
						if((*t12->center + u * t12->radius - *a3->position).norm() <= a3->radius + probe_radius
							&& (*t12->ccc1 + u * t12->ccr1 - *a3->position).norm() <= a3->radius + probe_radius
							&& (*t12->ccc2 + u * t12->ccr2 - *a3->position).norm() <= a3->radius + probe_radius){
							//*out << "Trous buried " << t12->a2->cindex << " " << t12->a2->cindex << endl;
							t12->is_free = false;
							t12->is_buried = true;
						}
	
						u = (*t23->center - *a1->position) - *t23->axis * (*t23->center - *a1->position).dot(*t23->axis);
						u.normalize();
						if((*t23->center + u * t23->radius - *a1->position).norm() <= a1->radius + probe_radius
							&& (*t23->ccc1 + u * t23->ccr1 - *a1->position).norm() <= a1->radius + probe_radius
							&& (*t23->ccc2 + u * t23->ccr2 - *a1->position).norm() <= a1->radius + probe_radius){
							//*out << "Trous buried " << t23->a1->cindex << " " << t23->a2->cindex << endl;
							t23->is_free = false;
							t23->is_buried = true;
						}
	
						u = (*t13->center - *a2->position) - *t13->axis * (*t13->center - *a2->position).dot(*t13->axis);
						u.normalize();
						if((*t13->center + u * t13->radius - *a2->position).norm() <= a2->radius + probe_radius
							&& (*t13->ccc1 + u * t13->ccr1 - *a2->position).norm() <= a2->radius + probe_radius
							&& (*t13->ccc2 + u * t13->ccr2 - *a2->position).norm() <= a2->radius + probe_radius){
							//*out << "Trous buried " << t13->a1->cindex << " " << t13->a2->cindex << endl;
							t13->is_free = false;
							t13->is_buried = true;
						}
					} // probe_height_squred
				} // triangle angle
			}
		}
	}
	*out << "#toruses examined " << count << " #toruses " << toruses.size() << " " << vtoruses.size() << endl; out->flush();

	//compute saddle faces
	//it is enough to look at the toruses of the given atom
	Atom *a = atom[atom_index];
	for(vector<Face*>::iterator titr = (a->toruses).begin(); titr != (a->toruses).end(); titr++){
		Torus* t = (Torus *) *titr;
		if(t->a1->cindex == atom_index || t->a2->cindex == atom_index){
			if(!t->is_buried)//t->trianglefaces.size() > 0)
				t->build_saddle_faces();
			//t->print_details();
			if(!t->is_free)
				t->is_buried = ((t->saddles).size() == 0);
		}
	}

	//compute convex faces of the target atom
	bool atom_buried = a->is_buried;
	if(!atom_buried){
		atom_buried = true;
	  	for(vector<Face*>::iterator titr = a->toruses.begin(); titr != a->toruses.end(); titr++){
    		Torus* t = (Torus *) *titr;
    		atom_buried &= (t->is_buried);
  		}
	}
   	if(!atom_buried)
		a->build_convex_faces();

	/**out << "convex faces on atom " << atom->cindex << endl;
	fitr = (atom->convexfaces).begin();
	while(fitr != (atom->convexfaces).end()){
		Convex_Face *cf = (Convex_Face *) *fitr;
		cf->print_details();
		fitr++;
	}*/
}

/*
 * For some weird reason VRML expects faces to be continuous, otherwise it jumps faces 
 * and creates weird triangles
 * 
 * the triangles on saddles and the triangles on convex faces are generated independently
 */
void print_triangles(vector<Triangle*> triangles, char* filename, short format){
	fstream fout(filename, fstream::out);
	
	hash_map<long,Vector*,hash<long>,eqlong> vertex;
	hash_map<const char*,long,hash<const char*>,eqstr> coord_vs_vertex;
	Graph* vertex_graph = new Graph();
	hash_map<long,GVertex*,hash<long>,eqlong> vertex_graph_node;
	Graph* triangle_graph = new Graph();
	hash_map<long,GVertex*,hash<long>,eqlong> triangle_graph_node;
	hash_map<long,hash_set<long,hash<long>,eqlong>,hash<long>,eqlong> vertex_vs_triangle;
	
	long num_vertices = 0;
	long num_faces = 0;
	
	switch(format){
		case VRML1:
			fout << "#VRML V1.0 ascii\nSeparator {\n";
		break;
		case VRML2:
			fout << "#VRML V2.0 utf8\n\n";
			//fout << "Shape {\n";
		break;
	}
	
	for(vector<Triangle *>::iterator itr = triangles.begin(); itr != triangles.end(); itr++){
		Triangle *t = *itr;
		long vi1,vi2,vi3;
		vi1 = vi2 = vi3 = -1;
		string vdetails;
		stringstream line1,line2,line3;
		Vector v = *(t->v1);
		line1 << "(" << v.x << "," << v.y << "," << v.z << ")";
		vdetails = line1.str();
		if(coord_vs_vertex.count(vdetails.c_str()) > 0){
			vi1 = coord_vs_vertex[vdetails.c_str()];
		} else {
			vi1 = num_vertices;
			vertex[num_vertices] = (t->v1);
			coord_vs_vertex[(new string(vdetails))->c_str()] = vi1;
			vertex_vs_triangle[vi1] = *(new hash_set<long,hash<long>,eqlong>);
			GVertex *gv = new GVertex(num_vertices);
			vertex_graph->insert_vertex(gv);
			vertex_graph_node[num_vertices] = gv;
			num_vertices++;
		}

		v = *(t->v2);
		line2 << "(" << v.x << "," << v.y << "," << v.z << ")";
		vdetails = line2.str();
		if(coord_vs_vertex.count(vdetails.c_str()) > 0){
			vi2 = coord_vs_vertex[vdetails.c_str()];
		} else {
			vi2 = num_vertices;
			vertex[num_vertices] = (t->v2);
			coord_vs_vertex[(new string(vdetails))->c_str()] = vi2;
			vertex_vs_triangle[vi2] = *(new hash_set<long,hash<long>,eqlong>);
			GVertex *gv = new GVertex(num_vertices);
			vertex_graph->insert_vertex(gv);
			vertex_graph_node[num_vertices] = gv;
			num_vertices++;
		}

		v = *(t->v3);
		line3 << "(" << v.x << "," << v.y << "," << v.z << ")";
		vdetails = line3.str();
		if(coord_vs_vertex.count(vdetails.c_str()) > 0){
			vi3 = coord_vs_vertex[vdetails.c_str()];
		} else {
			vi3 = num_vertices;
			vertex[num_vertices] = (t->v3);
			coord_vs_vertex[(new string(vdetails))->c_str()] = vi3;
			vertex_vs_triangle[vi3] = *(new hash_set<long,hash<long>,eqlong>);
			GVertex *gv = new GVertex(num_vertices);
			vertex_graph->insert_vertex(gv);
			vertex_graph_node[num_vertices] = gv;
			num_vertices++;
		}

		vertex_graph_node[vi1]->add_neighbor(vertex_graph_node[vi2]);
		vertex_graph_node[vi2]->add_neighbor(vertex_graph_node[vi1]);
		vertex_graph_node[vi2]->add_neighbor(vertex_graph_node[vi3]);
		vertex_graph_node[vi3]->add_neighbor(vertex_graph_node[vi2]);
		vertex_graph_node[vi3]->add_neighbor(vertex_graph_node[vi1]);
		vertex_graph_node[vi1]->add_neighbor(vertex_graph_node[vi3]);
		
		{
			GVertex *gv = new GVertex(num_faces);
			gv->data = t;
			triangle_graph->insert_vertex(gv);
			triangle_graph_node[num_faces] = gv;
		
			vertex_vs_triangle[vi1].insert(num_faces);
			vertex_vs_triangle[vi2].insert(num_faces);
			vertex_vs_triangle[vi3].insert(num_faces);
		
			num_faces++;
		}
		
		switch(format){
			case VRML1: case VRML2:
				fout << "Shape {\n";
				if(format == VRML2)
					fout << " appearance Appearance {\n\t material ";
				else
					fout << "\n\t";
				fout << "Material {\n\t\tdiffuseColor ";
				switch(t->orientation){
					case CONVEX:
						fout << "0.4 0.9 0.4";
					break;
					case PLANE:
						fout << "0.9 0.4 0.4";
					break;
					case CONCAVE:
						fout << "0.4 0.4 0.9";
					break;
				}
				/*{
					short color = (short) (10.0 * t->curvature);
					if(color > 4) color = 4;
					if(color < -4) color = -4;
					if(color > 0)
						fout << (0.9 - color*0.05) << " " << (0.4 + color*0.05) << " 0.4";
					else
						fout << (0.9 + color*0.05) << " 0.4 " << (0.4 - color*0.05);
				}*/
				//fout << "1 1 1";
				fout << "\n\t\ttransparency 1\n\t}\n";
				if(format == VRML2)
					fout << " }\n\tgeometry IndexedFaceSet {\n\t coord Coordinate ";
				else
					fout << "\tCoordinate3 ";
					
				fout << "{\n\t\tpoint [\n";
				v = *(t->v1);
				fout << "\t\t\t" << v.x << " " << v.y << " " << v.z << endl;
				v = *(t->v2);
				fout << "\t\t\t" << v.x << " " << v.y << " " << v.z << endl;
				v = *(t->v3);
				fout << "\t\t\t" << v.x << " " << v.y << " " << v.z << endl;
				fout <<"\t\t]\n\t}";
				fout << "\n\t";
				
				if(format == VRML1)
					fout << "IndexedFaceSet {\n\t\t";
				fout << "coordIndex [\n";
				fout << "\t\t\t 0,1,2,-1" << endl;
				fout << "\t\t]\n\t}\n}\n";
			break;
		}
	}

	// compute neighbors of each triangle
	for(hash_map<long,GVertex*,hash<long>,eqlong>::iterator itr = triangle_graph_node.begin(); itr != triangle_graph_node.end(); itr++){
		long tindex = itr->first;
		Triangle *t = ((Triangle*) ((GVertex*) itr->second)->data);
		long vi1,vi2,vi3;
		string vdetails;
		stringstream line1,line2,line3;
		
		Vector v = *(t->v1);
		line1 << "(" << v.x << "," << v.y << "," << v.z << ")";
		vdetails = line1.str();
		vi1 = coord_vs_vertex[vdetails.c_str()];

		v = *(t->v2);
		line2 << "(" << v.x << "," << v.y << "," << v.z << ")";
		vdetails = line2.str();
		vi2 = coord_vs_vertex[vdetails.c_str()];

		v = *(t->v3);
		line3 << "(" << v.x << "," << v.y << "," << v.z << ")";
		vdetails = line3.str();
		vi3 = coord_vs_vertex[vdetails.c_str()];
		
		//find common neighbors of vi1 and vi2
		hash_set<long,hash<long>,eqlong> ts1,ts2,ts3;
		ts1 = vertex_vs_triangle[vi1];
		ts2 = vertex_vs_triangle[vi2];
		ts3 = vertex_vs_triangle[vi3];
		for(hash_set<long,hash<long>,eqlong>::iterator titr = ts1.begin(); titr != ts1.end(); titr++){
			long ntindex = *titr;
			if(ntindex != tindex && (ts2.count(ntindex) >0 || ts3.count(ntindex) > 0))
				triangle_graph_node[tindex]->add_neighbor(triangle_graph_node[ntindex]);
		}
		for(hash_set<long,hash<long>,eqlong>::iterator titr = ts2.begin(); titr != ts2.end(); titr++){
			long ntindex = *titr;
			if(ntindex != tindex && ts3.count(ntindex) > 0)
				triangle_graph_node[tindex]->add_neighbor(triangle_graph_node[ntindex]);
		}
	}
	
	/*{
		if(format == VRML1 || format == VRML2){
			fout << "Shape {\n";
				if(format == VRML2)
					fout << "\tappearance Appearance {\n\t material ";
				fout << "\n\tMaterial {\n\t\tdiffuseColor 0.4 0.9 0.9";
				fout << "\n\t\ttransparency 0\n\t}\n";
				if(format == VRML2)
					fout << " }\n\tgeometry IndexedFaceSet {\n\t coord Coordinate ";
				else
					fout << "\tCoordinate3 ";	
				fout << "{\n\t\tpoint [\n";
		}
	
		list<long> triangle_trial;
		triangle_graph->compute_vertex_cycle(&triangle_trial);
		*out << "#triangles " << triangles.size() << " " << triangle_graph->num_vertices << " " << triangle_trial.size() << endl;
		
		for(int i = 0 ; i < num_vertices; i++){
			Vector v = *vertex[i];
			switch(format){
				case OBJ:
					fout << "v " << v.x << " " << v.y << " " << v.z << endl;
				break;
				case VRML1: case VRML2:
					fout << "\t\t\t" << v.x << " " << v.y << " " << v.z << "," << endl;
				break;
			}
		}
		
		if(format == VRML1 || format == VRML2){
			if(format == VRML1)
				fout << "IndexedFaceSet {\n\t\t";
			fout << "coordIndex [\n";
		}
			
		for(list<long>::iterator itr = triangle_trial.begin(); itr != triangle_trial.end(); itr++){
			Triangle *t = ((Triangle*) triangle_graph_node[*itr]->data);
			long vi1,vi2,vi3;
			string vdetails;
			stringstream line1,line2,line3;
			
			Vector v = *(t->v1);
			line1 << "(" << v.x << "," << v.y << "," << v.z << ")";
			vdetails = line1.str();
			vi1 = coord_vs_vertex[vdetails.c_str()];
	
			v = *(t->v2);
			line2 << "(" << v.x << "," << v.y << "," << v.z << ")";
			vdetails = line2.str();
			vi2 = coord_vs_vertex[vdetails.c_str()];
	
			v = *(t->v3);
			line3 << "(" << v.x << "," << v.y << "," << v.z << ")";
			vdetails = line3.str();
			vi3 = coord_vs_vertex[vdetails.c_str()];
			
			switch(format){
				case OBJ:
					fout << "f " << vi1 << " " << vi2 << " " << vi3 << endl;
				break;
				case VRML1: case VRML2:
					fout << "\t\t\t" << vi1 << "," << vi2 << "," << vi3 << ",-1" << endl;
				break;
			}
		}
		if(format == VRML1 || format == VRML2)
			fout << "\t\t]\n\t}\n" << endl;
	}*/
	
	if(format == VRML1)
		fout << "}\n";
	fout.close();
}

void print_spheres(float** sphere, unsigned int num_spheres, char *filename, short format){
	fstream fout(filename, fstream::out);
	
	switch(format){
		case VRML1:
			fout << "#VRML V1.0 ascii\nSeparator {\n";
		break;
	}
	
	for(unsigned int i = 0; i < num_spheres; i++){
		switch(format){
			case VRML1:
			fout << "Separator {\n";
			fout << "\tTransform { translation " << sphere[i][0] << " " << sphere[i][1] << " " << sphere[i][2] << " }\n";
			//fout << "\tTransform { scaleFactor 2 2 2 }\n";
			fout << "\tMaterial {\n\t\tdiffuseColor 0.9 0.6 0.4\n\t\ttransparency 0\n\t}\n";
			fout << "\tSphere {\n\t\tradius " << sphere[i][3] << "\n\t}\n";
			fout << "}\n";
			break;
		}
	}
	
	switch(format){
		case VRML1:
			fout << "}\n";
		break;
	}
}
		
void print_points(vector<Cluster*> *faces, char* filename, short format){
	fstream fout(filename, fstream::out);
	
	switch(format){
		case VRML1:
			fout << "#VRML V1.0 ascii\nSeparator {\n";
		break;
		case VRML2:
			fout << "#VRML V2.0 utf8\n\n";
			//fout << "Shape {\n";
		break;
	}
	
	// verifying the critical points and normals
	for(vector<Cluster*>::iterator fitr = faces->begin(); fitr != faces->end(); fitr++){
		Cluster *f = *fitr;
		if(f->point != NULL){
			fout << "Shape {\n";
			if(format == VRML2)
				fout << " appearance Appearance {\n\t material ";
			else
				fout << "\n\t";
			fout << "Material {\n\t\tdiffuseColor ";
			{
				short color = (short) (10.0 * f->curvature);
				if(color > 4) color = 4;
				if(color < -4) color = -4;
				if(color > 0)
					fout << (0.9 - color*0.05) << " " << (0.4 + color*0.05) << " 0.4";
				else
					fout << (0.9 + color*0.05) << " 0.4 " << (0.4 - color*0.05);
			}
			fout << "\n\t\ttransparency 0\n\t}\n";
			if(format == VRML2)
				fout << " }\n\tgeometry IndexedFaceSet {\n\t coord Coordinate ";
			else
				fout << "\tCoordinate3 ";
					
			fout << "{\n\t\tpoint [\n";
			Vector v = *(f->point) + *(f->normal);
			fout << "\t\t\t" << v.x << " " << v.y << " " << v.z << endl;
			
			Vector w = Vector(1,0,0);
			Vector ex = w - *(f->normal) * (w.dot(f->normal));
			if(ex == Vector(0,0,0)){
				w = Vector(0,1,0);
				ex = w - *(f->normal) * (w.dot(f->normal));
			}
			ex.normalize();
			Vector ey = Vector(f->normal);
			Vector ez = Vector(ex.cross(ey));
			
			v = *(f->point) + ex*(0.3);
			fout << "\t\t\t" << v.x << " " << v.y << " " << v.z << endl;
			
			v = *(f->point) - ex*(0.15) + ez*(sqrt(3)*0.15);
			fout << "\t\t\t" << v.x << " " << v.y << " " << v.z << endl;
			
			v = *(f->point) - ex*(0.15) - ez*(sqrt(3)*0.15);
			fout << "\t\t\t" << v.x << " " << v.y << " " << v.z << endl;
			
			fout <<"\t\t]\n\t}";
			fout << "\n\t";
				
			if(format == VRML1)
				fout << "IndexedFaceSet {\n\t\t";
			
			fout << "coordIndex [\n";
			fout << "\t\t\t 0,1,2,-1" << endl;
			fout << "\t\t\t 0,2,3,-1" << endl;
			fout << "\t\t\t 0,3,1,-1" << endl;
			fout << "\t\t\t 1,2,3,-1" << endl;
			fout << "\t\t]\n\t}\n}\n";
				
				/*case VRML2:
					fout << "Translation {\n\t\ttranslation " << v.x << " " << v.y << " " << v.z << "\n\t\t}\n";
					
					//fout << "\t\t\t" << v.x << " " << v.y << " " << v.z << endl;
					fout << "\tMaterial {\n\t\tdiffuseColor 0.9 0.6 0.4\n\t\ttransparency 0\n\t}\n";
					fout << "\tSphere {\n\t\tradius 1\n\t}\n";
					fout << "}\n";
				break;*/
		}
	}
	fout << "}\n";
	fout.close();
}

void Complex::triangulate_surface(){
	triangles.clear();

	int convex_face_triangles = 0;
	for(int i = 0 ; i < num_atoms; i++){
   		Atom* a = atom[i];
   	    //*out << " triangulating atom " << a->cindex << " buried " << a->is_buried << " ";
   	    if(!(a->is_buried)  && a->mass > 1){
   	    	a->triangulate();
		 	//*out << a->triangles.size() << endl;
			for(vector<Triangle *>::iterator itr = (a->triangles).begin(); itr != (a->triangles).end(); itr++){
				triangles.push_back(*itr);
			}
			convex_face_triangles += (a->triangles).size();
   	    }
	}
	
	for(hash_map<const long, Torus *, hash<long>,eqlong>::iterator titr = toruses.begin(); titr != toruses.end(); titr++){
		Torus* t = (Torus *) titr->second;
        if(!t->is_buried){
			t->triangulate();
			for(vector<Triangle *>::iterator itr = (t->triangles).begin(); itr != (t->triangles).end(); itr++){
				//triangles.push_back(*itr);
			}
        }
	}
	
	for(hash_map<const long, Triangle_Face *, hash<long>,eqlong>::iterator tritr = trianglefaces.begin(); tritr != trianglefaces.end(); tritr++){
		Triangle* t = (Triangle *) tritr->second;
		if(!((Triangle_Face*) t)->is_buried){
			//triangles.push_back(t);
		}
	}
	
	//eliminate_triangles_in_holes();
	
	// compute the critical point for each triangle
	for(vector<Triangle *>::iterator titr = triangles.begin(); titr != triangles.end(); titr++){
		Triangle* t = *titr;
		t->is_buried = false;
		{
			Vector centroid = (*(t->v1) + *(t->v2) + *(t->v3))*(1.0/3.0);
			triangle_centroids.push_back(new Vector(centroid));
		}
	}
	
	*out << "#triangles " << triangles.size() << endl;
	*out << "#convexface triangles " << convex_face_triangles << endl;
}

/*
 * Should triangulate saddle first, since convex face triangles stitch into saddle (especially in the case when the torus is free)
 */
void Complex::triangulate_surface(int atom_index){
	//*out << "c->triangulate surface " << atom_index << endl;
	int convex_face_triangles = 0;
	for(int i = 0 ; i < num_atoms ; i++){
   		Atom* a = atom[i];
   	    
		//*out << a->cindex << " " << (a->toruses).size() << endl;
       	for(vector<Face*>::iterator titr = a->toruses.begin(); titr != a->toruses.end(); titr++){
        	Torus* t = (Torus *) *titr;
        	bool flag = !t->is_buried && (((t->a1->cindex == atom_index) || (t->a2->cindex == atom_index)));
        	if(flag){
				t->triangulate();
				//*out << atom_index << " " << t->a1->cindex << " " << t->a2->cindex << " " << (t->triangles).size() << endl;
				for(vector<Triangle *>::iterator itr = (t->triangles).begin(); itr != (t->triangles).end(); itr++){
					triangles.push_back(*itr);
				}
        	}
		}
	}
}

void Complex::write_as_pdb(string filename){
	fstream pdbout;
	pdbout.open(filename.c_str(), fstream::out);
	for( int i = 0; i < num_atoms; i++){
		Atom *a = atom[i];
		Monomer *mr = monomer[a->monocindex];
		if((a->name).size() <= 3){
			//FORMAT('ATOM',I7,2X,A3,1X,A3,I6,4X,3F8.3,2(1X,F5.2))
			float coord[3];
			Vector v =*(a->position);
			coord[0] = v.x;
			coord[1] = v.y;
			coord[2] = v.z;
			int ci[3];
			char cf[3][7];
			for(int i=0;i<3;i++){
				ci[i] = (int) coord[i];
				float d = coord[i] - ci[i];
				sprintf(cf[i],"%6.3f",d);
			}
			
			char buf[80];
			sprintf(buf,"ATOM  %5d  %-3s %-4s%c%4s    %4d.%s%4d.%s%4d.%s  1.00  0.00",
			 (a->cindex)+1,(a->name).c_str(),(mr->name).c_str(),mr->chain,(mr->index).c_str(),
			 ci[0],&cf[0][3],ci[1],&cf[1][3],ci[2],&cf[2][3]); 
			pdbout << buf << endl;
		}
	}
	pdbout.close();
}

void Complex::print_details(ostream *out){
	for(hash_map<char, Molecule*, hash<char>, eqint>::iterator mitr = molecules.begin(); mitr != molecules.end(); mitr++){
		Molecule *m = mitr->second;
		m->print_details(out);
	}
	
	int num_clusters = 0;
	vector<Cluster*>::iterator cfitr;
	for(cfitr = clusters.begin(); cfitr != clusters.end(); cfitr++){
       	Cluster* c = (Cluster *) *cfitr;
		if(c->point != NULL)
			num_clusters++;
	}
	*out << num_clusters << endl;
	*out << "~ point normal coarser_normal normalstability curvature radius area #atoms atom_indices #residues residue_indices" << endl;	
	for(cfitr = clusters.begin(); cfitr != clusters.end(); cfitr++){
       	Cluster* c = (Cluster *) *cfitr;
       	if(c->point != NULL)
			c->print_details(out);
	}
	
	*out << triangle_centroids.size() << endl;
	*out << "~ centroid" << endl;
	for(vector<Vector*>::iterator tcitr = triangle_centroids.begin(); tcitr != triangle_centroids.end(); tcitr++){
		Vector v = *((Vector*) *tcitr);
		*out << v.x << " " << v.y << " " << v.z << endl;
	}
	
	/*hash_map<const long, Torus*, hash<long>, eqlong>::iterator titr = toruses.begin(), tend = toruses.end();
	while(titr != tend){
        Torus* t = (Torus *) titr->second;
		t->print_details();
		titr++;
	}

	hash_map<const long, Triangle_Face*, hash<long>, eqlong>::iterator tritr = trianglefaces.begin(), trend = trianglefaces.end();
    while(tritr != trend){
   	   	Triangle_Face* tr = (Triangle_Face *) tritr->second;
		tr->print_details();	
		tritr++;
	}
	
	vector<Triangle *>::iterator ritr = triangles.begin(), rend = triangles.end();
	while(ritr != rend){
		Triangle *t = *ritr;
		t->print_details();
		ritr++;
	}*/
}

/*
 * Output format
 * ~ at the beginning of a line indicates a comment
 * #atoms
 * atom_details with residue information
 * #residues
 * aminoacid details
 * #clustered points
 * clustered face details
 * #triangles
 * critical points of triangles
 */
void Protein::print_details(ostream *out){
	*out << aminoacid.size() << endl;
	*out << "~ index name sstructure pInterface" << endl;
	//*out << "~ index name type centroid amide_nitrogen.index carbonyl_oxygen.index" << endl;
	for(hash_map<const char*, Aminoacid*, hash<const char*>,eqstr>::iterator aaitr = aminoacid.begin(); aaitr != aminoacid.end(); aaitr++){
       	Aminoacid* aa = (Aminoacid *) aaitr->second;
		aa->print_details(out);
	}
		
	*out << atom.size() << endl;
	*out << "~ chain index name position aaindex is_buried" << endl;
	for(hash_map<unsigned short, Atom*, hash<unsigned short>, eqint>::iterator aitr = atom.begin(); aitr != atom.end(); aitr++){
       	Atom* a = (Atom *) aitr->second;
       	*out << chain << " ";
		a->print_details(out);
	}
}

void Triangle_Face::compute_point(){
	point = NULL;
	area = 0;
	Vector g = (*v1 + *v2 + *v3 ) * (1.0/3.0);

	//normal points out of the surface
	normal = new Vector(*center - g);
	if(normal->norm() > 0){
		normal->normalize();
		point = new Vector( *center - *normal * probe_radius);
		area = compute_area();
	}
}

float Triangle::compute_area(){
	/*float angle_sum = 0;
	if(orientation != PLANE){
		Vector t12 = (*v2 - *v1) - (*v1 - *center) * ( (*v2 - *v1).dot(*v1 - *center) /  (radius * radius) );
		Vector t13 = (*v3 - *v1) - (*v1 - *center) * ( (*v3 - *v1).dot(*v1 - *center) /  (radius * radius) ) ;
		float d = t12.dot(t13) / ( t12.norm() * t13.norm() ) ;
		angle_sum += acos(d);
		
		Vector t23 = (*v3 - *v2) - (*v2 - *center) * ( (*v3 - *v2).dot(*v2 - *center) / (radius * radius) );
		Vector t21 = (*v1 - *v2) - (*v2 - *center) * ( (*v1 - *v2).dot(*v2 - *center) / (radius * radius) );
		d = t23.dot(t21) / ( t23.norm() * t21.norm() ) ;
		angle_sum += acos(d);
	
		Vector t31 = (*v1 - *v3) - (*v3 - *center) * ( (*v1 - *v3).dot(*v3 - *center) / (radius * radius) );
		Vector t32 = (*v2 - *v3) - (*v3 - *center) * ( (*v2 - *v3).dot(*v3 - *center) / (radius * radius) );
		d = t31.dot(t32) / ( t31.norm() * t32.norm() ) ;
		angle_sum += acos(d);

		a = (angle_sum - PI)*radius*radius;

		*out << " angle_sum " << angle_sum << " radius " << radius << " " << (*v1 - *center).norm() << " " << (*v2 - *center).norm() << " " << (*v3 - *center).norm() << " " <<  " area " << a << " plane area " << 0.5 * ((*v2 - *v1).cross(*v3 - *v1)).norm() << endl;
	} else*/
		area = 0.5 * ((*v2 - *v1).cross(*v3 - *v1)).norm();
			
	return area;
}

/*
 * For plane triangles is the orientation correct?
 */
Vector* compute_gravitational_center(vector<Triangle*> triangles, float *area){
	Vector g = Vector(0,0,0);
	*area = 0;
	for(vector<Triangle *>::iterator titr = triangles.begin(); titr != triangles.end(); titr++){
		Triangle *tr = *titr;
		if(!tr->is_buried){
			float a = tr->compute_area();
			Vector n;
	
			if(tr->orientation == PLANE)
				n = (*(tr->v3) - *(tr->v2)).cross((*(tr->v2) - *(tr->v1)));
			else
				n = (*(tr->v1) + *(tr->v2) + *(tr->v3)) * (1.0/3.0) - *(tr->center);
			
			g = g + n * a;
			*area = *area + a;
		}
	}

	if(*area != 0)
		return new Vector(g * (1.0/(*area)));
	else
		return NULL;
}


void Saddle::compute_point(Vector *axis, Vector *center, float radius){
	point = NULL;
	Vector *g = compute_gravitational_center(triangles,&area);
	if(g != NULL){
		Vector n = (*center - *g);
		if(n.norm() > 0){
			n.normalize();
			normal = new Vector(n);
			point = new Vector(*center - *normal * radius);
		}
	}
}

void Atom::compute_point(){
	point = NULL;
	Vector *g = compute_gravitational_center(triangles,&area);
	if(g != NULL){
		Vector n = *g - *position;
		if(n.norm() > 0){
			n.normalize();
			normal = new Vector(n);
			point = new Vector( *(position) + *normal * radius);
			//*out << "#check " << radius << " " << Vector::distance(*point,*position) << " " << normal->norm() << endl;
		}
	}
}	

void Complex::compute_curvature(unsigned int n, Face **face){	
	// use a grid to compute neighbors
	{
		float minx, maxx, miny, maxy, minz, maxz;
		for(int i = 0 ; i < num_atoms ; i++){
			Atom *a = atom[i];
			Vector v = *(a->position);
			if(i == 0){
				minx = maxx = v.x;
				miny = maxy = v.y;
				minz = maxz = v.z;
			}
			
			maxx = (v.x > maxx)? v.x : maxx;
			minx = (v.x < minx)? v.x : minx;
			maxy = (v.y > maxy)? v.y : maxy;
			miny = (v.y < miny)? v.y : miny;
			maxz = (v.z > maxz)? v.z : maxz;
			minz = (v.z < minz)? v.z : minz;
		}
	
		//*out << "extent " << maxx - minx << " " << maxy - miny << " " << maxz - minz << endl; out->flush();
		Vector *grid_origin = new Vector(minx-NEIGHBOR_CUTOFF-1,miny-NEIGHBOR_CUTOFF-1,minz-NEIGHBOR_CUTOFF-1);
		hash_map<long,hash_set<int,hash<int>,eqint>,hash<long>,eqlong> grid;
		int NUM_DIVISIONS = 512;
		float INCREMENT = NEIGHBOR_CUTOFF;
		
		for(int i = 0 ; i < n; i++){
			Vector v = *(face[i]->point) - *grid_origin;
			int vx = (int) (v.x/INCREMENT);
			int vy = (int) (v.y/INCREMENT);
			int vz = (int) (v.z/INCREMENT);

			long index = (vx*NUM_DIVISIONS + vy)*NUM_DIVISIONS + vz;
			if(grid.count(index) == 0){
				grid[index] = *(new hash_set<int,hash<int>,eqint>);
			}
			grid[index].insert(i);
		}
		//*out << "grid size " << grid.size() << endl; out->flush();
		
		for(int pi = 0 ; pi < n; pi++){
			Vector n = *(face[pi]->normal);
			for(int k = 0; k < 1; k++){
				face[pi]->num_neighbors = 0;
				Vector v = *(face[pi]->point) + n*NEIGHBOR_CUTOFF - *grid_origin;
				int vx = (int) (v.x/INCREMENT);
				int vy = (int) (v.y/INCREMENT);
				int vz = (int) (v.z/INCREMENT);
				
				for(int i = vx - (int) (NEIGHBOR_CUTOFF/INCREMENT); i <= vx + (int) (NEIGHBOR_CUTOFF/INCREMENT); i++){
					for(int j = vy - (int) (NEIGHBOR_CUTOFF/INCREMENT); j <= vy + (int) (NEIGHBOR_CUTOFF/INCREMENT); j++){
						for(int k = vz - (int) (NEIGHBOR_CUTOFF/INCREMENT); k <= vz + (int) (NEIGHBOR_CUTOFF/INCREMENT); k++){
							long index = (i*NUM_DIVISIONS + j)*NUM_DIVISIONS + k;
							//*out << pi << " " << index << endl; out->flush();
							if(grid.count(index) > 0){
								hash_set<int,hash<int>,eqint> candidates = grid[index];
								for(hash_set<int,hash<int>,eqint>::iterator citr = candidates.begin(); citr != candidates.end(); citr++){
									int pj = *citr;
									if(Vector::distance(face[pi]->point,face[pj]->point) < NEIGHBOR_CUTOFF){
										face[pi]->num_neighbors++;
									}
								}
							}
						}
					}
				}
			}
		}
		
		// find the closest point on an atom
		for(int ai = 0 ; ai < num_atoms ; ai++){
			Atom *a = atom[ai];
			Vector v = *(a->position) - *grid_origin;
			int vx = (int) (v.x/INCREMENT);
			int vy = (int) (v.y/INCREMENT);
			int vz = (int) (v.z/INCREMENT);
			
			int num_points_close_to_atom=0;
			float closest_point_distance;
			for(int i = vx - (int) (NEIGHBOR_CUTOFF/INCREMENT); i <= vx + (int) (NEIGHBOR_CUTOFF/INCREMENT); i++){
				for(int j = vy - (int) (NEIGHBOR_CUTOFF/INCREMENT); j <= vy + (int) (NEIGHBOR_CUTOFF/INCREMENT); j++){
					for(int k = vz - (int) (NEIGHBOR_CUTOFF/INCREMENT); k <= vz + (int) (NEIGHBOR_CUTOFF/INCREMENT); k++){
						long index = (i*NUM_DIVISIONS + j)*NUM_DIVISIONS + k;
						if(grid.count(index) > 0){
							hash_set<int,hash<int>,eqint> candidates = grid[index];
							for(hash_set<int,hash<int>,eqint>::iterator citr = candidates.begin(); citr != candidates.end(); citr++){
								int pj = *citr;
								float d = Vector::distance(a->position,face[pj]->point);
								if(num_points_close_to_atom == 0 || d < closest_point_distance)
									closest_point_distance = d;
								num_points_close_to_atom++;
							}
						}
					}
				}
			}
			if(a->radius*MSMS_EXPANSION_FACTOR + 0.001 < closest_point_distance){
				a->is_buried = true;
				//*out << ai << " " << a->radius << " " << closest_point_distance << endl;
			}
		}
	}	

	float sum_degree = 0, sum_degree_squared = 0;
	for(int i = 0 ; i < n; i++){
		sum_degree += face[i]->area * face[i]->num_neighbors;
		sum_degree_squared += face[i]->area * face[i]->num_neighbors * face[i]->num_neighbors;
	}	
	float avg_degree = sum_degree/sa;
	float sigma_degree = sqrt((sum_degree_squared - avg_degree*avg_degree)/sa);
	*out << "avg neighbors " << avg_degree << " sigma " << sigma_degree << endl;
	
	for(int i = 0 ; i < n; i++)
		face[i]->curvature = - ((face[i]->num_neighbors - avg_degree)/sigma_degree);
}

void Complex::compute_sphere_based_curvature(Vector *grid_origin, int grid_num_xdivisions, int grid_num_ydivisions, int grid_num_zdivisions, GridCell **grid, 
 unsigned int n, Face **face){
 	short extent=ceil(CURVATURE_SPHERE_RADIUS/GRID_SPACING);
	for(int pi = 0 ; pi < n; pi++){
		face[pi]->num_neighbors = 0;
		Vector v = *(face[pi]->point) - *grid_origin;
		int vx = (int) (v.x/GRID_SPACING);
		int vy = (int) (v.y/GRID_SPACING);
		int vz = (int) (v.z/GRID_SPACING);
		
		for(int x = max(0,vx-extent); x <= min(grid_num_xdivisions-1,vx+extent); x++)
			for(int y = max(0,vy-extent); y <= min(grid_num_ydivisions-1,vy+extent); y++)
				for(int z = max(0,vz-extent); z <= min(grid_num_zdivisions-1,vz+extent); z++){
					Vector grid_center = *grid_origin + Vector( (((float) x) + 0.5)*GRID_SPACING, (((float) y) + 0.5)*GRID_SPACING, (((float) z) + 0.5)*GRID_SPACING);
					bool intersects_cell = (Vector::distance(*(face[pi]->point),grid_center) <= CURVATURE_SPHERE_RADIUS);
					if(intersects_cell){
						Vector vg = grid_center - *grid_origin;
						int vgx = (int) (vg.x/GRID_SPACING);
						int vgy = (int) (vg.y/GRID_SPACING);
						int vgz = (int) (vg.z/GRID_SPACING);
	
						if(vgx >= 0 && vgx < grid_num_xdivisions && vgy >= 0 && vgy < grid_num_ydivisions && vgz >= 0 && vgz < grid_num_zdivisions){
							unsigned int index = (vgx*grid_num_ydivisions + vgy)*grid_num_zdivisions + vgz;
		
							if(grid[index] != NULL){
								//*out << pi << "\t" << index << endl;
								GridCell *cell = grid[index];
								if(cell->type != EXTERIOR_ATOM_NEIGHBOR)
									face[pi]->num_neighbors++;
							}
						}
					}
				}
	}
	
	unsigned int max_neighbors=0;
	for(int x = -extent; x <= extent; x++)
		for(int y = -extent; y <= extent; y++)
			for(int z = -extent; z <= extent; z++){
				Vector grid_center = Vector( (((float) x) + 0.5)*GRID_SPACING, (((float) y) + 0.5)*GRID_SPACING, (((float) z) + 0.5)*GRID_SPACING);
				if(Vector::distance(grid_center,Vector(0,0,0)) <= CURVATURE_SPHERE_RADIUS)
					max_neighbors++;
			}
	
	for(int i = 0 ; i < n; i++){
		*out << "c\t" << face[i]->curvature << "\t" << face[i]->num_neighbors << "\t";
		face[i]->curvature = 2*(((float)face[i]->num_neighbors)/max_neighbors - 0.5);
		*out << face[i]->curvature << endl;
	}	
}

/* 
assume that the surface is triangulated, for each face, compute the critical point 
currently 70% of surface is in saddle faces
let clustering take care of optimization, triangles approx of equal area? 
kmeans works best when the points are uniformly distributed -- ENSURE trianglefaces have the same area as other triangles
Use the distance on the surface not the l2 norm
*/
void Complex::compute_points(int algorithm){
	vector<Face *> coarse_faces, coarse_convex_faces, coarse_concave_faces, faces, convex_faces, concave_faces;
	
	int num_tf_triangles = 0;
	for(hash_map<const long, Triangle_Face*, hash<long>, eqlong>::iterator tritr = trianglefaces.begin(); tritr != trianglefaces.end(); tritr++){
       	Triangle_Face* trf = (Triangle_Face *) tritr->second;
       	if(!trf->is_buried && !trf->a1->is_buried && !trf->a2->is_buried && !trf->a3->is_buried){
			trf->compute_point();
			if(trf->area > 0){
				trf->radius = - probe_radius;
				//*out << "trf " << trf->area << " " << trf->radius << endl; out->flush();
				coarse_concave_faces.push_back(trf);
				coarse_faces.push_back(trf);
				vector<Triangle*> triangles;
				trf->triangulate(&triangles);
				num_tf_triangles += triangles.size();
				if(triangles.size() > 1){
					for(vector<Triangle *>::iterator titr = triangles.begin(); titr != triangles.end(); titr++){
						Triangle *tr = *titr;
						tr->compute_area();
						if(tr->area > 0){
							Vector g = (*(tr->v1) + *(tr->v2) + *(tr->v3)) * (1.0/3.0);
							Vector n = *(tr->center) - g;
							if(n.norm() > 0){
								n.normalize();
								tr->normal = new Vector(n);
								tr->point = new Vector( *(tr->center) - *(tr->normal) * probe_radius);
								//tr->type = CONCAVE_FACE;
								tr->radius = - probe_radius;
								tr->face_atom_indices.push_back(trf->a1->cindex);
								tr->face_atom_indices.push_back(trf->a2->cindex);
								tr->face_atom_indices.push_back(trf->a3->cindex);
								concave_faces.push_back(tr);
								faces.push_back(tr);
							}
						}
					}
				} else {
					concave_faces.push_back(trf);
					faces.push_back(trf);
				}
			}
       	}
	}
	*out << "#trfaces " << trianglefaces.size() << " #not buried " << coarse_concave_faces.size() << " #triangles " << num_tf_triangles;
	*out << endl; out->flush();

	int num_saddle_triangles = 0, num_free_torus_triangles = 0, num_saddle_faces = 0, num_buried_toruses = 0, num_free_toruses = 0;
	for(hash_map<const long, Torus*, hash<long>, eqlong>::iterator titr = toruses.begin(); titr != toruses.end(); titr++){
        Torus* t = (Torus *) titr->second;
        //*out << t->is_free << " " << t->is_buried << " " << (t->saddles).size() << endl;
        if(!t->is_buried && !t->a1->is_buried && !t->a2->is_buried){
	        num_saddle_faces += (t->saddles).size();
			for(vector<Saddle *>::iterator sitr = (t->saddles).begin(); sitr != (t->saddles).end(); sitr++){
				Saddle *s = *sitr;
				if(!s->is_buried){
					s->compute_point(t->axis, t->center, t->radius);
					if(s->point != NULL){
						s->radius = 0;
						coarse_faces.push_back(s);	
						for(vector<Triangle *>::iterator titr = (s->triangles).begin(); titr != (s->triangles).end(); titr++){
							Triangle *tr = *titr;
							if(tr->area > 0){
								Vector g = (*(tr->v1) + *(tr->v2) + *(tr->v3)) * (1.0/3.0);
								Vector n = (g - *t->center);
								n.normalize();
								tr->normal = new Vector(n);
								tr->point = new Vector(g);
								tr->radius = 0;
								tr->face_atom_indices.push_back(t->a1->cindex);
								tr->face_atom_indices.push_back(t->a2->cindex);
								faces.push_back(tr);
							}
						}
						num_saddle_triangles += (s->triangles).size();
					}
				}
			}
        } else
        	num_buried_toruses++;
        if(t->is_free && !t->is_buried){
        	num_free_toruses++;
        	for(vector<Triangle *>::iterator titr = (t->triangles).begin(); titr != (t->triangles).end(); titr++){
				Triangle *tr = *titr;
				tr->compute_area();
				if(tr->area > 0){
					Vector g = (*(tr->v1) + *(tr->v2) + *(tr->v3)) * (1.0/3.0);
					Vector n = (g - *t->center);
					n.normalize();
					tr->normal = new Vector(n);
					tr->point = new Vector(g);
					tr->radius = 0;
					tr->face_atom_indices.push_back(t->a1->cindex);
					tr->face_atom_indices.push_back(t->a2->cindex);
					faces.push_back(tr);
				}
			}
			num_free_torus_triangles += (t->triangles).size();
        }
	}
	*out << "#toruses " << toruses.size() << " #buried " << num_buried_toruses << " #free " << num_free_toruses << endl;
	*out << "#saddle_faces " << num_saddle_faces << " #saddle triangles " << num_saddle_triangles << " #free_torus_triangles " << num_free_torus_triangles << endl; 
	
	for(int i = 0 ; i < num_atoms; i++){
		Atom* a = atom[i];
		if(!a->is_buried){
			a->compute_point();
			Vector *v = a->point;
			if(v != NULL){
				coarse_convex_faces.push_back(a);
				coarse_faces.push_back(a);
			}
				
			for(vector<Triangle *>::iterator titr = (a->triangles).begin(); titr != (a->triangles).end(); titr++){
				Triangle *tr = *titr;
				if(tr->area > 0 && !tr->is_buried){
					Vector g = (*(tr->v1) + *(tr->v2) + *(tr->v3)) * (1.0/3.0);
					Vector n = g - *(a->position);
					if(n.norm() > 0){
						n.normalize();
						tr->normal = new Vector(n);
						tr->point = new Vector( *(a->position) + *(tr->normal) * (a->radius));
						tr->orientation = CONVEX;
						tr->face_atom_indices.push_back(a->cindex);
						tr->radius = a->radius;
						convex_faces.push_back(tr);
						faces.push_back(tr);
					}
				}
			}
		}
	}
	
	*out << "#Cconvexfaces " << coarse_convex_faces.size() << " " << " #Cconcavefaces " << coarse_concave_faces.size() << " #Cfaces " << coarse_faces.size() << endl;
	*out << "#convexfaces " << convex_faces.size() << " " << " #concavefaces " << concave_faces.size() << " #faces " << faces.size() << endl;
	
	sa = 0;
	for(vector<Face*>::iterator itr = concave_faces.begin(); itr != concave_faces.end(); itr++){
		sa += ((Face*) *itr)->area;
	}
	*out << "concavefaces area " << sa << " <area> " << sa/concave_faces.size() << endl;
	
	sa = 0;
	for(vector<Face*>::iterator itr = convex_faces.begin(); itr != convex_faces.end(); itr++){
		sa += ((Face*) *itr)->area;
	}
	*out << "convexfaces area " << sa << " <area> " << sa/faces.size() << endl;
	
	sa = 0;
	for(vector<Face*>::iterator itr = faces.begin(); itr != faces.end(); itr++){
		sa += ((Face*) *itr)->area;
	}
	*out << "exposed area " << sa << " allfaces <area> " << sa/faces.size() << endl;
	
	int n = faces.size();
	Face* face[n + 1];
	int count = 0;
	for(vector<Face*>::iterator itr = faces.begin(); itr != faces.end(); itr++){
		face[count++] = (Face*) *itr;
		//*out << count << " " << ((Face*) *itr)->radius << endl;
	}
	
	float min_a, max_a;
	int a_divisions;
	float box_size = 0.1;
	{
		for(int i = 0 ; i < n; i++){
			float area = face[i]->area;
			if(i == 0)
				min_a= max_a= area;
			else{
				if(area > max_a)
					max_a = area;
				else if(area < min_a)
					min_a = area;
			}
		}
		*out << "a min " << min_a << " max " << max_a << endl;
		/*for(int i = 0 ; i < n; i++){
			if(face[i]->area == max_a)
				face[i]->print_details();
		}*/
		a_divisions = (int) (max_a/box_size)+1;
		int area_distribution[a_divisions];
		for(int i = 0 ; i < a_divisions; i++)
			area_distribution[i] = 0;
		for(int i = 0 ; i < n; i++){
			area_distribution[(int) ((face[i]->area)/box_size)]++;
		}
		for(int i = 0 ; i < a_divisions; i++)
			*out << area_distribution[i] << " ";
		*out << endl;
	}
	
	compute_curvature(n,face);
	
	// see what the assignment of curvature is
	//print_triangles(triangles, VRML1);
		
	/*int max_degree, min_degree;
	{
		for(int i = 0 ; i < n; i++){
			int num_neighbors = face[i]->num_neighbors;
			if(i == 0)
				max_degree = min_degree = num_neighbors;
			else{
				if(num_neighbors > max_degree)
					max_degree = num_neighbors;
				else if(num_neighbors < min_degree)
					min_degree = num_neighbors;
			}
		}
		*out << "neighbors min " << min_degree << " max " << max_degree << endl;
		
		int distribution[max_degree+1];
		for(int i = 0 ; i < max_degree+1; i++)
			distribution[i] = 0;
		for(int i = 0 ; i < n; i++)
			distribution[face[i]->num_neighbors]++;
		for(int i = min_degree ; i <= max_degree; i++)
			*out << distribution[i] << " ";
		*out << endl;
	}*/
	
	if(coarse_convex_faces.size() + coarse_concave_faces.size() > 0){
		int num_concave_faces = coarse_concave_faces.size();
		int num_convex_faces = coarse_convex_faces.size();
		int K = (int) (NUM_CLUSTERS(num_concave_faces,num_convex_faces));
		*out << "N " << faces.size() << " K " << K << endl;
		
		switch(algorithm){
			case KMEANS:
				cluster_points_KMedians(face,n,K);
			break;
			case DENSITY_CLUSTER:
				cluster_points_agglomerative(face,n,K);
			break;
		}
	}	
	
	int nc = clusters.size();
	Face* cluster[nc + 1];
	count = 0;
	
	int num_positive_curvature = 0, num_negative_curvature = 0;
	float sum_abs_curvature = 0;
	// compute the residues that each face is incident on
	// compute normals that use larger snapshots
	// assume that normals and curvatures were already calculated
	for(vector<Cluster*>::iterator itr = clusters.begin(); itr != clusters.end(); itr++){
		Cluster *cf = (Cluster*) *itr;
		//cout << "count " << count << endl;
		cf->compute_coarsenormal_radius(face,n);
		hash_set<int,hash<int>, eqint> merged_residue_indices;
		for(vector<unsigned short>::iterator aitr = cf->face_atom_indices.begin(); aitr != cf->face_atom_indices.end(); aitr++){
			int acindex = *aitr;
			Atom *a = atom[acindex];
			for(hash_map<char, Molecule*, hash<char>, eqint>::iterator mitr = molecules.begin(); mitr != molecules.end(); mitr++){
				Molecule *m = mitr->second;
				if(m->atom.count(a->index) > 0 && m->atom[a->index]->cindex == acindex){
					switch(m->type){
						case PROTEIN:
							int mrindex = ((Protein*) m)->aminoacid[a->monoindex]->cindex;
							merged_residue_indices.insert(mrindex);
							break;
					}
					break;
				}
			}
		}
		for(hash_set<int,hash<int>, eqint>::iterator mitr = merged_residue_indices.begin(); mitr != merged_residue_indices.end(); mitr++)
			cf->face_monomer_indices.push_back(*mitr);
		cluster[count++] = cf;
		if(cf->curvature < 0){
			num_negative_curvature++;
			sum_abs_curvature -= cf->curvature;
		} else if(cf->curvature > 0){
			num_positive_curvature++;
			sum_abs_curvature += cf->curvature;
		}
	}
	*out << "#clusters curvature +ve " << num_positive_curvature << " -ve " << num_negative_curvature << " avg || " << sum_abs_curvature/n << endl;

	/*{
		float cd[n+1][n+1];
		for(int i = 0 ; i < n; i++){
			float area = cluster[i]->area;
			if(i == 0)
				min_a= max_a= area;
			else{
				if(area > max_a)
					max_a = area;
				if(area < min_a)
					min_a = area;
			}
			for(int j = i+1; j < n; j++){
				float dis = cd[i][j] = Vector::distance(*(cluster[j]->point),*(cluster[i]->point)); 
				if(i == 0 && j == 1){
					min_distance = max_distance = dis;
				} else {
					if(dis < min_distance)
						min_distance = dis;
					if(dis > max_distance)
						max_distance = dis;
				}
			}
		}
		*out << "a min " << min_a << " max " << max_a << endl;
		*out << "d min " << min_distance << " max " << max_distance << endl;
		
		d_divisions = max_distance/box_size + 1;
		int cd_distribution[d_divisions];
		for(int i = 0 ; i < d_divisions; i++)
			cd_distribution[i] = 0;
		a_divisions = max_a/1.0;
		int carea_distribution[a_divisions];
		for(int i = 0 ; i < a_divisions; i++)
			carea_distribution[i] = 0;
			
		for(int i = 0 ; i < n; i++){
			carea_distribution[(int) ((cluster[i]->area)/1.0)]++;
			for(int j = i+1; j < n; j++)
				cd_distribution[(int) (cd[i][j]/box_size)]++;
		} 
	
		for(int i = 0 ; i < a_divisions; i++)
			*out << carea_distribution[i] << " ";
		*out << endl << "d: ";
		
		for(int i = 0 ; i < d_divisions; i++)
			*out << cd_distribution[i] << " ";
		*out << endl;
	}*/
}

/*
 * Use MSMS to triangulate the surface. Compute critical points from the triangulation
 */
void Complex::compute_points_msms(Vector *grid_origin, int grid_num_xdivisions, int grid_num_ydivisions, int grid_num_zdivisions, GridCell **grid, int algorithm){
	// write out the coordinates and the radii
	string filename = string(pdbcode) + "_" + chains + ".xyzr";
	fstream msmsinput;
	msmsinput.open(filename.c_str(), fstream::out);
	for(int i = 0 ; i < num_atoms; i++){
		Atom *a = atom[i];
		msmsinput << a->position->x << "\t" << a->position->y << "\t" << a->position->z << "\t" << a->radius*MSMS_EXPANSION_FACTOR << endl;
	}
	msmsinput.close();
	
	// run msms
	char msms_command[512];
	sprintf(msms_command, "%s/%s -if %s -of %s_%s -probe_radius %f -density %f > %s_%s.msmslog",piedock_home.c_str(),
		MSMS_EXECUTABLE,filename.c_str(),pdbcode.c_str(),chains.c_str(),probe_radius, msms_density,pdbcode.c_str(),chains.c_str());
	
	int iret = system(msms_command);
	*out << "msms successful? " << iret << endl;
	if(iret != 0){
		*out << "ERROR: MSMS unsuccessful" << endl;
		exit(-1);
	}

	//read the vertices and faces
	filename = string(pdbcode) + "_" + chains + ".vert";
	fstream vertin;
	vertin.open(filename.c_str(), fstream::in);
	char buf[8192];
	vertin.getline(buf,8192); vertin.getline(buf,8192); vertin.getline(buf,8192);
	int num_vertices;
	{
		stringstream ss(buf,stringstream::in);
		ss >> num_vertices;
	}
	Vector *vertex[num_vertices+1], *normal[num_vertices+1];
	for(int i = 0 ; i < num_vertices; i++){
		vertin.getline(buf,8192);
		stringstream ss(buf,stringstream::in);
		float x,y,z;
		ss >> x;
		ss >> y;
		ss >> z;
		vertex[i] = new Vector(x,y,z);
		ss >> x;
		ss >> y;
		ss >> z;
		normal[i] = new Vector(x,y,z);
	}
	vertin.close();
	
	filename = string(pdbcode) + "_" + chains + ".face";
	fstream trianglein;
	trianglein.open(filename.c_str(), fstream::in);
	trianglein.getline(buf,8192); trianglein.getline(buf,8192); trianglein.getline(buf,8192);
	int num_triangles;
	{
		stringstream ss(buf,stringstream::in);
		ss >> num_triangles;
	}
	vector<Face*> faces;
	for(int i = 0 ; i < num_triangles; i++){
		trianglein.getline(buf,8192);
		stringstream ss(buf,stringstream::in);
		int vi3,vi1,vi2;
		ss >> vi1;
		ss >> vi2;
		ss >> vi3;
		// indexing in the file starts at 1
		vi1--; vi2--; vi3--;
		
		Triangle *tr = new Triangle();
		tr->v1 = vertex[vi1];
		tr->v2 = vertex[vi2];
		tr->v3 = vertex[vi3];
		
		Vector *centroid = new Vector((*(tr->v1) + *(tr->v2) + *(tr->v3))*(1.0/3.0));
		tr->point = centroid;
		Vector n = (*(tr->v2) - *(tr->v1)).cross(*(tr->v3) - *(tr->v1));
		tr->area = 0.5 * n.norm();
		if(tr->area > 0){
			n.normalize();
			if(n.dot(*(normal[vi1]) + *(normal[vi2]) + *(normal[vi3])) < 0)
				*out << "opposite normals" << endl;
			tr->normal = new Vector(n);
			tr->radius = 0.0;
			triangles.push_back(tr);
			faces.push_back(tr);
			//*out << "(" << tr->point->x << "," << tr->point->y << "," << tr->point->z << ")("
			//	<< "(" << tr->normal->x << "," << tr->normal->y << "," << tr->normal->z << ") " << tr->area << endl;
		}
	}
	trianglein.close();
	
	eliminate_triangles_in_holes();
	
	Face* face[faces.size() + 1];
	int nfaces = 0;
	sa = 0;
	float saall=0;
	for(vector<Face*>::iterator itr = faces.begin(); itr != faces.end(); itr++){
		Face *f = (Face*) *itr;
		saall += f->area;
		if(!(f->is_buried)){
			sa += f->area;
			face[nfaces++] = f;
			triangle_centroids.push_back(((Triangle*) f)->point);
		}
	}
	*out << "exposed area " << sa << " <area> " << sa/faces.size() << " saall " << saall << endl;
	
	compute_curvature(nfaces,face);
	compute_sphere_based_curvature(grid_origin, grid_num_xdivisions, grid_num_ydivisions, grid_num_zdivisions, grid, nfaces,face);
	
	int K;
	K = sa/10.0;//((float) num_aminoacids)*4.5;
	//K = sa/12.0;
	*out << "sa " << sa << " K " << K << endl;
		
	switch(algorithm){
		case KMEANS:
			cluster_points_KMedians(face,nfaces,K);
		break;
		case DENSITY_CLUSTER:
			cluster_points_agglomerative(face,nfaces,K);
		break;
	}	
	
	int nc = clusters.size();
	Face* cluster[nc + 1];
	int count = 0;
	
	int num_positive_curvature = 0, num_negative_curvature = 0;
	float sum_abs_curvature = 0;
	// compute the residues that each face is incident on
	// compute normals that use larger snapshots
	// assume that normals and curvatures were already calculated
	for(vector<Cluster*>::iterator itr = clusters.begin(); itr != clusters.end(); itr++){
		Cluster *cf = (Cluster*) *itr;
		//cout << "count " << count << endl;
		cf->compute_coarsenormal_radius(face,nfaces);
		hash_set<int,hash<int>, eqint> merged_residue_indices;
		for(vector<unsigned short>::iterator aitr = cf->face_atom_indices.begin(); aitr != cf->face_atom_indices.end(); aitr++){
			int acindex = *aitr;
			Atom *a = atom[acindex];
			for(hash_map<char, Molecule*, hash<char>, eqint>::iterator mitr = molecules.begin(); mitr != molecules.end(); mitr++){
				Molecule *m = mitr->second;
				if(m->atom.count(a->index) > 0 && m->atom[a->index]->cindex == acindex){
					switch(m->type){
						case PROTEIN:
							int mrindex = ((Protein*) m)->aminoacid[a->monoindex]->cindex;
							merged_residue_indices.insert(mrindex);
						break;
					}
					break;
				}
			}
		}
		for(hash_set<int,hash<int>, eqint>::iterator mitr = merged_residue_indices.begin(); mitr != merged_residue_indices.end(); mitr++)
			cf->face_monomer_indices.push_back(*mitr);
		cluster[count++] = cf;
		if(cf->curvature < 0){
			num_negative_curvature++;
			sum_abs_curvature -= cf->curvature;
		} else if(cf->curvature > 0){
			num_positive_curvature++;
			sum_abs_curvature += cf->curvature;
		}
	}
	*out << "#clusters curvature +ve " << num_positive_curvature << " -ve " << num_negative_curvature << " avg || " << sum_abs_curvature/nfaces << endl;
}

Cluster::Cluster(){
	type = CLUSTER_FACE;
}

Cluster::Cluster(Face *face){
	faces.push_back(face);
	type = CLUSTER_FACE;
}

Cluster::Cluster(char *buf){
	type = CLUSTER_FACE;
	stringstream line(buf,stringstream::in);
	float x,y,z;
	line >> x;
	line >> y;
	line >> z;
	point = new Vector(x,y,z);
	line >> x;
	line >> y;
	line >> z;
	normal = new Vector(x,y,z);
	line >> x;
	line >> y;
	line >> z;
	coarser_normal = new Vector(x,y,z);
	line >> normalstability;
	line >> curvature;
	line >> radius;
	line >> area;
	
	int num_atoms_on_face;
	line >> num_atoms_on_face;
	for(int j = 0 ; j < num_atoms_on_face ; j++){
		int atom_index;
		line >> atom_index;
		face_atom_indices.push_back(atom_index);
	}
	int num_residues_on_face;
	line >> num_residues_on_face;
	for(int j = 0 ; j < num_residues_on_face ; j++){
		int residue_index;
		line >> residue_index;
		face_monomer_indices.push_back(residue_index);
	}
}

void Cluster::print_details(ostream *out){
	//*out << type << " ";
	Vector v = *point;
	*out << v.x << " " << v.y << " " << v.z << " ";
	v = *normal;
	*out << v.x << " " << v.y << " " << v.z << " ";
	v = *coarser_normal;
	*out << v.x << " " << v.y << " " << v.z << " ";
	*out << normalstability << " " << curvature << " " << radius << " " << area << " ";
	
	*out << face_atom_indices.size();
	for(vector<unsigned short>::iterator itr = face_atom_indices.begin(); itr != face_atom_indices.end(); itr++)
		*out << " " << *itr;
	
	*out << " " << face_monomer_indices.size();
	for(vector<unsigned short>::iterator itr = face_monomer_indices.begin(); itr != face_monomer_indices.end(); itr++)
		*out << " " << *itr;
	*out << endl;
}
	
/*bool Cluster::can_accomodate(Face *f){
	vector<Face *>::iterator itr = faces.begin(), end = faces.end();
	while(itr != end){
		float d = Vector::distance((*itr)->point, f->point);
		if(d > CONCAVE_FACE_MAX_CLUSTER_RADIUS)
			return false;
		itr++;
	}
	return true;
}*/

/* normal is no longer a unit normal
 * computing the atoms that are part of the cluster
*/
void Cluster::compute_area_point_normal_curvature_atoms(){
	int num_faces = faces.size();
	area = 0;
	curvature = 0;
	if(num_faces > 0){
		Vector c = Vector(0,0,0);
		Vector n = Vector(0,0,0);
		for(vector<Face *>::iterator itr = faces.begin(); itr != faces.end(); itr++){
			Face* f = *itr;
			c = c + *(f->point) * f->area;
			if(f->normal != NULL){
				n = n + *(f->normal) * f->area; //(f->curvature * f->area);
				//*out << "(" << f->normal->x << "," << f->normal->y << "," << f->normal->z << ")," << f->area
				//	<< "-(" << n.x << "," << n.y << "," << n.z << ")\t";
			}
			area += f->area;
		}
	
		Vector *mean_point = new Vector(c * (1.0/area));
		normal = new Vector(n);
		if(normal->norm() > 0)
			normal->normalize();
			
		float dmin = MAX_D_USED_FOR_NORMAL;
		median = NULL;
		for(vector<Face *>::iterator itr = faces.begin(); itr != faces.end(); itr++){
			Face* f = *itr;
			float d = Vector::distance(mean_point,f->point);
			if(d < dmin){
				dmin = d;
				median = f;
			}
		}
		if(median != NULL)
			point = median->point;
			
		// curvature is the negative of the zscore of the outdegree
		/*long out_degree = 0;
		for(vector<Face *>::iterator itr = faces.begin(); itr != faces.end(); itr++){
			Face* f = *itr;
			out_degree += f->num_neighbors;
			for(vector<Face *>::iterator itr2 = itr+1; itr2 != faces.end(); itr2++){
				Face* f2 = *itr2;
				if(Vector::distance(f->point,f2->point) < NEIGHBOR_CUTOFF)
					out_degree--;
			}
		}
		*/
		
		for(vector<Face *>::iterator itr = faces.begin(); itr != faces.end(); itr++){
			Face* f = *itr;
			curvature += (f->area/area)*f->curvature;
		}
		//*out << num_faces << " " << area << " (" << point->x  << "," << point->y << ","  << point->z << ") ("
		//	<< normal->x << "," << normal->y << "," << normal->z<< ") " << curvature << endl;
	}	
}

void Cluster::compute_coarsenormal_radius(Face **allfaces, int N){
	Vector n = Vector(0,0,0);
	float a = 0;
	for(int fi = 0; fi < N; fi++){
		Face *f = allfaces[fi];
		if(f->normal != NULL){
			float d = Vector::distance(point,f->point);
			if(d < MAX_D_USED_FOR_NORMAL){
				n = n + *(f->normal) * f->area; //(f->curvature * f->area);
				a += f->area;
			}
		}
	}
	coarser_normal = new Vector(n);
	normalstability = coarser_normal->norm()/a;
	coarser_normal->normalize();

	//*out << normalstability << " " << normal->dot(*coarser_normal) << endl;
	// use the coarser normal for the direction of the normal
	float A = 0, B = 0;
	for(vector<Face *>::iterator fitr = faces.begin(); fitr != faces.end(); fitr++){
		Face *f = *fitr;
		A += (*point - *(f->point)).dot(*point - *(f->point));
		B += (*point - *(f->point)).dot(*coarser_normal);
	}
	if(faces.size() > 1){
	 	if(A != 0 && B != 0)
			radius = -A/(2*B);
		else
			radius = 0;
	} else{
		vector<Face *>::iterator fitr = faces.begin();
		Face *f = *fitr;
		radius = f->radius;
	}
	//*out << faces.size() << " " << A << " " << B << " " << radius << endl;
	
	hash_set<int, hash<int>, eqint> merged_atom_indices;
	for(vector<Face *>::iterator itr = faces.begin(); itr != faces.end(); itr++){
		Face* f = *itr;
		for(vector<unsigned short>::iterator iitr = f->face_atom_indices.begin(); iitr != f->face_atom_indices.end(); iitr++)
			merged_atom_indices.insert(*iitr);
	}
	face_atom_indices.clear();
	for(hash_set<int, hash<int>, eqint>::iterator mitr = merged_atom_indices.begin(); mitr != merged_atom_indices.end(); mitr++)
		face_atom_indices.push_back(*mitr);
}

/*
 * Assumes that the centroid and the normal are updated
 * computes d(c,p_i)^2 + (area_i/area)*d(n,n_i)^2
 * computes d(c,p_i)^2 + (area_i/area)*|c - c_i|
 */
float Cluster::get_spread(){
	float spread = 0;
	for(vector<Face *>::iterator itr = faces.begin(); itr != faces.end(); itr++){
		Face* f = *itr;
		float dij = (*point - *(f->point)).norm_squared();
		if(f->normal != NULL){
			Vector n = *(f->normal);
			dij += DISTANCE_NORMAL_WEIGHT*(n - *normal).norm_squared();
		}
		
		float ccve = curvature;
		if(ccve > 1) ccve = 1;
		if(ccve < -1) ccve = -1;
		float absccve = ABS(ccve);
		float fac1 = (1 - absccve)*(1 - absccve);
		float fac = (fac1 > 0.5) ? fac1 : 0.5;
		
		float fcve = f->curvature;
		if(fcve > 1) fcve = 1;
		if(fcve < -1) fcve = -1;
		
		float delta_curvature = ccve - fcve;
		delta_curvature = ABS(delta_curvature);
		dij += DISTANCE_CURVATURE_WEIGHT*fac*delta_curvature;
		
		spread += (f->area)*dij;
	}
	//*out << spread << endl;
	return spread;
}

/*
 * Computes the (size of intersection)/(size of union)
 * */
float cluster_density(Cluster* c1, Cluster* c2){
	hash_set<int, hash<int>, eqint> r1;
	for(vector<unsigned short>::iterator itr = c1->face_monomer_indices.begin(); itr != c1->face_monomer_indices.end(); itr++)
		r1.insert(*itr);
	
	int usize = r1.size(), isize = 0;
	for(vector<unsigned short>::iterator itr = c2->face_monomer_indices.begin(); itr != c2->face_monomer_indices.end(); itr++){
		if(r1.count(*itr) == 0)
			usize++;
		else
			isize++;
	}
	
	return (isize/usize);
}

/*void Cluster::add(Cluster* c){
	int num_faces = faces.size();
	int c_num_faces = c->faces.size();
	*point = (*point) * area + *(c->point) * (c->area);
	*normal = (*normal) * area + *(c->normal) * (c->area);
	area += c->area;
	
	*point = *(point) * (1.0/area); 
	*normal = *(normal) * (1.0/area);
	normal->normalize();
	
	for(vector<Face *>::iterator itr = c->faces.begin(); itr != c->faces.end(); itr++)
		faces.push_back(*itr);
	
	hash_set<int, hash<int>, eqint> merged_atom_indices;
	for(vector<int>::iterator itr = face_atom_indices.begin(); itr != face_atom_indices.end(); itr++)
		merged_atom_indices.insert(*itr);
	for(vector<int>::iterator itr = c->face_atom_indices.begin(); itr != c->face_atom_indices.end(); itr++)
		merged_atom_indices.insert(*itr);			
	face_atom_indices.clear();
	
	for(hash_set<int, hash<int>, eqint>::iterator mitr = merged_atom_indices.begin(); mitr != merged_atom_indices.end(); mitr++)
		face_atom_indices.push_back(*mitr);
	
	hash_set<int, hash<int>, eqint> merged_residue_indices;
	for(vector<int>::iterator itr = face_residue_indices.begin(); itr != face_residue_indices.end(); itr++)
		merged_residue_indices.insert(*itr);
	for(vector<int>::iterator itr = c->face_residue_indices.begin(); itr != c->face_residue_indices.end(); itr++)
		merged_residue_indices.insert(*itr);
	face_residue_indices.clear();
	
	for(hash_set<int, hash<int>, eqint>::iterator mitr = merged_residue_indices.begin(); mitr != merged_residue_indices.end(); mitr++)
		face_residue_indices.push_back(*mitr);
}*/

/*
points and normals - KMedian, coarser normals - KMeans but 
*/
void Complex::cluster_points_KMedians(Face ** face, int n, int K){
	Cluster *cluster[KMEDIANS_TRIALS][K+1];
	
	int index = 0, best_trial;
	float min_spread;
	srand(time(0));
	for(int trial = 0 ; trial < KMEDIANS_TRIALS ; trial++){
		// set initial centroids
		index = 0;
		bool taken[n];
		for(int i = 0 ; i < n; i++)
			taken[i] = false;
	
		while(index < K){
			int p = rand()%(n-index);
			int k = 0;
			int permuted_index;
			for(int j = 0 ; j < n ; j++)
				if(!taken[j]){
					if(k++ == p){
						taken[j] = true;
						permuted_index = j;
						break;
					}
				}
						
			if(face[permuted_index]->normal != NULL){
				cluster[trial][index] = new Cluster(face[permuted_index]);
				cluster[trial][index]->compute_area_point_normal_curvature_atoms();
				index++;
			}
		}
		*out << "initialized means " << index << endl; out->flush();
		
		int num_iterations = 0;
		float spread, new_spread = 0;
		bool done = false;
		while(!done){
			done = false;
			spread = new_spread;
			new_spread = 0;
			
			//do not associate any points with the cluster
			for(int i = 0 ; i < K; i++){
				cluster[trial][i]->faces.clear();
			}
	
			// compute the closest centroid for each point and add it to the corresponding cluster
			float dij;
			for(int i = 0 ; i < n; i++){
				Face *f = face[i];
				//*out << i << " " << f->point << " " << f->normal << endl; out->flush();
				float closest_distance;
				int closest_cluster;
				for(int j = 0 ; j < K; j++){
					dij = (*(cluster[trial][j]->point) - *(f->point)).norm_squared();
					if(f->normal != NULL){
						Vector n = *(f->normal);
						dij += DISTANCE_NORMAL_WEIGHT*(n - *(cluster[trial][j]->normal)).norm_squared();
					}
					float ccve = cluster[trial][j]->curvature;
					if(ccve > 1) ccve = 1;
					if(ccve < -1) ccve = -1;
					float absccve = ABS(ccve);
					float fac1 = (1 - absccve)*(1 - absccve);
					float fac = (fac1 > 0.5) ? fac1 : 0.5;
					
					float fcve = f->curvature;
					if(fcve > 1) fcve = 1;
					if(fcve < -1) fcve = -1;
					
					float delta_curvature = ccve - fcve;
					delta_curvature = ABS(delta_curvature);
					dij += DISTANCE_CURVATURE_WEIGHT*fac*delta_curvature;

					if(j == 0){
						closest_distance = dij;
						closest_cluster = 0;
					} else if(dij < closest_distance){
						closest_distance = dij;
						closest_cluster = j;
					}
				}
				cluster[trial][closest_cluster]->faces.push_back(f);
				//*out << " done face " << i << " " << closest_distance << endl;
			}
			
			// update the cluster centroids
			for(int i = 0 ; i < K; i++){
				cluster[trial][i]->compute_area_point_normal_curvature_atoms();
				new_spread += cluster[trial][i]->get_spread();
				//*out << new_spread << endl;
			}
			
			num_iterations++;
			done = (num_iterations > 1 && (spread - new_spread)/spread <= .01);
			
			*out << "#iterations " << num_iterations << " spread " << new_spread << endl;
		}
		if(trial == 0 || new_spread < min_spread){
			best_trial = trial;
			min_spread = new_spread;
		}
	}

	for(int i = 0 ; i < K; i++)
		if(cluster[best_trial][i]->area > 0){
			clusters.push_back(cluster[best_trial][i]);
			//*out << i << " " << cluster[trial][i]->point->x << " " << cluster[trial][i]->point->y << " " << cluster[trial][i]->point->z << 
			//	" " << cluster[trial][i]->normal->x << " " << cluster[trial][i]->normal->y << " " << cluster[trial][i]->normal->z << " " << cluster[trial][i]->area << endl;
		}
		
	*out << "#clusters " << clusters.size() << endl;
}

/*
Ignoring saddle faces for now, 
Clustering algorithm - Aglomerative
*/
void Complex::cluster_points_agglomerative(Face ** faces, int N, int K){
/*	int num_concave_faces = concave_faces.size();
	int num_convex_faces = convex_faces.size();
	int K = (int) (NUM_CLUSTERS(num_concave_faces,num_convex_faces));
    *out << "K " << K << endl;
	
	int num_faces = num_concave_faces + num_convex_faces;
	Face* face[num_faces+1];
	
	int index = 0;
	vector<Face *>::iterator itr = convex_faces.begin(), end = convex_faces.end();
	while(itr != end){
		face[index++] = *itr;
		itr++;
	}
	itr = concave_faces.begin(); end = concave_faces.end();
	while(itr != end){
		face[index++] = *itr;
		itr++;
	}
	//*out << num_faces << " " << index << endl;
	
	Cluster *cluster[num_faces+1];
	for(int i = 0 ; i < num_faces; i++){
		cluster[i] = new Cluster(face[i]);
		cluster[i]->update();
	}
	//*out << "updated initial clusters " << num_faces << endl;
	
	float distance[num_faces+1][num_faces+1];
	//*out << "before initializing distances" << endl;
	//float density[num_faces+1][num_faces+1];
	
	for(int i = 0 ; i < num_faces; i++){
		for(int j = i+1 ; j < num_faces; j++){
			distance[i][j] = distance[j][i] = (*(cluster[j]->point) - *(cluster[i]->point)).norm() 
										+ DISTANCE_NORMAL_WEIGHT*(*(cluster[j]->normal) - *(cluster[i]->normal)).norm();
			//density[i][j] = density[j][i] = cluster_density(cluster[i],cluster[j]);
		}
	}
	//*out << "after initializing distances" << endl;
	
	for(int iteration = 0 ; iteration < num_faces-K;iteration++){
		// compute the pair of clusters that are closest to each other
		float closest_distance = -1;
		int c1,c2;
		for(int i = 0 ; i < num_faces; i++){
			if(cluster[i] != NULL){
				for(int j = i+1 ; j < num_faces; j++){
					if(cluster[j] != NULL){
						float d = distance[i][j] + DENSITY_WEIGHT*(1.0/cluster_density(cluster[i],cluster[j])-1);
						if(closest_distance == -1 || ((closest_distance != -1) && (d < closest_distance))){
							closest_distance = d;
							c1 = i; c2 = j;
						}
					}
				}
			}
		}
		
		// combine clusters
		int c1_size = cluster[c1]->faces.size();
		int c2_size = cluster[c2]->faces.size();
		cluster[c1]->add(cluster[c2]);
		cluster[c2] = NULL;
		
		for(int i = 0 ; i < num_faces; i++){
			if(cluster[i] != NULL && i != c1){
				distance[i][c1] = distance[c1][i] = (distance[i][c1]*c1_size + distance[i][c2]*c2_size)/(c1_size+c2_size); 
				//density[i][c1] = density[c1][i] = cluster_density(cluster[i],cluster[c1]);
			}
		}
		
		//*out << "#iterations " << iteration << endl;
	}

	for(int i = 0 ; i < num_faces; i++)
		if(cluster[i] != NULL && cluster[i]->area > 0){
			clusters.push_back(cluster[i]);
			//*out << i << " " << cluster[i]->point->x << " " << cluster[i]->point->y << " " << cluster[i]->point->z << 
			//	" " << cluster[i]->normal->x << " " << cluster[i]->normal->y << " " << cluster[i]->normal->z << " " << cluster[i]->area << endl;
		}
		
	*out << "#clusters " << clusters.size() << endl;*/
}

void Complex::cleanup_except_atom(int atom_index){
	int count = 0, i = 0;
	while(count < num_atoms){
		while(atom[i] == NULL) i++;
		
		Atom *a = atom[i];
		if(a->cindex != atom_index)
			delete a;
		count++; i++;
	}
	delete atom;
	
	for(hash_map<const long, Torus*, hash<long>, eqlong>::iterator titr = toruses.begin(); titr != toruses.end(); titr++){
        Torus* t = (Torus *) titr->second;
		delete t;
	}

	for(hash_map<const long, Triangle_Face*, hash<long>, eqlong>::iterator tritr = trianglefaces.begin(); tritr != trianglefaces.end(); tritr++){
       	Triangle_Face* tr = (Triangle_Face *) tritr->second;
		delete tr;
	}
}

GridCell::GridCell(unsigned int indx){
	index = indx;
	points = NULL;
	atom_contact_neighbors = NULL;
	atom_overlap_neighbors = NULL;
	num_points = 0;
	num_atom_contact_neighbors = 0;
	num_atom_overlap_neighbors = 0;
}

GridCell::GridCell(unsigned int indx, short type){
	index = indx;
	this->type = type;
	points = NULL;
	atom_contact_neighbors = NULL;
	atom_overlap_neighbors = NULL;
	num_points = 0;
	num_atom_contact_neighbors = 0;
	num_atom_overlap_neighbors = 0;
}

void GridCell::build_point_array(vector<int> *points_l){
	if(points_l != NULL){
		num_points = points_l->size();
		if(num_points > 0)
			points = (unsigned int*) malloc(sizeof(unsigned int)*num_points);
		int pi = 0;
		/*Lnode *currentnode = points_l->head->next;
		while(currentnode != points_l->head) {
			points[pi] = currentnode->data;
			currentnode = currentnode->next;
			pi++;
		}*/
		for(vector<int>::iterator itr = points_l->begin(); itr != points_l->end(); itr++)
			points[pi++] = *itr;
	}
}

void GridCell::build_atom_neighbor_array(vector<int>*atom_overlap_neighbors_l, vector<int>*atom_contact_neighbors_l){
	if(atom_overlap_neighbors_l != NULL){
		num_atom_overlap_neighbors = atom_overlap_neighbors_l->size();
		if(num_atom_overlap_neighbors > 0)
			atom_overlap_neighbors = (unsigned short*) malloc(sizeof(unsigned short)*num_atom_overlap_neighbors);
		int ani = 0;
		/*Lnode *currentnode = atom_neighbors_l->head->next;
		while(currentnode != atom_neighbors_l->head) {
			atom_neighbors[ani] = currentnode->data;
			currentnode = currentnode->next;
			ani++;
		}*/
		for(vector<int>::iterator itr = atom_overlap_neighbors_l->begin(); itr != atom_overlap_neighbors_l->end(); itr++)
			atom_overlap_neighbors[ani++] = *itr;
	}
	if(atom_contact_neighbors_l != NULL){
		num_atom_contact_neighbors = atom_contact_neighbors_l->size();
		if(num_atom_contact_neighbors > 0)
			atom_contact_neighbors = (unsigned short*) malloc(sizeof(unsigned short)*num_atom_contact_neighbors);
		int ani = 0;
		for(vector<int>::iterator itr = atom_contact_neighbors_l->begin(); itr != atom_contact_neighbors_l->end(); itr++)
			atom_contact_neighbors[ani++] = *itr;
	}
}

void get_next_line(fstream* fin, char* buf, int buf_size){
	fin->getline(buf,buf_size);
	while(buf[0] == '~' || buf[0] == '\r' || buf[0] == '\n' || buf[0] == '\0')
		fin->getline(buf,buf_size);
	//*out << "buf " << buf << " " << " " << (string(buf)).size() << endl;
}

float max_atom_radius;

float Atom::get_max_radius(){
	return max_atom_radius;
}

void read_molecule_config(short type){
	char buf[8192];
	
	piedock_home = string(getenv(string("PIEDOCK_HOME").c_str()));

	// read particle properties from moil
	{
		string filename;
		switch(type){
			case PROTEIN:
				filename = piedock_home + "/" + string(MOP_DIR) + "MOIL_JOIN_DARS.PROP";
			break;
			case NA : case DNA : case RNA :
				filename = piedock_home + "/" + string(MOP_DIR) + "OPLSUA_ALL.PROP";
				//filename = piedock_home + "/" + string(MOP_DIR) + "OPLSAA_ALL.PROP";
			break;
		}	
		//string filename = piedock_home + "/" + string(MOP_DIR) + "ALL.PROP";
		//string filename = piedock_home + "/" + string(MOP_DIR) + "AMBER.PROP";
		//string filename = piedock_home + "/" + string(MOP_DIR) + "T46.PROP";
	    fstream fprop(filename.c_str(), ios::in);
		//*out << filename << " " << fprop.is_open() << endl; out->flush();
	    
		// the first line is a dummy
		get_next_line(&fprop,buf,8192);
		
		get_next_line(&fprop,buf,8192);
		// read in the radius of the atom
		string s = buf;
		int atom_type_index = 0;
		while( s.find("DONE") == string::npos ){
			if(buf[0] != '~'){
				int start = s.find("PNAM=(") + 6;
				int end = s.find(")");
				string *patom = new string(s.substr(start,end - start));
		
				start = s.find("PSGM=") + 5;
				string sigma_s = s.substr(start, s.length());
				const char* sigma_c = sigma_s.c_str();
				int i;
				for(i = 0 ; i < sigma_s.length(); i++){
					char c = sigma_c[i];
					if( c != '.' && (c < '0' || c > '9') )
						break;
				}
				sigma_s = sigma_s.substr(0,i);
				sigma_c = sigma_s.c_str();
		
			    float sigma = particle_sigma[patom->c_str()] = atof(sigma_c);
			    if(particle_sigma.size() == 0 || sigma > max_atom_radius)
			    	max_atom_radius = sigma;
			    
			    start = s.find("PCHG=") + 5;
				string charge_s = s.substr(start, s.find("PEPS"));
				const char* charge_c = charge_s.c_str();
		
			    particle_charge[patom->c_str()] = atof(charge_c);
				
				start = s.find("PEPS=") + 5;
				string eps_s = s.substr(start, s.find("PSGM"));
				const char* eps_c = eps_s.c_str();
		
			    particle_eps[patom->c_str()] = atof(eps_c);
			    
				start = s.find("PMAS=") + 5;
				string mass_s = s.substr(start, s.find("PCHG"));
				const char* mass_c = mass_s.c_str();
		
			    particle_mass[patom->c_str()] = atof(mass_c);
			    
				//*out << *patom << " " << particle_sigma[patom->c_str()] << " " << particle_charge[patom->c_str()] << " " << particle_type.size() << " " << particle_charge.size() << endl; out->flush();
				
				/*if(atom_type_index < NUM_ATOM_TYPES-1 && */ if(patom->c_str()[0] != 'H'){
					particle_type[patom->c_str()] = atom_type_index++;
				}
			}				
			
			get_next_line(&fprop,buf,8192);
			s = buf;
		}
		fprop.close();
		
		//*out << "done with all.prop\n"; out->flush();
		
		switch(type){
			case PROTEIN:
				filename = piedock_home + "/" + string(MOP_DIR) + "MOIL_JOIN_DARS.MONO";
			break;
			case NA : case DNA : case RNA :
				filename = piedock_home + "/" + string(MOP_DIR) + "BQCAO.MONO";
				//filename = piedock_home + "/" + string(MOP_DIR) + "OPLSUA_ALL.MONO";
				//filename = piedock_home + "/" + string(MOP_DIR) + "OPLSAA_ALL.MONO";
			break;
		}
		//filename = piedock_home + "/" + string(MOP_DIR) + "Capri.MONO";
		//filename = piedock_home + "/" + string(MOP_DIR) + "ALL.MONO";
		//filename = piedock_home + "/" + string(MOP_DIR) + "JIAN.MONO";
		//filename = piedock_home + "/" + string(MOP_DIR) + "T46.MONO";
	    fstream fmono(filename.c_str(), ios::in);
		
		// the first line is a dummy
		get_next_line(&fmono,buf,8192);
		
		get_next_line(&fmono,buf,8192);
	
		// find the conversion for the amino acid
		s = buf;
		while( s.find("*EOD") == string::npos ){
			//*out << "s " << s << endl;
			int start = s.find("(") + 1;
			int end = s.find(")");
			string *aminoacid = new string(s.substr(start,end - start));
			hash_map<const char*, const char*, hash<const char*>, eqstr> aacid_atom_propname;
	
			//*out << *aminoacid << " " << amino_acid.count((string("TIP3")).c_str()) << endl;
		
			get_next_line(&fmono,buf,8192);
			
			s = buf;
			while( s.find("DONE") == string::npos ){
				start = s.find("(") + 1;
				end = s.find(")");
				string *atom = new string(s.substr(start,end - start));
	
				start = s.find("PRTC=(") + 6;
				end = s.rfind(")");
				string patom = s.substr(start,end - start);
				
				aacid_atom_propname[atom->c_str()] = (new string(patom.c_str()))->c_str();
				
				//*out << *aminoacid << " " << *atom << " " << patom <<  endl; out->flush();
	
				get_next_line(&fmono,buf,8192);
				s = buf;
			}
			
			monomer_atom_vs_particle[aminoacid->c_str()] = aacid_atom_propname;
			
			//*out << *aminoacid << endl;
			
			get_next_line(&fmono,buf,8192);
			s = buf;
	
			while( s.find("DONE") == string::npos ){
				get_next_line(&fmono,buf,8192);
				//*out << "looking for DONE " << string(buf) << endl; out->flush();
				s = buf;
			}
	
			get_next_line(&fmono,buf,8192);
			s = buf;
		}
		fmono.close();
	}
	
	monomer_name_vs_type["ALA"] = ALA;
	monomer_name_vs_type["ARG"] = ARG;
	monomer_name_vs_type["ASN"] = ASN;
	monomer_name_vs_type["ASP"] = ASP;
	monomer_name_vs_type["CYS"] = CYS;
	monomer_name_vs_type["GLN"] = GLN;
	monomer_name_vs_type["GLU"] = GLU;
	monomer_name_vs_type["GLY"] = GLY;
	monomer_name_vs_type["HIS"] = HIS;
	monomer_name_vs_type["ILE"] = ILE;
	monomer_name_vs_type["LEU"] = LEU;
	monomer_name_vs_type["LYS"] = LYS;
	monomer_name_vs_type["MET"] = MET;
	monomer_name_vs_type["MSE"] = MET;
	monomer_name_vs_type["PHE"] = PHE;
	monomer_name_vs_type["PRO"] = PRO;
	monomer_name_vs_type["SER"] = SER;
	monomer_name_vs_type["THR"] = THR;
	monomer_name_vs_type["TRP"] = TRP;
	monomer_name_vs_type["TYR"] = TYR;
	monomer_name_vs_type["VAL"] = VAL;
	//monomer_name_vs_type["N"] = AMIDE_N;
	//monomer_name_vs_type["O"] = CARBONYL_O;
	
	monomer_name_vs_type["DA"] = DA;
	monomer_name_vs_type["DT"] = DT;
	monomer_name_vs_type["DC"] = DC;
	monomer_name_vs_type["DG"] = DG;
	monomer_name_vs_type["A"] = RA;
	monomer_name_vs_type["U"] = RU;
	monomer_name_vs_type["C"] = C;
	monomer_name_vs_type["G"] = G;
	monomer_name_vs_type["RA"] = RA;
	monomer_name_vs_type["RU"] = RU;
	monomer_name_vs_type["RC"] = C;
	monomer_name_vs_type["RG"] = G;
	
	/*{
		//fstream sizout("delphi.siz", fstream::out), crgout("delphi.crg", fstream::out);
		//sizout << "atom__res_radius_" << endl;
		//crgout << "atom__resnumbc_charge_" << endl;
		hash_map<const char*, hash_map<const char*, const char*, hash<const char*>, eqstr>, hash<const char*>, eqstr>::iterator aaitr;
		hash_map<const char*, const char*, hash<const char*>, eqstr>::iterator itr;
		for(aaitr = monomer_atom_vs_particle.begin(); aaitr != monomer_atom_vs_particle.end(); aaitr++){
			string aaname = string((char*) (aaitr->first));
			if(aaname.size() <= 3){
				float net_charge=0;
				hash_map<const char*, const char*, hash<const char*>, eqstr> details = (hash_map<const char*, const char*, hash<const char*>, eqstr>) (aaitr->second);
				for(itr = details.begin(); itr != details.end(); itr++){
					char buf[128];
					sprintf(buf,"%-6s%-4s%6.2f\n", (char*) (itr->first) , aaname.c_str(),particle_radius[(char*) itr->second] );
					//sizout << buf;
					sprintf(buf,"ATOM  %-4s%4.2f %s\n", (char*) (itr->first) , particle_radius[(char*) itr->second], aaname.c_str());
					cout << buf;
				
					sprintf(buf,"%-6s%-9s%5.3f\n", (char*) (itr->first) , aaname.c_str(), (float) (itr->second));
					//crgout << buf;
					net_charge += (float) (itr->second);
				}
				cout << aaname << " charge " << net_charge << endl;
			}
		}
		//sizout.close(); crgout.close();
	}*/
	
	//output monomers containing a particular atom type 
	/*cout << "particle_type_size " << particle_type.size() << endl;
	for(hash_map<const char*, float, hash<const char*>, eqstr>::iterator aitr = particle_radius.begin(); aitr != particle_radius.end(); aitr++){
		const char* aname = (const char*) aitr->first;
		cout << string(aname) << " - ";
		for(hash_map<const char*, hash_map<const char*, const char*, hash<const char*>, eqstr>, hash<const char*>, eqstr>::iterator aaitr = monomer_atom_vs_particle.begin(); 
		 aaitr != monomer_atom_vs_particle.end(); aaitr++){
			hash_map<const char*, const char*, hash<const char*>, eqstr> details = (hash_map<const char*, const char*, hash<const char*>, eqstr>) (aaitr->second);
			for(hash_map<const char*, const char*, hash<const char*>, eqstr>::iterator nitr = details.begin(); nitr != details.end(); nitr++){
				if(string(aname) == string((const char*) nitr->second))
					cout << string((const char*) aaitr->first) << "|" << string((const char*) nitr->first) << " ";
			}
		}
		cout << endl;
	}*/
	
	switch(NUM_ATOM_TYPES){
		case 18 :
			coarsen_atomtypesto18();
			break;
		case 20 :
			coarsen_atomtypesto20();
			break;
		case 32 :
			coarsen_atomtypesto32();
			break;
		case 46 :
			coarsen_atomtypesto46();
			break;
	}
};

void read_molecule_config(){
	read_molecule_config(PROTEIN);
}

void coarsen_atomtypesto18(){
	for(hash_map<const char*, short, hash<const char*>, eqstr>::iterator itr = particle_type.begin(); itr != particle_type.end(); itr++){
		itr->second = -1;
	}
	// cout << "considering 18 atom types\n";
	
	particle_type["NH"]=0;
	particle_type["CAH"]=1;
	particle_type["CO"]=2;
	particle_type["OC"]=3;
	particle_type["CAG"]=4;
	particle_type["CH2B"]=5;
	particle_type["CH2C"]=5;
	particle_type["CH2D"]=5;
	particle_type["CHT"]=5;
	particle_type["CBH"]=5;
	particle_type["CH2"]=5;
	particle_type["CBA"]=5;
	particle_type["NX"]=6;
	particle_type["CH2K"]=6;
	particle_type["CDK"]=7;
	particle_type["CX1"]=8;
	particle_type["OX1"]=8;
	particle_type["CSX1"]=8;
	particle_type["OSX1"]=8;
	particle_type["CR3"]=9;
	particle_type["NR2"]=9;
	particle_type["NAS"]=10;
	particle_type["CSX"]=10;
	particle_type["OSX"]=10;
	particle_type["CGN"]=10;
	particle_type["OD1N"]=10;
	particle_type["CR2"]=11;
	particle_type["NR1"]=11;
	particle_type["CH2S"]=12;
	particle_type["OH"]=12;
	particle_type["CGHS"]=13;
	particle_type["CGHP"]=13;
	particle_type["NDHP"]=13;
	particle_type["CHDP"]=13;
	particle_type["CHEP"]=13;
	particle_type["NEHP"]=13;
	particle_type["NDHS"]=13;
	particle_type["CHDH"]=13;
	particle_type["CHEH"]=13;
	particle_type["NEHS"]=13;
	particle_type["CZ"]=14;
	particle_type["CEY"]=14;
	particle_type["CR1"]=15;
	particle_type["CH2A"]=15;
	particle_type["SM"]=15;
	particle_type["CG"]=15;
	particle_type["CH2M"]=15;
	particle_type["CFH"]=15;
	particle_type["CGHT"]=15;
	particle_type["CTR"]=15;
	particle_type["CHTR"]=15;
	particle_type["CGTR"]=15;
	particle_type["CGL"]=15;
	particle_type["CGQ"]=15;
	particle_type["CG2T"]=15;
	particle_type["CH3"]=16;
	particle_type["CH3D"]=16;
	particle_type["CH3M"]=16;
	particle_type["CH3G"]=16;
	particle_type["SH"]=17;
	
	particle_type["CHPR"]=1;
	particle_type["SC"]=17;
	
	atom18_potential = (float**) malloc(18*sizeof(float*));
	for(int j = 0 ; j < 18; j++)
		atom18_potential[j] = (float*) malloc(18*sizeof(float));
		
	char buf[8192*32];
	string filename = piedock_home + "/" + string(DOCK_DIR) + "atom18.mat";
	filename = piedock_home +  "/" + string(DOCK_DIR) + "iterative_learn/novdw_atom18_contacts1/atom18.mat";
	fstream fatom18(filename.c_str(), ios::in);
	int line = 0;
	while(fatom18.good()){
    	fatom18.getline(buf,8192);
    	if(fatom18.gcount() > 0){
	    	stringstream ss (stringstream::in | stringstream::out);
	    	ss << buf;
	    	//cout << line << " " << string(buf) << endl; cout.flush();
	    	if(line == 0){
	    		string tag;	ss >> tag;
	    		ss >> vdw_weight;
	    		//cout << "VDW weight " << vdw_weight << endl;
	    	} else if ( line == 19) {
	    		string tag; ss >> tag;
	    		ss >> res_bkbn_vs_atomp_weight;
	    		cout << "res_bkbn_vs_atomp_weight " << res_bkbn_vs_atomp_weight << endl;
	    	} else {
	    		int type;
	    		ss >> type;
	    		for(int j = 0 ; j < 18; j++){
					ss >> atom18_potential[type][j];
//					/cout << type << " " << j << " pot " << atom18_potential[type][j] << endl; cout.flush();
	    		}
	    	}
	    	line++;
    	}
    }
    fatom18.close();	
};

void coarsen_atomtypesto20(){
	for(hash_map<const char*, short, hash<const char*>, eqstr>::iterator itr = particle_type.begin(); itr != particle_type.end(); itr++){
		itr->second = -1;
	}
	// cout << "considering 20 atom types\n"; cout.flush();
	
	particle_type["NH"]=0;
	particle_type["CAH"]=1;
	particle_type["CO"]=2;
	particle_type["OC"]=3;
	particle_type["CAG"]=4;
	particle_type["CH2B"]=5;
	particle_type["CH2C"]=18;
	particle_type["CH2D"]=5;
	particle_type["CHT"]=5;
	particle_type["CBH"]=5;
	particle_type["CH2"]=5;
	particle_type["CBA"]=5;
	particle_type["NX"]=6;
	particle_type["CH2K"]=6;
	particle_type["CDK"]=7;
	particle_type["CX1"]=8;
	particle_type["OX1"]=8;
	particle_type["CSX1"]=8;
	particle_type["OSX1"]=8;
	particle_type["CR3"]=9;
	particle_type["NR2"]=9;
	particle_type["NAS"]=10;
	particle_type["CSX"]=10;
	particle_type["OSX"]=10;
	particle_type["CGN"]=10;
	particle_type["OD1N"]=10;
	particle_type["CR2"]=11;
	particle_type["NR1"]=11;
	particle_type["CH2S"]=12;
	particle_type["OH"]=12;
	particle_type["CGHS"]=13;
	particle_type["CGHP"]=13;
	particle_type["NDHP"]=13;
	particle_type["CHDP"]=13;
	particle_type["CHEP"]=13;
	particle_type["NEHP"]=13;
	particle_type["NDHS"]=13;
	particle_type["CHDH"]=13;
	particle_type["CHEH"]=13;
	particle_type["NEHS"]=13;
	particle_type["CZ"]=14;
	particle_type["CEY"]=14;
	particle_type["CR1"]=15;
	particle_type["CH2A"]=15;
	particle_type["SM"]=19;
	particle_type["CG"]=15;
	particle_type["CH2M"]=15;
	particle_type["CFH"]=15;
	particle_type["CGHT"]=15;
	particle_type["CTR"]=15;
	particle_type["CHTR"]=15;
	particle_type["CGTR"]=15;
	particle_type["CGL"]=15;
	particle_type["CGQ"]=15;
	particle_type["CG2T"]=15;
	particle_type["CH3"]=16;
	particle_type["CH3D"]=16;
	particle_type["CH3M"]=16;
	particle_type["CH3G"]=16;
	particle_type["SH"]=17;
	
	particle_type["CHPR"]=1;
	particle_type["SC"]=17;
	
	atom20_potential = (float**) malloc(21*sizeof(float*));
	for(int j = 0 ; j < 21; j++)
		atom20_potential[j] = (float*) malloc(20*sizeof(float));
		
	char buf[8192*32];
	string filename = piedock_home + "/" + string(DOCK_DIR) + "atom20.mat";
	filename = piedock_home +  "/" + string(DOCK_DIR) + "/iterative_learn/vdw_atom20_contacts2/atom20.mat";
#ifdef DEV_VERSION 
#endif	
	fstream fatom20(filename.c_str(), ios::in);
	int line = 0;
	while(fatom20.good()){
    	fatom20.getline(buf,8192);
    	if(fatom20.gcount() > 0){
	    	stringstream ss (stringstream::in | stringstream::out);
	    	ss << buf;
	    	if(line == 0){
	    		string tag;	ss >> tag;
	    		ss >> vdw_weight;
	    	} else {
	    		int type;
	    		ss >> type;
	    		for(int j = 0 ; j < 20; j++)
					ss >> atom20_potential[type][j];
    		}
    		line++;
    	}
    }
    fatom20.close();	
};

void coarsen_atomtypesto32(){
	for(hash_map<const char*, short, hash<const char*>, eqstr>::iterator itr = particle_type.begin(); itr != particle_type.end(); itr++){
		itr->second = -1;
	}
	// cout << "considering 32 atom types\n";
	
	particle_type["NX"]=0;
	particle_type["NX1"]=0;
	particle_type["NH"]=1;
	particle_type["CO"]=2;
	particle_type["CGN"]=2;
	particle_type["CSX"]=2;
	particle_type["OC"]=3;
	particle_type["OD1N"]=3;
	particle_type["OSX"]=3;
	particle_type["CAH"]=4;
	particle_type["CACX"]=4;
	particle_type["CAG"]=4;
	particle_type["CH3"]=5;
	particle_type["CBA"]=5;
	particle_type["CH3G"]=5;
	particle_type["CG2T"]=5;
	particle_type["CH3D"]=5;
	particle_type["CH2"]=6;
	particle_type["CGQ"]=6;
	particle_type["CDK"]=6;
	particle_type["CBH"]=6;
	particle_type["CGL"]=6;
	particle_type["CR1"]=6;
	particle_type["CG"]=7;
	particle_type["CFH"]=7;
	particle_type["CEY"]=7;
	particle_type["CZ"]=8;
	particle_type["OH"]=9;
	particle_type["CGTR"]=10;
	particle_type["CTR"]=10;
	particle_type["CHTR"]=11;
	particle_type["CGHT"]=11;
	particle_type["NDHS"]=12;
	particle_type["CH2M"]=13;
	particle_type["SM"]=14;
	particle_type["SC"]=14;
	particle_type["CH3M"]=15;
	particle_type["CH2K"]=16;
	particle_type["CH2S"]=17;
	particle_type["CHT"]=17;
	particle_type["CH2D"]=18;
	particle_type["CHPR"]=18;
	particle_type["CH2C"]=19;
	particle_type["CCH2"]=19;
	particle_type["SH"]=20;
	particle_type["CGHP"]=21;
	particle_type["CGHS"]=21;
	particle_type["CHDP"]=21;
	particle_type["CHDH"]=21;
	particle_type["NDHP"]=22;
	particle_type["NEHP"]=22;
	particle_type["NEHS"]=22;
	particle_type["CHEP"]=23;
	particle_type["CHEH"]=23;
	particle_type["CR2"]=24;
	particle_type["NR1"]=25;
	particle_type["CR3"]=26;
	particle_type["NR2"]=27;
	particle_type["NAS"]=28;
	particle_type["CH2B"]=29;
	particle_type["CH2A"]=29;
	particle_type["CX1"]=30;
	particle_type["CSX1"]=30;
	particle_type["OX1"]=31;
	particle_type["OSX1"]=31;
	
	
	atom32_dpotential = (float***) malloc(32*sizeof(float**));
	for(int j = 0 ; j < 32; j++){
		atom32_dpotential[j] = (float**) malloc(32*sizeof(float*));
		for(int k = 0 ; k < 32; k++)
			atom32_dpotential[j][k] = (float*) malloc(3*sizeof(float));
	}
	
	char buf[8192*32];
	string filename = piedock_home + "/" + string(DOCK_DIR) + "T32S3.pot";
	float in_jian_format[16*33*3];
	int index=0;
	 
	FILE *POT=fopen(filename.c_str(), "r");
	for (int i=0; i<16*33*3; i++)
    	fscanf(POT, "%f", &in_jian_format[i]);
	fclose(POT);
	
	/*	fstream fatom32(filename.c_str());
	while(fatom32.good()){
    	fatom32 >> in_jian_format[index++];
	    cout << index-1 << in_jian_format[index] << endl;
    }
    fatom32.close();*/
    
    for(int j = 0 ; j < 32; j++)
		for(int k = j ; k < 32; k++)
			for(int l=0; l < 3; l++){
				int index = k*(k+1)/2+j + l*16*33;
				atom32_dpotential[j][k][l] = atom32_dpotential[k][j][l] = (0-in_jian_format[index]);  	
				//cout << j << " " << k << " " << l << " " << atom32_dpotential[j][k][l] << endl;
			}
}

void coarsen_atomtypesto46(){
	for(hash_map<const char*, short, hash<const char*>, eqstr>::iterator itr = particle_type.begin(); itr != particle_type.end(); itr++){
		itr->second = -1;
	}
	// cout << "considering 46 atom types\n";
	
	particle_type["NX"]=0;
	particle_type["NX1"]=0;
	particle_type["NH"]=1;
	particle_type["CO"]=2;
	particle_type["CGN"]=2;
	particle_type["CSX"]=2;
	particle_type["OC"]=3;
	particle_type["OD1N"]=3;
	particle_type["OSX"]=3;
	particle_type["CAH"]=4;
	particle_type["CACX"]=4;
	particle_type["CAG"]=5;
	particle_type["CH3"]=6;
	particle_type["CBA"]=6;
	particle_type["CH3G"]=6;
	particle_type["CG2T"]=6;
	particle_type["CH3D"]=7;
	particle_type["CH2"]=8;
	particle_type["CGQ"]=8;
	particle_type["CDK"]=8;
	particle_type["CBH"]=9;
	particle_type["CGL"]=9;
	particle_type["CR1"]=10;
	particle_type["CG"]=11;
	particle_type["CFH"]=12;
	particle_type["CEY"]=12;
	particle_type["CZ"]=13;
	particle_type["OH"]=14;
	particle_type["CGTR"]=15;
	particle_type["CTR"]=16;
	particle_type["CHTR"]=17;
	particle_type["CGHT"]=18;
	particle_type["NDHS"]=19;
	particle_type["CH2M"]=20;
	particle_type["SM"]=21;
	particle_type["SC"]=21;
	particle_type["CH3M"]=22;
	particle_type["CH2K"]=23;
	particle_type["CH2S"]=24;
	particle_type["CHT"]=25;
	particle_type["CH2D"]=26;
	particle_type["CHPR"]=27;
	particle_type["CH2C"]=28;
	particle_type["CCH2"]=28;
	particle_type["SH"]=29;
	particle_type["CGHP"]=30;
	particle_type["CGHS"]=30;
	particle_type["CHDP"]=31;
	particle_type["CHDH"]=31;
	particle_type["NDHP"]=32;
	particle_type["NEHP"]=33;
	particle_type["NEHS"]=33;
	particle_type["CHEP"]=34;
	particle_type["CHEH"]=34;
	particle_type["CR2"]=35;
	particle_type["NR1"]=36;
	particle_type["CR3"]=37;
	particle_type["NR2"]=38;
	particle_type["NAS"]=39;
	particle_type["CH2B"]=40;
	particle_type["CH2A"]=41;
	particle_type["CX1"]=42;
	particle_type["CSX1"]=43;
	particle_type["OX1"]=44;
	particle_type["OSX1"]=45;
}


short get_aatype(char symbol){
	short aatype=-1;
	switch(symbol){
		case'A':
			aatype = ALA;
		break;
		case 'R':
			aatype = ARG; 
		break;
		case 'N':
			aatype = ASN; 
		break; 
		case 'D':
			aatype = ASP; 
		break; 
		case 'C':
			aatype = CYS; 
		break; 
		case 'Q':
			aatype = GLN; 
		break; 
		case 'E':
			aatype = GLU; 
		break;
		case 'G':
			aatype = GLY; 
		break; 
		case 'H':
			aatype = HIS; 
		break; 
		case 'I':
			aatype = ILE; 
		break; 
		case 'L':
			aatype = LEU; 
		break; 
		case 'K':
			aatype = LYS; 
		break; 
		case 'M':
			aatype = MET; 
		break; 
		case 'F':
			aatype = PHE; 
		break; 
		case 'P':
			aatype = PRO; 
		break; 
		case 'S':
			aatype = SER; 
		break; 
		case 'T':
			aatype = THR; 
		break; 
		case 'W':
			aatype = TRP; 
		break; 
		case 'Y':
			aatype = TYR; 
		break; 
		case 'V':
			aatype = VAL; 
		break;
	}
	return aatype;
}

float compute_rmsd(int num_points, Vector* reference_points, Vector* model_points){
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
	
	cubicsolve_rmsd(ref_xlist, mov_xlist, num_points,&rmsd);
	/*cout << " fast_rmsd " << rmsd << "\t";
	
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
	cout << " regular " << rmsd << endl;*/
	return rmsd;
}

float compute_lrmsd(int num_points, int ligand_start,Vector* reference_points, Vector* model_points, float *rmsd){
	float ref_xlist[num_points+1][3], mov_xlist[num_points+1][3];
	float mov_com[3], mov_to_ref[3], U[3][3];
	int num_ligand_points = num_points - ligand_start;	

	for(int i = 0; i < num_points; i++){
		Vector v = reference_points[i];
		ref_xlist[i][0] = v.x;
		ref_xlist[i][1] = v.y;
		ref_xlist[i][2] = v.z;
		
		v = model_points[i];
		mov_xlist[i][0] = v.x;
		mov_xlist[i][1] = v.y;
		mov_xlist[i][2] = v.z;
	
  		//printf("%d\t%f %f %f\t%f %f %f\n",i,ref_xlist[i][0],ref_xlist[i][1],ref_xlist[i][2],mov_xlist[i][0],mov_xlist[i][1],mov_xlist[i][2]);
	}
		
	// function writes on the arrays passed
	int rank;
	calculate_rotation_rmsd(ref_xlist, mov_xlist, num_points, mov_com, mov_to_ref,&rank,U,rmsd);
	
	// compute the rmsd between ligand points
	float lrmsd = 0;
	if(rank == 3){
		Vector ux = Vector(U[0][0],U[0][1],U[0][2]);
		ux.normalize();
		Vector uy = Vector(U[1][0],U[1][1],U[1][2]);
		uy.normalize();
		Vector uz = Vector(U[2][0],U[2][1],U[2][2]);
		uz.normalize();
		for(int i = 0 ; i < num_ligand_points; i++){
			Vector v = model_points[i+ligand_start];
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
		lrmsd = lrmsd/num_points;
		lrmsd = sqrt(lrmsd);
	}
	return lrmsd;
};

