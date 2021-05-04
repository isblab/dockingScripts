/* PIEDOCK                                                            *
 * Authors: D. V. S. Ravikant                                         *
 * (C) 2011 D. V. S. Ravikant and Ron Elber                           */

#include "utilities.hpp"

Vector::Vector(){
	};

Vector::Vector(float x, float y,float z){
		this->x = x;
		this->y = y;
		this->z = z;
	};

Vector::Vector(Vector* v){
		x = v->x;
		y = v->y;
		z = v->z;
	};

void Vector::toarray(float array[3]){
	array[0] = x;
	array[1] = y;
	array[2] = z;
};

float Vector::norm_squared(){
		return (x*x + y*y + z*z);
	};

float Vector::norm(){
		return sqrt(x*x + y*y + z*z);
	};

float Vector::distance(Vector* v1, Vector* v2){
		return distance(*v1,*v2);
	}
	
float Vector::distance(Vector v1, Vector v2){
		return (v1 - v2).norm();
	}	

float Vector::distance_squared(Vector* v1, Vector* v2){
		return distance_squared(*v1,*v2);
	}
	
float Vector::distance_squared(Vector v1, Vector v2){
		return (v1 - v2).norm_squared();
	}	

void Vector::normalize(){
		float n = norm();
		x /= n;
		y /= n;
		z /= n;
	};

float Vector::dot(Vector* v){
		return dot(*v);
	}

float Vector::dot(Vector v){
		return (x*v.x + y*v.y + z*v.z);
	};

Vector Vector::cross(Vector* v){
		return cross(*v);
	}
	
Vector Vector::cross(Vector v){
		float X,Y,Z;
		X = y*v.z - z*v.y;
		Y = z*v.x - x*v.z;
		Z = x*v.y - y*v.x;	
		return ( Vector(X,Y,Z) );
	};

Vector Vector::operator+(const Vector &v){
	        return (Vector(x + v.x, y + v.y, z + v.z));
	};

Vector Vector::operator-(const Vector &v){
	        return (Vector(x - v.x, y - v.y, z - v.z));
	};

Vector Vector::operator*(const float s){
	        return (Vector(x*s, y*s, z*s));
	};

bool Vector::operator==(const Vector &v){
		return ((*this - v).norm_squared() < EPSILON);
	};

/*
 * assumes that v2 is counterclockwise with respect to v1 when viewed along the axis
 * returns a value in [0,2PI]
 */
float Vector::angle(Vector v1, Vector v2, Vector axis){
	v1.normalize();	
	v2.normalize();	

	Vector v = v1.cross(v2) ;
	float angle = asin(v.norm());

	if( v1.dot(v2) < 0)
		angle = PI - angle;

	if(v.dot(axis) > 0 && angle > EPSILON)
		angle = 2 * PI - angle;
	
	return angle;
};

Vector* Vector::find_closest(vector<Vector*> *set,Vector* point){
		find_closest(set,NULL,point);
	};
	
Vector* Vector::find_closest(vector<Vector*> *set,vector<Vector*> *exclude_set, Vector* point){
	if(set->size() == 0)
		return NULL;
	else {
		hash_map<const long, bool, hash<long>,eqlong> excluded_vertices;
		if(exclude_set != NULL){
			for(vector<Vector*>::iterator itr = exclude_set->begin(); itr != exclude_set->end(); itr++){
				excluded_vertices[(long) *itr] = true;
			}
		}
		
		Vector* result = NULL;
		float min_distance = -1;

		for(vector<Vector*>::iterator itr = set->begin(); itr != set->end(); itr++){
			Vector *v = *itr;
			if(excluded_vertices.count((long) v) == 0){
				float d = distance(v,point);
				if(min_distance == -1 || d < min_distance){
					min_distance = d;
					result = v;
				}
			}
		}
		
		return result;
	}
}

Lnode::Lnode(){
	//cout << "nlnode" << endl;
	next = NULL;
}

Lnode::~Lnode(){
	//cout << "dlnode" << endl;
	next = NULL;
}

NodeList::NodeList(){
	head = new Lnode(); 
	head->next = head;
};

NodeList::~NodeList(){
	Lnode *currentnode = head->next;
	while(currentnode != head) {
		Lnode* nextnode = currentnode->next;
		//cout << currentnode << " ";
		delete currentnode;
		currentnode = nextnode;
	}
	//cout << head << endl;
	delete head;
}

//template< class T >
void NodeList::add( int data ){
	Lnode *p = new Lnode();
	p->data = data;
	p->next = head->next;
	head->next = p;
};

//template< class T >
bool NodeList::empty(){
	return (head->next == head);
};

int NodeList::size(){
	int size = 0;
	Lnode *currentnode = head->next;
	while(currentnode != head) {
		size++;
		currentnode = currentnode->next;
	}
	return size;
}

/*
template< class T >
T List< T >::remove(){
	T data;
	Lnode< T > *node;

	if (sentinel->next == sentinel) {
		cerr << "ERROR: `remove' called with empty list.\n";
		exit(1);
	}

	node = sentinel->next;
	data = node->data;

	sentinel->next = node->next;
	delete node;

	return data;
}*/

GVertex::GVertex(){
	static long count = 0;
	this->id = count++;
}

GVertex::GVertex(long id){
	this->id = id;
}

void GVertex::add_neighbor(GVertex* u){
	// do not add self loop	
	if(u->id != id){
		if(neighbors.count(u->id)==0){
			neighbors[u->id] = u;
			//neighbor_weights[u->id] = w;
		} /*else {
			neighbor_weights[u->id] += w;
		}*/
	}
}

bool GVertex::is_neighbor(GVertex* u){
	return (neighbors.count(u->id)!=0);
}


int GVertex::get_degree(){
	return (int) neighbors.size();
}

GEdge::GEdge(GVertex *s, GVertex *t){
	this->s = s;
	this->d = d;
}
	
Graph::Graph(){
	num_vertices = 0;
}

void Graph::insert_vertex(GVertex* v){
	vertex[v->id] = v;
	num_vertices++;
}
	

int Graph::bfs(){
	if(vertex.size() > 0){
		long start_vertex = (vertex.begin())->second->id;
		bool* visited = (bool *) malloc(sizeof(bool)*num_vertices);
		long* tree = (long *) malloc(sizeof(long)*num_vertices);
		for(int i = 0; i < num_vertices; i++){
			visited[i] = false;
		}
		
		int result = bfs(start_vertex,visited,tree);
		free((void*)visited);
		free((void*)tree);
		return result;
	} else
		return 0;
}

/*
 * returns the size of the component 
 */
int Graph::bfs(long start_vertex, bool* visited, long* parent){
	GVertex* current = vertex[start_vertex];

	queue<GVertex*> frontier;
	frontier.push(current);
	parent[current->id] = current->id;
	int component_size = 0;	

	while(!frontier.empty()){
		current = frontier.front();
		frontier.pop();
		if(!visited[current->id]){
			visited[current->id] = true;
			component_size++;
			for(hash_map<long, GVertex*, hash<long>, eqlong>::iterator itr = current->neighbors.begin(); itr != current->neighbors.end(); itr++){
				GVertex* child = itr->second;
				if(!visited[child->id]){
					parent[child->id] = current->id;
					//cout << "added to q " << child->id << endl;
					frontier.push(child);
				}
			}
		}
	}
	return component_size;
}


void Graph::compute_connected_component(list<long>*component, long start_vertex){
	bool* visited = (bool *) malloc(sizeof(bool)*num_vertices);
	long* tree = (long *) malloc(sizeof(long)*num_vertices);
	for(int i = 0; i < num_vertices; i++){
		visited[i] = false;
	}
	tree[start_vertex] = start_vertex;
	int result = dfs(start_vertex,visited,tree);
	
	for(int i = 0; i < num_vertices; i++){
		if(visited[i])
			component->push_back(i);
	}
	
	free((void*)visited);
	free((void*)tree);
}


/*
 * returns the size of the component 
 */
int Graph::dfs(){
	if(vertex.size() > 0){
		long start_vertex = (vertex.begin())->second->id;
		bool* visited = (bool *) malloc(sizeof(bool)*num_vertices);
		long* tree = (long *) malloc(sizeof(long)*num_vertices);
		for(int i = 0; i < num_vertices; i++){
			visited[i] = false;
			//distances[i] = -1;
		}
		tree[start_vertex] = start_vertex;
		int result = dfs(start_vertex,visited,tree);
		free((void*)visited);
		free((void*)tree);
		return result;
	} else
		return 0;
}


int Graph::dfs(long start_vertex, bool* visited, long* parent){
	GVertex* current = vertex[start_vertex];
	visited[current->id] = true;
	int component_size = 1;

	for(hash_map<long, GVertex*, hash<long>, eqlong>::iterator itr = current->neighbors.begin(); itr != current->neighbors.end(); itr++){
		GVertex* child = itr->second;
		if(!visited[child->id]){
			component_size += dfs(child->id, visited, parent);
			parent[child->id] = current->id;
		}

	}
	
	return component_size;
}

/*
 * computes trail starting and ending at the start vertex 
 * during a dfs traversal without visiting vertices in visited
 * root childresult root childresult ...
 */
void Graph::compute_vertex_cycle(list<long> *trail, long start_vertex, bool* visited){
	GVertex* current = vertex[start_vertex];
	visited[current->id] = true;
	trail->push_back(current->id);

	for(hash_map<long, GVertex*, hash<long>, eqlong>::iterator itr = current->neighbors.begin(); itr != current->neighbors.end(); itr++){
		GVertex* child = itr->second;
		if(!visited[child->id]){
			compute_vertex_cycle(trail, child->id, visited);
			//trail->push_back(current->id);
		}
	}
}

void Graph::compute_vertex_cycle(list<long> *trail, long start_vertex){
	bool* visited = (bool *) malloc(sizeof(bool)*num_vertices);
	for(int i = 0; i < num_vertices; i++){
		visited[i] = false;
	}
	
	compute_vertex_cycle(trail, start_vertex, visited);
	
	free((void*)visited);
}

void Graph::compute_vertex_cycle(list<long> *trail){
	bool* visited = (bool *) malloc(sizeof(bool)*num_vertices);
	for(int i = 0; i < num_vertices; i++){
		visited[i] = false;
	}
	
	long start_vertex = (vertex.begin())->second->id;
	compute_vertex_cycle(trail, start_vertex, visited);
	
	free((void*)visited);
};

void Graph::compute_common_neighbors(list<long> *vertices,list<long> *common_neighbors){
	common_neighbors->clear();
	int count = 0;
	for(list<long>::iterator vitr = vertices->begin(); vitr != vertices->end(); vitr++){
		long vindex = *vitr;
		if(vertex.count(vindex) > 0){
			GVertex *v = vertex[vindex];
			if(count == 0){
				for(hash_map<long, GVertex*, hash<long>, eqlong>::iterator nitr = v->neighbors.begin(); nitr != v->neighbors.end(); nitr++){
					long nid = ((GVertex*) nitr->second)->id;
					common_neighbors->push_back(nid);
				}
			} else {
				for(list<long>::iterator cnitr = common_neighbors->begin(); cnitr != common_neighbors->end(); cnitr++){
					long nid = *cnitr;
					if(v->neighbors.count(nid) == 0)
						common_neighbors->erase(cnitr);
				}
			}
			count++;
		}
	}
};

void Graph::print_details(){
	cout << "#vertices " << num_vertices << " edges (adjacency list)" << endl;
	for(hash_map<long, GVertex*, hash<long>, eqlong>::iterator vitr = vertex.begin(); vitr != vertex.end(); vitr++){
		GVertex* v = vitr->second;
		cout << v->id << " " << vitr->first << " - ";
		for(hash_map<long, GVertex*, hash<long>, eqlong>::iterator itr = v->neighbors.begin(); itr != v->neighbors.end(); itr++){
			GVertex* neighbor = itr->second;
			cout << neighbor->id << " ";
		}
		cout << endl;
	}
};

void permute(int n, int *pn){
	bool taken[n];
	for(int i = 0 ; i < n; i++)
		taken[i] = false;
	for(int i = 0 ; i < n; i++){
		int index = rand()%(n-i);
		int k = 0;
		for(int j = 0 ; j < n ; j++)
			if(!taken[j]){
				if(k++ == index){
					taken[j] = true;
					pn[i] = j;
					break;
				}
			}
	}
};

void cubic_solve(double  a,double  b,double  c, double  d,int *solutions, double *x){
	long double    a1 = b/a, a2 = c/a, a3 = d/a;
	long double    Q = (a1*a1 - 3.0*a2)/9.0;
	long double R = (2.0*a1*a1*a1 - 9.0*a1*a2 + 27.0*a3)/54.0;
	double    R2_Q3 = R*R - Q*Q*Q;

	double    theta;

	if (R2_Q3 <= 0){
		*solutions = 3;
		theta = acos(R/sqrt(Q*Q*Q));
		x[0] = -2.0*sqrt(Q)*cos(theta/3.0) - a1/3.0;
		x[1] = -2.0*sqrt(Q)*cos((theta+2.0*PI)/3.0) - a1/3.0;
		x[2] = -2.0*sqrt(Q)*cos((theta+4.0*PI)/3.0) - a1/3.0;
      }	else	{
		*solutions = 1;
		x[0] = pow(sqrt(R2_Q3)+fabs(R), 1/3.0);
		x[0] += Q/x[0];
		x[0] *= (R < 0.0) ? 1 : -1;
		x[0] -= a1/3.0;
	}
}

void solve(float a[3][3], double d[3], float v[3][3]){
	// compute the characteristic polynomial
	double c0,c1,c2,c3;
	c0=1;
	c1=0-(a[0][0] + a[1][1] + a[2][2]);
	c2=a[0][0]*a[1][1] + a[0][0]*a[2][2] + a[1][1]*a[2][2] - (a[2][1]*a[1][2] + a[0][1]*a[1][0] + a[0][2]*a[2][0]);
	c3=(0-a[0][0])*(a[1][1]*a[2][2] - a[2][1]*a[1][2]) + a[0][1]*(a[1][0]*a[2][2]-a[2][0]*a[1][2])
		+(0-a[0][2])*(a[1][0]*a[2][1]-a[2][0]*a[1][1]);
	
	int nroots;
	cubic_solve(c0,c1,c2,c3,&nroots,d);
}

/*
int* get_component_distribution(){
	long start_vertex = vertices[0]->id;
	int connected_component_sizes[num_vertices];
	int num_connected_components = 0;
	bool completed = false;
	bool visited[num_vertices];
	for(int i = 0 ; i < num_vertices; i++)
		visited[i] = false;

	int i;
	int distances[num_vertices];
	while(!completed){
		connected_component_sizes[num_connected_components++] = bfs(start_vertex, (bool*) visited, distances);
		for(i = 0 ; i < num_vertices; i++)
			if(!visited[i]){
				start_vertex = i;
				break;
			}
		if(i == num_vertices)
			completed = true;
	}

	cout << "num vertices " << num_vertices << endl;
	int max_size = 0, second_largest = 0;

	for(i = 0 ; i < num_connected_components; i++){
		if(connected_component_sizes[i] > max_size){
			second_largest = max_size;
			max_size = connected_component_sizes[i];
		} else if(connected_component_sizes[i] > second_largest)
			second_largest = connected_component_sizes[i];
		
	}

	cout << "max component size " << max_size << endl;
	cout << "2nd largest component size " << second_largest << endl;

	int component_size_distribution[max_size+1];
	for(i = 0 ; i <= max_size; i++)
		component_size_distribution[i] = 0;
	for(i = 0 ; i < num_connected_components; i++){
		component_size_distribution[connected_component_sizes[i]]++;
	}

	int check = 0;
            for(i = 0 ; i <= max_size; i++){
                    if(i > 0 && component_size_distribution[i] > 0)
			cout << i << " " << component_size_distribution[i] << endl;
		check += i * component_size_distribution[i];
	}

	cout << num_vertices << " " << check << endl;
	return component_size_distribution;
}*/
