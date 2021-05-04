/* PIEDOCK                                                            *
 * Authors: D. V. S. Ravikant                                         *
 * (C) 2011 D. V. S. Ravikant and Ron Elber                           */

#include <iostream>
#include <vector>
#include <math.h>
#include <ext/hash_map>
#include <queue>
#include <stack>
#include <list>
#include <string.h>

#ifndef EPSILON
#define EPSILON 0.0000000000000001
#endif
#ifndef PI
#define PI 3.14159265358
#endif

#define ABS(x) (((x) < 0) ? -(x) : (x))

#define minimum(x , y) (( x ) > ( y ) ? (y) : (x))
#define maximum(x , y) (( x ) < ( y ) ? (y) : (x))

#define stdext __gnu_cxx;
using namespace std;
using namespace stdext;

struct eqint
{
  bool operator()(const int i1, const int i2) const
  {
    return i1 == i2;
  }
};

struct eqlong
{
  bool operator()(const long l1, const long l2) const
  {
    return l1 == l2;
  }
};

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

class Vector{
	public:
		float x,y,z;
	
	Vector();

	Vector(float x, float y,float z);

	Vector(Vector* v);
	
	void toarray(float[3]);

	float norm_squared();
	
	float norm();

	static float distance(Vector* v1, Vector* v2);
	
	static float distance(Vector v1, Vector v2);

	static float distance_squared(Vector* v1, Vector* v2);
	
	static float distance_squared(Vector v1, Vector v2);
	
	void normalize();

	float dot(Vector* v);

	float dot(Vector v);

	Vector cross(Vector* v);
	
	Vector cross(Vector v);

	Vector operator+(const Vector &v);

	Vector operator-(const Vector &v);

	Vector operator*(const float s);
	
	bool operator==(const Vector &v);

	static float angle(Vector v1, Vector v2, Vector axis);
	
	static Vector* find_closest(vector<Vector*> *,Vector*);
	
	static Vector* find_closest(vector<Vector*> *,vector<Vector*> *,Vector*);
	
};

//template< class T >
class Lnode {
  public:
    int data;
    Lnode *next;
    
    Lnode();
    ~Lnode();
};

//template< class T >
class NodeList {
  public:
  	Lnode* head;
  
    NodeList();
    ~NodeList();
    void add( int );
    bool empty();
    int size();
    //T    remove();
};

class GVertex{
	public:
	long id;
	void *data;
	hash_map<long, GVertex*, hash<long>, eqlong> neighbors;
	//hash_map<long, float, hash<long>, eqlong> neighbor_weights;
	
	GVertex();
	
	GVertex(long id);
	
	void add_neighbor(GVertex* u);
	
	int get_degree();
	
	bool is_neighbor(GVertex* u);
};	

class GEdge{
	public:
	GVertex *s,*d;
	void *data;

	GEdge(GVertex *s, GVertex *t);
};

class Graph{
	// adjacency list representation
	public:
	long num_vertices;
	hash_map<long, GVertex*, hash<long>, eqlong> vertex;

	Graph();
	
	void insert_vertex(GVertex* v);
	
	//void get_degree_distribution();
	
	int bfs();
	
	int dfs();
	
	int bfs(long start_vertex, bool* visited, long* parent);
	
	int dfs(long start_vertex, bool* visited, long* parent);
	
	void compute_common_neighbors(list<long>*,list<long>*);
	
	void compute_vertex_cycle(list<long>*);
	
	void compute_vertex_cycle(list<long>*, long start_vertex);
	
	void compute_vertex_cycle(list<long> *trail, long start_vertex, bool* visited);
	
	void compute_connected_component(list<long>*, long start_vertex);
	
	void print_details();
};

void permute(int, int*);
