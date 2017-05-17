#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <sys/types.h>
#ifdef __APPLE__
#include <OpenCL/opencl.h>
#include <unistd.h>
#else
#include <CL/cl.h>
#endif
#include <time.h>
#include <stdlib.h>


#define LEAF_THRESHOLD pow(2, 14)
#define BODYCOUNT pow(2, 20)
#define BOUNDMAG 10
#define G 1.327 * __exp10(13) //kilometers, solar masses, (km/s)^2
#define SOFTENING 1000000
#define TIME_STEP 1

#define THETA 1.5 //once theta is 1 we always compare peers (neat!)

#define STARCOUNT pow(2, 11)
#define THREADCOUNT pow(2, 11)
#define GROUPSIZE 256
#define SOLAR_MASS 1
#define BIG_RADIUS 12
#define ANCHOR_MASS __exp10(8.5)
#define WINDOW_FACTOR __exp10(BIG_RADIUS - 2) * 4
#define WINDOW_DIM 1000 
#define FRAMECOUNT 240

# ifndef DEVICE
# define DEVICE CL_DEVICE_TYPE_DEFAULT
# endif

#define xmid c->bounds.xmax - (c->bounds.xmax - c->bounds.xmin) / 2
#define ymid c->bounds.ymax - (c->bounds.ymax - c->bounds.ymin) / 2
#define zmid c->bounds.zmax - (c->bounds.zmax - c->bounds.zmin) / 2


typedef struct s_context
{
    cl_device_id device_id;
    cl_context context;
    cl_command_queue commands;
}               t_context;

typedef struct s_sim
{
    cl_mem      d_p_start;
    cl_mem      d_a;
    cl_mem      d_v_start;
    cl_mem      d_v_end;
    cl_mem      d_p_end;
    cl_kernel   kernel;

    size_t global;
    size_t local;
    size_t soften;
    size_t timestep;
    float grav;
    size_t count;
}               t_sim;

typedef struct s_body
{
	cl_float4	position;
	cl_float4 	velocity;
}				t_body;

typedef struct s_bounds
{
	float xmin;
	float xmax;
	float ymin;
	float ymax;
	float zmin;
	float zmax;
}				t_bounds;

typedef struct s_cell
{
	t_body **bodies;
	struct s_cell *parent;
	struct s_cell **children;
	cl_float4 center;
	cl_float4 force_bias;
	t_bounds bounds;
}				t_cell;

typedef struct s_workunit
{
	int			id;
	int			localcount;
	int			neighborcount;
	t_body		*local_bodies;
	t_body		*neighborhood;
	cl_float4	force_bias;
}				t_workunit;

typedef struct s_resultunit
{
	int		identifier;
	int		localcount;
	t_body	*local_bodies;
}				t_resultunit;

typedef struct s_octree
{
	t_cell *root;
	t_body **bodies;
	size_t n_bodies;
	t_bounds bounds;
}				t_octree;

typedef struct s_ret
{
    cl_float4 *P;
    cl_float4 *V;
}               t_ret;


t_ret crunch_NxM(cl_float4 *N, cl_float4 *M, cl_float4 *V, size_t ncount, size_t mcount, cl_float4 force_bias);
t_ret gpu_magic(t_body **N0, t_body **M0, cl_float4 force_bias);

int count_bodies(t_body **bodies);


void print_bounds(t_bounds bounds);
void print_vec(cl_float4 v);
void print_body(t_body *b);
float rand_float(float max);
cl_float4 neg_vec(cl_float4 v);
cl_float4 vadd(cl_float4 a, cl_float4 b);
t_body *rand_body(int mag);

void pair_force_cell(t_cell *i, t_cell *j);
float multipole_acceptance_criterion(t_cell *us, t_cell *them);
t_cell **find_inners_do_outers(t_cell *cell, t_cell *root, t_octree *t);
t_body **bodies_from_cells(t_cell **cells);
t_workunit *make_workunit_for_cell(t_cell *cell, t_octree *t, int index);
void update(t_cell *c, t_resultunit *r);
t_workunit *new_workunit(t_body **local, t_body **neighborhood, cl_float4 force_bias, int index);

t_cell *init_cell(t_body **bodies, t_cell *parent, t_bounds bounds);
t_octree *init_tree(t_body **bodies, size_t n, t_bounds bounds);
void    paint_bodies_octants(t_body **bodies, t_cell *c);
t_body ***scoop_octants(t_body **bodies);
void divide_cell(t_cell *c);
void tree_it_up(t_cell *root);
t_cell **enumerate_leaves(t_cell *root);

int count_cell_array(t_cell **cells);
int count_leaves(t_cell *root);
cl_float4 center_of_mass(t_cell *c, t_body **bodies);
int count_bodies(t_body **bodies);
t_body **demo_bodies(size_t n, int mag);