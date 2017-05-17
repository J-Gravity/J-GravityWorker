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
#define BODYCOUNT pow(2, 22)
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

typedef struct s_env
{
    void    *mlx;
    void    *win;
    void    *img;

    t_context *cont;
    t_sim *sim;
    cl_event render;
    cl_event read;
    cl_float4 *stars;
    int flag;
}               t_env;

typedef struct s_vector
{
	float x;
	float y;
	float z;
	float w;
}				t_vector;

typedef struct s_body
{
	t_vector	position;
	t_vector 	velocity;
	t_vector	acceleration;
    int        oct;
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
	t_vector center;
	t_vector force_bias;
	t_bounds bounds;
}				t_cell;

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
t_ret gpu_magic(t_body **N0, t_body **M0, t_vector force_bias);

int count_bodies(t_body **bodies);


void print_bounds(t_bounds bounds);
void print_vec(t_vector v);
void print_body(t_body *b);
float rand_float(float max);
t_vector neg_vec(t_vector v);
t_vector vadd(t_vector a, t_vector b);
t_body *rand_body(int mag);

void pair_force_cell(t_cell *i, t_cell *j);
float multipole_acceptance_criterion(t_cell *us, t_cell *them);
t_cell **find_inners_do_outers(t_cell *cell, t_cell *root, t_octree *t);
t_body **bodies_from_cells(t_cell **cells);
t_ret compute_cell(t_cell *cell, t_octree *t);
void update(t_cell *c, t_ret r);

t_cell *init_cell(t_body **bodies, t_cell *parent, t_bounds bounds);
t_octree *init_tree(t_body **bodies, size_t n, t_bounds bounds);
void    paint_bodies_octants(t_body **bodies, t_cell *c);
t_body ***scoop_octants(t_body **bodies);
void divide_cell(t_cell *c);
void tree_it_up(t_cell *root);
t_cell **enumerate_leaves(t_cell *root);

int count_cell_array(t_cell **cells);
int count_leaves(t_cell *root);
t_vector center_of_mass(t_cell *c, t_body **bodies);
int count_bodies(t_body **bodies);
t_body **demo_bodies(size_t n, int mag);