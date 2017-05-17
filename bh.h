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


#define LEAF_THRESHOLD pow(2, 10)
#define BODYCOUNT pow(2, 15)
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