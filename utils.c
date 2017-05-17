#include "bh.h"

void print_bounds(t_bounds bounds)
{
	printf("xmin: %f xmax: %f\nymin: %f, ymax: %f\nzmin: %f, zmax: %f\n", \
		bounds.xmin, bounds.xmax, bounds.ymin, bounds.ymax, bounds.zmin, bounds.zmax);
}

void print_vec(cl_float4 v)
{
	printf("x: %f y: %f z: %f w: %f\n", v.x, v.y, v.z, v.w);
}

void print_body(t_body *b)
{
	printf("pos:\n");
	print_vec(b->position);
	printf("vel:\n");
	print_vec(b->velocity);
}

float rand_float(float max)
{
    //generate random float from 0..max
    float r;

    r = (float)rand() / (float)(RAND_MAX/max);
    return r;
}

cl_float4 neg_vec(cl_float4 v)
{
	return (cl_float4){-1 * v.x, -1 * v.y, -1 * v.z, -1 * v.w};
}

cl_float4 vadd(cl_float4 a, cl_float4 b)
{
	return (cl_float4){a.x + b.x, a.y + b.y, a.z + b.z, a.w + b.w};
}

t_body *rand_body(int mag)
{
	t_body *body;

	body = (t_body *)malloc(sizeof(t_body));
    double elevation = asin(rand_float(2) - 1);
    double azimuth = 2 * M_PI * rand_float(1);
    double radius = cbrt(rand_float(1)) * __exp10(mag);
    body->position = (cl_float4){radius * cos(elevation) * cos(azimuth), \
                                radius * cos(elevation) * sin(azimuth), \
                                0.3 * radius * sin(elevation), \
                            	1.0};
    body->velocity = (cl_float4){0, 0, 0, 0};
    return body;
}