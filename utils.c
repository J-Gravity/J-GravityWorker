#include "bh.h"

void print_bounds(t_bounds bounds)
{
	printf("xmin: %f xmax: %f\nymin: %f, ymax: %f\nzmin: %f, zmax: %f\n", \
		bounds.xmin, bounds.xmax, bounds.ymin, bounds.ymax, bounds.zmin, bounds.zmax);
}

void print_vec(t_vector v)
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

t_vector neg_vec(t_vector v)
{
	return (t_vector){-1 * v.x, -1 * v.y, -1 * v.z, -1 * v.w};
}

t_vector vadd(t_vector a, t_vector b)
{
	return (t_vector){a.x + b.x, a.y + b.y, a.z + b.z, a.w + b.w};
}

t_body *rand_body(int mag)
{
	t_body *body;

	body = (t_body *)malloc(sizeof(t_body));
    double elevation = asin(rand_float(2) - 1);
    double azimuth = 2 * M_PI * rand_float(1);
    double radius = cbrt(rand_float(1)) * __exp10(mag);
    body->position = (t_vector){radius * cos(elevation) * cos(azimuth), \
                                radius * cos(elevation) * sin(azimuth), \
                                0.3 * radius * sin(elevation), \
                            	1.0};
    body->velocity = (t_vector){0, 0, 0, 0};
    body->acceleration = (t_vector){0,0,0,0};
    return body;
}