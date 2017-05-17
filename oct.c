#include "bh.h"


#define xmid c->bounds.xmax - (c->bounds.xmax - c->bounds.xmin) / 2
#define ymid c->bounds.ymax - (c->bounds.ymax - c->bounds.ymin) / 2
#define zmid c->bounds.zmax - (c->bounds.zmax - c->bounds.zmin) / 2

//debug output stuff
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

t_body **demo_bodies(size_t n, int mag)
{
	t_body **bodies;

	bodies = (t_body **)malloc(sizeof(t_body *) * (n + 1));
	for (int i = 0; i < n; i++)
		bodies[i] = rand_body(mag);
	bodies[n] = NULL;
	return (bodies);
}

int count_bodies(t_body **bodies)
{
	int i;
	i = 0;
	while (bodies[i])
		i++;
	return i;
}

t_vector vadd(t_vector a, t_vector b)
{
	return (t_vector){a.x + b.x, a.y + b.y, a.z + b.z, a.w + b.w};
}

t_vector center_of_mass(t_cell *c, t_body **bodies)
{
	t_vector v;

	v = (t_vector){0,0,0, count_bodies(bodies)};
	if (v.w == 0)
	{
		return (t_vector){xmid, ymid, zmid, 0};
	}
	for (int i = 0; bodies[i]; i++)
		v = vadd(v, bodies[i]->position);
	return (t_vector){v.x / v.w, v.y / v.w, v.z / v.w, v.w};
}

t_cell *init_cell(t_body **bodies, t_cell *parent, t_bounds bounds)
{
	t_cell *c;

	c = (t_cell *)malloc(sizeof(t_cell));
	c->bodies = bodies;
	c->parent = parent;
	c->children = NULL;
	c->bounds = bounds;
	c->center = center_of_mass(c, bodies);
	c->force_bias = (t_vector){0, 0, 0, 0};
	return (c);
}

t_octree *init_tree(t_body **bodies, size_t n, t_bounds bounds)
{
	t_octree *t;

	t = (t_octree *)malloc(sizeof(t_octree));
	t->bodies = bodies;
	t->n_bodies = n;
	t->bounds = bounds;
	t->root = init_cell(bodies, NULL, bounds);
	return (t);
}

int in_bounds(t_body *body, t_bounds bounds)
{
	if (bounds.xmin <= body->position.x && body->position.x <= bounds.xmax)
		if (bounds.ymin <= body->position.y && body->position.y <= bounds.ymax)
			if (bounds.zmin <= body->position.z && body->position.z <= bounds.zmax)
				return (1);
	return (0);
}

void 	paint_bodies_octants(t_body **bodies, t_cell *c)
{
	for (int i = 0; bodies[i]; i++)
	{
		if (bodies[i]->position.x < xmid)
		{
			if (bodies[i]->position.y < ymid)
			{
				if (bodies[i]->position.z < zmid)
					bodies[i]->oct = 0;
				else
					bodies[i]->oct = 1;
			}
			else
			{
				if (bodies[i]->position.z < zmid)
					bodies[i]->oct = 2;
				else
					bodies[i]->oct = 3;
			}
		}
		else
		{
			if (bodies[i]->position.y < ymid)
			{
				if (bodies[i]->position.z < zmid)
					bodies[i]->oct = 6;
				else
					bodies[i]->oct = 7;
			}
			else
			{
				if (bodies[i]->position.z < zmid)
					bodies[i]->oct = 4;
				else
					bodies[i]->oct = 5;
			}
		}
	}
}

t_body **scoop_octant(t_body **bodies, int octant)
{
	int count = count_bodies(bodies);
	t_body **ret = (t_body **)malloc(sizeof(t_body *) * (count + 1));
	ret[count] = NULL;
	int j = 0;
	for (int i = 0; bodies[i]; i++)
		if (bodies[i]->oct == octant)
			ret[j++] = bodies[i];
	ret[j] = NULL;
	return (ret);
}

t_body ***scoop_octants(t_body **bodies)
{
	t_body ***ret = (t_body ***)malloc(sizeof(t_body **) * 9);
	ret[8] = NULL;
	int count = count_bodies(bodies);
	for (int i = 0; i < 8; i++)
	{
		ret[i] = (t_body **)malloc(sizeof(t_body *) * (count + 1));
		ret[i][count] = NULL;
	}
	int indices[8] = {0, 0, 0, 0, 0, 0, 0, 0};
	for (int i = 0; bodies[i]; i++)
	{
		ret[bodies[i]->oct][indices[bodies[i]->oct]] = bodies[i];
		indices[bodies[i]->oct] += 1;
	}
	int sum = 0;
	for (int i = 0; i < 8; i++)
	{
		ret[i][indices[i]] = NULL;
		sum += indices[i];
	}
	return (ret);
}

t_body **select_bodies_by_bounds(t_body **bodies, t_bounds bounds)
{
	int count;

	count = count_bodies(bodies);
	t_body **contained = (t_body **)malloc(sizeof(t_body *) * (count + 1));
	contained[count] = NULL;
	int j = 0;
	for (int i = 0; bodies[i]; i++)
	{
		if (in_bounds(bodies[i], bounds))
			contained[j++] = bodies[i];
	}
	contained[j] = NULL;
	return contained;
}

void divide_cell(t_cell *c)
{
	t_bounds subbounds[8];
	//min->max is left->right, bottom->top, near->far
	subbounds[0] = (t_bounds){c->bounds.xmin, xmid, \
								c->bounds.ymin, ymid, \
								c->bounds.zmin, zmid}; //bottom left near

	subbounds[6] = (t_bounds){xmid, c->bounds.xmax, \
								c->bounds.ymin, ymid, \
								c->bounds.zmin, zmid}; //bottom right near
	subbounds[2] = (t_bounds){c->bounds.xmin, xmid, \
								ymid, c->bounds.ymax, \
								c->bounds.zmin, zmid}; // top left near
	subbounds[1] = (t_bounds){c->bounds.xmin, xmid, \
								c->bounds.ymin, ymid, \
								zmid, c->bounds.zmax}; // bottom left far

	subbounds[3] = (t_bounds){c->bounds.xmin, xmid, \
								ymid, c->bounds.ymax, \
								zmid, c->bounds.zmax}; //top left far
	subbounds[7] = (t_bounds){xmid, c->bounds.xmax, \
								c->bounds.ymin, ymid, \
								zmid, c->bounds.zmax}; //bottom right far
	subbounds[4] = (t_bounds){xmid, c->bounds.xmax, \
								ymid, c->bounds.ymax, \
								c->bounds.zmin, zmid}; //right top near

	subbounds[5] = (t_bounds){xmid, c->bounds.xmax, \
								ymid, c->bounds.ymax, \
								zmid, c->bounds.zmax}; //right top far

	//////Numbering is a bit weird here ^ but the idea is that even indices are near, odd far
								//0..3 are left, 4..7 are right
								//0, 1, 6, 7 bottom, 2345 top
	t_cell **children = (t_cell **)malloc(sizeof(t_cell *) * 9);
	paint_bodies_octants(c->bodies, c);
	t_body ***kids = scoop_octants(c->bodies);
	for (int i = 0; i < 8; i++)
		children[i] = init_cell(kids[i], c, subbounds[i]);
	children[8] = NULL;
	c->children = children;
}

t_body *simple_body(float x, float y, float z)
{
	t_body *b;

	b = (t_body *)malloc(sizeof(t_body));
	b->position = (t_vector){x, y, z, 1};
	b->velocity = (t_vector){0, 0, 0, 0};
	return b;
}

void tree_it_up(t_cell *root)
{
	if (!root)
		return ;
	if (count_bodies(root->bodies) < LEAF_THRESHOLD)
		return ;
	divide_cell(root);
	for (int i = 0; i < 8; i++)
		tree_it_up(root->children[i]);
}

int count_leaves(t_cell *root)
{
	int sum;

	sum = 0;
	if (root->children == NULL)
		return (1);
	for (int i = 0; i < 8; i++)
		sum += count_leaves(root->children[i]);
	return sum;
}

int count_cell_array(t_cell **cells)
{
	//returns the length of a null terminated array of cell pointers
	int i = 0;
	if (!cells)
		return 0;
	while (cells[i])
		i++;
	return i;
}

t_cell **enumerate_leaves(t_cell *root)
{
	//goal is to return a linear t_cell** that's all the leaf nodes

	t_cell **ret;

	if (!root->children)
	{
		ret = (t_cell **)malloc(sizeof(t_cell *) * 2);
		ret[0] = root;
		ret[1] = NULL;
		return ret;
	}
	t_cell ***returned = (t_cell ***)malloc(sizeof(t_cell **) * 8);
	int total = 0;
	for (int i = 0; i < 8; i++)
	{
		returned[i] = enumerate_leaves(root->children[i]);
		total += count_cell_array(returned[i]);
	}
	ret = (t_cell **)malloc(sizeof(t_cell *) * (total + 1));
	for (int i = 0; i < total;)
	{
		for (int j = 0; j < 8; j++)
		{
			for (int k = 0; returned[j][k]; k++, i++)
			{
				ret[i] = returned[j][k];
			}
			free(returned[j]);
		}
		free(returned);
	}
	ret[total] = NULL;
	return (ret);
}

t_vector neg_vec(t_vector v)
{
	return (t_vector){-1 * v.x, -1 * v.y, -1 * v.z, -1 * v.w};
}

void pair_force_cell(t_cell *i, t_cell *j)
{
	t_vector r;
//needs negatives for janus
	r.x = j->center.x - i->center.x;
	r.y = j->center.y - i->center.y;
	r.z = j->center.z - i->center.z;

	float distSq = r.x * r.x + r.y * r.y + r.z * r.z + SOFTENING;
	float invDist = 1.0 / sqrt(distSq);
	float invDistCube = invDist * invDist * invDist;
    float f = j->center.w * invDistCube;
    i->force_bias = vadd(i->force_bias, (t_vector){r.x * f, r.y * f, r.z * f});
}

//above is oversimplified, you have to compare directly with neighbors, what i'm calling N x (N + M)
//need function to get neighborhood
//apparently you can just compare to the peers (ie same parent).

float multipole_acceptance_criterion(t_cell *us, t_cell *them)
{
	float s;
	float d;
	t_vector r;

	s = them->bounds.xmax - them->bounds.xmin;
	r.x = them->center.x - us->center.x;
	r.y = them->center.y - us->center.y;
	r.z = them->center.z - us->center.z;
	d = sqrt(r.x * r.x + r.y * r.y + r.z * r.z);
	if (d == 0)
		return (0);
	return (s/d);
}

t_cell **find_inners_do_outers(t_cell *cell, t_cell *root, t_octree *t)
{
	t_cell **ret;
	t_cell ***returned;

	if (root != t->root && multipole_acceptance_criterion(cell, root) < THETA)
	{
		//printf("cell is far\n");
		pair_force_cell(cell, root);
		return (NULL);
	}
	else if (!(root->children))
	{
		ret = (t_cell **)malloc(sizeof(t_cell *) * 2);
		ret[0] = root;
		ret[1] = NULL;
		//printf("childless cell that's near\n");
		return (ret);
	}
	else
	{
		returned = (t_cell ***)malloc(sizeof(t_cell **) * 9);
		returned[8] = NULL;
		int total = 0;
		for (int i = 0; i < 8; i++)
		{
			returned[i] = find_inners_do_outers(cell, root->children[i], t);
			total += count_cell_array(returned[i]);
		}
		ret = (t_cell **)malloc(sizeof(t_cell *) * (total + 1));
		for (int i = 0; i < total;)
		{
			for (int j = 0; j < 8; j++)
			{
				if (!returned[j])
					continue ;
				for (int k = 0; returned[j][k]; k++, i++)
				{
					ret[i] = returned[j][k];
				}
				free(returned[j]);
			}
			free(returned);
		}
		ret[total] = NULL;
		//printf("returning the findings\n");
		return (ret);
	}
}

t_body **bodies_from_cells(t_cell **cells)
{
	int count;
	t_body **bodies;

	//printf("lining up bodies\n");
	count = 0;
	for (int i = 0; cells[i]; i++)
		count += count_bodies(cells[i]->bodies);
	bodies = (t_body **)malloc(sizeof(t_body *) * (count + 1));
	bodies[count] = NULL;
	int k = 0;
	for (int i = 0; cells[i]; i++)
	{
		for (int j = 0; cells[i]->bodies[j]; j++, k++)
			bodies[k] = cells[i]->bodies[j];
	}
	//printf("returning bodies\n");
	return (bodies);
}

t_ret compute_cell(t_cell *cell, t_octree *t)
{
	t_cell **inners;
	t_body **direct_bodies;

	if (count_bodies(cell->bodies) == 0)
		return (t_ret){NULL, NULL};
	inners = find_inners_do_outers(cell, t->root, t);
	//printf("found inners\n");
	direct_bodies = bodies_from_cells(inners);
	//printf("lined em up\n");
	//GPU will do cell->bodies X direct_bodies with a bias of cell->force_bias, then integrate
	return(gpu_magic(cell->bodies, direct_bodies, cell->force_bias));
}

void update(t_cell *c, t_ret r)
{
	int count = count_bodies(c->bodies);
	for (int i = 0; i < count; i++)
	{
		// printf("body %d was at %f %f %f and is now at %f %f %f, %f\nmessage %f\n", i, \
		// 	c->bodies[i]->position.x, c->bodies[i]->position.y, c->bodies[i]->position.z, \
		// 	r.P[i].x, r.P[i].y, r.P[i].z, r.P[i].w, r.V[i].w);
		c->bodies[i]->position.x = r.P[i].x;
		c->bodies[i]->position.y = r.P[i].y;
		c->bodies[i]->position.z = r.P[i].z;
		c->bodies[i]->velocity.x = r.V[i].x;
		c->bodies[i]->velocity.y = r.V[i].y;
		c->bodies[i]->velocity.z = r.V[i].z;
		c->force_bias = (t_vector){0, 0, 0, 0};
	}
}

int main(void)
{
	//make some dummy bodies
	t_body **bodies;
	clock_t start = clock();

	bodies = demo_bodies(BODYCOUNT, BOUNDMAG);
	clock_t bods = clock();
	t_bounds bounds = (t_bounds){-1 * __exp10(BOUNDMAG), __exp10(BOUNDMAG), \
								-1 * __exp10(BOUNDMAG), __exp10(BOUNDMAG), \
								-1 * __exp10(BOUNDMAG), __exp10(BOUNDMAG)};
	t_octree *t = init_tree(bodies, BODYCOUNT, bounds);
	tree_it_up(t->root);
	clock_t tree = clock();
	t_cell **leaves = enumerate_leaves(t->root);
	clock_t enumerate = clock();
	printf("making bodies took %lu, tree-ing was %lu, enumerate was %lu, total %lu\n", \
			bods - start, tree - bods, enumerate - tree, enumerate - start);
	t_ret *rets = (t_ret *)malloc(sizeof(t_ret) * (count_cell_array(leaves) + 1));
	//for (int i = 0; leaves[i]; i++)
		rets[8] = compute_cell(leaves[8], t);
	clock_t crunch = clock();
	//for (int i = 0; leaves[i]; i++)
		update(leaves[8], rets[8]);
	clock_t integrate = clock();
	printf("gpu phase took %lu, setting new values took %lu total %lu\n", crunch - enumerate, integrate - crunch, integrate - start);
}


/*
to do list:

add visualizer
keep fingers crossed
nXm probably needs to become mXn for higher threadcount ? it's fast already, it seems?
mcount is still really high really frequently.
for theta = 2, M seems reasonably sized.

need better ways to guarantee good chunk sizes, I don't like bothering to send 100 X 18000 to the gpu
even when leaf threshold is 32768, most cells are like 5000-6000

some exciting and complex work to be done condensing better workunits out of nXm sets
for two nearby cells n1 n2, their ms likely overlap substantially, and they guaranteed
will contain each other's bodies (ie n2 is in m1, and n1 is in m2)
comparing all the ms is ~(average mlen)^2, but i think they can be compared smartly (eg with peers)

should probably free things ever
next step is verification. torn between just plugging in display stuff and hoping, vs doing step by step verification first

after that
gpu scheduling (put all batches on, then clfinish, then read)
investigation - how many comparisons are we doing now? can measure algorithmic improvement and gpu utilization
multithreading (save for actual server code maybe)
gpu tree generation, workset enumeration, etc

48 billion comparisons for 2^20 (~1 trillion brute force) at theta = 2
84 billion for same at theta = 1.5. total runtime did not double, gpu likely staying busier
so we're really underutilizing the gpu currently (~50%), but also any worry about redundant calculations is silly for now
coalescing nXm calls does sound awesome though so i'll still look into it

*/