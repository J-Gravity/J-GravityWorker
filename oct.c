#include "bh.h"

t_body **demo_bodies(size_t n, int mag)
{
	//whips up a bunch of randomized stars evenly distributed
	//over a sphere of radius 10^mag
	t_body **bodies;

	bodies = (t_body **)malloc(sizeof(t_body *) * (n + 1));
	for (int i = 0; i < n; i++)
		bodies[i] = rand_body(mag);
	bodies[n] = NULL;
	return (bodies);
}

int count_bodies(t_body **bodies)
{
	//returns the number of bodies in a null terminated array of bodies.
	int i;
	i = 0;
	while (bodies[i])
		i++;
	return i;
}

t_vector center_of_mass(t_cell *c, t_body **bodies)
{
	//center of mass (sum of all positions divided by total mass)
	//NB this is currently assuming all masses are 1.
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

int count_leaves(t_cell *root)
{
	//returns number of leaf nodes in the tree
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

int main(void)
{
	//make some dummy bodies and the bounds of the simulation
	t_body **bodies = demo_bodies(BODYCOUNT, BOUNDMAG);
	t_bounds bounds = (t_bounds){-1 * __exp10(BOUNDMAG), __exp10(BOUNDMAG), \
								-1 * __exp10(BOUNDMAG), __exp10(BOUNDMAG), \
								-1 * __exp10(BOUNDMAG), __exp10(BOUNDMAG)};
	//initialize the tree
	t_octree *t = init_tree(bodies, BODYCOUNT, bounds);
	//split the tree
	tree_it_up(t->root);
	//scan the tree for leaf nodes and return them in a list
	t_cell **leaves = enumerate_leaves(t->root);
	//set up space for return values
	t_ret *rets = (t_ret *)malloc(sizeof(t_ret) * (count_cell_array(leaves) + 1));
	//for each leaf node, compute the forces according to Barnes-Hut
	for (int i = 0; leaves[i]; i++)
		rets[i] = compute_cell(leaves[i], t);
	//replace current values with results (ie, step forward)
	for (int i = 0; leaves[i]; i++)
		update(leaves[i], rets[i]);
}