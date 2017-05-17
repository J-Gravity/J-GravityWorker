#include "bh.h"

t_cell *init_cell(t_body **bodies, t_cell *parent, t_bounds bounds)
{
	//allocates and sets up a cell struct.
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
	//allocates and sets up a tree and creates its root cell
	t_octree *t;

	t = (t_octree *)malloc(sizeof(t_octree));
	t->bodies = bodies;
	t->n_bodies = n;
	t->bounds = bounds;
	t->root = init_cell(bodies, NULL, bounds);
	return (t);
}

void 	paint_bodies_octants(t_body **bodies, t_cell *c)
{
	//mark each body with a value for which octant it will be in after cell is split
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

t_body ***scoop_octants(t_body **bodies)
{
	//return 8 arrays of bodies, one for each octant (used after paint_octants)
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

void divide_cell(t_cell *c)
{
	//split a cell into 8 octants.

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

void tree_it_up(t_cell *root)
{
	//recursively flesh out the barnes-hut tree from the root node.
	if (!root)
		return ;
	if (count_bodies(root->bodies) < LEAF_THRESHOLD)
		return ;
	divide_cell(root);
	for (int i = 0; i < 8; i++)
		tree_it_up(root->children[i]);
}

t_cell **enumerate_leaves(t_cell *root)
{
	//goal is to return a linear t_cell** that's all the leaf nodes in the tree

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