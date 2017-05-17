#include "bh.h"

void pair_force_cell(t_cell *i, t_cell *j)
{
	//compute the force between two distant cells, treating them as single particles
	t_vector r;

	r.x = j->center.x - i->center.x;
	r.y = j->center.y - i->center.y;
	r.z = j->center.z - i->center.z;

	float distSq = r.x * r.x + r.y * r.y + r.z * r.z + SOFTENING;
	float invDist = 1.0 / sqrt(distSq);
	float invDistCube = invDist * invDist * invDist;
    float f = j->center.w * invDistCube * i->center.w > 0 ? 1 : -1;
    i->force_bias = vadd(i->force_bias, (t_vector){r.x * f, r.y * f, r.z * f});
}

float multipole_acceptance_criterion(t_cell *us, t_cell *them)
{
	//assess whether a cell is "near" or "far" for the sake of barnes-hut
	//if the value returned is less than THETA, that cell is far
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

	/*
		recursively flow through the tree, determining if cells are near or far from
		the cell we're currently considering. We skip the root.
		
		If the cell is far away (m_a_c < THETA), that cell is far enough away to treat as 1 particle.
		we compute that force and add it to the total force on our cell (force_bias).
		
		if the cell is nearby and childless (ie leaf), it is near enough that direct calculation is necessary,
		so it returns a null-terminated array just containing a pointer to the cell.

		if the cell is nearby and has children, we recurse down to its children.
		we make space for the 8 arrays that will be returned (some might be null)
		then we copy them into one final array and return it.

		in this way, we end up with all necessary distant calculations done
		(and the net resulting force stored in cell->force_bias)
		and we have a nice handy list of all the cells whose bodies we'll need to compare against.
	*/

	if (root != t->root && multipole_acceptance_criterion(cell, root) < THETA)
	{
		pair_force_cell(cell, root);
		return (NULL);
	}
	else if (!(root->children))
	{
		ret = (t_cell **)malloc(sizeof(t_cell *) * 2);
		ret[0] = root;
		ret[1] = NULL;
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
		return (ret);
	}
}

t_body **bodies_from_cells(t_cell **cells)
{
	//given a null terminated list of cells, make an array of all the bodies in those cells.
	int count;
	t_body **bodies;

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
	return (bodies);
}

t_ret compute_cell(t_cell *cell, t_octree *t)
{
	t_cell **inners;
	t_body **direct_bodies;

	//skip empty cells (yes there are empty cells)
	if (count_bodies(cell->bodies) == 0)
		return (t_ret){NULL, NULL};
	//traverse tree doing faraway calculations and enumerating nearby cells (see above)
	inners = find_inners_do_outers(cell, t->root, t);
	//use the result to make a list of particles we need to direct compare against
	direct_bodies = bodies_from_cells(inners);
	//GPU will do cell->bodies X direct_bodies with a bias of cell->force_bias, then integrate
	return(gpu_magic(cell->bodies, direct_bodies, cell->force_bias));
}

void update(t_cell *c, t_ret r)
{
	//update our bodies with the results of the gpu calculation
	int count = count_bodies(c->bodies);
	for (int i = 0; i < count; i++)
	{
		c->bodies[i]->position.x = r.P[i].x;
		c->bodies[i]->position.y = r.P[i].y;
		c->bodies[i]->position.z = r.P[i].z;
		c->bodies[i]->velocity.x = r.V[i].x;
		c->bodies[i]->velocity.y = r.V[i].y;
		c->bodies[i]->velocity.z = r.V[i].z;
		c->force_bias = (t_vector){0, 0, 0, 0};
	}
}