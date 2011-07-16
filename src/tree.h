#ifndef _TREE_H
#define _TREE_H
struct cell;

struct cell {
	double m;
	double mx[3];
	double x[3];
	double w;
	struct cell *oct[8];
	int pt;							// has double usages
} cell;

extern struct cell* root;

struct cell *add_particle(struct cell *, int, struct cell *,  int);

int get_octant(int n, struct cell *node);

int isInside(struct cell *);

void init_tree();

struct cell *update_tree(struct cell *node);

//void finalize_tree();

#endif
