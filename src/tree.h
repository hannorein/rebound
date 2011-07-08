#ifndef _TREE_H
#define _TREE_H
struct cell;

struct cell {
	double m;
	double mx;
	double my;
	double mz;
	double x;
	double y;
	double z;
	double w;
//	particle* 
	struct cell *oct[8];
//	int	nchildren;
	int pt;							//number of particles
//	int lv;
	int num;
} cell;

extern struct cell* root;

struct cell *add_leaf(struct cell *, int, struct cell *,  int);

int get_octant(int n, struct cell *node);


void init_tree();

//void refresh_tree();

//void finalize_tree();

#endif
