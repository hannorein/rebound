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
	struct cell *oct[8];
	int pt;							// has double usages
} cell;

extern struct cell** root;
extern int root_nx;
extern int root_ny;
extern int root_nz;

struct cell *add_particle(struct cell *, int, struct cell *,  int);

int get_octant(int n, struct cell *node);

int isInside(struct cell *);

void tree_init();
void tree_update();

struct cell *update_tree(struct cell *node);

//void finalize_tree();

#endif
