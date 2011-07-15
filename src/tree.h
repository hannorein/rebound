#ifndef _TREE_H
#define _TREE_H
struct cell;
struct particle;

struct cell {
	double x;
	double y;
	double z;
	double w;
	struct cell **oct;
	struct particle *pt;
} cell;

extern struct cell** root;
extern int root_nx;
extern int root_ny;
extern int root_nz;
void tree_check_moved_particles();
void tree_update_cell(struct cell *node);
void tree_update();
void tree_init();

//void finalize_tree();

#endif
