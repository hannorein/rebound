#ifndef _TREE_H
#define _TREE_H
struct cell;
struct particle;

struct cell {
	double x[3];
	double w;
	struct cell **oct;
	struct particle *pt;
} cell;

extern struct cell* root;
void tree_check_moved_particles();
void tree_update(struct cell *node);
void tree_init();

//void finalize_tree();

#endif
