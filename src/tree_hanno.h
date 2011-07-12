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
void update_tree1();
void update_tree(struct cell *node);
void init_tree();

//void finalize_tree();

#endif
