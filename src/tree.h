#ifndef _TREE_H
#define _TREE_H
#if defined(GRAVITY_TREE) || defined(COLLISIONS_TREE)
#define TREE

struct cell;

struct cell {
#ifdef GRAVITY_TREE
	double m;
	double mx;
	double my;
	double mz;
#ifdef QUADRUPOLE
	double mxx;
	double mxy;
	double mxz;
	double myy;
	double myz;
	double mzz;
#endif // QUADRUPOLE
#endif // GRAVITY_TREE
	double x;
	double y;
	double z;
	double w;
	struct cell *oct[8];
	int pt;						// has double usages
};

extern struct cell** tree_root;

void tree_init();
void tree_update();
#ifdef GRAVITY_TREE
void tree_update_gravity_data();
#endif // GRAVITY_TREE
void tree_add_particle_to_tree(int pt);
#ifdef MPI
void tree_add_essential_node(struct cell* node);
void tree_prepare_essential_tree();
#endif // MPI


#endif // TREE
#endif // _TREE_H
