#ifndef _TREE_H
#define _TREE_H
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
	#endif
#endif
	double x;
	double y;
	double z;
	double w;
	struct cell *oct[8];
	int pt;						// has double usages
} cell;

extern struct cell** root;
extern int root_nx;
extern int root_ny;
extern int root_nz;

void tree_init();
void tree_update();

#endif
