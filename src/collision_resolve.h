#ifndef _COLLISION_RESOLVE_H
#define _COLLISION_RESOLVE_H
#include "boundaries.h"

struct collision{
	int p1;
	int p2;
	struct ghostbox gb;
	double time;
	int crossing;
	int ri;	 /**< Index of rootcell */
} collision;

extern double coefficient_of_restitution;
extern double minimum_collision_velocity;
extern double (*coefficient_of_restitution_for_velocity) (double);


void collision_resolve_single(struct collision c);

#endif // _COLLISION_RESOLVE_H

