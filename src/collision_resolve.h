#ifndef _COLLISION_RESOLVE_H
#define _COLLISION_RESOLVE_H
#include "boundaries.h"

struct collision{
	int p1;
	int p2;
	struct ghostbox gb;
	double time;
	int crossing;
} collision;

extern double coefficient_of_restitution;
extern double minimum_collision_velocity;

void collisions_resolve_single(struct collision c);

#endif // _COLLISION_RESOLVE_H

