#ifndef _PARTICLE_H
#define _PARTICLE_H

struct cell;

struct particle {
	double x;
	double y;
	double z;
	double vx;
	double vy;
	double vz;
	double ovx;
	double ovy;
	double ovz;
	double ax;
	double ay;
	double az;
	double m;
#ifndef COLLISIONS_NONE
	double r; 
	double lastcollision;
#endif
#if defined(GRAVITY_TREE) || defined(COLLISIONS_TREE)
	struct cell* c;
#endif
} particle;

extern struct particle* restrict particles;

void init_particles();

#endif
