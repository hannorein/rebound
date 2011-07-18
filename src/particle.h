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
	struct cell* c;
#ifndef COLLISIONS_NONE
	double r; 
	double lastcollision;
#endif
} particle;

extern struct particle* restrict particles;

void init_particles();

#endif
