#ifndef _PARTICLE_H
#define _PARTICLE_H

struct particle {
	double x;
	double y;
	double z;
	double vx;
	double vy;
	double vz;
	double ax;
	double ay;
	double az;
	double m;
} particle;

extern struct particle* particles;

void init_particles();

#endif
