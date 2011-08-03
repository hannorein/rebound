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
} ;

extern struct particle* particles;
void particles_add(struct particle pt);
int particles_get_rootbox_for_particle(struct particle pt);

#endif
