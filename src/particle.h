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

extern struct particle* particles;
extern int root_nx;
extern int root_ny;
extern int root_nz;
extern double boxsize;
extern double boxsize_x;
extern double boxsize_y;
extern double boxsize_z;
extern double boxsize_max;

void init_particles();
void init_box();
int get_rootbox_for_particle(struct particle pt);
int get_rootbox_for_particle_int(int pt);

#endif
