/**
 * Colliding and merging planets
 * 
 * This example integrates a densely packed planetary system 
 * which becomes unstable on a timescale of only a few orbits. The IAS15 
 * integrator with adaptive timestepping is used. The bodies have a finite
 * size and merge if they collide. Note that the size is unphysically large
 * in this example. 
 */
#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <math.h>
#include "rebound.h"

void collision_resolve_merger(struct reb_simulation* r, struct reb_collision c);
void heartbeat(struct reb_simulation* r);

int main(int argc, char* argv[]){
	struct reb_simulation* r = reb_create_simulation();
	r->dt 			= 0.01*2.*M_PI;				// initial timestep
	r->integrator 		= REB_INTEGRATOR_IAS15;
	r->collision		= REB_COLLISION_DIRECT;
	r->collision_resolve 	= collision_resolve_merger;		// Setup our own collision routine.
	r->heartbeat		= heartbeat;

	struct reb_particle star = {0};
	star.m = 1;
	star.r = 0.1;
	reb_add(r, star);
	
	// Add planets
	int N_planets = 7;
	for (int i=0;i<N_planets;i++){
		double a = 1.+(double)i/(double)(N_planets-1);		// semi major axis in AU
		double v = sqrt(1./a); 					// velocity (circular orbit)
		struct reb_particle planet = {0};
		planet.m = 1e-4; 
		planet.r = 4e-2; 					// radius in AU (it is unphysically large in this example)
		planet.lastcollision = 0;				// The first time particles can collide with each other
		planet.x = a; 
		planet.vy = v;
		reb_add(r, planet); 
	}
	reb_move_to_com(r);				// This makes sure the planetary systems stays within the computational domain and doesn't drift.

	reb_integrate(r, INFINITY);
}

void collision_resolve_merger(struct reb_simulation* r, struct reb_collision c){
	struct reb_particle* particles = r->particles;
	struct reb_particle p1 = particles[c.p1];
	struct reb_particle p2 = particles[c.p2];
	double x21  = p1.x  - p2.x; 
	double y21  = p1.y  - p2.y; 
	double z21  = p1.z  - p2.z; 
	double rp   = p1.r+p2.r;
	if (rp*rp < x21*x21 + y21*y21 + z21*z21) return; // not overlapping
	double vx21 = p1.vx - p2.vx; 
	double vy21 = p1.vy - p2.vy; 
	double vz21 = p1.vz - p2.vz; 
	if (vx21*x21 + vy21*y21 + vz21*z21 >0) return; // not approaching
	if (p1.lastcollision>=r->t || p2.lastcollision>=r->t) return; // already collided
	printf("\nCollision detected. Merging.\n");
	particles[c.p2].lastcollision = r->t;
	particles[c.p1].lastcollision = r->t;
	// Note: We assume only one collision per timestep. 
	// Setup new particle (in position of particle p1. Particle p2 will be discarded.
	struct reb_particle cm = reb_get_com_of_pair(p1, p2);
	particles[c.p1].x = cm.x;
	particles[c.p1].y = cm.y;
	particles[c.p1].z = cm.z;
	particles[c.p1].vx = cm.vx;
	particles[c.p1].vy = cm.vy;
	particles[c.p1].vz = cm.vz;
	particles[c.p1].r = p1.r*pow(cm.m/p1.m,1./3.);	// Assume a constant density
	particles[c.p1].m = cm.m;
	// Remove one particle.
	r->N--;
	particles[c.p2] = particles[r->N];
	// Make sure we don't drift, so let's go back to the center of momentum
	reb_move_to_com(r);	
}

void heartbeat(struct reb_simulation* r){
	if (reb_output_check(r, 10.*2.*M_PI)){  
		reb_output_timing(r, 0);
	}
}
