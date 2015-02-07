/**
 * @file 	problem.c
 * @brief 	Example problem: Mergers.
 * @author 	Hanno Rein <hanno@hanno-rein.de>
 * @detail	This example integrates a densly packed planetary system 
 * which becomes unstable on a timescale of only a few orbits. The IAS15 
 * integrator with adaptive timestepping is used. The bodies have a finite
 * size and merge if they collide. Note that the size is unphysically large
 * in this example. 
 * 
 * @section 	LICENSE
 * Copyright (c) 2014 Hanno Rein, Dave Spiegel
 *
 * This file is part of rebound.
 *
 * rebound is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * rebound is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with rebound.  If not, see <http://www.gnu.org/licenses/>.
 *
 */
#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <math.h>
#include <time.h>
#include "main.h"
#include "tools.h"
#include "output.h"
#include "particle.h"
#include "integrator.h"
#include "collision_resolve.h"

#ifdef OPENGL
extern int display_wire;
#endif // OPENGL

void collision_resolve_merger(struct collision c);

void problem_init(int argc, char* argv[]){
	dt = 0.01*2.*M_PI;						// initial timestep

#ifdef OPENGL
	display_wire	= 1;						// show instantaneous orbits
#endif // OPENGL
	collision_resolve = collision_resolve_merger;			// Setup our own collision routine.
	init_boxwidth(10); 					

	struct particle star;
	star.m = 1;
	star.r = 0.1;
	star.x = 0; 	star.y = 0; 	star.z = 0;
	star.vx = 0; 	star.vy = 0; 	star.vz = 0;
	particles_add(star);
	
	// Add planets
	int N_planets = 7;
	for (int i=0;i<N_planets;i++){
		double a = 1.+(double)i/(double)(N_planets-1);		// semi major axis in AU
		double v = sqrt(1./a); 					// velocity (circular orbit)
		struct particle planet;
		planet.m = 1e-4; 
		planet.r = 4e-2; 					// radius in AU (it is unphysically large in this example)
		planet.lastcollision = 0;				// The first time particles can collide with each other
		planet.x = a; 	planet.y = 0; 	planet.z = 0;
		planet.vx = 0; 	planet.vy = v; 	planet.vz = 0;
		particles_add(planet); 
	}
	tools_move_to_center_of_momentum();				// This makes sure the planetary systems stays within the computational domain and doesn't drift.
}

void collision_resolve_merger(struct collision c){
	struct particle p1 = particles[c.p1];
	struct particle p2 = particles[c.p2];
	double x21  = p1.x  - p2.x; 
	double y21  = p1.y  - p2.y; 
	double z21  = p1.z  - p2.z; 
	double rp   = p1.r+p2.r;
	printf("collision 1\n");
	if (rp*rp < x21*x21 + y21*y21 + z21*z21) return; // not overlapping
	double vx21 = p1.vx - p2.vx; 
	double vy21 = p1.vy - p2.vy; 
	double vz21 = p1.vz - p2.vz; 
	printf("collision 2\n");
	if (vx21*x21 + vy21*y21 + vz21*z21 >0) return; // not approaching
	printf("collision 3 %f %f %f\n", p1.lastcollision, p2.lastcollision, t);
	if (p1.lastcollision>=t || p2.lastcollision>=t) return; // already collided
	printf("collision 4\n");
	particles[c.p2].lastcollision = t;
	particles[c.p1].lastcollision = t;
	// Note: We assume only one collision per timestep. 
	// Setup new particle (in position of particle p1. Particle p2 will be discarded.
	struct particle cm = tools_get_center_of_mass(p1, p2);
	particles[c.p1].x = cm.x;
	particles[c.p1].y = cm.y;
	particles[c.p1].z = cm.z;
	particles[c.p1].vx = cm.vx;
	particles[c.p1].vy = cm.vy;
	particles[c.p1].vz = cm.vz;
	particles[c.p1].r = p1.r*pow(cm.m/p1.m,1./3.);	// Assume a constant density
	particles[c.p1].m = cm.m;
	// Remove one particle.
	N--;
	particles[c.p2] = particles[N];
	// Make sure we don't drift, so let's go back to the center of momentum
	tools_move_to_center_of_momentum();	
}

void problem_output(){
	if (output_check(10.*2.*M_PI)){  
		output_timing();
	}
}

void problem_finish(){
}
